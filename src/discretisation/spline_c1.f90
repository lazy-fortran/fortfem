module fortfem_spline_c1
    !! C1 extraction of a multipatch spline space on bilinear blocks: the
    !! Argyris-type isogeometric space of Kapl, Sangalli, Takacs
    !! (arXiv:1711.05161, CAGD 2019; survey arXiv:1812.09088, Sec. 4.2),
    !! implemented from the papers. The space is the direct sum of
    !! patch-interior, edge and vertex functions. Each function is a sparse
    !! combination of patch B-splines and is returned as rows of E
    !! (doc/design/multipatch_spline_spaces.md, eqs. 11-13).
    !!
    !! Every side carries a frame xi = o + P (sigma, tau): sigma is the
    !! distance from the side, tau runs along it, and eps = det P. Gluing data
    !! are alpha = eps det(d_sigma F, d_tau F) > 0 and, on an interface with
    !! the eps = +1 side L and the eps = -1 side R,
    !! beta = det(d_sigma F_R, d_sigma F_L) = alpha_L beta_R + alpha_R beta_L
    !! with linear beta_L, beta_R (minimal norm split).
    !!
    !! Interior sides flagged side_c0 are coupled in C0: the trace functions
    !! (kind 0) are shared, the transversal functions (kind 1) belong to each
    !! side, and both sides use the boundary gluing data. At a vertex the
    !! blocks joined by C1 sides form sectors. Each sector carries C2 data
    !! (value, gradient, Hessian), and across a C0 side the value and the
    !! first and second derivatives along the side agree, so the vertex
    !! functions span the null space of these conditions.
    use fortfem_kinds, only: dp
    use fortfem_spline_layout, only: block_layout_t, side_start, side_end, &
        opposite_side
    use fortfem_spline_extraction_1d, only: build_open_knot_vector, &
        build_bezier_extraction_1d
    use fortfem_interval_bernstein_jets, only: interval_bernstein_jet
    implicit none
    private

    type, public :: c1_block_knots_t
        real(dp), allocatable :: knots1(:), knots2(:)
    end type c1_block_knots_t

    type :: uspace_t
        integer :: degree = 0
        integer :: n = 0
        real(dp), allocatable :: knots(:), breaks(:), extraction(:, :, :)
        integer, allocatable :: first(:)
    end type uspace_t

    type :: gluing_t
        real(dp) :: alpha(2) = [1.0_dp, 0.0_dp] ! alpha(tau) = alpha(1) + alpha(2) tau
        real(dp) :: beta(2) = 0.0_dp
    end type gluing_t

    type :: triplets_t
        integer, allocatable :: dof(:), broken(:)
        real(dp), allocatable :: weight(:)
        integer :: count = 0
    end type triplets_t

    public :: build_c1_extraction

    interface
        subroutine dgesv(n, nrhs, a, lda, ipiv, b, ldb, info)
            import :: dp
            integer, intent(in) :: n, nrhs, lda, ldb
            real(dp), intent(inout) :: a(lda, *), b(ldb, *)
            integer, intent(out) :: ipiv(*), info
        end subroutine dgesv
        subroutine dsyev(jobz, uplo, n, a, lda, w, work, lwork, info)
            import :: dp
            character, intent(in) :: jobz, uplo
            integer, intent(in) :: n, lda, lwork
            real(dp), intent(inout) :: a(lda, *)
            real(dp), intent(out) :: w(*), work(*)
            integer, intent(out) :: info
        end subroutine dsyev
    end interface

contains

    subroutine build_c1_extraction(layout, degree, blocks, offsets, dof, broken, &
            weight, dof_count, status, dof_vertex)
        !! Triplets (dof, broken, weight) of E. Status 1: degree < 3, a curved
        !! side that is not a boundary side or whose adjacent sides are C1
        !! interfaces, or non-uniform interior multiplicity below 2; 2: an
        !! interface without a linear beta split (not AS-G1). dof_vertex is
        !! the vertex of a vertex function and 0 for the other functions.
        type(block_layout_t), intent(in) :: layout
        integer, intent(in) :: degree, offsets(:)
        type(c1_block_knots_t), intent(in) :: blocks(:)
        integer, allocatable, intent(out) :: dof(:), broken(:)
        real(dp), allocatable, intent(out) :: weight(:)
        integer, intent(out) :: dof_count, status
        integer, allocatable, intent(out), optional :: dof_vertex(:)

        type(gluing_t), allocatable :: glue(:, :)
        type(triplets_t) :: e
        integer, allocatable :: vertex_of(:)
        integer :: nb, k, s, a, interior_count

        status = 1
        nb = layout%block_count
        if (degree < 3 .or. any(layout%coons)) return
        do k = 1, nb
            s = layout%arc_side(k)
            if (layout%periodic(k)) return
            if (s == 0) cycle
            if (layout%neighbour(s, k) /= 0) return
            do a = 1, 4
                if (a == s .or. a == opposite_side(s)) cycle
                if (layout%neighbour(a, k) /= 0 .and. .not. layout%side_c0(a, k)) &
                    return
            end do
        end do
        do k = 1, nb
            if (.not. uniform_multiplicity(blocks(k)%knots1, degree)) return
            if (.not. uniform_multiplicity(blocks(k)%knots2, degree)) return
        end do
        allocate(glue(4, nb))
        call build_gluing(layout, glue, status)
        if (status /= 0) return
        allocate(e%dof(1024), e%broken(1024), e%weight(1024))
        dof_count = 0
        call add_interior(blocks, offsets, degree, e, dof_count)
        call add_edges(layout, blocks, offsets, degree, glue, e, dof_count, status)
        if (status /= 0) return
        interior_count = dof_count
        call add_vertices(layout, blocks, offsets, degree, glue, e, dof_count, &
            vertex_of, status)
        if (status /= 0) return
        if (present(dof_vertex)) then
            allocate(dof_vertex(dof_count))
            dof_vertex = 0
            dof_vertex(interior_count + 1:) = vertex_of
        end if
        dof = e%dof(:e%count)
        broken = e%broken(:e%count)
        weight = e%weight(:e%count)
        status = 0
    end subroutine build_c1_extraction

    logical function uniform_multiplicity(knots, degree) result(ok)
        real(dp), intent(in) :: knots(:)
        integer, intent(in) :: degree

        integer :: k, run, m

        ok = .false.
        m = 0
        run = 0
        do k = degree + 2, size(knots) - degree - 1
            if (knots(k) == knots(k - 1)) then
                run = run + 1
            else
                if (run > 0) then
                    if (m == 0) m = run
                    if (run /= m) return
                end if
                run = 1
            end if
        end do
        if (run > 0) then
            if (m == 0) m = run
            if (run /= m) return
        end if
        ok = m == 0 .or. m >= 2
    end function uniform_multiplicity

    ! ----------------------------------------------------------------
    ! Frames, bilinear geometry and gluing data
    ! ----------------------------------------------------------------

    pure subroutine side_frame(side, from_end, o, p, eps)
        integer, intent(in) :: side
        logical, intent(in) :: from_end
        real(dp), intent(out) :: o(2), p(2, 2)
        integer, intent(out) :: eps

        select case (side)
        case (1)
            o = [0.0_dp, 0.0_dp]
            p = reshape([0.0_dp, 1.0_dp, 1.0_dp, 0.0_dp], [2, 2])
        case (2)
            o = [1.0_dp, 0.0_dp]
            p = reshape([-1.0_dp, 0.0_dp, 0.0_dp, 1.0_dp], [2, 2])
        case (3)
            o = [0.0_dp, 1.0_dp]
            p = reshape([0.0_dp, -1.0_dp, 1.0_dp, 0.0_dp], [2, 2])
        case default
            o = [0.0_dp, 0.0_dp]
            p = reshape([1.0_dp, 0.0_dp, 0.0_dp, 1.0_dp], [2, 2])
        end select
        if (from_end) then
            o = o + p(:, 2)
            p(:, 2) = -p(:, 2)
        end if
        eps = nint(p(1, 1)*p(2, 2) - p(1, 2)*p(2, 1))
    end subroutine side_frame

    pure subroutine bilinear_jacobian(corner, xi, df)
        real(dp), intent(in) :: corner(2, 4), xi(2)
        real(dp), intent(out) :: df(2, 2)

        df(:, 1) = (1.0_dp - xi(2))*(corner(:, 2) - corner(:, 1)) + &
            xi(2)*(corner(:, 3) - corner(:, 4))
        df(:, 2) = (1.0_dp - xi(1))*(corner(:, 4) - corner(:, 1)) + &
            xi(1)*(corner(:, 3) - corner(:, 2))
    end subroutine bilinear_jacobian

    pure function mixed(corner, a, b) result(m)
        !! d2F/(da db) of the bilinear map along axis directions a and b.
        real(dp), intent(in) :: corner(2, 4), a(2), b(2)
        real(dp) :: m(2)

        m = (corner(:, 1) - corner(:, 2) + corner(:, 3) - corner(:, 4))* &
            (a(1)*b(2) + a(2)*b(1))
    end function mixed

    pure real(dp) function det2(a, b)
        real(dp), intent(in) :: a(2), b(2)
        det2 = a(1)*b(2) - a(2)*b(1)
    end function det2

    subroutine build_gluing(layout, glue, status)
        type(block_layout_t), intent(in) :: layout
        type(gluing_t), intent(out) :: glue(:, :)
        integer, intent(out) :: status

        real(dp) :: o1(2), p1(2, 2), o2(2), p2(2, 2), df(2, 2), c1(2, 4), c2(2, 4)
        real(dp) :: a1(0:2), a2(0:2), bq(0:2), tau, sl(2), sr(2), scale
        real(dp) :: m(3, 4), rhs(3), x(4), mm(3, 3), y(3), al(2), ar(2)
        integer :: k, s, k2, s2, e1, e2, q, pivot(3), info
        logical :: rev

        status = 0
        do k = 1, layout%block_count
            do s = 1, 4
                k2 = layout%neighbour(s, k)
                if (k2 == 0) cycle
                if (layout%side_c0(s, k)) cycle
                s2 = layout%neighbour_side(s, k)
                if (k2 < k .or. (k2 == k .and. s2 < s)) cycle
                rev = layout%reversed(s, k)
                call side_frame(s, .false., o1, p1, e1)
                call side_frame(s2, rev, o2, p2, e2)
                status = 2
                if (e1 == e2) return
                c1 = layout%vertices(:, layout%corners(:, k))
                c2 = layout%vertices(:, layout%corners(:, k2))
                do q = 0, 2
                    tau = 0.5_dp*q
                    call bilinear_jacobian(c1, o1 + tau*p1(:, 2), df)
                    sl = matmul(df, p1(:, 1))
                    a1(q) = e1*det2(sl, matmul(df, p1(:, 2)))
                    call bilinear_jacobian(c2, o2 + tau*p2(:, 2), df)
                    sr = matmul(df, p2(:, 1))
                    a2(q) = e2*det2(sr, matmul(df, p2(:, 2)))
                    if (e1 == 1) then
                        bq(q) = det2(sr, sl)
                    else
                        bq(q) = det2(sl, sr)
                    end if
                end do
                if (min(minval(a1), minval(a2)) <= 0.0_dp) return
                scale = 2.0_dp/(a1(1) + a2(1))
                a1 = scale*a1
                a2 = scale*a2
                bq = scale*bq
                if (e1 == 1) then
                    al = [a1(0), a1(2) - a1(0)]
                    ar = [a2(0), a2(2) - a2(0)]
                else
                    al = [a2(0), a2(2) - a2(0)]
                    ar = [a1(0), a1(2) - a1(0)]
                end if
                rhs = [bq(0), -3.0_dp*bq(0) + 4.0_dp*bq(1) - bq(2), &
                    2.0_dp*bq(0) - 4.0_dp*bq(1) + 2.0_dp*bq(2)]
                ! unknowns x = (bL0, bL1, bR0, bR1)
                m(1, :) = [ar(1), 0.0_dp, al(1), 0.0_dp]
                m(2, :) = [ar(2), ar(1), al(2), al(1)]
                m(3, :) = [0.0_dp, ar(2), 0.0_dp, al(2)]
                mm = matmul(m, transpose(m))
                y = rhs
                call dgesv(3, 1, mm, 3, pivot, y, 3, info)
                if (info /= 0) then
                    call minimal_norm_rank2(m, rhs, x, info)
                    if (info /= 0) return
                else
                    x = matmul(transpose(m), y)
                end if
                if (maxval(abs(matmul(m, x) - rhs)) > 1.0e-10_dp* &
                    max(1.0_dp, maxval(abs(rhs)))) return
                if (e1 == 1) then
                    glue(s, k)%alpha = al
                    glue(s, k)%beta = x(1:2)
                    call store_reversed(glue(s2, k2), ar, x(3:4), rev)
                else
                    glue(s, k)%alpha = ar
                    glue(s, k)%beta = x(3:4)
                    call store_reversed(glue(s2, k2), al, x(1:2), rev)
                end if
                status = 0
            end do
        end do
    end subroutine build_gluing

    subroutine minimal_norm_rank2(m, rhs, x, info)
        !! Proportional alpha_L, alpha_R make the split matrix rank 2; solve
        !! the normal equations on the row space by a pseudo-inverse.
        real(dp), intent(in) :: m(3, 4), rhs(3)
        real(dp), intent(out) :: x(4)
        integer, intent(out) :: info

        real(dp) :: a(4, 4), b(4), lambda
        integer :: pivot(4)

        lambda = 1.0e-14_dp*max(1.0_dp, maxval(abs(m)))**2
        a = matmul(transpose(m), m)
        a(1, 1) = a(1, 1) + lambda
        a(2, 2) = a(2, 2) + lambda
        a(3, 3) = a(3, 3) + lambda
        a(4, 4) = a(4, 4) + lambda
        b = matmul(transpose(m), rhs)
        call dgesv(4, 1, a, 4, pivot, b, 4, info)
        x = b
    end subroutine minimal_norm_rank2

    pure subroutine store_reversed(g, alpha, beta, rev)
        !! Store data given in the neighbour's common tau in its own
        !! orientation: tau -> 1 - tau flips the sign of beta.
        type(gluing_t), intent(out) :: g
        real(dp), intent(in) :: alpha(2), beta(2)
        logical, intent(in) :: rev

        if (rev) then
            g%alpha = [alpha(1) + alpha(2), -alpha(2)]
            g%beta = [-(beta(1) + beta(2)), beta(2)]
        else
            g%alpha = alpha
            g%beta = beta
        end if
    end subroutine store_reversed

    pure subroutine oriented_gluing(g, from_end, alpha, beta)
        type(gluing_t), intent(in) :: g
        logical, intent(in) :: from_end
        real(dp), intent(out) :: alpha(2), beta(2)

        if (from_end) then
            alpha = [g%alpha(1) + g%alpha(2), -g%alpha(2)]
            beta = [-(g%beta(1) + g%beta(2)), g%beta(2)]
        else
            alpha = g%alpha
            beta = g%beta
        end if
    end subroutine oriented_gluing

    ! ----------------------------------------------------------------
    ! Univariate spaces
    ! ----------------------------------------------------------------

    subroutine uspace_init(u, knots, degree)
        type(uspace_t), intent(out) :: u
        real(dp), intent(in) :: knots(:)
        integer, intent(in) :: degree

        integer :: status

        u%degree = degree
        u%knots = knots
        u%n = size(knots) - degree - 1
        call build_bezier_extraction_1d(knots, degree, u%extraction, u%first, &
            u%breaks, status)
        if (status /= 0) error stop "spline_c1: invalid knot vector"
    end subroutine uspace_init

    subroutine uspace_derived(u, tau_reversed, degree, multiplicity, v)
        !! Space of the given degree and interior multiplicity on the
        !! breakpoints of u, in tau = xi or tau = 1 - xi.
        type(uspace_t), intent(in) :: u
        logical, intent(in) :: tau_reversed
        integer, intent(in) :: degree, multiplicity
        type(uspace_t), intent(out) :: v

        real(dp), allocatable :: breaks(:), knots(:)
        integer, allocatable :: mult(:)
        integer :: status

        breaks = u%breaks
        if (tau_reversed) breaks = 1.0_dp - breaks(size(breaks):1:-1)
        allocate(mult(size(breaks)))
        mult = multiplicity
        call build_open_knot_vector(degree, breaks, mult, knots, status)
        if (status /= 0) error stop "spline_c1: derived knot vector"
        call uspace_init(v, knots, degree)
    end subroutine uspace_derived

    subroutine uspace_jets(u, x, jets)
        !! jets(i, d) = d^d N_i/dx^d at x, all functions.
        type(uspace_t), intent(in) :: u
        real(dp), intent(in) :: x
        real(dp), intent(out) :: jets(:, 0:)

        real(dp) :: b(u%degree + 1, 0:2), h, t
        integer :: e, a, d, status

        e = size(u%breaks) - 1
        do a = 1, size(u%breaks) - 2
            if (x < u%breaks(a + 1)) then
                e = a
                exit
            end if
        end do
        h = u%breaks(e + 1) - u%breaks(e)
        t = min(max((x - u%breaks(e))/h, 0.0_dp), 1.0_dp)
        call interval_bernstein_jet(u%degree, t, b(:, 0), b(:, 1), b(:, 2), status)
        jets = 0.0_dp
        do d = 0, min(2, ubound(jets, 2))
            do a = 1, u%degree + 1
                jets(u%first(e) + a - 1, d) = &
                    dot_product(u%extraction(a, :, e), b(:, d))/h**d
            end do
        end do
    end subroutine uspace_jets

    pure integer function interior_multiplicity(u) result(m)
        type(uspace_t), intent(in) :: u

        m = 1
        if (size(u%breaks) > 2) m = count(u%knots == u%breaks(2))
    end function interior_multiplicity

    subroutine greville_interpolate(u, values_at, coef)
        !! Coefficients in u of the spline with the given values at the
        !! Greville points of u (values_at(i) is the value at tau(x_i)).
        type(uspace_t), intent(in) :: u
        real(dp), intent(in) :: values_at(:)
        real(dp), intent(out) :: coef(:)

        real(dp) :: a(u%n, u%n), b(u%n, 1), jets(u%n, 0:0)
        integer :: i, pivot(u%n), info

        do i = 1, u%n
            call uspace_jets(u, greville(u, i), jets)
            a(i, :) = jets(:, 0)
        end do
        b(:, 1) = values_at
        call dgesv(u%n, 1, a, u%n, pivot, b, u%n, info)
        if (info /= 0) error stop "spline_c1: Greville collocation"
        coef = b(:, 1)
    end subroutine greville_interpolate

    pure real(dp) function greville(u, i)
        type(uspace_t), intent(in) :: u
        integer, intent(in) :: i
        greville = sum(u%knots(i + 1:i + u%degree))/real(u%degree, dp)
    end function greville

    ! ----------------------------------------------------------------
    ! Function assembly
    ! ----------------------------------------------------------------

    subroutine add_interior(blocks, offsets, degree, e, dof_count)
        type(c1_block_knots_t), intent(in) :: blocks(:)
        integer, intent(in) :: offsets(:), degree
        type(triplets_t), intent(inout) :: e
        integer, intent(inout) :: dof_count

        integer :: k, i, j, n1, n2

        do k = 1, size(blocks)
            n1 = size(blocks(k)%knots1) - degree - 1
            n2 = size(blocks(k)%knots2) - degree - 1
            do j = 3, n2 - 2
                do i = 3, n1 - 2
                    dof_count = dof_count + 1
                    call push(e, dof_count, offsets(k) + i + n1*(j - 1), 1.0_dp)
                end do
            end do
        end do
    end subroutine add_interior

    subroutine side_spaces(blocks, degree, k, side, from_end, along, across, &
            along_dir, across_dir, along_rev, across_rev)
        !! Univariate block spaces along and across a side seen from its
        !! start (or end) vertex.
        type(c1_block_knots_t), intent(in) :: blocks(:)
        integer, intent(in) :: degree, k, side
        logical, intent(in) :: from_end
        type(uspace_t), intent(out) :: along, across
        integer, intent(out) :: along_dir, across_dir
        logical, intent(out) :: along_rev, across_rev

        real(dp) :: o(2), p(2, 2)
        integer :: eps

        call side_frame(side, from_end, o, p, eps)
        if (p(1, 2) /= 0.0_dp) then
            along_dir = 1
            across_dir = 2
            call uspace_init(along, blocks(k)%knots1, degree)
            call uspace_init(across, blocks(k)%knots2, degree)
        else
            along_dir = 2
            across_dir = 1
            call uspace_init(along, blocks(k)%knots2, degree)
            call uspace_init(across, blocks(k)%knots1, degree)
        end if
        along_rev = sum(p(:, 2)) < 0.0_dp
        across_rev = sum(p(:, 1)) < 0.0_dp
    end subroutine side_spaces

    subroutine transversal_hermite(across, across_rev, c0, c1)
        !! Block coefficients of c0 = b0 + b1 and c1 = b1/b1'(0) in sigma,
        !! where sigma = xi or 1 - xi.
        type(uspace_t), intent(in) :: across
        logical, intent(in) :: across_rev
        real(dp), intent(out) :: c0(:), c1(:)

        real(dp) :: jets(across%n, 0:1)
        integer :: n

        n = across%n
        c0 = 0.0_dp
        c1 = 0.0_dp
        if (across_rev) then
            call uspace_jets(across, 1.0_dp, jets)
            c0(n - 1:n) = 1.0_dp
            c1(n - 1) = -1.0_dp/jets(n - 1, 1)
        else
            call uspace_jets(across, 0.0_dp, jets)
            c0(1:2) = 1.0_dp
            c1(2) = 1.0_dp/jets(2, 1)
        end if
    end subroutine transversal_hermite

    subroutine add_edges(layout, blocks, offsets, degree, glue, e, dof_count, status)
        type(block_layout_t), intent(in) :: layout
        type(c1_block_knots_t), intent(in) :: blocks(:)
        integer, intent(in) :: offsets(:), degree
        type(gluing_t), intent(in) :: glue(:, :)
        type(triplets_t), intent(inout) :: e
        integer, intent(inout) :: dof_count
        integer, intent(out) :: status

        integer :: k, s, k2, s2, side_count, j, which, kind, n0, n1, m
        integer :: ks(2), ss(2)
        logical :: revs(2), c0
        type(uspace_t) :: along, plus, minus

        status = 0
        do k = 1, layout%block_count
            do s = 1, 4
                k2 = layout%neighbour(s, k)
                s2 = layout%neighbour_side(s, k)
                if (k2 /= 0) then
                    if (k2 < k .or. (k2 == k .and. s2 < s)) cycle
                    side_count = 2
                    ks = [k, k2]
                    ss = [s, s2]
                    revs = [.false., layout%reversed(s, k)]
                else
                    side_count = 1
                    ks = [k, 0]
                    ss = [s, 0]
                    revs = .false.
                end if
                c0 = .false.
                if (k2 /= 0) c0 = layout%side_c0(s, k)
                call side_tau_spaces(blocks, degree, k, s, .false., along, plus, minus, m)
                n0 = plus%n
                n1 = minus%n
                do kind = 0, 1
                    do j = merge(4, 3, kind == 0), merge(n0 - 3, n1 - 2, kind == 0)
                        if (.not. c0) dof_count = dof_count + 1
                        do which = 1, side_count
                            if (c0 .and. (kind == 1 .or. which == 1)) &
                                dof_count = dof_count + 1
                            call edge_function(blocks, offsets, degree, &
                                glue(ss(which), ks(which)), ks(which), ss(which), &
                                revs(which), kind, j, e, dof_count)
                        end do
                    end do
                end do
            end do
        end do
    end subroutine add_edges

    subroutine side_tau_spaces(blocks, degree, k, side, from_end, along, plus, &
            minus, m)
        type(c1_block_knots_t), intent(in) :: blocks(:)
        integer, intent(in) :: degree, k, side
        logical, intent(in) :: from_end
        type(uspace_t), intent(out) :: along, plus, minus
        integer, intent(out) :: m

        type(uspace_t) :: across
        integer :: along_dir, across_dir
        logical :: along_rev, across_rev

        call side_spaces(blocks, degree, k, side, from_end, along, across, &
            along_dir, across_dir, along_rev, across_rev)
        m = interior_multiplicity(along)
        call uspace_derived(along, along_rev, degree, max(m - 1, 1), plus)
        call uspace_derived(along, along_rev, degree - 1, max(m - 1, 1), minus)
    end subroutine side_tau_spaces

    subroutine edge_function(blocks, offsets, degree, g, k, side, from_end, kind, &
            j, e, dof)
        !! kind 0: b+_j(tau) c0(sigma) - beta(tau) b+_j'(tau) c1(sigma);
        !! kind 1: eps alpha(tau) b-_j(tau) c1(sigma) (1-based j).
        type(c1_block_knots_t), intent(in) :: blocks(:)
        integer, intent(in) :: offsets(:), degree, k, side, kind, j, dof
        type(gluing_t), intent(in) :: g
        logical, intent(in) :: from_end
        type(triplets_t), intent(inout) :: e

        type(uspace_t) :: along, across, plus, minus
        real(dp) :: o(2), p(2, 2), alpha(2), beta(2), tau
        real(dp), allocatable :: v0(:), v1(:), t0(:), t1(:), c0(:), c1(:)
        real(dp), allocatable :: jp(:, :), jm(:, :)
        integer :: eps, along_dir, across_dir, i, m
        logical :: along_rev, across_rev

        call side_frame(side, from_end, o, p, eps)
        call side_spaces(blocks, degree, k, side, from_end, along, across, &
            along_dir, across_dir, along_rev, across_rev)
        m = interior_multiplicity(along)
        call uspace_derived(along, along_rev, degree, max(m - 1, 1), plus)
        call uspace_derived(along, along_rev, degree - 1, max(m - 1, 1), minus)
        call oriented_gluing(g, from_end, alpha, beta)
        allocate(v0(along%n), v1(along%n), t0(along%n), t1(along%n), &
            c0(across%n), c1(across%n), jp(plus%n, 0:1), jm(minus%n, 0:0))
        do i = 1, along%n
            tau = greville(along, i)
            if (along_rev) tau = 1.0_dp - tau
            call uspace_jets(plus, tau, jp)
            call uspace_jets(minus, tau, jm)
            if (kind == 0) then
                v0(i) = jp(j, 0)
                v1(i) = -(beta(1) + beta(2)*tau)*jp(j, 1)
            else
                v0(i) = 0.0_dp
                v1(i) = eps*(alpha(1) + alpha(2)*tau)*jm(j, 0)
            end if
        end do
        call greville_interpolate(along, v0, t0)
        call greville_interpolate(along, v1, t1)
        call transversal_hermite(across, across_rev, c0, c1)
        call push_tensor(e, dof, offsets(k), along_dir, t0, c0, along%n, across%n)
        call push_tensor(e, dof, offsets(k), along_dir, t1, c1, along%n, across%n)
    end subroutine edge_function

    subroutine add_vertices(layout, blocks, offsets, degree, glue, e, dof_count, &
            vertex_of, status)
        type(block_layout_t), intent(in) :: layout
        type(c1_block_knots_t), intent(in) :: blocks(:)
        integer, intent(in) :: offsets(:), degree
        type(gluing_t), intent(in) :: glue(:, :)
        type(triplets_t), intent(inout) :: e
        integer, intent(inout) :: dof_count
        integer, allocatable, intent(out) :: vertex_of(:)
        integer, intent(out) :: status

        real(dp), allocatable :: basis(:, :)
        integer, allocatable :: ck(:, :), sector(:)
        real(dp) :: sigma
        integer :: v, k, c, f, valence, nsector, q

        status = 0
        allocate(vertex_of(0))
        do v = 1, layout%vertex_count
            call vertex_corners(layout, v, ck)
            valence = size(ck, 2)
            if (valence == 0) cycle
            sigma = 0.0_dp
            do q = 1, valence
                sigma = sigma + corner_scale(layout, blocks, degree, ck(1, q), ck(2, q))
            end do
            sigma = real(valence, dp)/sigma
            call vertex_sectors(layout, ck, sector, nsector)
            call vertex_basis(layout, ck, sector, nsector, sigma, basis, status)
            if (status /= 0) return
            do f = 1, size(basis, 2)
                dof_count = dof_count + 1
                vertex_of = [vertex_of, v]
                do q = 1, valence
                    k = ck(1, q)
                    c = ck(2, q)
                    call vertex_function(layout, blocks, offsets, degree, glue, &
                        k, c, basis(6*(sector(q) - 1) + 1:6*sector(q), f), e, dof_count)
                end do
            end do
        end do
    end subroutine add_vertices

    subroutine vertex_corners(layout, v, ck)
        !! (block, corner) pairs at vertex v.
        type(block_layout_t), intent(in) :: layout
        integer, intent(in) :: v
        integer, allocatable, intent(out) :: ck(:, :)

        integer :: k, c, n

        n = count(layout%corners == v)
        allocate(ck(2, n))
        n = 0
        do k = 1, layout%block_count
            do c = 1, 4
                if (layout%corners(c, k) /= v) cycle
                n = n + 1
                ck(:, n) = [k, c]
            end do
        end do
    end subroutine vertex_corners

    subroutine vertex_sectors(layout, ck, sector, nsector)
        !! Sectors: corners joined through C1 sides at the vertex.
        type(block_layout_t), intent(in) :: layout
        integer, intent(in) :: ck(:, :)
        integer, allocatable, intent(out) :: sector(:)
        integer, intent(out) :: nsector

        integer, allocatable :: parent(:)
        integer :: q, s, k2, r, a, b, label

        allocate(parent(size(ck, 2)), sector(size(ck, 2)))
        parent = [(q, q = 1, size(ck, 2))]
        do q = 1, size(ck, 2)
            do s = 1, 4
                if (side_start(s) /= ck(2, q) .and. side_end(s) /= ck(2, q)) cycle
                k2 = layout%neighbour(s, ck(1, q))
                if (k2 == 0 .or. layout%side_c0(s, ck(1, q))) cycle
                do r = 1, size(ck, 2)
                    if (ck(1, r) /= k2) cycle
                    a = root(q)
                    b = root(r)
                    parent(max(a, b)) = min(a, b)
                end do
            end do
        end do
        nsector = 0
        sector = 0
        do q = 1, size(ck, 2)
            label = root(q)
            if (sector(label) == 0) then
                nsector = nsector + 1
                sector(label) = nsector
            end if
            sector(q) = sector(label)
        end do

    contains

        integer function root(i) result(r0)
            integer, intent(in) :: i

            r0 = i
            do while (parent(r0) /= r0)
                r0 = parent(r0)
            end do
        end function root

    end subroutine vertex_sectors

    subroutine vertex_basis(layout, ck, sector, nsector, sigma, basis, status)
        !! Columns: C2 data (value, grad, H11, H12, H22) per sector of the
        !! vertex functions. One sector: the scaled unit vectors. Several:
        !! an orthonormal null space of the C0 matching conditions in the
        !! scaled data.
        type(block_layout_t), intent(in) :: layout
        integer, intent(in) :: ck(:, :), sector(:), nsector
        real(dp), intent(in) :: sigma
        real(dp), allocatable, intent(out) :: basis(:, :)
        integer, intent(out) :: status

        real(dp), allocatable :: rows(:, :), gram(:, :), lambda(:), work(:)
        real(dp) :: scale(6), corner(2, 4), df(2, 2), o(2), p(2, 2), t(2)
        integer :: exponent(6), q, s, k2, r, m, n, i, eps, a, b, info, kept

        status = 0
        exponent = [0, 1, 1, 2, 2, 2]
        scale = sigma**exponent
        n = 6*nsector
        if (nsector == 1) then
            allocate(basis(6, 6))
            basis = 0.0_dp
            do i = 1, 6
                basis(i, i) = scale(i)
            end do
            return
        end if
        allocate(rows(n, 0))
        do q = 1, size(ck, 2)
            do s = 1, 4
                if (side_start(s) /= ck(2, q) .and. side_end(s) /= ck(2, q)) cycle
                k2 = layout%neighbour(s, ck(1, q))
                if (k2 == 0 .or. .not. layout%side_c0(s, ck(1, q))) cycle
                do r = 1, size(ck, 2)
                    if (ck(1, r) /= k2 .or. sector(r) == sector(q)) cycle
                    corner = layout%vertices(:, layout%corners(:, ck(1, q)))
                    call side_frame(s, side_end(s) == ck(2, q), o, p, eps)
                    call bilinear_jacobian(corner, o, df)
                    t = matmul(df, p(:, 2))
                    a = 6*(sector(q) - 1)
                    b = 6*(sector(r) - 1)
                    call add_row([1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp])
                    call add_row([0.0_dp, t(1), t(2), 0.0_dp, 0.0_dp, 0.0_dp])
                    call add_row([0.0_dp, 0.0_dp, 0.0_dp, t(1)**2, 2.0_dp*t(1)*t(2), &
                        t(2)**2])
                end do
            end do
        end do
        m = size(rows, 2)
        allocate(gram(n, n), lambda(n), work(10*n))
        gram = matmul(rows, transpose(rows))
        call dsyev('V', 'U', n, gram, n, lambda, work, 10*n, info)
        status = 1
        if (info /= 0) return
        kept = count(lambda <= 1.0e-10_dp*max(maxval(lambda), 1.0_dp))
        allocate(basis(n, kept))
        basis = gram(:, :kept)
        do i = 1, nsector
            do q = 1, kept
                basis(6*(i - 1) + 1:6*i, q) = scale*basis(6*(i - 1) + 1:6*i, q)
            end do
        end do
        status = 0

    contains

        subroutine add_row(w)
            real(dp), intent(in) :: w(6)
            real(dp), allocatable :: grown(:, :)
            real(dp) :: row(n)

            row = 0.0_dp
            row(a + 1:a + 6) = scale*w
            row(b + 1:b + 6) = -scale*w
            if (maxval(abs(row)) == 0.0_dp) return
            row = row/maxval(abs(row))
            allocate(grown(n, size(rows, 2) + 1))
            grown(:, :size(rows, 2)) = rows
            grown(:, size(rows, 2) + 1) = row
            call move_alloc(grown, rows)
        end subroutine add_row

    end subroutine vertex_basis

    real(dp) function corner_scale(layout, blocks, degree, k, c) result(scale)
        type(block_layout_t), intent(in) :: layout
        type(c1_block_knots_t), intent(in) :: blocks(:)
        integer, intent(in) :: degree, k, c

        real(dp) :: df(2, 2), corner(2, 4), xi(2), h

        select case (c)
        case (1)
            xi = [0.0_dp, 0.0_dp]
        case (2)
            xi = [1.0_dp, 0.0_dp]
        case (3)
            xi = [1.0_dp, 1.0_dp]
        case default
            xi = [0.0_dp, 1.0_dp]
        end select
        corner = layout%vertices(:, layout%corners(:, k))
        call bilinear_jacobian(corner, xi, df)
        h = 0.5_dp*(blocks(k)%knots1(degree + 2) - blocks(k)%knots1(degree + 1) + &
            blocks(k)%knots2(degree + 2) - blocks(k)%knots2(degree + 1))
        scale = h/real(degree, dp)*norm2(df)
    end function corner_scale

    subroutine vertex_function(layout, blocks, offsets, degree, glue, k, c, phi, &
            e, dof)
        !! Vertex function on block k at corner c for C2 data
        !! phi = (value, grad, H11, H12, H22): the two edge parts minus the
        !! corner part (survey arXiv:1812.09088, p. 12-13).
        type(block_layout_t), intent(in) :: layout
        type(c1_block_knots_t), intent(in) :: blocks(:)
        integer, intent(in) :: offsets(:), degree, k, c, dof
        type(gluing_t), intent(in) :: glue(:, :)
        real(dp), intent(in) :: phi(6)
        type(triplets_t), intent(inout) :: e

        real(dp) :: corner(2, 4), grad(2), hess(2, 2), tangent(2, 2)
        real(dp) :: dcorner(2), coef_u(2), coef_v(2)
        real(dp), allocatable :: cu0(:), cu1(:), cv0(:), cv1(:)
        integer :: sides(2), which, s, eps(2), dirs(2), n(2)
        logical :: ends(2)
        type(uspace_t) :: along, across
        integer :: along_dir, across_dir, a_side, b_side
        logical :: along_rev, across_rev
        real(dp) :: o(2), p(2, 2), axis(2, 2)

        corner = layout%vertices(:, layout%corners(:, k))
        grad = phi(2:3)
        hess = reshape([phi(4), phi(5), phi(5), phi(6)], [2, 2])
        which = 0
        do s = 1, 4
            if (side_start(s) == c .or. side_end(s) == c) then
                which = which + 1
                sides(which) = s
                ends(which) = side_end(s) == c
                call side_frame(s, ends(which), o, p, eps(which))
                axis(:, which) = p(:, 2)
            end if
        end do
        do which = 1, 2
            call edge_part(layout, blocks, offsets, degree, &
                glue(sides(which), k), k, sides(which), ends(which), grad, hess, &
                phi(1), e, dof)
        end do
        ! corner part: u along the eps = -1 side A, v along the eps = +1 side B
        a_side = merge(1, 2, eps(1) == -1)
        b_side = 3 - a_side
        call side_frame(sides(a_side), ends(a_side), o, p, eps(1))
        call bilinear_jacobian(corner, o, tangent)
        tangent = matmul(tangent, axis(:, [a_side, b_side]))
        dcorner = mixed(corner, axis(:, a_side), axis(:, b_side))
        call side_spaces(blocks, degree, k, sides(b_side), ends(b_side), along, &
            across, along_dir, across_dir, along_rev, across_rev)
        n(1) = across%n
        allocate(cu0(n(1)), cu1(n(1)))
        call transversal_hermite(across, across_rev, cu0, cu1)
        dirs(1) = across_dir
        call side_spaces(blocks, degree, k, sides(a_side), ends(a_side), along, &
            across, along_dir, across_dir, along_rev, across_rev)
        n(2) = across%n
        allocate(cv0(n(2)), cv1(n(2)))
        call transversal_hermite(across, across_rev, cv0, cv1)
        dirs(2) = across_dir
        coef_u = [phi(1), dot_product(grad, tangent(:, 1))]
        coef_v = [dot_product(grad, tangent(:, 2)), &
            dot_product(tangent(:, 1), matmul(hess, tangent(:, 2))) + &
            dot_product(grad, dcorner)]
        ! -(d00 c0 c0 + d10 c1(u) c0(v) + d01 c0(u) c1(v) + d11 c1(u) c1(v))
        call push_corner(e, dof, offsets(k), dirs(1), -coef_u(1), cu0, cv0, n)
        call push_corner(e, dof, offsets(k), dirs(1), -coef_u(2), cu1, cv0, n)
        call push_corner(e, dof, offsets(k), dirs(1), -coef_v(1), cu0, cv1, n)
        call push_corner(e, dof, offsets(k), dirs(1), -coef_v(2), cu1, cv1, n)
    end subroutine vertex_function

    subroutine edge_part(layout, blocks, offsets, degree, g, k, side, from_end, &
            grad, hess, value, e, dof)
        type(block_layout_t), intent(in) :: layout
        type(c1_block_knots_t), intent(in) :: blocks(:)
        integer, intent(in) :: offsets(:), degree, k, side, dof
        type(gluing_t), intent(in) :: g
        logical, intent(in) :: from_end
        real(dp), intent(in) :: grad(2), hess(2, 2), value
        type(triplets_t), intent(inout) :: e

        type(uspace_t) :: along, across, plus, minus
        real(dp) :: o(2), p(2, 2), corner(2, 4), df(2, 2), alpha(2), beta(2)
        real(dp) :: t(2), ds(2), dvec(2), dprime(2), dts(2), d0(0:2), d1(0:1)
        real(dp) :: hp(3, 3), hm(2, 2), cp(3, 3), cm(2, 2), tau
        real(dp), allocatable :: v0(:), v1(:), t0(:), t1(:), c0(:), c1(:)
        real(dp), allocatable :: jp(:, :), jm(:, :)
        integer :: eps, along_dir, across_dir, i, j, m, pivot(3), info
        logical :: along_rev, across_rev

        call side_frame(side, from_end, o, p, eps)
        corner = layout%vertices(:, layout%corners(:, k))
        call oriented_gluing(g, from_end, alpha, beta)
        call bilinear_jacobian(corner, o, df)
        t = matmul(df, p(:, 2))
        ds = matmul(df, p(:, 1))
        dts = mixed(corner, p(:, 1), p(:, 2))
        dvec = eps*(ds + beta(1)*t)/alpha(1)
        dprime = eps*((dts + beta(2)*t)*alpha(1) - (ds + beta(1)*t)*alpha(2))/ &
            alpha(1)**2
        d0 = [value, dot_product(grad, t), dot_product(t, matmul(hess, t))]
        d1 = [dot_product(grad, dvec), dot_product(t, matmul(hess, dvec)) + &
            dot_product(grad, dprime)]
        call side_spaces(blocks, degree, k, side, from_end, along, across, &
            along_dir, across_dir, along_rev, across_rev)
        m = interior_multiplicity(along)
        call uspace_derived(along, along_rev, degree, max(m - 1, 1), plus)
        call uspace_derived(along, along_rev, degree - 1, max(m - 1, 1), minus)
        ! Hermite bases: d^i c+_j(0) = delta_ij (i, j = 0..2), c-_j likewise
        allocate(jp(plus%n, 0:2), jm(minus%n, 0:2))
        call uspace_jets(plus, 0.0_dp, jp)
        call uspace_jets(minus, 0.0_dp, jm)
        hp = transpose(jp(1:3, 0:2))
        cp = reshape([1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, &
            0.0_dp, 1.0_dp], [3, 3])
        call dgesv(3, 3, hp, 3, pivot, cp, 3, info)
        if (info /= 0) error stop "spline_c1: Hermite basis c+"
        hm = transpose(jm(1:2, 0:1))
        cm = reshape([1.0_dp, 0.0_dp, 0.0_dp, 1.0_dp], [2, 2])
        call dgesv(2, 2, hm, 2, pivot, cm, 2, info)
        if (info /= 0) error stop "spline_c1: Hermite basis c-"
        allocate(v0(along%n), v1(along%n), t0(along%n), t1(along%n), &
            c0(across%n), c1(across%n))
        do i = 1, along%n
            tau = greville(along, i)
            if (along_rev) tau = 1.0_dp - tau
            call uspace_jets(plus, tau, jp)
            call uspace_jets(minus, tau, jm)
            v0(i) = 0.0_dp
            v1(i) = 0.0_dp
            do j = 0, 2
                v0(i) = v0(i) + d0(j)*dot_product(cp(:, j + 1), jp(1:3, 0))
                v1(i) = v1(i) - d0(j)*(beta(1) + beta(2)*tau)* &
                    dot_product(cp(:, j + 1), jp(1:3, 1))
            end do
            do j = 0, 1
                v1(i) = v1(i) + eps*d1(j)*(alpha(1) + alpha(2)*tau)* &
                    dot_product(cm(:, j + 1), jm(1:2, 0))
            end do
        end do
        call greville_interpolate(along, v0, t0)
        call greville_interpolate(along, v1, t1)
        call transversal_hermite(across, across_rev, c0, c1)
        call push_tensor(e, dof, offsets(k), along_dir, t0, c0, along%n, across%n)
        call push_tensor(e, dof, offsets(k), along_dir, t1, c1, along%n, across%n)
    end subroutine edge_part

    subroutine push_tensor(e, dof, offset, along_dir, t, c, n_along, n_across)
        type(triplets_t), intent(inout) :: e
        integer, intent(in) :: dof, offset, along_dir, n_along, n_across
        real(dp), intent(in) :: t(:), c(:)

        integer :: i, j

        do j = 1, n_across
            if (c(j) == 0.0_dp) cycle
            do i = 1, n_along
                if (abs(t(i)*c(j)) < 1.0e-14_dp) cycle
                if (along_dir == 1) then
                    call push(e, dof, offset + i + n_along*(j - 1), t(i)*c(j))
                else
                    call push(e, dof, offset + j + n_across*(i - 1), t(i)*c(j))
                end if
            end do
        end do
    end subroutine push_tensor

    subroutine push_corner(e, dof, offset, u_dir, factor, cu, cv, n)
        !! factor * cu(u) * cv(v); u runs along block direction u_dir.
        type(triplets_t), intent(inout) :: e
        integer, intent(in) :: dof, offset, u_dir, n(2)
        real(dp), intent(in) :: factor, cu(:), cv(:)

        integer :: i, j

        if (factor == 0.0_dp) return
        do j = 1, n(2)
            if (cv(j) == 0.0_dp) cycle
            do i = 1, n(1)
                if (cu(i) == 0.0_dp) cycle
                if (u_dir == 1) then
                    call push(e, dof, offset + i + n(1)*(j - 1), factor*cu(i)*cv(j))
                else
                    call push(e, dof, offset + j + n(2)*(i - 1), factor*cu(i)*cv(j))
                end if
            end do
        end do
    end subroutine push_corner

    subroutine push(e, dof, broken, weight)
        type(triplets_t), intent(inout) :: e
        integer, intent(in) :: dof, broken
        real(dp), intent(in) :: weight

        integer, allocatable :: idof(:), ibroken(:)
        real(dp), allocatable :: w(:)

        if (e%count == size(e%dof)) then
            allocate(idof(2*e%count), ibroken(2*e%count), w(2*e%count))
            idof(:e%count) = e%dof
            ibroken(:e%count) = e%broken
            w(:e%count) = e%weight
            call move_alloc(idof, e%dof)
            call move_alloc(ibroken, e%broken)
            call move_alloc(w, e%weight)
        end if
        e%count = e%count + 1
        e%dof(e%count) = dof
        e%broken(e%count) = broken
        e%weight(e%count) = weight
    end subroutine push

end module fortfem_spline_c1
