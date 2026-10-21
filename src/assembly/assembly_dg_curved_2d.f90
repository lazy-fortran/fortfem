module fortfem_dg_curved_2d
    !! Discontinuous Galerkin advection-diffusion-reaction on curved 2D meshes.
    !!
    !! The space is the broken Lagrange space of degree k (0 to 3 tested to 3)
    !! on the curved isoparametric triangles of a meridian_mesh_t (geometry
    !! order 1 to 3).  The global index of local dof i of triangle t is
    !! (t - 1)*nl + i, as in fortfem_triangle_discontinuous_dof_map.
    !!
    !! The operator is
    !!   L u = div(b u) + c u - div(K grad u)        (DG_FORM_CONSERVATIVE)
    !!   L u = b . grad u + c u - div(K grad u)      (DG_FORM_ADVECTIVE)
    !! with a real velocity b(2), a complex reaction c and a complex, possibly
    !! anisotropic and nonsymmetric 2x2 tensor K, and Dirichlet data g on the
    !! whole boundary (the advective part uses g on the inflow part only).
    !! The bilinear form a(v, u) (rows v, columns u, no complex conjugation)
    !! uses the upwind flux for the advection and the symmetric interior
    !! penalty (SIPG) method for the diffusion:
    !!   - {K grad u . n}[v] - {K^T grad v . n}[u] + sigma [u][v],
    !!   sigma = eta (k+1)(k+2)/2 max_s |n . K_s n| max_s |F|/|T_s|,
    !! pointwise on every face, where n is the unit normal of the face, [.]
    !! the jump and {.} the average; on a boundary face [u] = u and {.} is
    !! the trace.  The penalty n . K n is the anisotropic scaling of Ern,
    !! Stephansen and Zunino: it is small across faces normal to the weak
    !! direction of K and coercive for positive definite K and eta large
    !! enough (default 4).  The conservative form is locally conservative:
    !! with c = 0 every column of an element without boundary faces sums to
    !! zero.  Measure weights (e.g. R for the meridian R dR dZ measure) are
    !! carried by the coefficients and the source.
    !!
    !! Coefficients and data are batched callbacks per triangle t at the
    !! physical points x(2, m); they are called concurrently from OpenMP
    !! threads.  Element and face loops are OpenMP parallel.
    use fortfem_kinds, only: dp
    use fortfem_meridian_curved_mesh, only: meridian_mesh_t
    use fortfem_triangle_duffy_quadrature, only: triangle_duffy_quadrature
    use fortfem_triangle_lagrange_arbitrary_order, only: &
        evaluate_triangle_lagrange_basis, initialize_triangle_lagrange_basis, &
        triangle_lagrange_basis_t
    use fortnum_quadrature, only: gauss_legendre_ab
    use fortsparse, only: csc_z_t
    implicit none
    private

    integer, parameter, public :: DG_FORM_CONSERVATIVE = 1
    integer, parameter, public :: DG_FORM_ADVECTIVE = 2
    real(dp), parameter, public :: DG_DEFAULT_PENALTY = 4.0_dp
    real(dp), parameter :: default_penalty = DG_DEFAULT_PENALTY
    integer, parameter :: par_min = 32

    type, public :: dg_curved_space_t
        type(meridian_mesh_t) :: mm
        integer :: k = 0, nl = 0, ng = 0, ndof = 0, nt = 0
        integer :: nq = 0, nfq = 0, nface = 0, nint = 0
        real(dp), allocatable :: qw(:), tl(:, :), tlg(:, :, :)
        real(dp), allocatable :: gvq(:, :), gg(:, :, :)
        real(dp), allocatable :: fw(:), ftl(:, :, :), ftlg(:, :, :, :)
        real(dp), allocatable :: fgv(:, :, :), fgg(:, :, :, :)
        real(dp) :: dref(2, 3) = 0.0_dp
        integer, allocatable :: face_t(:, :), face_e(:, :)
        integer, allocatable :: elem_face(:, :), elem_side(:, :)
        logical, allocatable :: face_rev(:)
        real(dp), allocatable :: area(:), flen(:)
    end type dg_curved_space_t

    type, abstract, public :: dg_coefficients_t
    contains
        procedure(dg_coefficients_i), deferred :: evaluate
    end type dg_coefficients_t

    type, abstract, public :: dg_data_t
    contains
        procedure(dg_field_i), deferred :: source
        procedure(dg_field_i), deferred :: boundary
    end type dg_data_t

    abstract interface
        subroutine dg_coefficients_i(self, t, x, kk, b, c, status)
            !! K(2, 2, m), b(2, m) and c(m) at the points x(2, m) of t.
            import :: dp, dg_coefficients_t
            class(dg_coefficients_t), intent(in) :: self
            integer, intent(in) :: t
            real(dp), intent(in) :: x(:, :)
            complex(dp), intent(out) :: kk(:, :, :)
            real(dp), intent(out) :: b(:, :)
            complex(dp), intent(out) :: c(:)
            integer, intent(out) :: status
        end subroutine dg_coefficients_i

        subroutine dg_field_i(self, t, x, f, status)
            !! Source f or boundary value g at the points x(2, m) of t.
            import :: dp, dg_data_t
            class(dg_data_t), intent(in) :: self
            integer, intent(in) :: t
            real(dp), intent(in) :: x(:, :)
            complex(dp), intent(out) :: f(:)
            integer, intent(out) :: status
        end subroutine dg_field_i
    end interface

    public :: dg_curved_space_init
    public :: dg_curved_assemble
    public :: dg_curved_rhs
    public :: dg_curved_qpoint
    public :: penalty_factor

contains

    subroutine dg_curved_space_init(mm, degree, space, status)
        !! Broken degree-k space on the curved mesh mm, its quadrature rules
        !! and the face connectivity (interior faces first).
        type(meridian_mesh_t), intent(in) :: mm
        integer, intent(in) :: degree
        type(dg_curved_space_t), intent(out) :: space
        integer, intent(out) :: status

        type(triangle_lagrange_basis_t) :: basis, gbasis
        real(dp), allocatable :: xi(:), eta(:), s(:)
        real(dp) :: vref(2, 3), p(2)
        integer :: q, e, nt

        status = 1
        if (degree < 0 .or. degree > 6) return
        if (.not. allocated(mm%nodes)) return
        nt = mm%topology%n_triangles
        if (nt < 1) return
        call initialize_triangle_lagrange_basis(degree, basis, status)
        if (status /= 0) return
        call initialize_triangle_lagrange_basis(mm%geo_order, gbasis, status)
        if (status /= 0) return
        space%mm = mm
        space%k = degree
        space%nt = nt
        space%nl = (degree + 1)*(degree + 2)/2
        space%ng = size(mm%nodes, 2)
        space%ndof = space%nl*nt

        call triangle_duffy_quadrature(2*degree + 3*mm%geo_order, xi, eta, &
            space%qw, status)
        if (status /= 0) return
        space%nq = size(space%qw)
        allocate (space%tl(space%nl, space%nq), space%tlg(2, space%nl, space%nq))
        allocate (space%gg(2, space%ng, space%nq), space%gvq(space%ng, space%nq))
        do q = 1, space%nq
            call evaluate_triangle_lagrange_basis(basis, xi(q), eta(q), &
                space%tl(:, q), space%tlg(:, :, q), status)
            if (status /= 0) return
            call evaluate_triangle_lagrange_basis(gbasis, xi(q), eta(q), &
                space%gvq(:, q), space%gg(:, :, q), status)
            if (status /= 0) return
        end do

        space%nfq = degree + mm%geo_order + 2
        allocate (s(space%nfq), space%fw(space%nfq))
        call gauss_legendre_ab(space%nfq, 0.0_dp, 1.0_dp, s, space%fw)
        allocate (space%ftl(space%nl, space%nfq, 3))
        allocate (space%ftlg(2, space%nl, space%nfq, 3))
        allocate (space%fgv(space%ng, space%nfq, 3))
        allocate (space%fgg(2, space%ng, space%nfq, 3))
        vref = reshape([0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, 1.0_dp], [2, 3])
        do e = 1, 3
            space%dref(:, e) = vref(:, mod(e, 3) + 1) - vref(:, e)
            do q = 1, space%nfq
                p = vref(:, e) + s(q)*space%dref(:, e)
                call evaluate_triangle_lagrange_basis(basis, p(1), p(2), &
                    space%ftl(:, q, e), space%ftlg(:, :, q, e), status)
                if (status /= 0) return
                call evaluate_triangle_lagrange_basis(gbasis, p(1), p(2), &
                    space%fgv(:, q, e), space%fgg(:, :, q, e), status)
                if (status /= 0) return
            end do
        end do
        call build_faces(space, status)
        if (status /= 0) return
        call measure_cells(space, status)
    end subroutine dg_curved_space_init

    subroutine dg_curved_assemble(space, coef, form, a, status, penalty)
        !! Sparse matrix a(v, u) of the DG operator (rows: test functions).
        !! The block pattern of the broken space is known, so the CSC arrays
        !! are filled directly, one element column block per OpenMP task:
        !! column block t couples to t and its face neighbours.  The face
        !! terms of an interior face are evaluated by both sides, each for
        !! its own columns, which keeps the element loop free of races.
        type(dg_curved_space_t), intent(in) :: space
        class(dg_coefficients_t), intent(in) :: coef
        integer, intent(in) :: form
        type(csc_z_t), intent(out) :: a
        integer, intent(out) :: status
        real(dp), intent(in), optional :: penalty

        integer, allocatable :: nblk(:), blk(:, :)
        complex(dp), allocatable :: ml(:, :)
        real(dp) :: eta
        integer :: nl, nt, t, e, f, s, r, i, j, jj, bad, st, p, col, nr

        status = 1
        if (form /= DG_FORM_CONSERVATIVE .and. form /= DG_FORM_ADVECTIVE) return
        eta = default_penalty
        if (present(penalty)) eta = penalty
        nl = space%nl
        nt = space%nt
        allocate (nblk(nt), blk(4, nt))
        do t = 1, nt
            call row_blocks(space, t, nblk(t), blk(:, t))
        end do
        a%nrow = space%ndof
        a%ncol = space%ndof
        allocate (a%col_ptr(space%ndof + 1))
        a%col_ptr(1) = 1
        do t = 1, nt
            do jj = 1, nl
                j = (t - 1)*nl + jj
                a%col_ptr(j + 1) = a%col_ptr(j) + nblk(t)*nl
            end do
        end do
        a%nnz = a%col_ptr(space%ndof + 1) - 1
        allocate (a%row_idx(a%nnz), a%val(a%nnz))
        bad = 0
        !$omp parallel if(nt >= par_min) default(shared) &
        !$omp private(t, e, f, s, r, i, j, jj, st, p, col, nr, ml) &
        !$omp reduction(max:bad)
        allocate (ml(2*nl, 2*nl))
        !$omp do schedule(dynamic, 8)
        do t = 1, nt
            do jj = 1, nl
                col = (t - 1)*nl + jj
                p = a%col_ptr(col)
                do r = 1, nblk(t)
                    do i = 1, nl
                        a%row_idx(p) = (blk(r, t) - 1)*nl + i
                        a%val(p) = (0.0_dp, 0.0_dp)
                        p = p + 1
                    end do
                end do
            end do
            call volume_block(space, coef, form, t, ml(:nl, :nl), st)
            if (st /= 0) then
                bad = 1
                cycle
            end if
            call add_block(t, t, ml(:nl, :nl))
            do e = 1, 3
                f = space%elem_face(e, t)
                s = space%elem_side(e, t)
                call face_block(space, coef, form, eta, f, ml, st, s)
                if (st /= 0) exit
                nr = 1
                if (f <= space%nint) nr = 2
                do r = 1, nr
                    call add_block(space%face_t(r, f), t, &
                        ml((r - 1)*nl + 1:r*nl, (s - 1)*nl + 1:s*nl))
                end do
            end do
            if (st /= 0) bad = 1
        end do
        !$omp end do
        !$omp end parallel
        if (bad == 0) status = 0

    contains

        subroutine add_block(rt, ct, mb)
            !! Add the (row element rt, column element ct) block mb.
            integer, intent(in) :: rt, ct
            complex(dp), intent(in) :: mb(:, :)

            integer :: k, q, base

            do k = 1, nblk(ct)
                if (blk(k, ct) == rt) exit
            end do
            do q = 1, nl
                base = a%col_ptr((ct - 1)*nl + q) + (k - 1)*nl - 1
                a%val(base + 1:base + nl) = a%val(base + 1:base + nl) + mb(:, q)
            end do
        end subroutine add_block

    end subroutine dg_curved_assemble

    subroutine row_blocks(space, t, n, blk)
        !! Sorted elements coupled to the column block of triangle t.
        type(dg_curved_space_t), intent(in) :: space
        integer, intent(in) :: t
        integer, intent(out) :: n, blk(4)

        integer :: e, f, i, v

        n = 1
        blk(1) = t
        do e = 1, 3
            f = space%elem_face(e, t)
            if (f > space%nint) cycle
            n = n + 1
            blk(n) = space%face_t(3 - space%elem_side(e, t), f)
        end do
        do i = 2, n
            v = blk(i)
            do e = i - 1, 1, -1
                if (blk(e) <= v) exit
                blk(e + 1) = blk(e)
            end do
            blk(e + 1) = v
        end do
    end subroutine row_blocks

    subroutine dg_curved_rhs(space, coef, data, form, rhs, status, penalty)
        !! Load vector: volume source and the boundary terms of the data g
        !! (SIPG Dirichlet terms and the upwind inflow value).
        type(dg_curved_space_t), intent(in) :: space
        class(dg_coefficients_t), intent(in) :: coef
        class(dg_data_t), intent(in) :: data
        integer, intent(in) :: form
        complex(dp), intent(out) :: rhs(:)
        integer, intent(out) :: status
        real(dp), intent(in), optional :: penalty

        real(dp), allocatable :: x(:, :), w(:), gr(:, :, :), n(:, :), ws(:)
        real(dp), allocatable :: val(:, :), b(:, :)
        complex(dp), allocatable :: f(:), kk(:, :, :), c(:), rl(:)
        complex(dp) :: kn(2), kgt
        real(dp) :: eta, beta, sig
        integer :: t, q, fc, i, j, bad, st, nl, e

        status = 1
        if (form /= DG_FORM_CONSERVATIVE .and. form /= DG_FORM_ADVECTIVE) return
        if (size(rhs) /= space%ndof) return
        eta = default_penalty
        if (present(penalty)) eta = penalty
        nl = space%nl
        rhs = (0.0_dp, 0.0_dp)
        bad = 0
        !$omp parallel if(space%nt >= par_min) default(shared) &
        !$omp private(x, w, gr, n, ws, val, b, f, kk, c, rl, sig, kn, kgt, beta, t, q, &
        !$omp fc, i, j, st, e) reduction(max:bad)
        allocate (x(2, max(space%nq, space%nfq)), w(space%nq))
        allocate (gr(2, nl, max(space%nq, space%nfq)), n(2, space%nfq))
        allocate (ws(space%nfq), val(nl, space%nfq), b(2, space%nfq))
        allocate (f(max(space%nq, space%nfq)), kk(2, 2, space%nfq), c(space%nfq))
        allocate (rl(nl))
        !$omp do schedule(dynamic, 8)
        do t = 1, space%nt
            call volume_points(space, t, x(:, :space%nq), w, gr(:, :, :space%nq), &
                st)
            if (st == 0) call data%source(t, x(:, :space%nq), f(:space%nq), st)
            if (st /= 0) then
                bad = 1
                cycle
            end if
            do q = 1, space%nq
                i = (t - 1)*nl
                rhs(i + 1:i + nl) = rhs(i + 1:i + nl) + w(q)*f(q)*space%tl(:, q)
            end do
        end do
        !$omp end do
        !$omp do schedule(dynamic, 8)
        do fc = space%nint + 1, space%nface
            t = space%face_t(1, fc)
            e = space%face_e(1, fc)
            call face_points(space, t, e, x(:, :space%nfq), n, ws, val, &
                gr(:, :, :space%nfq), st)
            if (st == 0) call coef%evaluate(t, x(:, :space%nfq), kk, b, c, st)
            if (st == 0) call data%boundary(t, x(:, :space%nfq), f(:space%nfq), st)
            if (st /= 0) then
                bad = 1
                cycle
            end if
            i = (t - 1)*nl
            rl = (0.0_dp, 0.0_dp)
            do q = 1, space%nfq
                beta = dot_product(b(:, q), n(:, q))
                kn(1) = kk(1, 1, q)*n(1, q) + kk(1, 2, q)*n(2, q)
                kn(2) = kk(2, 1, q)*n(1, q) + kk(2, 2, q)*n(2, q)
                sig = penalty_factor(space, eta, space%flen(fc)/space%area(t)) &
                    *abs(dot_product(n(:, q), kn))
                rl = rl + ws(q)*f(q)*sig*val(:, q)
                do j = 1, nl
                    kgt = gr(1, j, q)*kn(1) + gr(2, j, q)*kn(2)
                    rl(j) = rl(j) - ws(q)*f(q)*kgt
                end do
                if (beta < 0.0_dp) then
                    rl = rl - ws(q)*beta*f(q)*val(:, q)
                end if
            end do
            ! A corner triangle owns two boundary faces that different threads
            ! may process at the same time; serialize the update of its rows.
            !$omp critical (dg_curved_rhs_boundary)
            rhs(i + 1:i + nl) = rhs(i + 1:i + nl) + rl
            !$omp end critical (dg_curved_rhs_boundary)
        end do
        !$omp end do
        !$omp end parallel
        if (bad == 0) status = 0
    end subroutine dg_curved_rhs

    pure function penalty_factor(space, eta, ratio) result(p)
        type(dg_curved_space_t), intent(in) :: space
        real(dp), intent(in) :: eta, ratio
        real(dp) :: p

        p = eta*0.5_dp*real((space%k + 1)*(space%k + 2), dp)*ratio
    end function penalty_factor

    subroutine volume_block(space, coef, form, t, ml, status)
        !! Element matrix of the volume terms of triangle t.
        type(dg_curved_space_t), intent(in) :: space
        class(dg_coefficients_t), intent(in) :: coef
        integer, intent(in) :: form, t
        complex(dp), intent(out) :: ml(:, :)
        integer, intent(out) :: status

        real(dp) :: x(2, space%nq), w(space%nq), gr(2, space%nl, space%nq)
        real(dp) :: b(2, space%nq), bg(space%nl), phi(space%nl)
        complex(dp) :: kk(2, 2, space%nq), c(space%nq), kg(2, space%nl)
        integer :: q, j

        ml = (0.0_dp, 0.0_dp)
        call volume_points(space, t, x, w, gr, status)
        if (status /= 0) return
        call coef%evaluate(t, x, kk, b, c, status)
        if (status /= 0) return
        do q = 1, space%nq
            phi = space%tl(:, q)
            kg(1, :) = kk(1, 1, q)*gr(1, :, q) + kk(1, 2, q)*gr(2, :, q)
            kg(2, :) = kk(2, 1, q)*gr(1, :, q) + kk(2, 2, q)*gr(2, :, q)
            bg = b(1, q)*gr(1, :, q) + b(2, q)*gr(2, :, q)
            do j = 1, space%nl
                ml(:, j) = ml(:, j) + w(q)*(gr(1, :, q)*kg(1, j) &
                    + gr(2, :, q)*kg(2, j) + c(q)*phi*phi(j))
                if (form == DG_FORM_CONSERVATIVE) then
                    ml(:, j) = ml(:, j) - w(q)*bg*phi(j)
                else
                    ml(:, j) = ml(:, j) + w(q)*phi*bg(j)
                end if
            end do
        end do
    end subroutine volume_block

    subroutine face_block(space, coef, form, eta, f, ml, status, side)
        !! Face matrix of face f: (2 nl)^2 on an interior face, ordered (side
        !! 1, side 2), and the leading nl^2 on a boundary face.  Only the
        !! columns of the trial functions of the given side are computed.
        type(dg_curved_space_t), intent(in) :: space
        class(dg_coefficients_t), intent(in) :: coef
        integer, intent(in) :: form, f
        real(dp), intent(in) :: eta
        complex(dp), intent(out) :: ml(:, :)
        integer, intent(out) :: status
        integer, intent(in) :: side

        integer :: nl, nfq, m, q, q2, s, j, t(2), jlo, jhi
        real(dp) :: x(2, space%nfq), n(2, space%nfq), ws(space%nfq)
        real(dp) :: n2(2, space%nfq), ws2(space%nfq), x2(2, space%nfq)
        real(dp) :: val(space%nl, space%nfq, 2), gr(2, space%nl, space%nfq, 2)
        real(dp) :: val2(space%nl, space%nfq), gr2(2, space%nl, space%nfq)
        real(dp) :: b(2, space%nfq, 2), beta, ratio, half, sig
        complex(dp) :: kk(2, 2, space%nfq, 2), c(space%nfq, 2), kn(2, 2)
        complex(dp) :: jmp(2*space%nl), avg(2*space%nl), adj(2*space%nl)
        real(dp) :: up(2*space%nl), down(2*space%nl)

        nl = space%nl
        nfq = space%nfq
        ml = (0.0_dp, 0.0_dp)
        m = nl
        if (f <= space%nint) m = 2*nl
        jlo = (side - 1)*nl + 1
        jhi = side*nl
        t = space%face_t(:, f)
        call face_points(space, t(1), space%face_e(1, f), x, n, ws, val(:, :, 1), &
            gr(:, :, :, 1), status)
        if (status /= 0) return
        call coef%evaluate(t(1), x, kk(:, :, :, 1), b(:, :, 1), c(:, 1), status)
        if (status /= 0) return
        ratio = space%flen(f)/space%area(t(1))
        half = 1.0_dp
        if (m == 2*nl) then
            call face_points(space, t(2), space%face_e(2, f), x2, n2, ws2, &
                val2, gr2, status)
            if (status /= 0) return
            do q = 1, nfq
                q2 = q
                if (space%face_rev(f)) q2 = nfq + 1 - q
                val(:, q, 2) = val2(:, q2)
                gr(:, :, q, 2) = gr2(:, :, q2)
            end do
            call coef%evaluate(t(2), x, kk(:, :, :, 2), b(:, :, 2), c(:, 2), &
                status)
            if (status /= 0) return
            ratio = max(ratio, space%flen(f)/space%area(t(2)))
            half = 0.5_dp
        end if
        do q = 1, nfq
            sig = 0.0_dp
            beta = 0.0_dp
            do s = 1, m/nl
                q2 = (s - 1)*nl
                kn(1, s) = kk(1, 1, q, s)*n(1, q) + kk(1, 2, q, s)*n(2, q)
                kn(2, s) = kk(2, 1, q, s)*n(1, q) + kk(2, 2, q, s)*n(2, q)
                sig = max(sig, abs(dot_product(n(:, q), kn(:, s))))
                beta = beta + dot_product(b(:, q, s), n(:, q))/real(m/nl, dp)
                jmp(q2 + 1:q2 + nl) = merge(1.0_dp, -1.0_dp, s == 1)*val(:, q, s)
                ! {K grad u . n} and {K^T grad v . n}.
                avg(q2 + 1:q2 + nl) = half*(n(1, q)*(kk(1, 1, q, s)*gr(1, :, q, s) &
                    + kk(1, 2, q, s)*gr(2, :, q, s)) + n(2, q)*(kk(2, 1, q, s) &
                    *gr(1, :, q, s) + kk(2, 2, q, s)*gr(2, :, q, s)))
                adj(q2 + 1:q2 + nl) = half*(gr(1, :, q, s)*kn(1, s) &
                    + gr(2, :, q, s)*kn(2, s))
            end do
            sig = penalty_factor(space, eta, ratio)*sig
            up = 0.0_dp
            down = 0.0_dp
            if (m == 2*nl) then
                if (beta >= 0.0_dp) then
                    up(:nl) = val(:, q, 1)
                    down(nl + 1:) = val(:, q, 2)
                else
                    up(nl + 1:) = val(:, q, 2)
                    down(:nl) = val(:, q, 1)
                end if
            else if (beta >= 0.0_dp) then
                up(:nl) = val(:, q, 1)
            else
                down(:nl) = val(:, q, 1)
            end if
            do j = jlo, jhi
                ml(:m, j) = ml(:m, j) + ws(q)*(-jmp(:m)*avg(j) - adj(:m)*jmp(j) &
                    + sig*jmp(:m)*jmp(j))
                if (form == DG_FORM_CONSERVATIVE) then
                    ml(:m, j) = ml(:m, j) + ws(q)*beta*jmp(:m)*up(j)
                else
                    ml(:m, j) = ml(:m, j) - ws(q)*beta*down(:m)*jmp(j)
                end if
            end do
        end do
        status = 0
    end subroutine face_block

    subroutine build_faces(space, status)
        !! Match the local edges of all triangles through their vertex pairs.
        type(dg_curved_space_t), intent(inout) :: space
        integer, intent(out) :: status

        integer, allocatable :: start(:), fill(:), ent(:), vmax(:), mate(:)
        integer, allocatable :: tv(:, :)
        integer :: nv, nt, t, e, a, b, lo, i, j, slot, ni, nb, f, m

        status = 1
        nt = space%nt
        tv = space%mm%topology%triangles
        nv = maxval(tv)
        allocate (start(nv + 1), fill(nv), ent(3*nt), vmax(3*nt), mate(3*nt))
        start = 0
        do t = 1, nt
            do e = 1, 3
                lo = min(tv(e, t), tv(mod(e, 3) + 1, t))
                start(lo + 1) = start(lo + 1) + 1
            end do
        end do
        start(1) = 1
        do i = 1, nv
            start(i + 1) = start(i + 1) + start(i)
        end do
        fill = 0
        do t = 1, nt
            do e = 1, 3
                a = tv(e, t)
                b = tv(mod(e, 3) + 1, t)
                lo = min(a, b)
                slot = start(lo) + fill(lo)
                fill(lo) = fill(lo) + 1
                ent(slot) = 3*(t - 1) + e
                vmax(slot) = max(a, b)
            end do
        end do
        mate = 0
        do i = 1, nv
            do a = start(i), start(i + 1) - 1
                do b = a + 1, start(i + 1) - 1
                    if (vmax(a) /= vmax(b)) cycle
                    if (mate(a) /= 0) return
                    if (mate(b) /= 0) return
                    mate(a) = b
                    mate(b) = a
                end do
            end do
        end do
        ni = count(mate > 0)/2
        nb = count(mate == 0)
        space%nint = ni
        space%nface = ni + nb
        allocate (space%face_t(2, ni + nb), space%face_e(2, ni + nb))
        allocate (space%face_rev(ni + nb))
        space%face_t = 0
        space%face_e = 0
        space%face_rev = .false.
        f = 0
        m = ni
        do a = 1, 3*nt
            t = (ent(a) - 1)/3 + 1
            e = ent(a) - 3*(t - 1)
            if (mate(a) == 0) then
                m = m + 1
                space%face_t(1, m) = t
                space%face_e(1, m) = e
            else if (mate(a) > a) then
                f = f + 1
                space%face_t(1, f) = t
                space%face_e(1, f) = e
                b = mate(a)
                j = (ent(b) - 1)/3 + 1
                space%face_t(2, f) = j
                space%face_e(2, f) = ent(b) - 3*(j - 1)
                space%face_rev(f) = tv(e, t) /= tv(space%face_e(2, f), j)
            end if
        end do
        allocate (space%elem_face(3, nt), space%elem_side(3, nt))
        do f = 1, space%nface
            do i = 1, 2
                if (space%face_t(i, f) == 0) cycle
                space%elem_face(space%face_e(i, f), space%face_t(i, f)) = f
                space%elem_side(space%face_e(i, f), space%face_t(i, f)) = i
            end do
        end do
        status = 0
    end subroutine build_faces

    subroutine measure_cells(space, status)
        !! Areas of the curved triangles and lengths of the curved faces.
        type(dg_curved_space_t), intent(inout) :: space
        integer, intent(out) :: status

        real(dp), allocatable :: x(:, :), n(:, :), ws(:), val(:, :), gr(:, :, :)
        real(dp), allocatable :: w(:)
        integer :: t, f, st

        status = 1
        allocate (space%area(space%nt), space%flen(space%nface))
        allocate (x(2, space%nq), w(space%nq), gr(2, space%nl, space%nq))
        do t = 1, space%nt
            call volume_points(space, t, x, w, gr, st)
            if (st /= 0) return
            space%area(t) = sum(w)
        end do
        deallocate (x, gr)
        allocate (x(2, space%nfq), gr(2, space%nl, space%nfq))
        allocate (n(2, space%nfq), ws(space%nfq), val(space%nl, space%nfq))
        do f = 1, space%nface
            call face_points(space, space%face_t(1, f), space%face_e(1, f), x, &
                n, ws, val, gr, st)
            if (st /= 0) return
            space%flen(f) = sum(ws)
        end do
        status = 0
    end subroutine measure_cells

    subroutine dg_curved_qpoint(space, t, q, x, wq, val, grad, status)
        !! Volume quadrature point q of triangle t: physical point x, weight
        !! wq (dx dy), basis values val(nl) and physical gradients grad(2, nl).
        type(dg_curved_space_t), intent(in) :: space
        integer, intent(in) :: t, q
        real(dp), intent(out) :: x(2), wq, val(:), grad(:, :)
        integer, intent(out) :: status

        real(dp) :: jac(2, 2), jit(2, 2), det
        integer :: i

        status = 1
        x = 0.0_dp
        jac = 0.0_dp
        do i = 1, space%ng
            x = x + space%gvq(i, q)*space%mm%nodes(:, i, t)
            jac(:, 1) = jac(:, 1) + space%gg(1, i, q)*space%mm%nodes(:, i, t)
            jac(:, 2) = jac(:, 2) + space%gg(2, i, q)*space%mm%nodes(:, i, t)
        end do
        det = jac(1, 1)*jac(2, 2) - jac(1, 2)*jac(2, 1)
        if (det <= 0.0_dp) return
        jit(1, 1) = jac(2, 2)/det
        jit(2, 1) = -jac(1, 2)/det
        jit(1, 2) = -jac(2, 1)/det
        jit(2, 2) = jac(1, 1)/det
        val = space%tl(:, q)
        grad = matmul(jit, space%tlg(:, :, q))
        wq = space%qw(q)*det
        status = 0
    end subroutine dg_curved_qpoint

    subroutine volume_points(space, t, x, w, grad, status)
        !! All volume quadrature points of triangle t.
        type(dg_curved_space_t), intent(in) :: space
        integer, intent(in) :: t
        real(dp), intent(out), contiguous :: x(:, :), grad(:, :, :)
        real(dp), intent(out) :: w(:)
        integer, intent(out) :: status

        real(dp) :: val(space%nl)
        integer :: q

        do q = 1, space%nq
            call dg_curved_qpoint(space, t, q, x(:, q), w(q), val, grad(:, :, q), &
                status)
            if (status /= 0) return
        end do
    end subroutine volume_points

    subroutine face_points(space, t, e, x, n, ws, val, grad, status)
        !! Face quadrature points of local edge e of triangle t in the edge
        !! direction of t: physical points x, outward unit normals n, weights
        !! ws (arc length), basis values val(nl, nfq) and gradients.
        type(dg_curved_space_t), intent(in) :: space
        integer, intent(in) :: t, e
        real(dp), intent(out) :: x(:, :), n(:, :), ws(:), val(:, :)
        real(dp), intent(out) :: grad(:, :, :)
        integer, intent(out) :: status

        real(dp) :: jac(2, 2), jit(2, 2), det, tg(2), len
        integer :: q, i

        status = 1
        do q = 1, space%nfq
            x(:, q) = 0.0_dp
            jac = 0.0_dp
            do i = 1, space%ng
                x(:, q) = x(:, q) + space%fgv(i, q, e)*space%mm%nodes(:, i, t)
                jac(:, 1) = jac(:, 1) + space%fgg(1, i, q, e)*space%mm%nodes(:, i, t)
                jac(:, 2) = jac(:, 2) + space%fgg(2, i, q, e)*space%mm%nodes(:, i, t)
            end do
            det = jac(1, 1)*jac(2, 2) - jac(1, 2)*jac(2, 1)
            if (det <= 0.0_dp) return
            tg = matmul(jac, space%dref(:, e))
            len = norm2(tg)
            n(1, q) = tg(2)/len
            n(2, q) = -tg(1)/len
            ws(q) = space%fw(q)*len
            jit(1, 1) = jac(2, 2)/det
            jit(2, 1) = -jac(1, 2)/det
            jit(1, 2) = -jac(2, 1)/det
            jit(2, 2) = jac(1, 1)/det
            val(:, q) = space%ftl(:, q, e)
            grad(:, :, q) = matmul(jit, space%ftlg(:, :, q, e))
        end do
        status = 0
    end subroutine face_points

end module fortfem_dg_curved_2d
