module fortfem_spline_eval
    !! Quadrature and jets of multipatch spline fields
    !! (doc/design/multipatch_spline_spaces.md, eqs. 5-8).
    !!
    !! A jet holds, for each of nfield fields, the value and the derivatives
    !! with respect to the reference coordinates s up to the requested order.
    !! Component f + nfield*d with d = 0 value, d = 1, 2 for d/ds1, d/ds2 and,
    !! at order 2, d = 3, 4, 5 for (11, 21, 22). Coefficients are coef(nfield,
    !! dof_count). Cell points are tensor Gauss points, q = q1 + nq*(q2 - 1).
    use fortfem_kinds, only: dp
    use fortfem_gll_line, only: gauss_rule
    use fortfem_interval_bernstein_jets, only: interval_bernstein_jet
    use fortfem_spline_space, only: spline_space_t, spline_cell_broken
    use fortfem_spline_layout, only: side_start, side_end
    implicit none
    private

    integer, parameter, public :: SPLINE_JET_COUNT(0:2) = [1, 3, 6]

    type, public :: spline_quadrature_t
        integer :: points_per_direction = 0
        integer :: points_per_cell = 0
        integer :: point_count = 0
        real(dp), allocatable :: t(:) ! Gauss points on [0, 1]
        real(dp), allocatable :: t_weight(:)
        real(dp), allocatable :: bernstein(:, :, :) ! (p + 1, nq, 0:2)
        integer, allocatable :: point_offset(:) ! (cell): first point
        integer, allocatable :: block(:) ! (point)
        real(dp), allocatable :: xi(:, :) ! (2, point)
        real(dp), allocatable :: s(:, :) ! (2, point)
        real(dp), allocatable :: weight(:) ! (point), includes det DG
        real(dp), allocatable :: chain(:, :, :) ! (6, 6, point): s-jet = chain xi-jet
    end type spline_quadrature_t

    public :: spline_quadrature_build
    public :: spline_evaluate
    public :: spline_evaluate_at
    public :: spline_broken_coefficients
    public :: spline_cell_tables
    public :: spline_cell_xi_jets
    public :: spline_xi_to_s_chain
    public :: spline_geometry_at
    public :: spline_locate
    public :: spline_side_quadrature
    public :: spline_point_matrix
    public :: spline_boundary_nodes

    type, public :: spline_side_quadrature_t
        !! Gauss points on the block sides of one boundary tag, ordered by
        !! block and side. tangent is the unit ds/dt, counter-clockwise along
        !! the boundary, and weight = |ds/dt| dt.
        integer :: point_count = 0
        integer, allocatable :: block(:), side(:)
        real(dp), allocatable :: xi(:, :), s(:, :), tangent(:, :), weight(:)
    end type spline_side_quadrature_t

contains

    subroutine spline_quadrature_build(space, points_per_direction, quad, status)
        !! Status 1: invalid input, 2: a block map is singular at a point.
        type(spline_space_t), intent(in) :: space
        integer, intent(in) :: points_per_direction
        type(spline_quadrature_t), intent(out) :: quad
        integer, intent(out) :: status

        real(dp), allocatable :: x(:), w(:), u(:, :), jets(:, :, :)
        real(dp), allocatable :: table1(:, :, :), table2(:, :, :)
        real(dp) :: h(2), det
        integer :: nq, p, cell, q, q1, q2, first, k, d, ierr, local

        status = 1
        nq = points_per_direction
        p = space%degree
        if (nq < 1 .or. space%cell_count < 1) return
        call gauss_rule(nq, x, w, status)
        if (status /= 0) return
        quad%points_per_direction = nq
        quad%points_per_cell = nq*nq
        quad%point_count = nq*nq*space%cell_count
        quad%t = 0.5_dp*(x + 1.0_dp)
        quad%t_weight = 0.5_dp*w
        allocate(quad%bernstein(p + 1, nq, 0:2))
        do q = 1, nq
            call interval_bernstein_jet(p, quad%t(q), quad%bernstein(:, q, 0), &
                quad%bernstein(:, q, 1), quad%bernstein(:, q, 2), status)
            if (status /= 0) return
        end do
        allocate(quad%point_offset(space%cell_count), &
            quad%block(quad%point_count), quad%xi(2, quad%point_count), &
            quad%s(2, quad%point_count), quad%weight(quad%point_count), &
            quad%chain(6, 6, quad%point_count))
        allocate(u(2, (p + 1)**2), jets(2, 6, nq*nq))
        allocate(table1(p + 1, nq, 0:2), table2(p + 1, nq, 0:2))
        ierr = 0
        do cell = 1, space%cell_count
            first = (cell - 1)*nq*nq + 1
            quad%point_offset(cell) = first
            k = space%cell_block(cell)
            call spline_cell_tables(space, quad, cell, table1, table2, h)
            call geometry_coefficients(space, cell, u)
            call spline_cell_xi_jets(table1, table2, u, jets)
            do q2 = 1, nq
                do q1 = 1, nq
                    q = q1 + nq*(q2 - 1)
                    local = first + q - 1
                    quad%block(local) = space%patch_block(k)
                    associate (blk => space%blocks(k))
                        quad%xi(1, local) = blk%breaks1(space%cell_element(1, cell)) &
                            + h(1)*quad%t(q1)
                        quad%xi(2, local) = blk%breaks2(space%cell_element(2, cell)) &
                            + h(2)*quad%t(q2)
                    end associate
                    quad%s(:, local) = jets(:, 1, q)
                    det = jets(1, 2, q)*jets(2, 3, q) - jets(1, 3, q)*jets(2, 2, q)
                    if (det <= 0.0_dp) ierr = 1
                    quad%weight(local) = abs(det)*h(1)*h(2)* &
                        quad%t_weight(q1)*quad%t_weight(q2)
                    call spline_xi_to_s_chain(jets(:, :, q), quad%chain(:, :, local), d)
                    if (d /= 0) ierr = 1
                end do
            end do
        end do
        status = 0
        if (ierr /= 0) status = 2
    end subroutine spline_quadrature_build

    subroutine spline_broken_coefficients(space, coef, broken)
        !! broken = E^T c for coefficients coef(nfield, dof_count).
        type(spline_space_t), intent(in) :: space
        real(dp), intent(in) :: coef(:, :)
        real(dp), intent(out) :: broken(:, :)

        integer :: j, q

        !$omp parallel do default(none) shared(space, coef, broken) private(j, q)
        do j = 1, space%broken_count
            broken(:, j) = 0.0_dp
            do q = space%broken_ptr(j), space%broken_ptr(j + 1) - 1
                broken(:, j) = broken(:, j) + &
                    space%broken_weight(q)*coef(:, space%broken_dof(q))
            end do
        end do
        !$omp end parallel do
    end subroutine spline_broken_coefficients

    subroutine spline_evaluate(space, quad, coef, order, jet, status)
        !! jet(nfield*SPLINE_JET_COUNT(order), point_count) at all
        !! quadrature points.
        type(spline_space_t), intent(in) :: space
        type(spline_quadrature_t), intent(in) :: quad
        real(dp), intent(in) :: coef(:, :)
        integer, intent(in) :: order
        real(dp), intent(out) :: jet(:, :)
        integer, intent(out) :: status

        real(dp), allocatable :: broken(:, :), u(:, :), xi_jets(:, :, :)
        real(dp), allocatable :: table1(:, :, :), table2(:, :, :)
        integer, allocatable :: cell_broken(:)
        real(dp) :: h(2)
        integer :: nf, p, nq, cell, q, nd

        status = 1
        nf = size(coef, 1)
        p = space%degree
        nq = quad%points_per_direction
        if (order < 0 .or. order > 2 .or. size(coef, 2) /= space%dof_count) return
        nd = SPLINE_JET_COUNT(order)
        if (size(jet, 1) /= nf*nd .or. size(jet, 2) /= quad%point_count) return
        allocate(broken(nf, space%broken_count))
        call spline_broken_coefficients(space, coef, broken)
        !$omp parallel default(none) &
        !$omp shared(space, quad, broken, jet, nf, p, nq, nd) &
        !$omp private(u, xi_jets, table1, table2, cell_broken, h, cell, q)
        allocate(u(nf, (p + 1)**2), xi_jets(nf, 6, nq*nq), &
            table1(p + 1, nq, 0:2), table2(p + 1, nq, 0:2), &
            cell_broken((p + 1)**2))
        !$omp do
        do cell = 1, space%cell_count
            call spline_cell_broken(space, cell, cell_broken)
            u = broken(:, cell_broken)
            call spline_cell_tables(space, quad, cell, table1, table2, h)
            call spline_cell_xi_jets(table1, table2, u, xi_jets)
            do q = 1, nq*nq
                call apply_chain(quad%chain(:, :, quad%point_offset(cell) + q - 1), &
                    xi_jets(:, :, q), nd, jet(:, quad%point_offset(cell) + q - 1))
            end do
        end do
        !$omp end do
        !$omp end parallel
        status = 0
    end subroutine spline_evaluate

    subroutine spline_evaluate_at(space, coef, block, xi, order, jet, status)
        !! Jet at one block parameter xi in [0, 1]^2. Only the coefficients of
        !! the functions on the containing cell are gathered.
        type(spline_space_t), intent(in) :: space
        real(dp), intent(in) :: coef(:, :), xi(2)
        integer, intent(in) :: block, order
        real(dp), intent(out) :: jet(:)
        integer, intent(out) :: status

        real(dp), allocatable :: bern(:, :, :), u(:, :), g(:, :), xi_jets(:, :, :)
        real(dp), allocatable :: g_jets(:, :, :), table1(:, :, :), table2(:, :, :)
        integer, allocatable :: cell_broken(:)
        real(dp) :: chain(6, 6), h(2), t(2)
        integer :: nf, p, cell, k, a, q, ierr

        status = 1
        nf = size(coef, 1)
        p = space%degree
        if (block < 1 .or. block > space%layout_block_count .or. order < 0 .or. &
            order > 2) return
        if (any(xi < 0.0_dp) .or. any(xi > 1.0_dp)) return
        if (size(jet) /= nf*SPLINE_JET_COUNT(order)) return
        if (size(coef, 2) /= space%dof_count) return
        call point_cell(space, block, xi, cell, t)
        if (cell == 0) return
        allocate(bern(p + 1, 2, 0:2), table1(p + 1, 1, 0:2), table2(p + 1, 1, 0:2))
        do k = 1, 2
            call interval_bernstein_jet(p, t(k), bern(:, k, 0), bern(:, k, 1), &
                bern(:, k, 2), status)
            if (status /= 0) return
        end do
        call cell_tables_at(space, cell, bern(:, 1:1, :), bern(:, 2:2, :), table1, &
            table2, h)
        allocate(u(nf, (p + 1)**2), g(2, (p + 1)**2), xi_jets(nf, 6, 1), &
            g_jets(2, 6, 1), cell_broken((p + 1)**2))
        call spline_cell_broken(space, cell, cell_broken)
        u = 0.0_dp
        do a = 1, size(cell_broken)
            do q = space%broken_ptr(cell_broken(a)), &
                    space%broken_ptr(cell_broken(a) + 1) - 1
                u(:, a) = u(:, a) + space%broken_weight(q)*coef(:, space%broken_dof(q))
            end do
        end do
        call geometry_coefficients(space, cell, g)
        call spline_cell_xi_jets(table1, table2, u, xi_jets)
        call spline_cell_xi_jets(table1, table2, g, g_jets)
        call spline_xi_to_s_chain(g_jets(:, :, 1), chain, ierr)
        status = 2
        if (ierr /= 0) return
        call apply_chain(chain, xi_jets(:, :, 1), SPLINE_JET_COUNT(order), jet)
        status = 0
    end subroutine spline_evaluate_at

    subroutine spline_cell_tables(space, quad, cell, table1, table2, h)
        !! B-splines of the cell and their xi-derivatives at the Gauss points,
        !! table(a, q, d) = d^d N_a/dxi^d.
        type(spline_space_t), intent(in) :: space
        type(spline_quadrature_t), intent(in) :: quad
        integer, intent(in) :: cell
        real(dp), intent(out) :: table1(:, :, 0:), table2(:, :, 0:), h(2)

        call cell_tables_at(space, cell, quad%bernstein, quad%bernstein, &
            table1, table2, h)
    end subroutine spline_cell_tables

    pure subroutine spline_cell_xi_jets(table1, table2, u, jets)
        !! Sum-factorised xi-jets (value, 1, 2, 11, 12, 22) at the tensor
        !! points from cell coefficients u(nfield, (p + 1)^2).
        real(dp), intent(in) :: table1(:, :, 0:), table2(:, :, 0:), u(:, :)
        real(dp), intent(out) :: jets(:, :, :)

        integer, parameter :: d1(6) = [0, 1, 0, 2, 1, 0], d2(6) = [0, 0, 1, 0, 1, 2]
        real(dp) :: partial(size(u, 1), size(table1, 2), size(table2, 1), 0:2)
        integer :: nf, n1, nq1, nq2, a, b, q1, q2, d, c

        nf = size(u, 1)
        n1 = size(table1, 1)
        nq1 = size(table1, 2)
        nq2 = size(table2, 2)
        partial = 0.0_dp
        do d = 0, 2
            do b = 1, size(table2, 1)
                do a = 1, n1
                    c = a + n1*(b - 1)
                    do q1 = 1, nq1
                        partial(:, q1, b, d) = partial(:, q1, b, d) + &
                            table1(a, q1, d)*u(:, c)
                    end do
                end do
            end do
        end do
        jets = 0.0_dp
        do q2 = 1, nq2
            do q1 = 1, nq1
                do c = 1, 6
                    do b = 1, size(table2, 1)
                        jets(:, c, q1 + nq1*(q2 - 1)) = &
                            jets(:, c, q1 + nq1*(q2 - 1)) + &
                            table2(b, q2, d2(c))*partial(:, q1, b, d1(c))
                    end do
                end do
            end do
        end do
    end subroutine spline_cell_xi_jets

    pure subroutine spline_xi_to_s_chain(geometry_jet, chain, status)
        !! Linear map from the xi-jet (u, u_1, u_2, u_11, u_12, u_22) to the
        !! s-jet (u, u_s1, u_s2, u_s11, u_s21, u_s22) through the block map
        !! with xi-jet geometry_jet(2, 6), eqs. (7) and (8).
        real(dp), intent(in) :: geometry_jet(2, 6)
        real(dp), intent(out) :: chain(6, 6)
        integer, intent(out) :: status

        real(dp) :: dg(2, 2), inv(2, 2), hg(2, 2, 2), det, hxi(2, 2), grad(2)
        real(dp) :: hs(2, 2), unit(6)
        integer :: column, r

        status = 1
        chain = 0.0_dp
        dg(:, 1) = geometry_jet(:, 2)
        dg(:, 2) = geometry_jet(:, 3)
        det = dg(1, 1)*dg(2, 2) - dg(1, 2)*dg(2, 1)
        if (det == 0.0_dp) return
        inv(1, :) = [dg(2, 2), -dg(1, 2)]/det
        inv(2, :) = [-dg(2, 1), dg(1, 1)]/det
        do r = 1, 2
            hg(r, 1, 1) = geometry_jet(r, 4)
            hg(r, 1, 2) = geometry_jet(r, 5)
            hg(r, 2, 1) = geometry_jet(r, 5)
            hg(r, 2, 2) = geometry_jet(r, 6)
        end do
        do column = 1, 6
            unit = 0.0_dp
            unit(column) = 1.0_dp
            grad = matmul(transpose(inv), unit(2:3))
            hxi = reshape([unit(4), unit(5), unit(5), unit(6)], [2, 2])
            do r = 1, 2
                hxi = hxi - grad(r)*hg(r, :, :)
            end do
            hs = matmul(transpose(inv), matmul(hxi, inv))
            chain(:, column) = [unit(1), grad(1), grad(2), hs(1, 1), hs(2, 1), &
                hs(2, 2)]
        end do
        status = 0
    end subroutine spline_xi_to_s_chain

    pure subroutine apply_chain(chain, xi_jets, nd, jet)
        real(dp), intent(in) :: chain(6, 6), xi_jets(:, :)
        integer, intent(in) :: nd
        real(dp), intent(out) :: jet(:)

        integer :: nf, d, c

        nf = size(xi_jets, 1)
        do d = 1, nd
            jet(nf*(d - 1) + 1:nf*d) = 0.0_dp
            do c = 1, 6
                if (chain(d, c) == 0.0_dp) cycle
                jet(nf*(d - 1) + 1:nf*d) = jet(nf*(d - 1) + 1:nf*d) + &
                    chain(d, c)*xi_jets(:, c)
            end do
        end do
    end subroutine apply_chain

    subroutine spline_geometry_at(space, block, xi, jet, status)
        !! xi-jet (G, G_1, G_2, G_11, G_12, G_22) of the block map at xi.
        type(spline_space_t), intent(in) :: space
        integer, intent(in) :: block
        real(dp), intent(in) :: xi(2)
        real(dp), intent(out) :: jet(2, 6)
        integer, intent(out) :: status

        real(dp), allocatable :: bern(:, :, :), table1(:, :, :), table2(:, :, :)
        real(dp), allocatable :: g(:, :), g_jets(:, :, :)
        real(dp) :: h(2), t(2)
        integer :: p, cell, k

        status = 1
        p = space%degree
        if (block < 1 .or. block > space%layout_block_count) return
        call point_cell(space, block, min(max(xi, 0.0_dp), 1.0_dp), cell, t)
        if (cell == 0) return
        allocate(bern(p + 1, 2, 0:2), table1(p + 1, 1, 0:2), table2(p + 1, 1, 0:2), &
            g(2, (p + 1)**2), g_jets(2, 6, 1))
        do k = 1, 2
            call interval_bernstein_jet(p, t(k), bern(:, k, 0), bern(:, k, 1), &
                bern(:, k, 2), status)
            if (status /= 0) return
        end do
        call cell_tables_at(space, cell, bern(:, 1:1, :), bern(:, 2:2, :), table1, &
            table2, h)
        call geometry_coefficients(space, cell, g)
        call spline_cell_xi_jets(table1, table2, g, g_jets)
        jet = g_jets(:, :, 1)
        status = 0
    end subroutine spline_geometry_at

    subroutine spline_locate(space, block_hint, s, block, xi, status)
        !! Block and parameters of a reference point s, trying block_hint
        !! first. Projected Newton iteration on each block map (xi kept in
        !! [0, 1]^2, xi1 wrapped on periodic blocks), restarted from the
        !! nearest point of a 9 x 9 parameter grid when the first start
        !! fails. Points on block sides are found in one of the adjacent
        !! blocks. Status 1: no block contains s.
        type(spline_space_t), intent(in) :: space
        integer, intent(in) :: block_hint
        real(dp), intent(in) :: s(2)
        integer, intent(out) :: block
        real(dp), intent(out) :: xi(2)
        integer, intent(out) :: status

        real(dp) :: scale
        integer :: trial, k, start
        logical :: found

        status = 1
        block = 0
        scale = maxval(abs(space%geometry)) + 1.0_dp
        do start = 1, 2
            do trial = 0, space%layout_block_count
                if (trial == 0) then
                    k = block_hint
                    if (k < 1 .or. k > space%layout_block_count) cycle
                else
                    k = trial
                    if (k == block_hint) cycle
                end if
                if (start == 1) then
                    xi = 0.5_dp
                else
                    call grid_start(k, xi)
                end if
                call newton(k, xi, found)
                if (.not. found) cycle
                block = k
                status = 0
                return
            end do
        end do

    contains

        subroutine newton(k, xi, found)
            integer, intent(in) :: k
            real(dp), intent(inout) :: xi(2)
            logical, intent(out) :: found

            real(dp) :: jet(2, 6), r(2), det, step(2)
            integer :: it, ierr

            found = .false.
            do it = 1, 60
                call spline_geometry_at(space, k, xi, jet, ierr)
                if (ierr /= 0) return
                r = jet(:, 1) - s
                if (norm2(r) < 1.0e-14_dp*scale) exit
                det = jet(1, 2)*jet(2, 3) - jet(1, 3)*jet(2, 2)
                if (det == 0.0_dp) return
                step = [jet(2, 3)*r(1) - jet(1, 3)*r(2), &
                    -jet(2, 2)*r(1) + jet(1, 2)*r(2)]/det
                xi = project(k, xi - step)
            end do
            call spline_geometry_at(space, k, xi, jet, ierr)
            found = ierr == 0 .and. norm2(jet(:, 1) - s) <= 1.0e-12_dp*scale
        end subroutine newton

        pure function project(k, x) result(y)
            integer, intent(in) :: k
            real(dp), intent(in) :: x(2)
            real(dp) :: y(2)

            y = min(max(x, 0.0_dp), 1.0_dp)
            if (space%layout%periodic(k)) y(1) = modulo(x(1), 1.0_dp)
        end function project

        subroutine grid_start(k, xi)
            integer, intent(in) :: k
            real(dp), intent(out) :: xi(2)

            real(dp) :: jet(2, 6), best, d, x(2)
            integer :: a, b, ierr

            best = huge(1.0_dp)
            xi = 0.5_dp
            do b = 0, 8
                do a = 0, 8
                    x = [a, b]/8.0_dp
                    call spline_geometry_at(space, k, x, jet, ierr)
                    d = norm2(jet(:, 1) - s)
                    if (d < best) then
                        best = d
                        xi = x
                    end if
                end do
            end do
        end subroutine grid_start

    end subroutine spline_locate

    subroutine spline_side_quadrature(space, tag, points_per_cell, side_quad, status)
        !! Gauss points on every cell side that lies on a layout side with the
        !! given tag, ordered by cell.
        type(spline_space_t), intent(in) :: space
        integer, intent(in) :: tag, points_per_cell
        type(spline_side_quadrature_t), intent(out) :: side_quad
        integer, intent(out) :: status

        real(dp), allocatable :: x(:), w(:)
        real(dp) :: jet(2, 6), xi(2), tangent(2), speed, h, t, t0
        integer :: cell, patch, k, s, q, n, along, ierr, pass, e(2)
        real(dp), parameter :: orientation(4) = [1.0_dp, 1.0_dp, -1.0_dp, -1.0_dp]

        status = 1
        if (points_per_cell < 1) return
        call gauss_rule(points_per_cell, x, w, status)
        if (status /= 0) return
        do pass = 1, 2
            n = 0
            do cell = 1, space%cell_count
                patch = space%cell_block(cell)
                k = space%patch_block(patch)
                e = space%cell_element(:, cell)
                associate (blk => space%blocks(patch))
                    do s = 1, 4
                        if (space%layout%side_tag(s, k) /= tag) cycle
                        if (.not. on_side(s, e, blk%element_count)) cycle
                        if (pass == 1) then
                            n = n + points_per_cell
                            cycle
                        end if
                        along = merge(1, 2, s == 1 .or. s == 3)
                        if (along == 1) then
                            t0 = blk%breaks1(e(1))
                            h = blk%breaks1(e(1) + 1) - t0
                        else
                            t0 = blk%breaks2(e(2))
                            h = blk%breaks2(e(2) + 1) - t0
                        end if
                        do q = 1, points_per_cell
                            t = t0 + h*0.5_dp*(x(q) + 1.0_dp)
                            xi = side_xi(s, t)
                            call spline_geometry_at(space, k, xi, jet, ierr)
                            if (ierr /= 0) return
                            tangent = orientation(s)*jet(:, 1 + along)
                            speed = norm2(tangent)
                            n = n + 1
                            side_quad%block(n) = k
                            side_quad%side(n) = s
                            side_quad%xi(:, n) = xi
                            side_quad%s(:, n) = jet(:, 1)
                            side_quad%tangent(:, n) = tangent/speed
                            side_quad%weight(n) = speed*h*0.5_dp*w(q)
                        end do
                    end do
                end associate
            end do
            if (pass == 1) then
                side_quad%point_count = n
                allocate(side_quad%block(n), side_quad%side(n), side_quad%xi(2, n), &
                    side_quad%s(2, n), side_quad%tangent(2, n), side_quad%weight(n))
            end if
        end do
        status = 0
    end subroutine spline_side_quadrature

    subroutine spline_point_matrix(space, block, xi, component, rowptr, colidx, &
            val, status)
        !! Sparse rows M(q, I) = jet component of phi_I at (block(q), xi(:, q)),
        !! component 0 for the value, 1, 2 for d/ds1, d/ds2 and 3, 4, 5 for
        !! (11, 21, 22). Columns are sorted within each row.
        type(spline_space_t), intent(in) :: space
        integer, intent(in) :: block(:), component
        real(dp), intent(in) :: xi(:, :)
        integer, allocatable, intent(out) :: rowptr(:), colidx(:)
        real(dp), allocatable, intent(out) :: val(:)
        integer, intent(out) :: status

        real(dp), allocatable :: bern(:, :, :), table1(:, :, :), table2(:, :, :)
        real(dp), allocatable :: u(:, :), g(:, :), xi_jets(:, :, :), g_jets(:, :, :)
        real(dp), allocatable :: jet(:), local(:), row(:), grow_val(:)
        integer, allocatable :: cell_broken(:), cols(:), grow_col(:)
        logical, allocatable :: used(:)
        real(dp) :: chain(6, 6), h(2), t(2)
        integer :: np, nloc, p, q, k, a, j, cell, used_count, nnz, ierr

        status = 1
        np = size(block)
        p = space%degree
        nloc = (p + 1)**2
        if (component < 0 .or. component > 5 .or. size(xi, 1) /= 2 .or. &
            size(xi, 2) /= np) return
        if (any(block < 1) .or. any(block > space%layout_block_count)) return
        if (any(xi < 0.0_dp) .or. any(xi > 1.0_dp)) return
        allocate(bern(p + 1, 2, 0:2), table1(p + 1, 1, 0:2), table2(p + 1, 1, 0:2), &
            u(nloc, nloc), g(2, nloc), xi_jets(nloc, 6, 1), g_jets(2, 6, 1), &
            jet(nloc*(component + 1)), cell_broken(nloc), row(space%dof_count), &
            used(space%dof_count), cols(space%dof_count), rowptr(np + 1), &
            colidx(16*np), val(16*np))
        u = 0.0_dp
        do a = 1, nloc
            u(a, a) = 1.0_dp
        end do
        row = 0.0_dp
        used = .false.
        rowptr(1) = 1
        nnz = 0
        do q = 1, np
            call point_cell(space, block(q), xi(:, q), cell, t)
            if (cell == 0) return
            do k = 1, 2
                call interval_bernstein_jet(p, t(k), bern(:, k, 0), bern(:, k, 1), &
                    bern(:, k, 2), status)
                if (status /= 0) return
            end do
            call cell_tables_at(space, cell, bern(:, 1:1, :), bern(:, 2:2, :), &
                table1, table2, h)
            call geometry_coefficients(space, cell, g)
            call spline_cell_xi_jets(table1, table2, u, xi_jets)
            call spline_cell_xi_jets(table1, table2, g, g_jets)
            call spline_xi_to_s_chain(g_jets(:, :, 1), chain, ierr)
            status = 2
            if (ierr /= 0) return
            call apply_chain(chain, xi_jets(:, :, 1), component + 1, jet)
            local = jet(nloc*component + 1:nloc*(component + 1))
            call spline_cell_broken(space, cell, cell_broken)
            used_count = 0
            do a = 1, nloc
                do k = space%broken_ptr(cell_broken(a)), &
                        space%broken_ptr(cell_broken(a) + 1) - 1
                    j = space%broken_dof(k)
                    if (.not. used(j)) then
                        used(j) = .true.
                        used_count = used_count + 1
                        cols(used_count) = j
                    end if
                    row(j) = row(j) + space%broken_weight(k)*local(a)
                end do
            end do
            call sort_integers(cols(:used_count))
            if (nnz + used_count > size(colidx)) then
                allocate(grow_col(2*(nnz + used_count)), grow_val(2*(nnz + used_count)))
                grow_col(:nnz) = colidx(:nnz)
                grow_val(:nnz) = val(:nnz)
                call move_alloc(grow_col, colidx)
                call move_alloc(grow_val, val)
            end if
            do k = 1, used_count
                j = cols(k)
                nnz = nnz + 1
                colidx(nnz) = j
                val(nnz) = row(j)
                row(j) = 0.0_dp
                used(j) = .false.
            end do
            rowptr(q + 1) = nnz + 1
        end do
        colidx = colidx(:nnz)
        val = val(:nnz)
        status = 0
    end subroutine spline_point_matrix

    subroutine spline_boundary_nodes(space, tag, panels_per_cell, block, xi, s, status)
        !! Nodes of a closed counter-clockwise curve of quadratic panels on
        !! the sides with the given tag: each cell side is split into
        !! panels_per_cell equal parameter intervals, and panel p has the
        !! nodes 2p - 1 (start) and 2p (parameter mid point); its end is the
        !! start of panel p + 1. Status 1: invalid input, 2: the tagged sides
        !! do not form one closed loop.
        type(spline_space_t), intent(in) :: space
        integer, intent(in) :: tag, panels_per_cell
        integer, allocatable, intent(out) :: block(:)
        real(dp), allocatable, intent(out) :: xi(:, :), s(:, :)
        integer, intent(out) :: status

        integer, allocatable :: sides(:, :), order(:)
        real(dp), allocatable :: breaks(:)
        real(dp) :: jet(2, 6), t0, t1, h
        integer :: nside, k, sd, i, n, e, m, first, last, along, node, ierr, c
        logical :: forward

        status = 1
        if (panels_per_cell < 1) return
        nside = count(space%layout%side_tag == tag)
        if (nside == 0) return
        allocate(sides(2, nside), order(nside))
        i = 0
        do k = 1, space%layout_block_count
            do sd = 1, 4
                if (space%layout%side_tag(sd, k) /= tag) cycle
                i = i + 1
                sides(:, i) = [k, sd]
            end do
        end do
        status = 2
        order(1) = 1
        do i = 2, nside
            order(i) = 0
            do n = 1, nside
                if (any(order(:i - 1) == n)) cycle
                if (ccw_vertex(sides(:, n), .true.) == &
                    ccw_vertex(sides(:, order(i - 1)), .false.)) order(i) = n
            end do
            if (order(i) == 0) return
        end do
        if (ccw_vertex(sides(:, order(nside)), .false.) /= &
            ccw_vertex(sides(:, 1), .true.)) return
        n = 0
        do i = 1, nside
            call side_breaks(space, sides(1, i), sides(2, i), breaks)
            n = n + 2*panels_per_cell*(size(breaks) - 1)
        end do
        allocate(block(n), xi(2, n), s(2, n))
        node = 0
        do i = 1, nside
            k = sides(1, order(i))
            sd = sides(2, order(i))
            forward = sd <= 2
            call side_breaks(space, k, sd, breaks)
            m = size(breaks) - 1
            first = merge(1, m, forward)
            last = merge(m, 1, forward)
            do e = first, last, merge(1, -1, forward)
                t0 = merge(breaks(e), breaks(e + 1), forward)
                t1 = merge(breaks(e + 1), breaks(e), forward)
                h = (t1 - t0)/panels_per_cell
                do c = 0, 2*panels_per_cell - 1
                    node = node + 1
                    block(node) = k
                    xi(:, node) = side_xi(sd, t0 + 0.5_dp*c*h)
                    call spline_geometry_at(space, k, xi(:, node), jet, ierr)
                    if (ierr /= 0) return
                    s(:, node) = jet(:, 1)
                end do
            end do
        end do
        status = 0

    contains

        integer function ccw_vertex(side, start) result(v)
            integer, intent(in) :: side(2)
            logical, intent(in) :: start

            if (start .eqv. side(2) <= 2) then
                v = space%layout%corners(side_start(side(2)), side(1))
            else
                v = space%layout%corners(side_end(side(2)), side(1))
            end if
        end function ccw_vertex

    end subroutine spline_boundary_nodes

    pure function side_xi(side, t) result(xi)
        integer, intent(in) :: side
        real(dp), intent(in) :: t
        real(dp) :: xi(2)

        select case (side)
        case (1)
            xi = [t, 0.0_dp]
        case (2)
            xi = [1.0_dp, t]
        case (3)
            xi = [t, 1.0_dp]
        case default
            xi = [0.0_dp, t]
        end select
        xi = min(max(xi, 0.0_dp), 1.0_dp)
    end function side_xi

    pure subroutine sort_integers(a)
        integer, intent(inout) :: a(:)
        integer :: i, j, x

        do i = 2, size(a)
            x = a(i)
            j = i - 1
            do while (j >= 1)
                if (a(j) <= x) exit
                a(j + 1) = a(j)
                j = j - 1
            end do
            a(j + 1) = x
        end do
    end subroutine sort_integers

    subroutine point_cell(space, block, xi, cell, t)
        !! Active cell of layout block `block` containing xi, searched from the
        !! finest patch down, and the local cell coordinates t. cell = 0 if
        !! none is found.
        type(spline_space_t), intent(in) :: space
        integer, intent(in) :: block
        real(dp), intent(in) :: xi(2)
        integer, intent(out) :: cell
        real(dp), intent(out) :: t(2)

        integer :: e(2), patch

        cell = 0
        t = 0.0_dp
        do patch = size(space%patch_block), 1, -1
            if (space%patch_block(patch) /= block) cycle
            associate (blk => space%blocks(patch))
                e(1) = locate(blk%breaks1, xi(1))
                e(2) = locate(blk%breaks2, xi(2))
                cell = space%element_cell(space%patch_cell_offset(patch) + e(1) + &
                    blk%element_count(1)*(e(2) - 1))
                if (cell == 0) cycle
                t(1) = (xi(1) - blk%breaks1(e(1)))/(blk%breaks1(e(1) + 1) - &
                    blk%breaks1(e(1)))
                t(2) = (xi(2) - blk%breaks2(e(2)))/(blk%breaks2(e(2) + 1) - &
                    blk%breaks2(e(2)))
            end associate
            return
        end do
    end subroutine point_cell

    subroutine side_breaks(space, block, side, breaks)
        !! Sorted end points of the active cells along a side of a layout block.
        type(spline_space_t), intent(in) :: space
        integer, intent(in) :: block, side
        real(dp), allocatable, intent(out) :: breaks(:)

        real(dp), allocatable :: ends(:)
        integer :: cell, patch, e(2), i, n
        real(dp) :: x

        allocate(ends(0))
        do cell = 1, space%cell_count
            patch = space%cell_block(cell)
            if (space%patch_block(patch) /= block) cycle
            e = space%cell_element(:, cell)
            if (.not. on_side(side, e, space%blocks(patch)%element_count)) cycle
            if (side == 1 .or. side == 3) then
                ends = [ends, space%blocks(patch)%breaks1(e(1):e(1) + 1)]
            else
                ends = [ends, space%blocks(patch)%breaks2(e(2):e(2) + 1)]
            end if
        end do
        do i = 2, size(ends)
            x = ends(i)
            n = i - 1
            do while (n >= 1)
                if (ends(n) <= x) exit
                ends(n + 1) = ends(n)
                n = n - 1
            end do
            ends(n + 1) = x
        end do
        allocate(breaks(size(ends)))
        n = 0
        do i = 1, size(ends)
            if (n > 0) then
                if (ends(i) - breaks(n) <= 1.0e-14_dp) cycle
            end if
            n = n + 1
            breaks(n) = ends(i)
        end do
        breaks = breaks(:n)
    end subroutine side_breaks

    pure logical function on_side(side, element, element_count) result(on)
        integer, intent(in) :: side, element(2), element_count(2)

        select case (side)
        case (1)
            on = element(2) == 1
        case (2)
            on = element(1) == element_count(1)
        case (3)
            on = element(2) == element_count(2)
        case default
            on = element(1) == 1
        end select
    end function on_side

    subroutine cell_tables_at(space, cell, bernstein1, bernstein2, table1, table2, h)
        type(spline_space_t), intent(in) :: space
        integer, intent(in) :: cell
        real(dp), intent(in) :: bernstein1(:, :, 0:), bernstein2(:, :, 0:)
        real(dp), intent(out) :: table1(:, :, 0:), table2(:, :, 0:), h(2)

        integer :: k, e1, e2, d

        k = space%cell_block(cell)
        e1 = space%cell_element(1, cell)
        e2 = space%cell_element(2, cell)
        associate (blk => space%blocks(k))
            h(1) = blk%breaks1(e1 + 1) - blk%breaks1(e1)
            h(2) = blk%breaks2(e2 + 1) - blk%breaks2(e2)
            do d = 0, 2
                table1(:, :, d) = matmul(blk%extraction1(:, :, e1), &
                    bernstein1(:, :, d))/h(1)**d
                table2(:, :, d) = matmul(blk%extraction2(:, :, e2), &
                    bernstein2(:, :, d))/h(2)**d
            end do
        end associate
    end subroutine cell_tables_at

    subroutine geometry_coefficients(space, cell, u)
        type(spline_space_t), intent(in) :: space
        integer, intent(in) :: cell
        real(dp), intent(out) :: u(:, :)

        integer :: cell_broken(size(u, 2))

        call spline_cell_broken(space, cell, cell_broken)
        u = space%geometry(:, cell_broken)
    end subroutine geometry_coefficients

    pure integer function locate(breaks, x) result(e)
        real(dp), intent(in) :: breaks(:), x

        do e = 1, size(breaks) - 2
            if (x < breaks(e + 1)) return
        end do
        e = size(breaks) - 1
    end function locate

end module fortfem_spline_eval
