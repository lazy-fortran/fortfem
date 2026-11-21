program test_spline_space_c0
    !! C0 multipatch spline space on the 5-block O-grid disc: layout
    !! topology, DOF counts, disc area, exact jets of the reference
    !! coordinates, continuity across block sides, and second-order jets
    !! against finite differences in the block parameters.
    use check, only: check_condition, check_summary
    use fortfem_kinds, only: dp
    use fortfem_spline_layout, only: block_layout_t, layout_ogrid, LAYOUT_WALL
    use fortfem_spline_space, only: spline_space_t, spline_space_options_t, &
        spline_space_build, spline_boundary_dofs
    use fortfem_spline_eval, only: spline_quadrature_t, spline_quadrature_build, &
        spline_evaluate, spline_evaluate_at
    implicit none

    type(block_layout_t) :: layout
    type(spline_space_t) :: space
    type(spline_quadrature_t) :: quad
    real(dp), allocatable :: coef(:, :), jet(:, :)
    integer, allocatable :: wall(:)
    real(dp), parameter :: pi = acos(-1.0_dp)
    integer :: status, k, n

    call layout_ogrid([0.0_dp, 0.0_dp], 1.0_dp, 0.5_dp, layout, status)
    call check_condition(status == 0 .and. layout%block_count == 5 .and. &
        all(layout%valence(1:4) == 3) .and. all(layout%valence(5:8) == 2) .and. &
        all(layout%boundary_vertex(5:8)) .and. .not. any(layout%boundary_vertex(1:4)), &
        "O-grid has four valence-3 interior vertices and a wall of four arcs")

    call ogrid_space(3, 4, space, status)
    n = 7
    call check_condition(status == 0 .and. space%broken_count == 5*n*n .and. &
        space%dof_count == 5*(n - 2)**2 + 12*(n - 2) + 8, &
        "C0 O-grid space has the conforming DOF count")
    call spline_boundary_dofs(space, LAYOUT_WALL, wall)
    call check_condition(size(wall) == 4*(n - 2) + 4, &
        "wall DOFs are the traces on the four arcs")

    call spline_quadrature_build(space, 4, quad, status)
    call check_condition(status == 0 .and. area_rate() > 3.8_dp, &
        "disc area error of the arc blocks converges at order p + 1")

    allocate(coef(2, space%dof_count))
    call coordinate_coefficients(space, coef)
    allocate(jet(2*6, quad%point_count))
    call spline_evaluate(space, quad, coef, 2, jet, status)
    call check_condition(status == 0 .and. &
        maxval(abs(jet(1:2, :) - quad%s)) < 1.0e-13_dp .and. &
        maxval(abs(jet(3, :) - 1.0_dp)) < 1.0e-12_dp .and. &
        maxval(abs(jet(4, :))) < 1.0e-12_dp .and. &
        maxval(abs(jet(5, :))) < 1.0e-12_dp .and. &
        maxval(abs(jet(6, :) - 1.0_dp)) < 1.0e-12_dp .and. &
        maxval(abs(jet(7:12, :))) < 1.0e-10_dp, &
        "jets of the reference coordinates are exact through curved blocks")

    call random_coefficients(space%dof_count, coef)
    call check_condition(side_jump(space, layout, coef) < 1.0e-13_dp, &
        "random C0 fields are continuous across every block side")
    do k = 1, 5
        call check_condition(fd_error(space, coef, k) < 2.0e-5_dp, &
            "second-order jets match finite differences in block parameters")
    end do
    call check_summary("spline_space_c0")

contains

    subroutine ogrid_space(degree, cells, space, status)
        integer, intent(in) :: degree, cells
        type(spline_space_t), intent(out) :: space
        integer, intent(out) :: status

        type(spline_space_options_t) :: options

        options%degree = degree
        allocate(options%knots(2, 5))
        options%knots(:, :)%element_count = cells
        call spline_space_build(layout, options, space, status)
    end subroutine ogrid_space

    subroutine coordinate_coefficients(space, coef)
        !! Global coefficients of s: the block maps are continuous, so any
        !! broken representative of a global function carries its value.
        type(spline_space_t), intent(in) :: space
        real(dp), intent(out) :: coef(:, :)
        integer :: j

        do j = 1, space%broken_count
            coef(:, space%broken_dof(space%broken_ptr(j))) = space%geometry(:, j)
        end do
    end subroutine coordinate_coefficients

    subroutine random_coefficients(n, coef)
        integer, intent(in) :: n
        real(dp), intent(out) :: coef(:, :)
        integer :: seed_size
        integer, allocatable :: seed(:)

        call random_seed(size=seed_size)
        allocate(seed(seed_size))
        seed = 20261121
        call random_seed(put=seed)
        call random_number(coef(:, :n))
    end subroutine random_coefficients

    real(dp) function area_rate() result(rate)
        type(spline_space_t) :: coarse, fine
        type(spline_quadrature_t) :: qc, qf
        integer :: status

        call ogrid_space(3, 8, coarse, status)
        call ogrid_space(3, 16, fine, status)
        call spline_quadrature_build(coarse, 4, qc, status)
        call spline_quadrature_build(fine, 4, qf, status)
        rate = log(abs(sum(qc%weight) - pi)/abs(sum(qf%weight) - pi))/log(2.0_dp)
    end function area_rate

    real(dp) function side_jump(space, layout, coef) result(jump)
        type(spline_space_t), intent(in) :: space
        type(block_layout_t), intent(in) :: layout
        real(dp), intent(in) :: coef(:, :)

        real(dp) :: mine(2), theirs(2), t
        integer :: k, s, k2, s2, q, status

        jump = 0.0_dp
        do k = 1, layout%block_count
            do s = 1, 4
                k2 = layout%neighbour(s, k)
                if (k2 == 0) cycle
                s2 = layout%neighbour_side(s, k)
                do q = 1, 7
                    t = real(q, dp)/8.0_dp
                    call spline_evaluate_at(space, coef, k, side_xi(s, t), 0, &
                        mine, status)
                    if (layout%reversed(s, k)) then
                        call spline_evaluate_at(space, coef, k2, &
                            side_xi(s2, 1.0_dp - t), 0, theirs, status)
                    else
                        call spline_evaluate_at(space, coef, k2, side_xi(s2, t), &
                            0, theirs, status)
                    end if
                    jump = max(jump, maxval(abs(mine - theirs)))
                end do
            end do
        end do
    end function side_jump

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
    end function side_xi

    real(dp) function fd_error(space, coef, block) result(error)
        !! d2u/dxi_i dxi_j = DG_i^T H_s u DG_j + grad_s u . d2G/dxi_i dxi_j,
        !! with DG and d2G from differences of s and d2u from differences
        !! of u.
        type(spline_space_t), intent(in) :: space
        real(dp), intent(in) :: coef(:, :)
        integer, intent(in) :: block

        real(dp), parameter :: delta = 1.0e-4_dp
        real(dp) :: xi(2), jet(12), u(2, -1:1, -1:1), s(2, -1:1, -1:1)
        real(dp) :: scoef_jet(2), dg(2, 2), d2g(2, 2, 2), hs(2, 2), d2u, ref
        real(dp), allocatable :: scoef(:, :)
        integer :: i, j, f, a, b, status

        allocate(scoef(2, space%dof_count))
        call coordinate_coefficients(space, scoef)
        xi = [0.37_dp, 0.61_dp]
        call spline_evaluate_at(space, coef, block, xi, 2, jet, status)
        do b = -1, 1
            do a = -1, 1
                call spline_evaluate_at(space, coef, block, xi + delta*[a, b], 0, &
                    u(:, a, b), status)
                call spline_evaluate_at(space, scoef, block, xi + delta*[a, b], 0, &
                    scoef_jet, status)
                s(:, a, b) = scoef_jet
            end do
        end do
        dg(:, 1) = (s(:, 1, 0) - s(:, -1, 0))/(2*delta)
        dg(:, 2) = (s(:, 0, 1) - s(:, 0, -1))/(2*delta)
        d2g(:, 1, 1) = (s(:, 1, 0) - 2*s(:, 0, 0) + s(:, -1, 0))/delta**2
        d2g(:, 2, 2) = (s(:, 0, 1) - 2*s(:, 0, 0) + s(:, 0, -1))/delta**2
        d2g(:, 1, 2) = (s(:, 1, 1) - s(:, 1, -1) - s(:, -1, 1) + s(:, -1, -1))/ &
            (4*delta**2)
        d2g(:, 2, 1) = d2g(:, 1, 2)
        error = 0.0_dp
        do f = 1, 2
            hs = reshape([jet(f + 6), jet(f + 8), jet(f + 8), jet(f + 10)], [2, 2])
            do j = 1, 2
                do i = 1, 2
                    ref = dot_product(dg(:, i), matmul(hs, dg(:, j))) + &
                        jet(f + 2)*d2g(1, i, j) + jet(f + 4)*d2g(2, i, j)
                    if (i == 1 .and. j == 1) d2u = (u(f, 1, 0) - 2*u(f, 0, 0) + &
                        u(f, -1, 0))/delta**2
                    if (i == 2 .and. j == 2) d2u = (u(f, 0, 1) - 2*u(f, 0, 0) + &
                        u(f, 0, -1))/delta**2
                    if (i /= j) d2u = (u(f, 1, 1) - u(f, 1, -1) - u(f, -1, 1) + &
                        u(f, -1, -1))/(4*delta**2)
                    error = max(error, abs(d2u - ref))
                end do
            end do
        end do
    end function fd_error

end program test_spline_space_c0
