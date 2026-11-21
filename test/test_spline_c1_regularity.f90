module spline_c1_regularity_kernel
    !! W = |grad u|^2/2 + u^2/2 - f u on the disc |s - c| <= 1: u - lap u = f
    !! with zero normal derivative for u = cos(pi r) + (1 - r^2/3) (s1 - c1).
    use fortfem_kinds, only: dp
    use fortfem_spline_energy, only: pointwise_energy_t
    implicit none
    private

    public :: disc_kernel_t, exact_solution, source

    real(dp), parameter :: pi = acos(-1.0_dp)
    real(dp), parameter, public :: center(2) = [0.3_dp, -0.2_dp]

    type, extends(pointwise_energy_t) :: disc_kernel_t
        real(dp), allocatable :: f(:)
    contains
        procedure :: evaluate
    end type disc_kernel_t

contains

    subroutine evaluate(self, cell, first_point, npoint, jet, weight, need, w, &
            gradient, hessian, status)
        class(disc_kernel_t), intent(in) :: self
        integer, intent(in) :: cell, first_point, npoint, need
        real(dp), intent(in) :: jet(:, :), weight(:)
        real(dp), intent(out) :: w(:), gradient(:, :), hessian(:, :)
        integer, intent(out) :: status

        real(dp) :: f
        integer :: q

        associate (unused => cell, unused_weight => weight)
        end associate
        do q = 1, npoint
            f = self%f(first_point + q - 1)
            w(q) = 0.5_dp*(jet(2, q)**2 + jet(3, q)**2 + jet(1, q)**2) - f*jet(1, q)
            if (need >= 1) gradient(:, q) = [jet(1, q) - f, jet(2, q), jet(3, q)]
            if (need >= 2) hessian(:, q) = [1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, &
                0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp]
        end do
        status = 0
    end subroutine evaluate

    pure real(dp) function exact_solution(s) result(u)
        real(dp), intent(in) :: s(2)

        real(dp) :: d(2), r

        d = s - center
        r = norm2(d)
        u = cos(pi*r) + (1.0_dp - r**2/3.0_dp)*d(1)
    end function exact_solution

    pure real(dp) function source(s) result(f)
        real(dp), intent(in) :: s(2)

        real(dp) :: d(2), r, lap

        d = s - center
        r = norm2(d)
        if (r < 1.0e-12_dp) then
            lap = -2.0_dp*pi**2
        else
            lap = -pi**2*cos(pi*r) - pi*sin(pi*r)/r
        end if
        lap = lap - 8.0_dp*d(1)/3.0_dp
        f = exact_solution(s) - lap
    end function source

end module spline_c1_regularity_kernel

program test_spline_c1_regularity
    !! Per-edge regularity in C1 spaces: C0 sides next to C1 sides with mixed
    !! vertices (doc/design/multipatch_spline_spaces.md, C1 construction).
    !! Polygonal O-grid with C0 ring-ring sides: dimension, rank, C1 across
    !! the C1 sides and gradient jumps across the C0 sides. Curved O-grid
    !! and C1 core inside a periodic annulus: continuity and order p + 1.
    use check, only: check_condition, check_summary
    use fortfem_kinds, only: dp
    use fortfem_element_to_csr, only: csr_matrix_t
    use fortfem_spline_layout, only: block_layout_t, layout_ogrid, layout_set_c0, &
        layout_set_curve, LAYOUT_WALL, LAYOUT_ANNULUS_INTERFACE
    use fortfem_spline_space, only: spline_space_t, spline_space_options_t, &
        spline_space_build, spline_space_update_geometry, SPLINE_C1
    use fortfem_spline_eval, only: spline_quadrature_t, spline_quadrature_build, &
        spline_evaluate, spline_evaluate_at
    use fortfem_spline_energy, only: spline_energy, spline_energy_hessian_csr
    use spline_c1_regularity_kernel, only: disc_kernel_t, exact_solution, source, &
        center
    implicit none

    interface
        subroutine dsyev(jobz, uplo, n, a, lda, w, work, lwork, info)
            import :: dp
            character, intent(in) :: jobz, uplo
            integer, intent(in) :: n, lda, lwork
            real(dp), intent(inout) :: a(lda, *)
            real(dp), intent(out) :: w(*), work(*)
            integer, intent(out) :: info
        end subroutine dsyev
    end interface

    type(block_layout_t) :: layout
    type(spline_space_t) :: space
    real(dp) :: c1_jump, c0_jump, value_jump, errors(3), rate
    integer :: status, k, level, cells, kst

    call layout_ogrid([0.0_dp, 0.0_dp], sqrt(2.0_dp), 0.5_dp, layout, status, &
        polygon=.true.)
    do k = 2, 5
        call layout_set_c0(layout, k, 2, .true., status)
    end do
    cells = 4
    call build(3, cells, 5, space, status)
    kst = 5*(2*(cells - 1) + 4 - 4)**2 + 12*(2*((cells - 1) + 4) - 11) + 6*8
    call check_condition(status == 0 .and. space%dof_count == kst + &
        4*(cells + 3 - 5) + 4*3, &
        "C0 sides add their transversal functions and three per split vertex")
    call check_condition(extraction_rank_ok(space), "extraction has full row rank")
    call jumps(space, 5, c1_jump, c0_jump, value_jump)
    call check_condition(value_jump < 1.0e-11_dp .and. c1_jump < 1.0e-11_dp .and. &
        c0_jump > 1.0e-2_dp, &
        "values continuous, gradients continuous on C1 and free on C0 sides")
    call check_condition(size(space%c0_interfaces, 2) == 4 .and. &
        all(space%c0_interfaces(1, :) >= 2 .and. space%c0_interfaces(3, :) >= 2) &
        .and. all(mod(space%c0_interfaces(2, :), 2) == 0), &
        "the four ring-ring sides are listed")

    call layout_ogrid(center, 1.0_dp, 0.3_dp, layout, status)
    call build(3, 4, 5, space, status)
    call jumps(space, 5, c1_jump, c0_jump, value_jump)
    call check_condition(status == 0 .and. value_jump < 1.0e-11_dp .and. &
        c1_jump < 1.0e-11_dp, "curved O-grid: C1 except on the ring-ring sides")
    do level = 1, 3
        call solve(5, 2**(level + 1), errors(level))
    end do
    rate = log(errors(2)/errors(3))/log(2.0_dp)
    print '(a,3es11.3,a,f6.2)', "curved O-grid C1 L2 errors", errors, " rate", rate
    call check_condition(rate > 3.7_dp, "curved O-grid C1 solution converges")

    call layout_ogrid(center, 1.0_dp, 0.3_dp, layout, status, &
        interface_fraction=0.6_dp)
    call build(3, 4, 6, space, status)
    call jumps(space, 5, c1_jump, c0_jump, value_jump)
    call check_condition(status == 0 .and. value_jump < 1.0e-11_dp .and. &
        c1_jump < 1.0e-11_dp .and. interface_jump(space) < 1.0e-11_dp .and. &
        size(space%c0_interfaces, 2) == 8, &
        "C1 core in the annulus: C1 inside, C0 on ring-ring sides and interface")
    call check_condition(extraction_rank_ok(space), &
        "annulus coupling keeps full row rank")
    do level = 1, 3
        call solve(6, 2**(level + 1), errors(level))
    end do
    rate = log(errors(2)/errors(3))/log(2.0_dp)
    print '(a,3es11.3,a,f6.2)', "C1 core + annulus L2 errors", errors, " rate", rate
    call check_condition(rate > 3.7_dp, "C1 core with annulus converges")

    call layout_set_curve(layout, LAYOUT_WALL, d_shape(0.8_dp, 1.8_dp, 0.5_dp), &
        status)
    call layout_set_curve(layout, LAYOUT_ANNULUS_INTERFACE, &
        d_shape(0.5_dp, 1.8_dp, 0.5_dp), status)
    call build(3, 8, 6, space, status)
    call jumps(space, 5, c1_jump, c0_jump, value_jump)
    call check_condition(status == 0 .and. value_jump < 1.0e-11_dp .and. &
        c1_jump < 1.0e-10_dp .and. interface_jump(space) < 1.0e-11_dp, &
        "C1 core in a D-shaped annulus")
    call layout_set_curve(layout, LAYOUT_ANNULUS_INTERFACE, &
        d_shape(0.45_dp, 1.6_dp, 0.3_dp), status)
    call spline_space_update_geometry(space, layout, status)
    call jumps(space, 5, c1_jump, c0_jump, value_jump)
    call check_condition(status == 0 .and. value_jump < 1.0e-11_dp .and. &
        c1_jump < 1.0e-10_dp .and. interface_jump(space) < 1.0e-11_dp, &
        "C1 stays C1 after the interface moves")
    call check_summary("spline_c1_regularity")

contains

    subroutine build(degree, cells, nb, space, status)
        integer, intent(in) :: degree, cells, nb
        type(spline_space_t), intent(out) :: space
        integer, intent(out) :: status

        type(spline_space_options_t) :: options

        options%degree = degree
        options%regularity = SPLINE_C1
        allocate(options%knots(2, nb))
        options%knots(:, :)%element_count = cells
        call spline_space_build(layout, options, space, status)
    end subroutine build

    logical function extraction_rank_ok(space) result(ok)
        type(spline_space_t), intent(in) :: space

        real(dp), allocatable :: e(:, :), gram(:, :), w(:), work(:)
        integer :: j, q, n, info

        n = space%dof_count
        allocate(e(n, space%broken_count))
        e = 0.0_dp
        do j = 1, space%broken_count
            do q = space%broken_ptr(j), space%broken_ptr(j + 1) - 1
                e(space%broken_dof(q), j) = e(space%broken_dof(q), j) + &
                    space%broken_weight(q)
            end do
        end do
        gram = matmul(e, transpose(e))
        allocate(w(n), work(10*n))
        call dsyev('N', 'U', n, gram, n, w, work, 10*n, info)
        ok = info == 0 .and. w(1) > 1.0e-10_dp*w(n)
    end function extraction_rank_ok

    subroutine jumps(space, core_blocks, c1_jump, c0_jump, value_jump)
        !! Largest jumps of value and gradient of a random field across the
        !! interior sides of the core blocks, split by the side coupling.
        type(spline_space_t), intent(in) :: space
        integer, intent(in) :: core_blocks
        real(dp), intent(out) :: c1_jump, c0_jump, value_jump

        real(dp), allocatable :: coef(:, :)
        real(dp) :: mine(3), theirs(3), t
        integer :: k, s, k2, s2, q, status

        allocate(coef(1, space%dof_count))
        call random_number(coef)
        c1_jump = 0.0_dp
        c0_jump = 0.0_dp
        value_jump = 0.0_dp
        do k = 1, core_blocks
            do s = 1, 4
                k2 = layout%neighbour(s, k)
                if (k2 == 0) cycle
                s2 = layout%neighbour_side(s, k)
                do q = 0, 8
                    t = real(q, dp)/8.0_dp
                    call spline_evaluate_at(space, coef, k, side_xi(s, t), 1, mine, &
                        status)
                    if (layout%reversed(s, k)) t = 1.0_dp - t
                    call spline_evaluate_at(space, coef, k2, side_xi(s2, t), 1, &
                        theirs, status)
                    value_jump = max(value_jump, abs(mine(1) - theirs(1)))
                    if (layout%side_c0(s, k)) then
                        c0_jump = max(c0_jump, maxval(abs(mine(2:3) - theirs(2:3))))
                    else
                        c1_jump = max(c1_jump, maxval(abs(mine(2:3) - theirs(2:3)))/ &
                            max(1.0_dp, maxval(abs(mine(2:3)))))
                    end if
                end do
            end do
        end do
    end subroutine jumps

    real(dp) function interface_jump(space) result(jump)
        type(spline_space_t), intent(in) :: space

        real(dp), allocatable :: coef(:, :)
        real(dp) :: core(1), ring(1), t
        integer :: k, q, status

        allocate(coef(1, space%dof_count))
        call random_number(coef)
        jump = 0.0_dp
        do k = 2, 5
            do q = 0, 10
                t = real(q, dp)/10.0_dp
                call spline_evaluate_at(space, coef, k, [t, 0.0_dp], 0, core, status)
                call spline_evaluate_at(space, coef, 6, &
                    [modulo(0.25_dp*(k - 2 + t), 1.0_dp), 1.0_dp], 0, ring, status)
                jump = max(jump, abs(core(1) - ring(1)))
            end do
        end do
    end function interface_jump

    subroutine solve(nb, cells, error)
        integer, intent(in) :: nb, cells
        real(dp), intent(out) :: error

        type(spline_space_t) :: space
        type(spline_quadrature_t) :: quad
        type(disc_kernel_t) :: kernel
        type(csr_matrix_t) :: matrix
        real(dp), allocatable :: c(:, :), g(:, :), jet(:, :)
        real(dp) :: energy
        integer :: status, q

        call build(3, cells, nb, space, status)
        if (status /= 0) error stop "space"
        call spline_quadrature_build(space, 5, quad, status)
        allocate(kernel%f(quad%point_count))
        do q = 1, quad%point_count
            kernel%f(q) = source(quad%s(:, q))
        end do
        allocate(c(1, space%dof_count), g(1, space%dof_count))
        c = 0.0_dp
        call spline_energy(space, quad, kernel, c, energy, status, g)
        call spline_energy_hessian_csr(space, quad, kernel, c, matrix, status)
        call pcg(matrix, -g(1, :), c(1, :))
        allocate(jet(1, quad%point_count))
        call spline_evaluate(space, quad, c, 0, jet, status)
        error = 0.0_dp
        do q = 1, quad%point_count
            error = error + quad%weight(q)*(jet(1, q) - exact_solution(quad%s(:, q)))**2
        end do
        error = sqrt(error)
    end subroutine solve

    pure function d_shape(minor, elongation, triangularity) result(samples)
        !! D shape about the layout centre, 200 samples from the outboard
        !! midplane, counter-clockwise.
        real(dp), intent(in) :: minor, elongation, triangularity
        real(dp) :: samples(2, 200), theta
        integer :: q

        do q = 1, 200
            theta = 2.0_dp*acos(-1.0_dp)*(q - 1)/200
            samples(:, q) = center + [minor*cos(theta + triangularity*sin(theta)), &
                elongation*minor*sin(theta)]
        end do
    end function d_shape

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

    subroutine pcg(matrix, b, x)
        type(csr_matrix_t), intent(in) :: matrix
        real(dp), intent(in) :: b(:)
        real(dp), intent(out) :: x(:)

        real(dp), allocatable :: r(:), z(:), p(:), ap(:), diag(:)
        real(dp) :: rz, rz_new, alpha, norm_b
        integer :: it, i, k

        allocate(r(size(b)), z(size(b)), p(size(b)), ap(size(b)), diag(size(b)))
        do i = 1, matrix%n
            diag(i) = 1.0_dp
            do k = matrix%rowptr(i), matrix%rowptr(i + 1) - 1
                if (matrix%colidx(k) == i) diag(i) = matrix%val(k)
            end do
        end do
        x = 0.0_dp
        r = b
        norm_b = norm2(r)
        z = r/diag
        p = z
        rz = dot_product(r, z)
        do it = 1, 20*size(b)
            call matvec(matrix, p, ap)
            alpha = rz/dot_product(p, ap)
            x = x + alpha*p
            r = r - alpha*ap
            if (norm2(r) < 1.0e-13_dp*norm_b) exit
            z = r/diag
            rz_new = dot_product(r, z)
            p = z + (rz_new/rz)*p
            rz = rz_new
        end do
    end subroutine pcg

    subroutine matvec(matrix, v, y)
        type(csr_matrix_t), intent(in) :: matrix
        real(dp), intent(in) :: v(:)
        real(dp), intent(out) :: y(:)
        integer :: row, k

        do row = 1, matrix%n
            y(row) = 0.0_dp
            do k = matrix%rowptr(row), matrix%rowptr(row + 1) - 1
                y(row) = y(row) + matrix%val(k)*v(matrix%colidx(k))
            end do
        end do
    end subroutine matvec

end program test_spline_c1_regularity
