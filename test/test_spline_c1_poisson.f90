module spline_c1_poisson_kernel
    !! W = |grad u|^2/2 + u^2/2 - f u on [-1, 1]^2, so the minimiser solves
    !! u - lap u = f with zero normal derivative. Exact solution
    !! u = cos(pi s1) cos(pi s2), f = (1 + 2 pi^2) u.
    use fortfem_kinds, only: dp
    use fortfem_spline_energy, only: pointwise_energy_t
    implicit none
    private

    public :: c1_kernel_t, exact_solution

    real(dp), parameter :: pi = acos(-1.0_dp)

    type, extends(pointwise_energy_t) :: c1_kernel_t
        real(dp), allocatable :: s(:, :)
    contains
        procedure :: evaluate
    end type c1_kernel_t

contains

    subroutine evaluate(self, cell, first_point, npoint, jet, weight, need, w, &
            gradient, hessian, status)
        class(c1_kernel_t), intent(in) :: self
        integer, intent(in) :: cell, first_point, npoint, need
        real(dp), intent(in) :: jet(:, :), weight(:)
        real(dp), intent(out) :: w(:), gradient(:, :), hessian(:, :)
        integer, intent(out) :: status

        real(dp) :: f
        integer :: q

        associate (unused => cell, unused_weight => weight)
        end associate
        do q = 1, npoint
            f = (1.0_dp + 2.0_dp*pi**2)*exact_solution(self%s(:, first_point + q - 1))
            w(q) = 0.5_dp*(jet(2, q)**2 + jet(3, q)**2 + jet(1, q)**2) - f*jet(1, q)
            if (need >= 1) gradient(:, q) = [jet(1, q) - f, jet(2, q), jet(3, q)]
            if (need >= 2) hessian(:, q) = [1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, &
                0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp]
        end do
        status = 0
    end subroutine evaluate

    pure real(dp) function exact_solution(s) result(u)
        real(dp), intent(in) :: s(2)

        u = cos(pi*s(1))*cos(pi*s(2))
    end function exact_solution

end module spline_c1_poisson_kernel

program test_spline_c1_poisson
    !! Gate 1 of doc/design/multipatch_spline_spaces.md: the C1 (AS-G1)
    !! space on the polygonal 5-block O-grid with four valence-3 vertices.
    !! The assembled E H E^T Hessian equals the matrix-free action, the
    !! discrete solution is C1 across every block side, and the L2 error
    !! converges at order p + 1 for p = 3, 4, 5 with r = p - 2.
    use check, only: check_condition, check_summary
    use fortfem_kinds, only: dp
    use fortfem_element_to_csr, only: csr_matrix_t
    use fortfem_spline_layout, only: block_layout_t, layout_ogrid
    use fortfem_spline_space, only: spline_space_t, spline_space_options_t, &
        spline_space_build, SPLINE_C1
    use fortfem_spline_eval, only: spline_quadrature_t, spline_quadrature_build, &
        spline_evaluate, spline_evaluate_at
    use fortfem_spline_energy, only: spline_energy, spline_energy_hessian_csr, &
        spline_energy_hessian_apply
    use spline_c1_poisson_kernel, only: c1_kernel_t, exact_solution
    implicit none

    type(block_layout_t) :: layout
    real(dp) :: errors(3), jumps(3), rate
    integer :: status, degree, level

    call layout_ogrid([0.0_dp, 0.0_dp], sqrt(2.0_dp), 0.5_dp, layout, status, &
        polygon=.true.)
    call check_condition(hessians_agree(), &
        "assembled E H E^T and matrix-free C1 Hessians agree")
    do degree = 3, 5
        do level = 1, 3
            call solve(degree, 2**(level + 1), errors(level), jumps(level))
        end do
        rate = log(errors(2)/errors(3))/log(2.0_dp)
        print '(a,i0,a,3es11.3,a,f6.2,a,es9.2)', "p = ", degree, " L2 errors", &
            errors, " rate", rate, " max C1 jump", maxval(jumps)
        call check_condition(rate > real(degree, dp) + 0.7_dp, &
            "C1 O-grid solution converges at order p + 1 in L2")
        call check_condition(maxval(jumps) < 1.0e-9_dp, &
            "C1 O-grid solution has continuous gradients across all sides")
    end do
    call check_summary("spline_c1_poisson")

contains

    subroutine setup(degree, cells, space, quad, kernel)
        integer, intent(in) :: degree, cells
        type(spline_space_t), intent(out) :: space
        type(spline_quadrature_t), intent(out) :: quad
        type(c1_kernel_t), intent(out) :: kernel

        type(spline_space_options_t) :: options
        integer :: status

        options%degree = degree
        options%regularity = SPLINE_C1
        allocate(options%knots(2, 5))
        options%knots(:, :)%element_count = cells
        call spline_space_build(layout, options, space, status)
        if (status /= 0) error stop "space"
        call spline_quadrature_build(space, degree + 1, quad, status)
        if (status /= 0) error stop "quadrature"
        kernel%s = quad%s
    end subroutine setup

    logical function hessians_agree() result(ok)
        type(spline_space_t) :: space
        type(spline_quadrature_t) :: quad
        type(c1_kernel_t) :: kernel
        type(csr_matrix_t) :: matrix
        real(dp), allocatable :: c(:, :), v(:, :), hv(:, :), cv(:)
        integer :: status, i, k

        call setup(3, 4, space, quad, kernel)
        allocate(c(1, space%dof_count), v(1, space%dof_count), &
            hv(1, space%dof_count), cv(space%dof_count))
        call random_number(c)
        call random_number(v)
        call spline_energy_hessian_csr(space, quad, kernel, c, matrix, status)
        ok = status == 0
        call spline_energy_hessian_apply(space, quad, kernel, c, v, hv, status)
        ok = ok .and. status == 0
        do i = 1, matrix%n
            cv(i) = 0.0_dp
            do k = matrix%rowptr(i), matrix%rowptr(i + 1) - 1
                cv(i) = cv(i) + matrix%val(k)*v(1, matrix%colidx(k))
            end do
        end do
        ok = ok .and. maxval(abs(cv - hv(1, :))) < 1.0e-11_dp*maxval(abs(cv))
    end function hessians_agree

    subroutine solve(degree, cells, error, jump)
        integer, intent(in) :: degree, cells
        real(dp), intent(out) :: error, jump

        type(spline_space_t) :: space
        type(spline_quadrature_t) :: quad
        type(c1_kernel_t) :: kernel
        type(csr_matrix_t) :: matrix
        real(dp), allocatable :: c(:, :), g(:, :), jet(:, :)
        real(dp) :: energy
        integer :: status, q

        call setup(degree, cells, space, quad, kernel)
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
        jump = gradient_jump(space, c)
    end subroutine solve

    real(dp) function gradient_jump(space, coef) result(jump)
        type(spline_space_t), intent(in) :: space
        real(dp), intent(in) :: coef(:, :)

        real(dp) :: mine(3), theirs(3), t
        integer :: k, s, k2, s2, q, status

        jump = 0.0_dp
        do k = 1, layout%block_count
            do s = 1, 4
                k2 = layout%neighbour(s, k)
                if (k2 == 0) cycle
                s2 = layout%neighbour_side(s, k)
                do q = 0, 16
                    t = real(q, dp)/16.0_dp
                    call spline_evaluate_at(space, coef, k, side_xi(s, t), 1, &
                        mine, status)
                    if (layout%reversed(s, k)) t = 1.0_dp - t
                    call spline_evaluate_at(space, coef, k2, side_xi(s2, t), 1, &
                        theirs, status)
                    jump = max(jump, maxval(abs(mine(2:3) - theirs(2:3))))
                end do
            end do
        end do
    end function gradient_jump

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

end program test_spline_c1_poisson
