module spline_gs_bem_field
    !! Exact free-boundary field: w = S(t) G(x, x0) with the Grad-Shafranov
    !! Green function G, t = (|x - x0| - a)/(b - a) and the C5 smoothstep S.
    !! w is smooth, equals G outside the disc |x - x0| < b, and so decays
    !! and solves L w = -div((1/R) grad w) = 0 there. The source is f = L w.
    !! The energy W = |grad psi|^2/(2R) - f psi on Omega plus the exterior
    !! energy (1/2) <S_h u, u> of the trace u has the minimiser psi = w.
    use fortfem_kinds, only: dp
    use fortfem_axisymmetric_gs_kernel, only: gs_green_hessian
    use fortfem_spline_energy, only: pointwise_energy_t
    implicit none
    private

    public :: gs_kernel_t, field, source

    real(dp), parameter, public :: x0(2) = [3.1_dp, 0.05_dp]
    real(dp), parameter :: a = 0.15_dp, b = 0.6_dp

    type, extends(pointwise_energy_t) :: gs_kernel_t
        real(dp), allocatable :: s(:, :), f(:)
    contains
        procedure :: evaluate
    end type gs_kernel_t

contains

    subroutine evaluate(self, cell, first_point, npoint, jet, weight, need, w, &
            gradient, hessian, status)
        class(gs_kernel_t), intent(in) :: self
        integer, intent(in) :: cell, first_point, npoint, need
        real(dp), intent(in) :: jet(:, :), weight(:)
        real(dp), intent(out) :: w(:), gradient(:, :), hessian(:, :)
        integer, intent(out) :: status

        real(dp) :: r, f
        integer :: q

        associate (unused => cell, unused_weight => weight)
        end associate
        do q = 1, npoint
            r = self%s(1, first_point + q - 1)
            f = self%f(first_point + q - 1)
            w(q) = 0.5_dp*(jet(2, q)**2 + jet(3, q)**2)/r - f*jet(1, q)
            if (need >= 1) gradient(:, q) = [-f, jet(2, q)/r, jet(3, q)/r]
            if (need >= 2) hessian(:, q) = [0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
                1.0_dp/r, 0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp/r]
        end do
        status = 0
    end subroutine evaluate

    pure subroutine field(x, w, grad, hess)
        real(dp), intent(in) :: x(2)
        real(dp), intent(out) :: w, grad(2), hess(2, 2)

        real(dp) :: g, g4(4), h4(4, 4), d(2), r, t, st, dst, d2st, dt(2)
        integer :: i, j

        d = x - x0
        r = norm2(d)
        t = (r - a)/(b - a)
        if (t <= 0.0_dp) then
            w = 0.0_dp
            grad = 0.0_dp
            hess = 0.0_dp
            return
        end if
        call gs_green_hessian(x(1), x(2), x0(1), x0(2), g, g4, h4)
        if (t >= 1.0_dp) then
            w = g
            grad = g4(1:2)
            hess = h4(1:2, 1:2)
            return
        end if
        st = t**6*(462.0_dp + t*(-1980.0_dp + t*(3465.0_dp + t*(-3080.0_dp + &
            t*(1386.0_dp - 252.0_dp*t)))))
        dst = 2772.0_dp*(t*(1.0_dp - t))**5
        d2st = 13860.0_dp*(t*(1.0_dp - t))**4*(1.0_dp - 2.0_dp*t)
        dt = d/(r*(b - a))
        w = st*g
        grad = dst*dt*g + st*g4(1:2)
        do j = 1, 2
            do i = 1, 2
                hess(i, j) = d2st*dt(i)*dt(j)*g + dst*g*(merge(1.0_dp, 0.0_dp, &
                    i == j) - d(i)*d(j)/r**2)/(r*(b - a)) + dst*(dt(i)*g4(j) + &
                    dt(j)*g4(i)) + st*h4(i, j)
            end do
        end do
    end subroutine field

    pure real(dp) function source(x) result(f)
        real(dp), intent(in) :: x(2)

        real(dp) :: w, grad(2), hess(2, 2)

        call field(x, w, grad, hess)
        f = -(hess(1, 1) + hess(2, 2))/x(1) + grad(1)/x(1)**2
    end function source

end module spline_gs_bem_field

program test_spline_gs_bem_coupling
    !! Spline-BEM coupling of the axisymmetric Grad-Shafranov operator on the
    !! O-grid disc. The trace matrix T interpolates the spline field at the
    !! nodes of the quadratic boundary panels, the coupled Hessian is
    !! A + T^T S_h T with the symmetric exterior Steklov-Poincare matrix S_h,
    !! and the discrete free-boundary solution converges to an exact field
    !! whose exterior part is the decaying Green function of a filament.
    use check, only: check_condition, check_summary
    use fortfem_kinds, only: dp
    use fortfem_element_to_csr, only: csr_matrix_t
    use fortfem_spline_layout, only: block_layout_t, layout_ogrid, LAYOUT_WALL
    use fortfem_spline_space, only: spline_space_t, spline_space_options_t, &
        spline_space_build
    use fortfem_spline_eval, only: spline_quadrature_t, spline_quadrature_build, &
        spline_evaluate, spline_point_matrix, spline_boundary_nodes, &
        spline_evaluate_at
    use fortfem_spline_energy, only: spline_energy, spline_energy_hessian_csr
    use fortfem_axisymmetric_gs_boundary_operators_2d, only: gs_bem_curve_t, &
        gs_bem_operators_t, gs_bem_curve_from_nodes, gs_bem_assemble, &
        gs_bem_steklov_symmetric
    use spline_gs_bem_field, only: gs_kernel_t, field, source
    implicit none

    interface
        subroutine dgesv(n, nrhs, a, lda, ipiv, b, ldb, info)
            import :: dp
            integer, intent(in) :: n, nrhs, lda, ldb
            real(dp), intent(inout) :: a(lda, *), b(ldb, *)
            integer, intent(out) :: ipiv(*), info
        end subroutine dgesv
    end interface

    type(block_layout_t) :: layout
    real(dp) :: errors(3), trace_errors(3), rate
    integer :: status, level

    call layout_ogrid([3.0_dp, 0.0_dp], 1.0_dp, 0.5_dp, layout, status)
    call check_condition(status == 0, "O-grid layout")
    call check_condition(trace_matrix_ok(), &
        "point matrix reproduces spline_evaluate_at values and gradients")
    do level = 1, 3
        call solve(2**(level + 1), errors(level), trace_errors(level))
    end do
    rate = log(errors(2)/errors(3))/log(2.0_dp)
    print "(a,3es11.3,a,f6.2,a,3es11.3)", "coupled GS L2 errors", errors, &
        " rate", rate, " max trace errors", trace_errors
    call check_condition(rate > 3.5_dp .and. errors(3) < 1.0e-4_dp, &
        "coupled spline-BEM GS solution converges to the free-boundary field")
    call check_summary("spline_gs_bem_coupling")

contains

    subroutine build(cells, space, quad)
        integer, intent(in) :: cells
        type(spline_space_t), intent(out) :: space
        type(spline_quadrature_t), intent(out) :: quad

        type(spline_space_options_t) :: options
        integer :: status

        options%degree = 3
        allocate(options%knots(2, 5))
        options%knots(:, :)%element_count = cells
        call spline_space_build(layout, options, space, status)
        if (status /= 0) error stop "space"
        call spline_quadrature_build(space, 5, quad, status)
        if (status /= 0) error stop "quadrature"
    end subroutine build

    logical function trace_matrix_ok() result(ok)
        type(spline_space_t) :: space
        type(spline_quadrature_t) :: quad
        real(dp), allocatable :: coef(:, :), xi(:, :), s(:, :), val(:)
        integer, allocatable :: block(:), rowptr(:), colidx(:)
        real(dp) :: jet(3), worst
        integer :: status, component, q

        call build(4, space, quad)
        allocate(coef(1, space%dof_count))
        call random_number(coef)
        call spline_boundary_nodes(space, LAYOUT_WALL, 2, block, xi, s, status)
        ok = status == 0 .and. size(block) == 2*2*4*4
        worst = 0.0_dp
        do component = 0, 2
            call spline_point_matrix(space, block, xi, component, rowptr, colidx, &
                val, status)
            ok = ok .and. status == 0
            do q = 1, size(block)
                call spline_evaluate_at(space, coef, block(q), xi(:, q), 1, jet, status)
                worst = max(worst, abs(jet(component + 1) - dot_product( &
                    val(rowptr(q):rowptr(q + 1) - 1), &
                    coef(1, colidx(rowptr(q):rowptr(q + 1) - 1)))))
            end do
        end do
        ok = ok .and. worst < 1.0e-12_dp*maxval(abs(coef))*10.0_dp
    end function trace_matrix_ok

    subroutine solve(cells, error, trace_error)
        integer, intent(in) :: cells
        real(dp), intent(out) :: error, trace_error

        type(spline_space_t) :: space
        type(spline_quadrature_t) :: quad
        type(gs_kernel_t) :: kernel
        type(csr_matrix_t) :: matrix
        type(gs_bem_curve_t) :: curve
        type(gs_bem_operators_t) :: ops
        real(dp), allocatable :: c(:, :), g(:, :), jet(:, :), xi(:, :), s(:, :)
        real(dp), allocatable :: val(:), steklov(:, :), t(:, :), k(:, :), rhs(:)
        integer, allocatable :: block(:), rowptr(:), colidx(:), pivot(:)
        real(dp) :: energy, w, grad(2), hess(2, 2)
        integer :: status, q, n, i, j, info

        call build(cells, space, quad)
        n = space%dof_count
        kernel%s = quad%s
        allocate(kernel%f(quad%point_count))
        do q = 1, quad%point_count
            kernel%f(q) = source(quad%s(:, q))
        end do
        allocate(c(1, n), g(1, n))
        c = 0.0_dp
        call spline_energy(space, quad, kernel, c, energy, status, g)
        call spline_energy_hessian_csr(space, quad, kernel, c, matrix, status)
        if (status /= 0) error stop "hessian"
        call spline_boundary_nodes(space, LAYOUT_WALL, 2, block, xi, s, status)
        if (status /= 0) error stop "nodes"
        call spline_point_matrix(space, block, xi, 0, rowptr, colidx, val, status)
        if (status /= 0) error stop "trace matrix"
        call gs_bem_curve_from_nodes(s, curve, status)
        if (status /= 0) error stop "curve"
        call gs_bem_assemble(curve, 1, ops, status)
        if (status /= 0) error stop "bem"
        call gs_bem_steklov_symmetric(ops, steklov, status)
        if (status /= 0) error stop "steklov"
        allocate(t(size(block), n), k(n, n), rhs(n), pivot(n))
        t = 0.0_dp
        do q = 1, size(block)
            t(q, colidx(rowptr(q):rowptr(q + 1) - 1)) = val(rowptr(q):rowptr(q + 1) - 1)
        end do
        k = matmul(transpose(t), matmul(steklov, t))
        do i = 1, n
            do j = matrix%rowptr(i), matrix%rowptr(i + 1) - 1
                k(i, matrix%colidx(j)) = k(i, matrix%colidx(j)) + matrix%val(j)
            end do
        end do
        rhs = -g(1, :)
        call dgesv(n, 1, k, n, pivot, rhs, n, info)
        if (info /= 0) error stop "solve"
        c(1, :) = rhs
        allocate(jet(1, quad%point_count))
        call spline_evaluate(space, quad, c, 0, jet, status)
        error = 0.0_dp
        do q = 1, quad%point_count
            call field(quad%s(:, q), w, grad, hess)
            error = error + quad%weight(q)*(jet(1, q) - w)**2
        end do
        error = sqrt(error)
        trace_error = 0.0_dp
        do q = 1, size(block)
            call field(s(:, q), w, grad, hess)
            trace_error = max(trace_error, abs(dot_product(t(q, :), c(1, :)) - w))
        end do
    end subroutine solve

end program test_spline_gs_bem_coupling
