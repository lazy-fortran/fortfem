module spline_toroidal_problem
    !! u - lap u = f in a tube about a figure-8 guide curve, in the chart
    !! coordinates (s1, s2, t): W = (|J^-T grad u|^2/2 + u^2/2 - f u) |det J|.
    use fortfem_kinds, only: dp
    use fortfem_spline_toroidal, only: pointwise_energy3_t
    implicit none
    private

    public :: tube_kernel_t, exact, source, guide_curve

    type, extends(pointwise_energy3_t) :: tube_kernel_t
        real(dp), allocatable :: jdet(:, :), ginv(:, :, :, :), f(:, :) ! (point, plane)
    contains
        procedure :: evaluate
    end type tube_kernel_t

contains

    subroutine evaluate(self, cell, plane, first_point, npoint, jet, weight, need, w, &
            gradient, hessian, status)
        class(tube_kernel_t), intent(in) :: self
        integer, intent(in) :: cell, plane, first_point, npoint, need
        real(dp), intent(in) :: jet(:, :), weight(:)
        real(dp), intent(out) :: w(:), gradient(:, :), hessian(:, :)
        integer, intent(out) :: status

        real(dp) :: a, gi(3, 3), u, gr(3), h(4, 4)
        integer :: q, p

        associate (unused => cell, unused_weight => weight)
        end associate
        do q = 1, npoint
            p = first_point + q - 1
            a = self%jdet(p, plane)
            gi = self%ginv(:, :, p, plane)
            u = jet(1, q)
            gr = jet(2:4, q)
            w(q) = (0.5_dp*dot_product(gr, matmul(gi, gr)) + 0.5_dp*u**2 - &
                self%f(p, plane)*u)*a
            if (need >= 1) gradient(:, q) = [(u - self%f(p, plane))*a, &
                matmul(gi, gr)*a]
            if (need >= 2) then
                h = 0.0_dp
                h(1, 1) = a
                h(2:4, 2:4) = gi*a
                hessian(:, q) = reshape(h, [16])
            end if
        end do
        status = 0
    end subroutine evaluate

    pure subroutine exact(x, problem, u, grad)
        !! 1: x1 + x2^2 + x1 x3 (full torus), 2: x1^2 - x2^2 + x3 (period pi),
        !! 3: sin(x1/2) cos(x2/2) + 0.3 x1 x3 (figure-8 tube).
        real(dp), intent(in) :: x(3)
        integer, intent(in) :: problem
        real(dp), intent(out) :: u, grad(3)

        select case (problem)
        case (1)
            u = x(1) + x(2)**2 + x(1)*x(3)
            grad = [1.0_dp + x(3), 2.0_dp*x(2), x(1)]
        case (2)
            u = x(1)**2 - x(2)**2 + x(3)
            grad = [2.0_dp*x(1), -2.0_dp*x(2), 1.0_dp]
        case default
            u = sin(0.5_dp*x(1))*cos(0.5_dp*x(2)) + 0.3_dp*x(1)*x(3)
            grad = [0.5_dp*cos(0.5_dp*x(1))*cos(0.5_dp*x(2)) + 0.3_dp*x(3), &
                -0.5_dp*sin(0.5_dp*x(1))*sin(0.5_dp*x(2)), 0.3_dp*x(1)]
        end select
    end subroutine exact

    pure real(dp) function source(x, problem) result(f)
        !! f = u - lap u.
        real(dp), intent(in) :: x(3)
        integer, intent(in) :: problem
        real(dp) :: u, grad(3)

        call exact(x, problem, u, grad)
        select case (problem)
        case (1)
            f = u - 2.0_dp
        case (2)
            f = u
        case default
            f = u + 0.5_dp*sin(0.5_dp*x(1))*cos(0.5_dp*x(2))
        end select
    end function source

    pure function guide_curve(n) result(points)
        !! Figure-8: lemniscate of Gerono in (x, y), lifted by z = 1.2 sin.
        integer, intent(in) :: n
        real(dp) :: points(3, n), theta
        integer :: k

        do k = 1, n
            theta = 2.0_dp*acos(-1.0_dp)*(k - 1)/n
            points(:, k) = [2.0_dp*cos(theta), sin(2.0_dp*theta), 1.2_dp*sin(theta)]
        end do
    end function guide_curve

end module spline_toroidal_problem

program test_spline_toroidal
    !! Gate 4 of doc/design/multipatch_spline_spaces.md: O-grid sections
    !! extended with Fourier modes along a cylindrical torus (one field
    !! period of nfp = 1 and 2) and along a figure-8 guide curve.
    !! u - lap u = f with exact Neumann data: order p + 1 in the section when
    !! the modes resolve u, geometric decay in the mode number on the
    !! figure-8 tube; the assembled Hessian equals the matrix-free one and
    !! the gradient matches finite differences.
    use check, only: check_condition, check_summary
    use fortfem_kinds, only: dp
    use fortfem_element_to_csr, only: csr_matrix_t
    use fortfem_spline_layout, only: block_layout_t, layout_ogrid, LAYOUT_WALL
    use fortfem_spline_space, only: spline_space_t, spline_space_options_t, &
        spline_space_build
    use fortfem_spline_eval, only: spline_quadrature_t, spline_quadrature_build, &
        spline_side_quadrature_t, spline_side_quadrature, spline_point_matrix
    use fortfem_spline_toroidal, only: spline_toroidal_t, spline_toroidal_fourier, &
        spline_energy3, spline_energy3_hessian_csr, spline_energy3_hessian_apply, &
        spline_evaluate3, spline_toroidal_periodic
    use fortfem_curve_chart, only: curve_chart_t, curve_chart_from_samples, &
        curve_chart_map, curve_chart_cylindrical
    use spline_toroidal_problem, only: tube_kernel_t, exact, source, guide_curve
    implicit none

    type(block_layout_t) :: layout
    type(curve_chart_t) :: chart
    type(spline_toroidal_t) :: tor
    real(dp) :: errors(3), rate
    integer :: status, level, problem, nfp

    do nfp = 1, 2
        problem = nfp
        call layout_ogrid([2.0_dp, 0.0_dp], 0.5_dp, 0.5_dp, layout, status)
        call curve_chart_cylindrical(chart, nfp)
        call spline_toroidal_fourier(2, 12, tor, status)
        do level = 1, 3
            call solve(2**level, errors(level))
        end do
        rate = log(errors(2)/errors(3))/log(2.0_dp)
        print '(a,i0,a,3es11.3,a,f6.2)', "torus nfp = ", nfp, " L2 errors", errors, &
            " rate", rate
        call check_condition(rate > 3.5_dp, &
            "torus field period: order p + 1 in the section with 5 modes")
    end do

    problem = 1
    call layout_ogrid([2.0_dp, 0.0_dp], 0.5_dp, 0.5_dp, layout, status)
    call curve_chart_cylindrical(chart, 1)
    do level = 1, 3
        call spline_toroidal_periodic(3, 2**(level + 2), 4, tor, status)
        call solve(4, errors(level))
    end do
    rate = log(errors(2)/errors(3))/log(2.0_dp)
    print '(a,3es11.3,a,f6.2)', "torus, periodic splines in t, L2 errors", errors, &
        " rate", rate
    call check_condition(rate > 3.5_dp, &
        "periodic splines along t converge at order q + 1")

    problem = 3
    call layout_ogrid([0.0_dp, 0.0_dp], 0.3_dp, 0.5_dp, layout, status)
    call curve_chart_from_samples(guide_curve(129), chart, status)
    call check_condition(status == 0, "rotation-minimising chart of a figure-8 curve")
    call spline_toroidal_fourier(6, 40, tor, status)
    call check_condition(consistent(), &
        "assembled and matrix-free Hessians agree, gradient matches differences")
    do level = 1, 3
        call spline_toroidal_fourier(2 + 3*level, 8*level + 16, tor, status)
        call solve(2, errors(level))
    end do
    print '(a,3es11.3)', "figure-8 tube, 11/17/23 modes, L2 errors", errors
    call check_condition(errors(2) < errors(1)/2.5_dp .and. &
        errors(3) < errors(2)/2.5_dp, &
        "figure-8 tube: the error decays geometrically in the mode number")
    call check_summary("spline_toroidal")

contains

    subroutine setup(cells, space, quad, kernel)
        integer, intent(in) :: cells
        type(spline_space_t), intent(out) :: space
        type(spline_quadrature_t), intent(out) :: quad
        type(tube_kernel_t), intent(out) :: kernel

        type(spline_space_options_t) :: options
        real(dp) :: x(3), jac(3, 3), g(3, 3), det
        integer :: status, q, j

        options%degree = 3
        allocate(options%knots(2, 5))
        options%knots(:, :)%element_count = cells
        call spline_space_build(layout, options, space, status)
        if (status /= 0) error stop "space"
        call spline_quadrature_build(space, 4, quad, status)
        allocate(kernel%jdet(quad%point_count, tor%plane_count), &
            kernel%ginv(3, 3, quad%point_count, tor%plane_count), &
            kernel%f(quad%point_count, tor%plane_count))
        do j = 1, tor%plane_count
            do q = 1, quad%point_count
                call curve_chart_map(chart, quad%s(:, q), tor%t(j), x, jac)
                det = jac(1, 1)*(jac(2, 2)*jac(3, 3) - jac(3, 2)*jac(2, 3)) - &
                    jac(1, 2)*(jac(2, 1)*jac(3, 3) - jac(3, 1)*jac(2, 3)) + &
                    jac(1, 3)*(jac(2, 1)*jac(3, 2) - jac(3, 1)*jac(2, 2))
                g = matmul(transpose(jac), jac)
                kernel%jdet(q, j) = abs(det)
                kernel%ginv(:, :, q, j) = inverse3(g)
                kernel%f(q, j) = source(x, problem)
            end do
        end do
    end subroutine setup

    logical function consistent() result(ok)
        type(spline_space_t) :: space
        type(spline_quadrature_t) :: quad
        type(tube_kernel_t) :: kernel
        type(csr_matrix_t) :: matrix
        real(dp), allocatable :: c(:, :, :), v(:, :, :), hv(:, :, :), g(:, :, :)
        real(dp), allocatable :: cv(:), flat(:)
        real(dp) :: e1, e2, step, directional
        integer :: status, i, k

        call setup(2, space, quad, kernel)
        allocate(c(1, tor%mode_count, space%dof_count), &
            v(1, tor%mode_count, space%dof_count), hv(1, tor%mode_count, space%dof_count), &
            g(1, tor%mode_count, space%dof_count))
        call random_number(c)
        call random_number(v)
        call spline_energy3_hessian_csr(space, quad, tor, kernel, c, matrix, status)
        ok = status == 0
        call spline_energy3_hessian_apply(space, quad, tor, kernel, c, v, hv, status)
        ok = ok .and. status == 0
        flat = reshape(v, [size(v)])
        allocate(cv(matrix%n))
        do i = 1, matrix%n
            cv(i) = 0.0_dp
            do k = matrix%rowptr(i), matrix%rowptr(i + 1) - 1
                cv(i) = cv(i) + matrix%val(k)*flat(matrix%colidx(k))
            end do
        end do
        print *, "csr vs apply", maxval(abs(cv - reshape(hv, [size(hv)]))), maxval(abs(cv))
        ok = ok .and. maxval(abs(cv - reshape(hv, [size(hv)]))) < &
            1.0e-11_dp*maxval(abs(cv))
        call spline_energy3(space, quad, tor, kernel, c, e1, status, g)
        step = 1.0e-6_dp
        call spline_energy3(space, quad, tor, kernel, c + step*v, e2, status)
        directional = sum(g*v)
        call spline_energy3(space, quad, tor, kernel, c - step*v, e1, status)
        print *, "fd", (e2 - e1)/(2.0_dp*step), directional
        ok = ok .and. abs((e2 - e1)/(2.0_dp*step) - directional) < &
            1.0e-6_dp*abs(directional)
    end function consistent

    subroutine solve(cells, error)
        integer, intent(in) :: cells
        real(dp), intent(out) :: error

        type(spline_space_t) :: space
        type(spline_quadrature_t) :: quad
        type(spline_side_quadrature_t) :: side
        type(tube_kernel_t) :: kernel
        type(csr_matrix_t) :: matrix
        real(dp), allocatable :: c(:, :, :), g(:, :, :), b(:, :, :), jet(:, :, :)
        real(dp), allocatable :: val(:), x(:)
        integer, allocatable :: rowptr(:), colidx(:)
        real(dp) :: energy, xx(3), jac(3, 3), xl(3), xn(3), nrm(3), u, gr(3)
        integer :: status, q, j, k

        call setup(cells, space, quad, kernel)
        allocate(c(1, tor%mode_count, space%dof_count), &
            g(1, tor%mode_count, space%dof_count), b(1, tor%mode_count, space%dof_count))
        c = 0.0_dp
        call spline_energy3(space, quad, tor, kernel, c, energy, status, g)
        call spline_energy3_hessian_csr(space, quad, tor, kernel, c, matrix, status)
        b = -g
        call spline_side_quadrature(space, LAYOUT_WALL, 5, side, status)
        call spline_point_matrix(space, side%block, side%xi, 0, rowptr, colidx, val, &
            status)
        do q = 1, side%point_count
            do j = 1, tor%plane_count
                call curve_chart_map(chart, side%s(:, q), tor%t(j), xx, jac)
                xl = jac(:, 1)*side%tangent(1, q) + jac(:, 2)*side%tangent(2, q)
                xn = jac(:, 1)*side%tangent(2, q) - jac(:, 2)*side%tangent(1, q)
                nrm = [xl(2)*jac(3, 3) - xl(3)*jac(2, 3), xl(3)*jac(1, 3) - &
                    xl(1)*jac(3, 3), xl(1)*jac(2, 3) - xl(2)*jac(1, 3)]
                if (dot_product(nrm, xn) < 0.0_dp) nrm = -nrm
                call exact(xx, problem, u, gr)
                do k = rowptr(q), rowptr(q + 1) - 1
                    b(1, :, colidx(k)) = b(1, :, colidx(k)) + side%weight(q)* &
                        tor%weight(j)*dot_product(gr, nrm)*tor%table(:, j, 0)*val(k)
                end do
            end do
        end do
        allocate(x(matrix%n))
        call pcg(matrix, reshape(b, [size(b)]), x)
        c = reshape(x, shape(c))
        allocate(jet(1, quad%point_count, tor%plane_count))
        call spline_evaluate3(space, quad, tor, c, 0, jet, status)
        error = 0.0_dp
        do j = 1, tor%plane_count
            do q = 1, quad%point_count
                call curve_chart_map(chart, quad%s(:, q), tor%t(j), xx, jac)
                call exact(xx, problem, u, gr)
                error = error + quad%weight(q)*tor%weight(j)*kernel%jdet(q, j)* &
                    (jet(1, q, j) - u)**2
            end do
        end do
        error = sqrt(error)
    end subroutine solve

    pure function inverse3(a) result(b)
        real(dp), intent(in) :: a(3, 3)
        real(dp) :: b(3, 3), det

        b(1, 1) = a(2, 2)*a(3, 3) - a(2, 3)*a(3, 2)
        b(1, 2) = a(1, 3)*a(3, 2) - a(1, 2)*a(3, 3)
        b(1, 3) = a(1, 2)*a(2, 3) - a(1, 3)*a(2, 2)
        b(2, 1) = a(2, 3)*a(3, 1) - a(2, 1)*a(3, 3)
        b(2, 2) = a(1, 1)*a(3, 3) - a(1, 3)*a(3, 1)
        b(2, 3) = a(1, 3)*a(2, 1) - a(1, 1)*a(2, 3)
        b(3, 1) = a(2, 1)*a(3, 2) - a(2, 2)*a(3, 1)
        b(3, 2) = a(1, 2)*a(3, 1) - a(1, 1)*a(3, 2)
        b(3, 3) = a(1, 1)*a(2, 2) - a(1, 2)*a(2, 1)
        det = a(1, 1)*b(1, 1) + a(1, 2)*b(2, 1) + a(1, 3)*b(3, 1)
        b = b/det
    end function inverse3

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
            ap = 0.0_dp
            do i = 1, matrix%n
                do k = matrix%rowptr(i), matrix%rowptr(i + 1) - 1
                    ap(i) = ap(i) + matrix%val(k)*p(matrix%colidx(k))
                end do
            end do
            alpha = rz/dot_product(p, ap)
            x = x + alpha*p
            r = r - alpha*ap
            if (norm2(r) < 1.0e-12_dp*norm_b) exit
            z = r/diag
            rz_new = dot_product(r, z)
            p = z + (rz_new/rz)*p
            rz = rz_new
        end do
    end subroutine pcg

end program test_spline_toroidal
