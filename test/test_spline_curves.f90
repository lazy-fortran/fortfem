module spline_curves_problem
    !! W = |grad u|^2/2 + u^2/2 - f u on a D-shaped domain with the boundary
    !! load g = du/dn, so the minimiser of the energy minus int g u ds solves
    !! u - lap u = f with du/dn = g. Exact u = cos(0.7 R) sin(0.9 Z) + R Z.
    use fortfem_kinds, only: dp
    use fortfem_spline_energy, only: pointwise_energy_t
    implicit none
    private

    public :: problem_kernel_t, exact_solution, exact_gradient, source, d_shape

    type, extends(pointwise_energy_t) :: problem_kernel_t
        real(dp), allocatable :: f(:)
    contains
        procedure :: evaluate
    end type problem_kernel_t

contains

    subroutine evaluate(self, cell, first_point, npoint, jet, weight, need, w, &
            gradient, hessian, status)
        class(problem_kernel_t), intent(in) :: self
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

        u = cos(0.7_dp*s(1))*sin(0.9_dp*s(2)) + s(1)*s(2)
    end function exact_solution

    pure function exact_gradient(s) result(g)
        real(dp), intent(in) :: s(2)
        real(dp) :: g(2)

        g = [-0.7_dp*sin(0.7_dp*s(1))*sin(0.9_dp*s(2)) + s(2), &
            0.9_dp*cos(0.7_dp*s(1))*cos(0.9_dp*s(2)) + s(1)]
    end function exact_gradient

    pure real(dp) function source(s) result(f)
        real(dp), intent(in) :: s(2)

        f = exact_solution(s) + (0.49_dp + 0.81_dp)*cos(0.7_dp*s(1))*sin(0.9_dp*s(2))
    end function source

    pure function d_shape(n, minor, elongation, triangularity) result(samples)
        !! R = 1.7 + a cos(theta + delta sin theta), Z = kappa a sin theta at
        !! n equally spaced theta from the outboard midplane.
        integer, intent(in) :: n
        real(dp), intent(in) :: minor, elongation, triangularity
        real(dp) :: samples(2, n)

        real(dp) :: theta
        integer :: q

        do q = 1, n
            theta = 2.0_dp*acos(-1.0_dp)*(q - 1)/n
            samples(:, q) = [1.7_dp + minor*cos(theta + triangularity*sin(theta)), &
                elongation*minor*sin(theta)]
        end do
    end function d_shape

end module spline_curves_problem

program test_spline_curves
    !! Closed sampled curves as wall and annulus interface: a D-shaped
    !! interface (kappa = 1.8, delta = 0.5) inside a D-shaped wall around the
    !! O-grid core. The disc area and the solution of a Neumann problem
    !! converge at order p + 1, core and annulus maps agree exactly on the
    !! interface, and a new interface curve changes only the geometry.
    use check, only: check_condition, check_summary
    use fortfem_kinds, only: dp
    use fortfem_element_to_csr, only: csr_matrix_t
    use fortfem_spline_layout, only: block_layout_t, layout_ogrid, &
        layout_set_curve, LAYOUT_WALL, LAYOUT_ANNULUS_INTERFACE
    use fortfem_spline_space, only: spline_space_t, spline_space_options_t, &
        spline_space_build, spline_space_update_geometry
    use fortfem_spline_eval, only: spline_quadrature_t, spline_quadrature_build, &
        spline_evaluate, spline_evaluate_at, spline_geometry_at, &
        spline_side_quadrature_t, spline_side_quadrature, spline_point_matrix
    use fortfem_spline_energy, only: spline_energy, spline_energy_hessian_csr
    use spline_curves_problem, only: problem_kernel_t, exact_solution, &
        exact_gradient, source, d_shape
    implicit none

    type(block_layout_t) :: layout
    type(spline_space_t) :: space
    real(dp) :: errors(3), areas(3), exact_area, rate, area_rate, error, area
    real(dp) :: samples(2, 200), old_corner(2), corner_jet(2, 6)
    integer :: status, level, dofs

    call make_layout(1.8_dp, 0.5_dp, layout, status)
    call check_condition(status == 0, "D-shaped wall and interface curves")
    exact_area = d_area(0.8_dp, 1.8_dp, 0.5_dp)
    do level = 1, 3
        call build(2**(level + 1), space, status)
        if (status /= 0) error stop "space"
        call solve(space, errors(level), areas(level))
    end do
    ! mean rates over n = 4 .. 16: at n = 8 the D tips are not yet resolved
    rate = log(errors(1)/errors(3))/log(4.0_dp)
    area_rate = log(abs(areas(1) - exact_area)/abs(areas(3) - exact_area))/ &
        log(4.0_dp)
    print '(a,3es11.3,a,f6.2,a,3es11.3,a,f6.2)', "p = 3 L2 errors", errors, &
        " rate", rate, " area errors", abs(areas - exact_area), " rate", area_rate
    call check_condition(rate > 3.5_dp, &
        "Neumann solution on the D-shaped domain converges at order p + 1")
    call check_condition(area_rate > 3.5_dp, &
        "area of the D-shaped domain converges at order p + 1")
    call check_condition(interface_jump(space) < 1.0e-12_dp, &
        "core and annulus agree in values and maps on the curved interface")

    dofs = space%dof_count
    old_corner = layout%vertices(:, 5)
    call layout_set_curve(layout, LAYOUT_ANNULUS_INTERFACE, &
        d_shape(200, 0.5_dp, 1.6_dp, 0.3_dp), status)
    call spline_space_update_geometry(space, layout, status)
    call spline_geometry_at(space, 2, [0.0_dp, 0.0_dp], corner_jet, status)
    call check_condition(status == 0 .and. space%dof_count == dofs .and. &
        maxval(abs(corner_jet(:, 1) - layout%vertices(:, 5))) < 1.0e-12_dp .and. &
        norm2(layout%vertices(:, 5) - old_corner) > 1.0e-2_dp, &
        "a new interface curve moves the maps on the same space")
    call solve(space, error, area)
    call check_condition(interface_jump(space) < 1.0e-12_dp .and. &
        abs(area - areas(3)) < 1.0e-10_dp .and. error < 2.0_dp*errors(3), &
        "after the update the wall is unchanged and the solution still converges")
    samples = d_shape(200, 0.8_dp, 1.8_dp, 0.5_dp)
    call layout_set_curve(layout, LAYOUT_WALL, samples(:, 200:1:-1), status)
    call check_condition(status == 2, "clockwise samples are rejected")
    call check_summary("spline_curves")

contains

    subroutine make_layout(elongation, triangularity, layout, status)
        real(dp), intent(in) :: elongation, triangularity
        type(block_layout_t), intent(out) :: layout
        integer, intent(out) :: status

        call layout_ogrid([1.7_dp, 0.0_dp], 0.8_dp, 0.3_dp, layout, status, &
            interface_fraction=0.6_dp)
        if (status /= 0) return
        call layout_set_curve(layout, LAYOUT_WALL, d_shape(256, 0.8_dp, &
            elongation, triangularity), status)
        if (status /= 0) return
        call layout_set_curve(layout, LAYOUT_ANNULUS_INTERFACE, d_shape(256, &
            0.5_dp, elongation, triangularity), status)
    end subroutine make_layout

    real(dp) function d_area(minor, elongation, triangularity) result(area)
        !! Area int R dZ of the D shape (trapezoidal rule, spectrally exact).
        real(dp), intent(in) :: minor, elongation, triangularity
        real(dp) :: theta
        integer :: q
        integer, parameter :: n = 4096

        area = 0.0_dp
        do q = 0, n - 1
            theta = 2.0_dp*acos(-1.0_dp)*q/n
            area = area + (1.7_dp + minor*cos(theta + triangularity*sin(theta)))* &
                elongation*minor*cos(theta)
        end do
        area = area*2.0_dp*acos(-1.0_dp)/n
    end function d_area

    subroutine build(cells, space, status)
        integer, intent(in) :: cells
        type(spline_space_t), intent(out) :: space
        integer, intent(out) :: status

        type(spline_space_options_t) :: options

        options%degree = 3
        allocate(options%knots(2, 6))
        options%knots(:, :)%element_count = cells
        call spline_space_build(layout, options, space, status)
    end subroutine build

    real(dp) function interface_jump(space) result(jump)
        type(spline_space_t), intent(in) :: space

        real(dp), allocatable :: coef(:, :)
        real(dp) :: core(1), ring(1), gcore(2, 6), gring(2, 6), t, xi1
        integer :: k, q, status

        allocate(coef(1, space%dof_count))
        call random_number(coef)
        jump = 0.0_dp
        do k = 2, 5
            do q = 0, 10
                t = real(q, dp)/10.0_dp
                xi1 = modulo(0.25_dp*(k - 2 + t), 1.0_dp)
                call spline_evaluate_at(space, coef, k, [t, 0.0_dp], 0, core, status)
                call spline_evaluate_at(space, coef, 6, [xi1, 1.0_dp], 0, ring, status)
                call spline_geometry_at(space, k, [t, 0.0_dp], gcore, status)
                call spline_geometry_at(space, 6, [xi1, 1.0_dp], gring, status)
                jump = max(jump, abs(core(1) - ring(1)), &
                    maxval(abs(gcore(:, 1) - gring(:, 1))))
            end do
        end do
    end function interface_jump

    subroutine solve(space, error, area)
        type(spline_space_t), intent(in) :: space
        real(dp), intent(out) :: error, area

        type(spline_quadrature_t) :: quad
        type(spline_side_quadrature_t) :: side
        type(problem_kernel_t) :: kernel
        type(csr_matrix_t) :: matrix
        real(dp), allocatable :: c(:, :), g(:, :), jet(:, :), val(:), b(:)
        integer, allocatable :: rowptr(:), colidx(:)
        real(dp) :: energy, normal(2)
        integer :: status, q, k

        call spline_quadrature_build(space, 5, quad, status)
        if (status /= 0) error stop "quadrature"
        allocate(kernel%f(quad%point_count))
        do q = 1, quad%point_count
            kernel%f(q) = source(quad%s(:, q))
        end do
        allocate(c(1, space%dof_count), g(1, space%dof_count))
        c = 0.0_dp
        call spline_energy(space, quad, kernel, c, energy, status, g)
        call spline_energy_hessian_csr(space, quad, kernel, c, matrix, status)
        b = -g(1, :)
        call spline_side_quadrature(space, LAYOUT_WALL, 5, side, status)
        call spline_point_matrix(space, side%block, side%xi, 0, rowptr, colidx, val, &
            status)
        if (status /= 0) error stop "boundary"
        do q = 1, side%point_count
            normal = [side%tangent(2, q), -side%tangent(1, q)]
            do k = rowptr(q), rowptr(q + 1) - 1
                b(colidx(k)) = b(colidx(k)) + side%weight(q)*val(k)* &
                    dot_product(exact_gradient(side%s(:, q)), normal)
            end do
        end do
        call pcg(matrix, b, c(1, :))
        allocate(jet(1, quad%point_count))
        call spline_evaluate(space, quad, c, 0, jet, status)
        error = 0.0_dp
        do q = 1, quad%point_count
            error = error + quad%weight(q)*(jet(1, q) - exact_solution(quad%s(:, q)))**2
        end do
        error = sqrt(error)
        area = sum(quad%weight)
    end subroutine solve

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

end program test_spline_curves
