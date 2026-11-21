module spline_xpoint_thb_geometry
    !! u - lap u = f with exact Neumann data, u = cos(0.7 R) sin(0.9 Z) + R Z,
    !! and a double-null geometry: Gamma_p made of four quadratic Bezier arcs
    !! with 90 degree corners at the X-points (R0, -h), (R0, h), straight legs
    !! with a slight bend, and an elliptic wall.
    use fortfem_kinds, only: dp
    use fortfem_spline_energy, only: pointwise_energy_t
    use fortfem_spline_xpoint, only: xpoint_curves_t
    implicit none
    private

    public :: problem_kernel_t, exact_solution, exact_gradient, source, &
        double_null, core_area

    real(dp), parameter, public :: r0 = 1.7_dp
    real(dp), parameter, public :: wall_a = 1.05_dp, wall_b = 1.6_dp

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

        f = exact_solution(s) + 1.3_dp*cos(0.7_dp*s(1))*sin(0.9_dp*s(2))
    end function source

    pure function bezier_controls(a, h, shift) result(c)
        !! (2, 3, 4): control points of the four Gamma_p arcs.
        real(dp), intent(in) :: a, h, shift
        real(dp) :: c(2, 3, 4)
        real(dp) :: k(2, 4)

        k(:, 1) = [r0 + shift, -h]
        k(:, 2) = [r0 + a, 0.0_dp]
        k(:, 3) = [r0 - shift, h]
        k(:, 4) = [r0 - a, 0.0_dp]
        c(:, :, 1) = reshape([k(:, 1), k(:, 1) + [a - shift, a - shift], k(:, 2)], [2, 3])
        c(:, :, 2) = reshape([k(:, 2), k(:, 3) + [a + shift, -(a + shift)], k(:, 3)], &
            [2, 3])
        c(:, :, 3) = reshape([k(:, 3), k(:, 3) + [-(a - shift), -(a - shift)], &
            k(:, 4)], [2, 3])
        c(:, :, 4) = reshape([k(:, 4), k(:, 1) + [-(a + shift), a + shift], k(:, 1)], &
            [2, 3])
    end function bezier_controls

    function double_null(a, h, shift, legs) result(curves)
        !! legs(q): corner q is an X-point (only corners 1 and 3 may be).
        real(dp), intent(in) :: a, h, shift
        logical, intent(in) :: legs(4)
        type(xpoint_curves_t) :: curves

        real(dp) :: c(2, 3, 4), t, theta
        integer :: q, i
        integer, parameter :: n = 40

        c = bezier_controls(a, h, shift)
        do q = 1, 4
            allocate(curves%segment(q)%points(2, n))
            do i = 1, n
                t = real(i - 1, dp)/(n - 1)
                curves%segment(q)%points(:, i) = (1.0_dp - t)**2*c(:, 1, q) + &
                    2.0_dp*t*(1.0_dp - t)*c(:, 2, q) + t**2*c(:, 3, q)
            end do
        end do
        if (legs(1)) then
            curves%leg_after(1)%points = leg(c(:, 1, 1), [1.0_dp, -1.0_dp])
            curves%leg_before(1)%points = leg(c(:, 1, 1), [-1.0_dp, -1.0_dp])
        end if
        if (legs(3)) then
            curves%leg_after(3)%points = leg(c(:, 1, 3), [-1.0_dp, 1.0_dp])
            curves%leg_before(3)%points = leg(c(:, 1, 3), [1.0_dp, 1.0_dp])
        end if
        allocate(curves%wall(2, 300))
        do i = 1, 300
            theta = 2.0_dp*acos(-1.0_dp)*(i - 1)/300
            curves%wall(:, i) = [r0 + wall_a*cos(theta), wall_b*sin(theta)]
        end do
    end function double_null

    function leg(x, d) result(points)
        !! From x along d to the wall, bent by 3 % of its length.
        real(dp), intent(in) :: x(2), d(2)
        real(dp), allocatable :: points(:, :)

        real(dp) :: e(2), lo, hi, mid, s(2), nrm(2), t
        integer :: i, it

        e = d/norm2(d)
        lo = 0.0_dp
        hi = 5.0_dp
        do it = 1, 80
            mid = 0.5_dp*(lo + hi)
            s = x + mid*e
            if (((s(1) - r0)/wall_a)**2 + (s(2)/wall_b)**2 < 1.0_dp) then
                lo = mid
            else
                hi = mid
            end if
        end do
        nrm = [-e(2), e(1)]
        allocate(points(2, 25))
        do i = 1, 25
            t = real(i - 1, dp)/24
            points(:, i) = x + t*lo*e + 0.12_dp*lo*t*(1.0_dp - t)*nrm
        end do
    end function leg

    real(dp) function core_area(a, h, shift) result(area)
        !! int R dZ over Gamma_p (fine trapezoidal rule on the arcs).
        real(dp), intent(in) :: a, h, shift
        real(dp) :: c(2, 3, 4), t, x(2), dx(2)
        integer :: q, i
        integer, parameter :: n = 20000

        c = bezier_controls(a, h, shift)
        area = 0.0_dp
        do q = 1, 4
            do i = 0, n - 1
                t = (i + 0.5_dp)/n
                x = (1.0_dp - t)**2*c(:, 1, q) + 2.0_dp*t*(1.0_dp - t)*c(:, 2, q) + &
                    t**2*c(:, 3, q)
                dx = 2.0_dp*(1.0_dp - t)*(c(:, 2, q) - c(:, 1, q)) + &
                    2.0_dp*t*(c(:, 3, q) - c(:, 2, q))
                area = area + x(1)*dx(2)/n
            end do
        end do
    end function core_area

end module spline_xpoint_thb_geometry

module spline_xpoint_thb_problem
    !! Narrow bump at the lower X-point: u = exp(-|x - X|^2/w^2), u - lap u = f
    !! with exact Neumann data.
    use fortfem_kinds, only: dp
    use fortfem_spline_energy, only: pointwise_energy_t
    implicit none
    private

    public :: bump_kernel_t, bump, bump_source

    real(dp), parameter, public :: xpt(2) = [1.7_dp, -1.0_dp]
    real(dp), parameter :: width = 0.08_dp

    type, extends(pointwise_energy_t) :: bump_kernel_t
        real(dp), allocatable :: f(:)
    contains
        procedure :: evaluate
    end type bump_kernel_t

contains

    subroutine evaluate(self, cell, first_point, npoint, jet, weight, need, w, &
            gradient, hessian, status)
        class(bump_kernel_t), intent(in) :: self
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

    pure subroutine bump(s, u, grad)
        real(dp), intent(in) :: s(2)
        real(dp), intent(out) :: u, grad(2)

        u = exp(-sum((s - xpt)**2)/width**2)
        grad = -2.0_dp*(s - xpt)/width**2*u
    end subroutine bump

    pure real(dp) function bump_source(s) result(f)
        real(dp), intent(in) :: s(2)
        real(dp) :: u, grad(2), r2

        call bump(s, u, grad)
        r2 = sum((s - xpt)**2)
        f = u - (4.0_dp*r2/width**4 - 4.0_dp/width**2)*u
    end function bump_source

end module spline_xpoint_thb_problem

program test_spline_xpoint_thb
    !! THB refinement on the double-null layout: radial levels refine the
    !! core ring and scrape-off blocks across Gamma_p and the private flux
    !! blocks along the legs. Partition of unity, exact projection, at most
    !! two levels per function, and adaptive refinement of a bump at the
    !! lower X-point beating uniform knots per DOF.
    use check, only: check_condition, check_summary
    use fortfem_kinds, only: dp
    use fortfem_element_to_csr, only: csr_matrix_t
    use fortfem_spline_layout, only: block_layout_t, LAYOUT_WALL
    use fortfem_spline_xpoint, only: layout_xpoint
    use fortfem_spline_space, only: spline_space_t, spline_space_options_t, &
        spline_space_build
    use fortfem_spline_eval, only: spline_quadrature_t, spline_quadrature_build, &
        spline_evaluate, spline_evaluate_at, spline_side_quadrature_t, &
        spline_side_quadrature, spline_point_matrix, spline_geometry_at
    use fortfem_spline_energy, only: spline_energy, spline_energy_hessian_csr
    use fortfem_spline_thb, only: spline_thb_t, spline_thb_init, spline_thb_refine, &
        spline_thb_mark, spline_project, THB_ISOTROPIC, THB_RADIAL
    use spline_xpoint_thb_geometry, only: double_null, r0
    use spline_xpoint_thb_problem, only: bump_kernel_t, bump, bump_source, xpt
    implicit none

    type(block_layout_t) :: layout
    type(spline_thb_t) :: thb
    type(spline_space_t) :: coarse, uniform
    type(spline_quadrature_t) :: quad
    real(dp), allocatable :: c0(:, :), c1(:, :), u(:, :), eta(:), jet(:, :)
    real(dp) :: uniform_error(3), error, worst, v(1)
    integer, allocatable :: marked(:)
    integer :: status, step, level, uniform_dofs(3), direction, q

    call layout_xpoint([r0, 0.0_dp], 0.4_dp, double_null(0.6_dp, 1.0_dp, 0.0_dp, &
        [.true., .false., .true., .false.]), layout, status)
    call spline_thb_init(layout, options(4), thb, status)
    call check_condition(status == 0, "THB hierarchy on the double null")
    coarse = thb%space
    allocate(c0(1, coarse%dof_count))
    call random_number(c0)
    call spline_thb_refine(thb, near(thb%space, 0.5_dp), status, THB_RADIAL)
    call spline_thb_refine(thb, near(thb%space, 0.3_dp), status, THB_RADIAL)
    call spline_thb_refine(thb, near(thb%space, 0.2_dp), status, THB_ISOTROPIC)
    call check_condition(status == 0 .and. size(thb%level) == 4, &
        "radial, radial and isotropic levels around the X-point")
    call spline_quadrature_build(thb%space, 4, quad, status)
    allocate(u(1, thb%space%dof_count), c1(1, thb%space%dof_count), &
        jet(1, quad%point_count))
    u = 1.0_dp
    call spline_evaluate(thb%space, quad, u, 0, jet, status)
    call check_condition(maxval(abs(jet(1, :) - 1.0_dp)) < 1.0e-12_dp, &
        "partition of unity")
    call spline_project(coarse, c0, thb%space, quad, c1, status)
    call spline_evaluate(thb%space, quad, c1, 0, jet, status)
    worst = 0.0_dp
    do q = 1, quad%point_count
        call spline_evaluate_at(coarse, c0, quad%block(q), quad%xi(:, q), 0, v, status)
        worst = max(worst, abs(v(1) - jet(1, q)))
    end do
    call check_condition(worst < 1.0e-10_dp .and. max_levels(thb%space) <= 2, &
        "exact projection and at most two levels per function")

    do level = 1, 3
        call spline_space_build(layout, options(2**(level + 1)), uniform, status)
        call solve(uniform, error, eta)
        uniform_error(level) = error
        uniform_dofs(level) = uniform%dof_count
    end do
    call spline_thb_init(layout, options(4), thb, status)
    do step = 1, 9
        call solve(thb%space, error, eta)
        print '(a,i2,a,i6,a,es10.3)', "THB step", step, " dofs", &
            thb%space%dof_count, " H1 error", error
        if (step == 9) exit
        call spline_thb_mark(eta, 0.5_dp, marked)
        direction = merge(THB_RADIAL, THB_ISOTROPIC, size(thb%level) <= 2)
        call spline_thb_refine(thb, marked, status, direction)
        if (status /= 0) error stop "refine"
    end do
    print '(a,3i7,a,3es10.3)', "uniform dofs", uniform_dofs, " H1 errors", &
        uniform_error
    call check_condition(thb%space%dof_count < uniform_dofs(3) .and. &
        error < uniform_error(3), &
        "adaptive THB on the double null beats uniform knots per DOF")
    call check_summary("spline_xpoint_thb")

contains

    function options(cells) result(opt)
        integer, intent(in) :: cells
        type(spline_space_options_t) :: opt

        opt%degree = 3
        allocate(opt%knots(2, layout%block_count))
        opt%knots(:, :)%element_count = cells
    end function options

    function near(space, radius) result(cells)
        type(spline_space_t), intent(in) :: space
        real(dp), intent(in) :: radius
        integer, allocatable :: cells(:)

        real(dp) :: g(2, 6), xi(2)
        integer :: cell, patch, e(2), status

        allocate(cells(0))
        do cell = 1, space%cell_count
            patch = space%cell_block(cell)
            e = space%cell_element(:, cell)
            xi(1) = 0.5_dp*sum(space%blocks(patch)%breaks1(e(1):e(1) + 1))
            xi(2) = 0.5_dp*sum(space%blocks(patch)%breaks2(e(2):e(2) + 1))
            call spline_geometry_at(space, space%patch_block(patch), xi, g, status)
            if (norm2(g(:, 1) - xpt) < radius) cells = [cells, cell]
        end do
    end function near

    integer function max_levels(space) result(span)
        type(spline_space_t), intent(in) :: space

        integer, allocatable :: lo(:), hi(:)
        integer :: j, q, patch, level

        allocate(lo(space%dof_count), hi(space%dof_count))
        lo = huge(1)
        hi = -1
        do patch = 1, space%block_count
            level = space%patch_level(patch)
            do j = space%blocks(patch)%offset + 1, space%blocks(patch)%offset + &
                    product(space%blocks(patch)%n)
                do q = space%broken_ptr(j), space%broken_ptr(j + 1) - 1
                    lo(space%broken_dof(q)) = min(lo(space%broken_dof(q)), level)
                    hi(space%broken_dof(q)) = max(hi(space%broken_dof(q)), level)
                end do
            end do
        end do
        span = maxval(hi - lo) + 1
    end function max_levels

    subroutine solve(space, error, eta)
        type(spline_space_t), intent(in) :: space
        real(dp), intent(out) :: error
        real(dp), allocatable, intent(out) :: eta(:)

        type(spline_quadrature_t) :: quad
        type(spline_side_quadrature_t) :: side
        type(bump_kernel_t) :: kernel
        type(csr_matrix_t) :: matrix
        real(dp), allocatable :: c(:, :), g(:, :), jet(:, :), val(:), b(:)
        integer, allocatable :: rowptr(:), colidx(:)
        real(dp) :: energy, normal(2), ue, ge(2)
        integer :: status, q, k, cell

        call spline_quadrature_build(space, 5, quad, status)
        allocate(kernel%f(quad%point_count))
        do q = 1, quad%point_count
            kernel%f(q) = bump_source(quad%s(:, q))
        end do
        allocate(c(1, space%dof_count), g(1, space%dof_count))
        c = 0.0_dp
        call spline_energy(space, quad, kernel, c, energy, status, g)
        call spline_energy_hessian_csr(space, quad, kernel, c, matrix, status)
        b = -g(1, :)
        call spline_side_quadrature(space, LAYOUT_WALL, 5, side, status)
        call spline_point_matrix(space, side%block, side%xi, 0, rowptr, colidx, val, &
            status)
        do q = 1, side%point_count
            normal = [side%tangent(2, q), -side%tangent(1, q)]
            call bump(side%s(:, q), ue, ge)
            do k = rowptr(q), rowptr(q + 1) - 1
                b(colidx(k)) = b(colidx(k)) + side%weight(q)*val(k)* &
                    dot_product(ge, normal)
            end do
        end do
        call pcg(matrix, b, c(1, :))
        allocate(jet(3, quad%point_count), eta(space%cell_count))
        call spline_evaluate(space, quad, c, 1, jet, status)
        eta = 0.0_dp
        do cell = 1, space%cell_count
            do q = quad%point_offset(cell), quad%point_offset(cell) + &
                    quad%points_per_cell - 1
                call bump(quad%s(:, q), ue, ge)
                eta(cell) = eta(cell) + quad%weight(q)*((jet(1, q) - ue)**2 + &
                    sum((jet(2:3, q) - ge)**2))
            end do
        end do
        error = sqrt(sum(eta))
        eta = sqrt(eta)
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

end program test_spline_xpoint_thb
