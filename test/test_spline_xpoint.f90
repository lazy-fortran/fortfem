module spline_xpoint_problem
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

end module spline_xpoint_problem

program test_spline_xpoint
    !! Double-null layout (doc/design/multipatch_spline_spaces.md, X-point
    !! layouts): 19 blocks, order p + 1 for a Neumann problem on the whole
    !! domain, Gamma_p as a tagged interior curve (core area by Green, closed
    !! node loop for the boundary element coupling), a moved shape on the
    !! same space, and a topology change reported by the geometry update.
    use check, only: check_condition, check_summary
    use fortfem_kinds, only: dp
    use fortfem_element_to_csr, only: csr_matrix_t
    use fortfem_spline_layout, only: block_layout_t, LAYOUT_WALL, &
        LAYOUT_PLASMA_BOUNDARY
    use fortfem_spline_xpoint, only: layout_xpoint, xpoint_curves_t
    use fortfem_spline_space, only: spline_space_t, spline_space_options_t, &
        spline_space_build, spline_space_update_geometry
    use fortfem_spline_eval, only: spline_quadrature_t, spline_quadrature_build, &
        spline_evaluate, spline_side_quadrature_t, spline_side_quadrature, &
        spline_point_matrix, spline_boundary_nodes
    use fortfem_spline_energy, only: spline_energy, spline_energy_hessian_csr
    use spline_xpoint_problem, only: problem_kernel_t, exact_solution, &
        exact_gradient, source, double_null, core_area, r0, wall_a, wall_b
    implicit none

    real(dp), parameter :: pi = acos(-1.0_dp)
    type(block_layout_t) :: layout, moved, single
    type(spline_space_t) :: space
    real(dp) :: errors(3), areas(3), green, rate, area_rate, error, area
    integer, allocatable :: block(:)
    real(dp), allocatable :: xi(:, :), s(:, :)
    integer :: status, level, dofs
    logical, parameter :: dn(4) = [.true., .false., .true., .false.]

    call layout_xpoint([r0, 0.0_dp], 0.4_dp, double_null(0.6_dp, 1.0_dp, 0.0_dp, &
        dn), layout, status)
    call check_condition(status == 0 .and. layout%block_count == 19, &
        "double-null layout with 19 blocks")
    do level = 1, 3
        call build(2**(level + 1), space, status)
        if (status /= 0) error stop "space"
        call solve(space, errors(level), areas(level), green)
    end do
    rate = log(errors(2)/errors(3))/log(2.0_dp)
    area_rate = log(abs(areas(2) - pi*wall_a*wall_b)/abs(areas(3) - &
        pi*wall_a*wall_b))/log(2.0_dp)
    print '(a,3es11.3,a,f6.2,a,3es11.3,a,f6.2)', "L2 errors", errors, " rate", &
        rate, " area errors", abs(areas - pi*wall_a*wall_b), " rate", area_rate
    call check_condition(rate > 3.5_dp .and. area_rate > 3.5_dp, &
        "solution and area converge at order p + 1 through the X-points")
    call check_condition(abs(green - core_area(0.6_dp, 1.0_dp, 0.0_dp)) < 2.0e-6_dp, &
        "the tagged Gamma_p encloses the core area (Green)")
    call spline_boundary_nodes(space, LAYOUT_PLASMA_BOUNDARY, 2, block, xi, s, status)
    call check_condition(status == 0 .and. size(block) == 2*2*4*16, &
        "Gamma_p is a closed loop of boundary element nodes")

    dofs = space%dof_count
    call layout_xpoint([r0, 0.0_dp], 0.4_dp, double_null(0.55_dp, 1.1_dp, 0.05_dp, &
        dn), moved, status)
    call spline_space_update_geometry(space, moved, status)
    call check_condition(status == 0 .and. space%dof_count == dofs, &
        "a moved double null keeps the space")
    call solve(space, error, area, green)
    call check_condition(error < 2.0_dp*errors(3) .and. &
        abs(green - core_area(0.55_dp, 1.1_dp, 0.05_dp)) < 2.0e-6_dp, &
        "after the move Gamma_p follows the new shape and the error stays small")
    call layout_xpoint([r0, 0.0_dp], 0.4_dp, double_null(0.6_dp, 1.0_dp, 0.0_dp, &
        [.true., .false., .false., .false.]), single, status)
    call check_condition(status == 0 .and. single%block_count == 16, &
        "single-null layout with 16 blocks")
    call spline_space_update_geometry(space, single, status)
    call check_condition(status == 1, &
        "a topology change is reported (status 1: re-project)")
    call near_double_null()
    call check_summary("spline_xpoint")

contains

    subroutine near_double_null()
        !! Lower X-point on Gamma_p, upper X-point 15 cm above the smooth top
        !! corner, inside the scrape-off layer.
        type(xpoint_curves_t) :: curves
        type(block_layout_t) :: near
        type(spline_space_options_t) :: options
        type(spline_space_t) :: nspace
        type(spline_quadrature_t) :: quad
        real(dp), parameter :: x2(2) = [r0 + 0.02_dp, 1.15_dp]
        integer :: status, v

        curves = double_null(0.6_dp, 1.0_dp, 0.0_dp, [.true., .false., .false., &
            .false.])
        curves%has_secondary(3) = .true.
        curves%secondary(:, 3) = x2
        call layout_xpoint([r0, 0.0_dp], 0.4_dp, curves, near, status)
        v = minloc(norm2(near%vertices - spread(x2, 2, near%vertex_count), dim=1), 1)
        options%degree = 3
        allocate(options%knots(2, near%block_count))
        options%knots(:, :)%element_count = 8
        call spline_space_build(near, options, nspace, status)
        call spline_quadrature_build(nspace, 5, quad, status)
        call check_condition(status == 0 .and. near%block_count == 16 .and. &
            norm2(near%vertices(:, v) - x2) < 1.0e-14_dp .and. &
            abs(sum(quad%weight) - pi*wall_a*wall_b) < 1.0e-4_dp, &
            "near double null: secondary X-point is a vertex of the SOL layers")
    end subroutine near_double_null

    subroutine build(cells, space, status)
        integer, intent(in) :: cells
        type(spline_space_t), intent(out) :: space
        integer, intent(out) :: status

        type(spline_space_options_t) :: options

        options%degree = 3
        allocate(options%knots(2, layout%block_count))
        options%knots(:, :)%element_count = cells
        call spline_space_build(layout, options, space, status)
    end subroutine build

    subroutine solve(space, error, area, green)
        type(spline_space_t), intent(in) :: space
        real(dp), intent(out) :: error, area, green

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
        call spline_side_quadrature(space, LAYOUT_PLASMA_BOUNDARY, 5, side, status)
        green = sum(side%weight*side%s(1, :)*side%tangent(2, :))
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

end program test_spline_xpoint
