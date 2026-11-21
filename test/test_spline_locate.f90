module spline_locate_geometry
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

end module spline_locate_geometry

program test_spline_locate
    !! Point location on every block side: points on all sides and random
    !! points within 1e-3 of every interior side, annulus interface and
    !! periodic seam (away from the wall) are found, with hints 1 and the
    !! own block, and the map reproduces them to 1e-12. Layouts: D-shaped O-grid with periodic
    !! annulus (including the interface samples) and the double null.
    use check, only: check_condition, check_summary
    use fortfem_kinds, only: dp
    use fortfem_spline_layout, only: block_layout_t, layout_ogrid, &
        layout_set_curve, LAYOUT_WALL, LAYOUT_ANNULUS_INTERFACE, side_start, side_end
    use fortfem_spline_xpoint, only: layout_xpoint
    use fortfem_spline_space, only: spline_space_t, spline_space_options_t, &
        spline_space_build
    use fortfem_spline_eval, only: spline_locate, spline_geometry_at
    use spline_locate_geometry, only: double_null, r0
    implicit none

    type(block_layout_t) :: layout
    type(spline_space_t) :: space
    real(dp) :: samples(2, 256)
    integer :: status, misses, tried, q

    call layout_ogrid([1.7_dp, 0.0_dp], 0.8_dp, 0.3_dp, layout, status, &
        interface_fraction=0.6_dp)
    call layout_set_curve(layout, LAYOUT_WALL, d_shape(0.8_dp), status)
    samples = d_shape(0.5_dp)
    call layout_set_curve(layout, LAYOUT_ANNULUS_INTERFACE, samples, status)
    call build(layout, 8, space)
    call sweep(space, misses, tried)
    do q = 1, 256
        call try(space, samples(:, q), 1, misses, tried)
        call try(space, samples(:, q), 6, misses, tried)
    end do
    print '(a,i0,a,i0)', "annulus layout: misses ", misses, " of ", tried
    call check_condition(misses == 0, &
        "D-shaped O-grid with annulus: sides, interface samples and seam located")

    call layout_xpoint([r0, 0.0_dp], 0.4_dp, double_null(0.6_dp, 1.0_dp, 0.0_dp, &
        [.true., .false., .true., .false.]), layout, status)
    call build(layout, 6, space)
    call sweep(space, misses, tried)
    print '(a,i0,a,i0)', "double null: misses ", misses, " of ", tried
    call check_condition(misses == 0, "double null: all sides and X-points located")
    call check_summary("spline_locate")

contains

    pure function d_shape(minor) result(points)
        real(dp), intent(in) :: minor
        real(dp) :: points(2, 256), theta
        integer :: q

        do q = 1, 256
            theta = 2.0_dp*acos(-1.0_dp)*(q - 1)/256
            points(:, q) = [1.7_dp + minor*cos(theta + 0.5_dp*sin(theta)), &
                1.8_dp*minor*sin(theta)]
        end do
    end function d_shape

    subroutine build(layout, cells, space)
        type(block_layout_t), intent(in) :: layout
        integer, intent(in) :: cells
        type(spline_space_t), intent(out) :: space

        type(spline_space_options_t) :: options
        integer :: status

        options%degree = 3
        allocate(options%knots(2, layout%block_count))
        options%knots(:, :)%element_count = cells
        call spline_space_build(layout, options, space, status)
        if (status /= 0) error stop "space"
    end subroutine build

    subroutine sweep(space, misses, tried)
        type(spline_space_t), intent(in) :: space
        integer, intent(out) :: misses, tried

        real(dp) :: jet(2, 6), t, off(2)
        integer :: k, s, q, r, status
        logical :: interior

        misses = 0
        tried = 0
        do k = 1, space%layout%block_count
            do s = 1, 4
                interior = space%layout%neighbour(s, k) /= 0 .or. &
                    space%layout%side_tag(s, k) == LAYOUT_ANNULUS_INTERFACE .or. &
                    (space%layout%periodic(k) .and. (s == 2 .or. s == 4))
                do q = 0, 40
                    t = real(q, dp)/40.0_dp
                    call spline_geometry_at(space, k, side_xi(s, t), jet, status)
                    call try(space, jet(:, 1), 1, misses, tried)
                    call try(space, jet(:, 1), k, misses, tried)
                    ! perturbed points near the wall may leave the domain
                    if (.not. interior) cycle
                    if (t < 0.05_dp .and. space%layout%boundary_vertex( &
                        space%layout%corners(side_start(s), k))) cycle
                    if (t > 0.95_dp .and. space%layout%boundary_vertex( &
                        space%layout%corners(side_end(s), k))) cycle
                    do r = 1, 3
                        call random_number(off)
                        call try(space, jet(:, 1) + 2.0e-3_dp*(off - 0.5_dp), k, &
                            misses, tried)
                    end do
                end do
            end do
        end do
    end subroutine sweep

    subroutine try(space, s, hint, misses, tried)
        type(spline_space_t), intent(in) :: space
        real(dp), intent(in) :: s(2)
        integer, intent(in) :: hint
        integer, intent(inout) :: misses, tried

        real(dp) :: xi(2), jet(2, 6)
        integer :: block, status

        tried = tried + 1
        call spline_locate(space, hint, s, block, xi, status)
        if (status == 0) call spline_geometry_at(space, block, xi, jet, status)
        if (status /= 0) then
            misses = misses + 1
        else if (norm2(jet(:, 1) - s) > 1.0e-12_dp*(1.0_dp + norm2(s))) then
            misses = misses + 1
        end if
    end subroutine try

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

end program test_spline_locate
