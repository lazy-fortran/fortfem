program test_harmonic_disc
    !! Harmonic map onto the unit disc and the pullback O-grid: the circle
    !! (exact map), an ellipse (symmetry, round trip, convergence in the
    !! samples), a nonconvex bean and a rounded triangle (valid meshes with
    !! certificate, falling wall distance and area error, preserved area),
    !! and the curvature-graded boundary correspondence (the circle is
    !! unchanged by it; on the tipped shapes the coarse mesh wall is closer
    !! to the curve and stays certified).
    use check, only: check_condition, check_summary
    use fortfem_kinds, only: dp
    use fortfem_spline_layout, only: block_layout_t, layout_ogrid
    use fortfem_spline_space, only: spline_space_t, spline_space_options_t, &
        spline_space_build, SPLINE_C0
    use fortfem_spline_eval, only: spline_geometry_at
    use fortfem_harmonic_disc, only: harmonic_disc_t, harmonic_disc_build, &
        harmonic_disc_forward, harmonic_disc_inverse, harmonic_disc_boundary, &
        harmonic_disc_pull_back, harmonic_disc_wall_error, &
        spline_jacobian_certificate
    implicit none

    real(dp), parameter :: pi = acos(-1.0_dp)
    complex(dp), parameter :: iu = (0.0_dp, 1.0_dp)
    real(dp), parameter :: GRADING_GAIN = 2.5_dp ! measured 3.9 (bean), 3.0 (triangle)

    call circle_is_exact()
    call ellipse_properties()
    call shaped_mesh("bean", 1)
    call shaped_mesh("rounded triangle", 2)
    call graded_boundary()
    call check_summary("harmonic disc map")

contains

    subroutine curve(kind, n, samples)
        !! 0: circle, 1: bean (concave at t = pi), 2: rounded triangle,
        !! 3: ellipse; equally spaced in the parameter t.
        integer, intent(in) :: kind, n
        real(dp), allocatable, intent(out) :: samples(:, :)

        real(dp) :: t
        complex(dp) :: z
        integer :: i

        allocate(samples(2, n))
        do i = 1, n
            t = 2.0_dp*pi*real(i - 1, dp)/real(n, dp)
            select case (kind)
            case (0)
                z = cmplx(0.3_dp, -0.2_dp, dp) + 1.7_dp*exp(iu*t)
            case (1)
                z = cmplx(cos(t) + 0.45_dp*cos(2.0_dp*t) - 0.45_dp, &
                    1.6_dp*sin(t), dp)
            case (2)
                z = exp(iu*t) + 0.2_dp*exp(-2.0_dp*iu*t)
            case default
                z = cmplx(2.0_dp*cos(t), sin(t), dp)
            end select
            samples(:, i) = [real(z, dp), aimag(z)]
        end do
    end subroutine curve

    subroutine circle_is_exact()
        type(harmonic_disc_t) :: map
        real(dp), allocatable :: samples(:, :)
        complex(dp) :: z, g, dg(2), c, exact
        real(dp) :: err, derr, ierr, radius(4)
        integer :: status, i, j

        call curve(0, 64, samples)
        call harmonic_disc_build(map, samples, status)
        call check_condition(status == 0 .and. map%residual < 1.0e-12_dp, &
            "circle: the map is built and its boundary values are exact")
        c = cmplx(0.3_dp, -0.2_dp, dp)
        radius = [0.0_dp, 0.5_dp, 0.99_dp, 0.99999_dp]
        err = 0.0_dp
        derr = 0.0_dp
        ierr = 0.0_dp
        do j = 1, 4
            do i = 1, 7
                exact = radius(j)*exp(iu*0.9_dp*real(i, dp))
                z = c + 1.7_dp*exact
                call harmonic_disc_forward(map, z, g, dg)
                err = max(err, abs(g - exact))
                derr = max(derr, abs(dg(1) - 1.0_dp/1.7_dp), abs(dg(2) + iu/1.7_dp))
                call harmonic_disc_inverse(map, exact, g, status)
                ierr = max(ierr, abs(g - z))
            end do
        end do
        call check_condition(err < 1.0e-12_dp, &
            "circle: G is the affine map up to the wall")
        call check_condition(derr < 1.0e-9_dp, "circle: its derivative is exact")
        call check_condition(ierr < 1.0e-11_dp, "circle: the inverse is exact")
    end subroutine circle_is_exact

    subroutine ellipse_properties()
        type(harmonic_disc_t) :: map, fine
        real(dp), allocatable :: samples(:, :)
        complex(dp) :: z, g, g2, dg(2), back
        real(dp) :: sym, trip, conv, det, detmin
        integer :: status, i, j

        call curve(3, 128, samples)
        call harmonic_disc_build(map, samples, status)
        call check_condition(status == 0 .and. map%residual < 1.0e-10_dp, &
            "ellipse: built with boundary residual below 1e-10")
        call curve(3, 256, samples)
        call harmonic_disc_build(fine, samples, status)
        call harmonic_disc_boundary(map, pi, z)
        call check_condition(abs(z + 2.0_dp) < 1.0e-11_dp, &
            "ellipse: half the arc length ends at the opposite vertex")
        sym = 0.0_dp
        trip = 0.0_dp
        conv = 0.0_dp
        detmin = huge(1.0_dp)
        do j = 1, 5
            do i = 1, 9
                z = 0.19_dp*real(j, dp)*cmplx(2.0_dp*cos(0.7_dp*real(i, dp)), &
                    sin(0.7_dp*real(i, dp)), dp)
                call harmonic_disc_forward(map, z, g, dg)
                call harmonic_disc_forward(map, conjg(z), g2, dg)
                sym = max(sym, abs(g2 - conjg(g)))
                call harmonic_disc_inverse(map, g, back, status)
                trip = max(trip, abs(back - z))
                call harmonic_disc_forward(fine, z, g2, dg)
                conv = max(conv, abs(g2 - g))
                det = -real(dg(1), dp)*aimag(dg(2)) + aimag(dg(1))*real(dg(2), dp)
                detmin = min(detmin, det)
            end do
        end do
        call check_condition(sym < 1.0e-11_dp, "ellipse: G commutes with the mirror")
        call check_condition(trip < 1.0e-10_dp, "ellipse: G^-1 o G is the identity")
        call check_condition(conv < 1.0e-10_dp, &
            "ellipse: 128 and 256 samples give the same map")
        call check_condition(detmin > 0.0_dp, "ellipse: det grad G is positive")
    end subroutine ellipse_properties

    subroutine graded_boundary()
        !! grading > 0 weights the boundary correspondence with the
        !! curvature. A circle has constant curvature, so its map does not
        !! change. On the bean and the rounded triangle the wall of a coarse
        !! mesh (4 cells per block side, degree 2) comes closer to the curve
        !! and the mesh stays certified, also with the preserved area.
        type(harmonic_disc_t) :: map, plain
        type(spline_space_t) :: space
        real(dp), allocatable :: samples(:, :)
        complex(dp) :: z, g, dg(2), c
        real(dp) :: distance(0:1), area, ratio
        logical :: certified
        integer :: status, kind, graded
        character(16) :: name

        call curve(0, 64, samples)
        call harmonic_disc_build(map, samples, status, grading=2.0_dp)
        call check_condition(status == 0, "graded circle: built")
        c = cmplx(0.3_dp, -0.2_dp, dp)
        z = c + 1.7_dp*cmplx(0.3_dp, 0.5_dp, dp)
        call harmonic_disc_forward(map, z, g, dg)
        call check_condition(abs(g - (z - c)/1.7_dp) < 1.0e-12_dp .and. &
            map%residual < 1.0e-12_dp, "graded circle: the map is unchanged")
        call harmonic_disc_build(map, samples, status, grading=-1.0_dp)
        call check_condition(status == 1, "negative grading is refused")
        call curve(2, 256, samples)
        call harmonic_disc_build(plain, samples, status)
        call harmonic_disc_build(map, samples, status, grading=0.0_dp)
        call harmonic_disc_boundary(plain, 1.0_dp, z)
        call harmonic_disc_boundary(map, 1.0_dp, g)
        call check_condition(z == g, "grading 0 is the arc-length correspondence")
        do kind = 1, 2
            name = merge("bean            ", "rounded triangle", kind == 1)
            do graded = 0, 1
                call mesh(kind, 2, 4, .true., space, map, status, real(graded, dp))
                call check_condition(status == 0, trim(name)//": graded mesh built")
                call harmonic_disc_wall_error(space, map, distance(graded), area, status)
                call spline_jacobian_certificate(space, certified, ratio, status)
                print '(a,a,i2,a,es9.2,a,es9.2,a,es9.2,a,es9.2)', trim(name), ": grading", &
                    graded, " wall distance", distance(graded), " area error", area, &
                    " Bernstein ratio", ratio, " map residual", map%residual
                call check_condition(certified .and. abs(area) < 1.0e-11_dp .and. &
                    map%residual < 1.0e-8_dp, &
                    trim(name)//": graded mesh certified, area kept, samples resolve the data")
            end do
            call check_condition(distance(1) < distance(0)/GRADING_GAIN, &
                trim(name)//": grading brings the coarse wall closer to the curve")
        end do
    end subroutine graded_boundary

    subroutine mesh(kind, degree, cells, keep_area, space, map, status, grading)
        integer, intent(in) :: kind, degree, cells
        logical, intent(in) :: keep_area
        type(spline_space_t), intent(out) :: space
        type(harmonic_disc_t), intent(out) :: map
        integer, intent(out) :: status
        real(dp), intent(in), optional :: grading

        type(block_layout_t) :: layout
        type(spline_space_options_t) :: options
        real(dp), allocatable :: samples(:, :)

        call curve(kind, 256, samples)
        call harmonic_disc_build(map, samples, status, grading)
        if (status /= 0) return
        call layout_ogrid([0.0_dp, 0.0_dp], 1.0_dp, 0.5_dp, layout, status)
        if (status /= 0) return
        options%degree = degree
        options%regularity = SPLINE_C0
        allocate(options%knots(2, layout%block_count))
        options%knots(:, :)%element_count = cells
        call spline_space_build(layout, options, space, status)
        if (status /= 0) return
        call harmonic_disc_pull_back(space, map, status, preserve_area=keep_area)
    end subroutine mesh

    subroutine shaped_mesh(name, kind)
        character(*), intent(in) :: name
        integer, intent(in) :: kind

        type(spline_space_t) :: space
        type(harmonic_disc_t) :: map
        real(dp) :: distance(2), area(2), ratio, jet(2, 6), kept
        complex(dp) :: g, dg(2)
        logical :: certified, all_certified
        integer :: status, degree, level, failed

        do degree = 2, 3
            all_certified = .true.
            do level = 1, 2
                call mesh(kind, degree, 8*level, .false., space, map, status)
                call check_condition(status == 0, name//": the pullback mesh is built")
                call spline_jacobian_certificate(space, certified, ratio, status, &
                    failed)
                all_certified = all_certified .and. certified .and. status == 0
                call harmonic_disc_wall_error(space, map, distance(level), &
                    area(level), status)
                print '(a,a,i2,a,i3,a,es9.2,a,es9.2,a,es10.2,a,i0)', name, &
                    ": degree", degree, " cells", 8*level, " Bernstein ratio", ratio, &
                    " wall distance", distance(level), " area error", area(level), &
                    " failed cells ", failed
            end do
            call check_condition(all_certified, &
                name//": every cell has positive Bernstein coefficients of det")
            ! the tips are not resolved at 8 cells, where the order of the
            ! splines is not reached yet (measured 3.3 to 12.6 per halving)
            call check_condition(distance(2) < distance(1)/3.0_dp .and. &
                abs(area(2)) < abs(area(1))/3.0_dp, &
                name//": wall distance and area error fall by more than 3 per halving")
            call mesh(kind, degree, 4, .true., space, map, status)
            call harmonic_disc_wall_error(space, map, distance(1), kept, status)
            call spline_jacobian_certificate(space, certified, ratio, status)
            print '(a,a,es10.2)', name, ": area error, area-preserving wall", kept
            call check_condition(abs(kept) < 1.0e-11_dp .and. certified, &
                name//": the area-preserving wall keeps the area and validity")
        end do
        ! the mesh is the pullback: G of the image of the core centre is 0
        call mesh(kind, 3, 8, .false., space, map, status)
        call spline_geometry_at(space, 1, [0.5_dp, 0.5_dp], jet, status)
        call harmonic_disc_forward(map, cmplx(jet(1, 1), jet(2, 1), dp), g, dg)
        call check_condition(abs(g) < 1.0e-4_dp, &
            name//": the centre of the core block is the preimage of 0")
    end subroutine shaped_mesh

end program test_harmonic_disc
