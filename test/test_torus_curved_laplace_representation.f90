program test_torus_curved_laplace_representation
    use check, only: check_condition, check_summary
    use fortfem_core, only: cartesian_to_toroidal, generate_torus_surface_mesh, &
        toroidal_point_to_cartesian
    use fortfem_boundary, only: evaluate_laplace_representation_torus_curved_3d
    use fortfem_fourier, only: evaluate_toroidal_harmonic_p, &
        toroidal_poisson_exterior_dtn_p
    use fortfem_kinds, only: dp
    implicit none

    real(dp), parameter :: major_radius = 2.0_dp, minor_radius = 0.6_dp
    real(dp), parameter :: scale = &
        sqrt(major_radius**2 - minor_radius**2)
    real(dp), parameter :: boundary_eta = acosh(major_radius/minor_radius)
    integer, parameter :: counts(2) = [48, 72]
    integer, allocatable :: triangles(:, :)
    real(dp), allocatable :: flux(:), parameters(:, :), trace(:), vertices(:, :)
    real(dp) :: dtn_value, errors(2), eta, exact, numerical
    real(dp) :: normal_derivative, phi, target(3), theta, value
    integer :: level, status, vertex
    logical :: all_passed

    all_passed = .true.
    call toroidal_point_to_cartesian(scale, 0.35_dp, 0.8_dp, 0.3_dp, target)
    call evaluate_toroidal_harmonic_p( &
        2, 1, 0.35_dp, 0.8_dp, 0.3_dp, exact, status)
    do level = 1, 2
        call generate_torus_surface_mesh( &
            major_radius, minor_radius, counts(level), counts(level) + 2, &
            vertices, triangles, parameters)
        allocate(trace(size(vertices, 2)), flux(size(vertices, 2)))
        do vertex = 1, size(vertices, 2)
            call cartesian_to_toroidal( &
                vertices(:, vertex), scale, eta, theta, phi)
            call evaluate_toroidal_harmonic_p( &
                2, 1, boundary_eta, theta, phi, trace(vertex), status)
            if (status /= 0) error stop "curved torus trace failed"
            call toroidal_poisson_exterior_dtn_p( &
                2, 1, scale, boundary_eta, theta, phi, value, &
                normal_derivative, dtn_value, status)
            if (status /= 0) error stop "curved torus flux failed"
            flux(vertex) = normal_derivative
        end do
        call evaluate_laplace_representation_torus_curved_3d( &
            parameters, triangles, major_radius, minor_radius, trace, flux, &
            target, 10, numerical, status)
        errors(level) = abs(numerical - exact)/abs(exact)
        deallocate(vertices, triangles, parameters, trace, flux)
    end do

    call record_condition(errors(2) < 0.65_dp*errors(1), &
        "curved torus BEM representation converges under refinement")
    if (errors(2) >= 2.5e-2_dp) then
        write(*, '(a,2es14.5)') "curved torus representation errors ", errors
    end if
    call record_condition(errors(2) < 2.5e-2_dp, &
        "curved torus BEM matches the half-integer Legendre harmonic")
    call check_summary("Curved torus Laplace representation")
    if (.not. all_passed) error stop 1

contains

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_torus_curved_laplace_representation
