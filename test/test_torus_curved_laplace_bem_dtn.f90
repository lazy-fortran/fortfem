program test_torus_curved_laplace_bem_dtn
    use check, only: check_condition, check_summary
    use fortfem_core, only: &
        cartesian_to_toroidal, evaluate_torus_curved_panel, &
        generate_torus_surface_mesh
    use fortfem_boundary, only: solve_laplace_bem_dtn_torus_curved_3d
    use fortfem_fourier, only: &
        evaluate_toroidal_harmonic_p, &
        toroidal_poisson_exterior_dtn_p
    use fortfem_kinds, only: dp
    implicit none

    real(dp), parameter :: major_radius = 2.0_dp, minor_radius = 0.6_dp
    real(dp), parameter :: scale = &
        sqrt(major_radius**2 - minor_radius**2)
    real(dp), parameter :: boundary_eta = acosh(major_radius/minor_radius)
    integer, parameter :: counts(2) = [7, 9]
    integer, allocatable :: triangles(:, :)
    real(dp), allocatable :: exact_flux(:), numerical_flux(:)
    real(dp), allocatable :: parameters(:, :), trace(:), vertices(:, :)
    real(dp) :: denominator, dtn_value, errors(2), eta, jacobian
    real(dp) :: normal_derivative, phi, point(3), tangent_eta(3)
    real(dp) :: tangent_xi(3), theta, value
    integer :: level, panel, status, vertex
    logical :: all_passed

    all_passed = .true.
    do level = 1, 2
        call generate_torus_surface_mesh( &
            major_radius, minor_radius, counts(level), counts(level) + 2, &
            vertices, triangles, parameters)
        allocate(trace(size(vertices, 2)), exact_flux(size(triangles, 2)))
        do vertex = 1, size(vertices, 2)
            call cartesian_to_toroidal( &
                vertices(:, vertex), scale, eta, theta, phi)
            call evaluate_toroidal_harmonic_p( &
                2, 1, boundary_eta, theta, phi, trace(vertex), status)
            if (status /= 0) error stop "toroidal BEM trace failed"
        end do
        do panel = 1, size(triangles, 2)
            call evaluate_torus_curved_panel( &
                parameters(:, triangles(:, panel)), major_radius, minor_radius, &
                1.0_dp/3.0_dp, 1.0_dp/3.0_dp, point, tangent_xi, tangent_eta, &
                jacobian, status)
            call cartesian_to_toroidal(point, scale, eta, theta, phi)
            call toroidal_poisson_exterior_dtn_p( &
                2, 1, scale, boundary_eta, theta, phi, value, &
                normal_derivative, dtn_value, status)
            exact_flux(panel) = normal_derivative
        end do
        call solve_laplace_bem_dtn_torus_curved_3d( &
            parameters, triangles, major_radius, minor_radius, trace, 5, &
            numerical_flux, status)
        if (status /= 0) error stop "curved torus BEM DtN solve failed"
        denominator = max(norm2(exact_flux), tiny(1.0_dp))
        errors(level) = norm2(numerical_flux - exact_flux)/denominator
        deallocate(vertices, triangles, parameters, trace, exact_flux)
        deallocate(numerical_flux)
    end do

    write (*, '(a,2(es12.4,1x))') "relative flux errors: ", errors
    call record_condition(errors(2) < 0.8_dp*errors(1), &
        "curved torus BEM DtN flux converges under refinement")
    call record_condition(errors(2) < 3.5e-1_dp, &
        "curved torus BEM DtN matches the toroidal-harmonic flux")
    call check_summary("Curved torus Laplace BEM DtN")
    if (.not. all_passed) error stop 1

contains

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_torus_curved_laplace_bem_dtn
