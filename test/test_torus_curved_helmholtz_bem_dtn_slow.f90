program test_torus_curved_helmholtz_bem_dtn_slow
    use check, only: check_condition, check_summary
    use fortfem_boundary, only: &
        evaluate_torus_curved_panel, &
        solve_helmholtz_bem_dtn_torus_curved_3d
    use fortfem_core, only: generate_torus_surface_mesh
    use fortfem_kinds, only: dp
    implicit none

    real(dp), parameter :: major_radius = 2.0_dp, minor_radius = 0.6_dp
    real(dp), parameter :: wave_number = 0.7_dp
    integer, parameter :: counts(2) = [7, 9]
    integer, allocatable :: triangles(:, :)
    real(dp), allocatable :: parameters(:, :), vertices(:, :)
    complex(dp), allocatable :: exact_flux(:), numerical_flux(:), trace(:)
    complex(dp) :: value
    real(dp) :: denominator, displacement(3), errors(2), jacobian
    real(dp) :: normal(3), point(3), radius, source(3), tangent_eta(3)
    real(dp) :: tangent_xi(3)
    integer :: level, panel, status, vertex
    logical :: all_passed

    all_passed = .true.
    source = [major_radius, 0.0_dp, 0.0_dp]
    do level = 1, 2
        call generate_torus_surface_mesh( &
            major_radius, minor_radius, counts(level), counts(level) + 2, &
            vertices, triangles, parameters)
        allocate(trace(size(vertices, 2)), exact_flux(size(triangles, 2)))
        do vertex = 1, size(vertices, 2)
            displacement = vertices(:, vertex) - source
            radius = norm2(displacement)
            trace(vertex) = &
                exp(cmplx(0.0_dp, wave_number*radius, dp))/radius
        end do
        do panel = 1, size(triangles, 2)
            call evaluate_torus_curved_panel( &
                parameters(:, triangles(:, panel)), major_radius, minor_radius, &
                1.0_dp/3.0_dp, 1.0_dp/3.0_dp, point, tangent_xi, tangent_eta, &
                jacobian, status)
            if (status /= 0) error stop "curved torus panel evaluation failed"
            normal = cross_product(tangent_xi, tangent_eta)/jacobian
            displacement = point - source
            radius = norm2(displacement)
            value = exp(cmplx(0.0_dp, wave_number*radius, dp))/radius
            exact_flux(panel) = value* &
                cmplx(-1.0_dp/radius, wave_number, dp)* &
                dot_product(displacement/radius, normal)
        end do
        call solve_helmholtz_bem_dtn_torus_curved_3d( &
            parameters, triangles, major_radius, minor_radius, wave_number, &
            trace, 7, numerical_flux, status)
        if (status /= 0) error stop "curved torus Helmholtz DtN solve failed"
        denominator = max(sqrt(sum(abs(exact_flux)**2)), tiny(1.0_dp))
        errors(level) = &
            sqrt(sum(abs(numerical_flux - exact_flux)**2))/denominator
        deallocate(vertices, triangles, parameters, trace, exact_flux)
        deallocate(numerical_flux)
    end do

    write (*, '(a,2(es12.4,1x))') "relative flux errors: ", errors
    call record_condition(errors(2) < 0.8_dp*errors(1), &
        "curved torus Helmholtz BEM DtN flux converges under refinement")
    call record_condition(errors(2) < 3.5e-1_dp, &
        "curved torus Helmholtz BEM DtN matches the outgoing-source flux")
    call check_summary("Curved torus Helmholtz BEM DtN")
    if (.not. all_passed) error stop 1

contains

    pure function cross_product(first, second) result(product)
        real(dp), intent(in) :: first(3), second(3)
        real(dp) :: product(3)

        product = [ &
            first(2)*second(3) - first(3)*second(2), &
            first(3)*second(1) - first(1)*second(3), &
            first(1)*second(2) - first(2)*second(1)]
    end function cross_product

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_torus_curved_helmholtz_bem_dtn_slow
