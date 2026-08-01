program test_torus_curved_laplace_representation_gradient
    use check, only: check_condition, check_summary
    use fortfem_api, only: &
        evaluate_laplace_representation_torus_curved_3d, &
        generate_torus_surface_mesh
    use fortfem_kinds, only: dp
    implicit none

    real(dp), parameter :: pi = acos(-1.0_dp)
    real(dp), parameter :: major_radius = 2.0_dp, minor_radius = 0.6_dp
    integer, parameter :: counts(2) = [8, 12]
    integer, allocatable :: triangles(:, :)
    real(dp), allocatable :: parameters(:, :), vertices(:, :)
    real(dp), allocatable :: flux(:), trace(:)
    real(dp) :: displacement(3), errors(2), exact, numerical
    real(dp) :: exact_gradient(3), gradient(3), normal(3), radius
    real(dp) :: source(3), target(3), theta, phi
    integer :: level, status, vertex
    logical :: all_passed

    all_passed = .true.
    source = [major_radius, 0.0_dp, 0.0_dp]
    target = [0.0_dp, 0.0_dp, 0.0_dp]
    displacement = target - source
    radius = norm2(displacement)
    exact = 1.0_dp/(4.0_dp*pi*radius)
    exact_gradient = -displacement/(4.0_dp*pi*radius**3)

    do level = 1, 2
        call generate_torus_surface_mesh( &
            major_radius, minor_radius, counts(level), counts(level) + 2, &
            vertices, triangles, parameters)
        allocate(trace(size(vertices, 2)), flux(size(vertices, 2)))
        do vertex = 1, size(vertices, 2)
            displacement = vertices(:, vertex) - source
            radius = norm2(displacement)
            trace(vertex) = 1.0_dp/(4.0_dp*pi*radius)
            theta = parameters(1, vertex)
            phi = parameters(2, vertex)
            normal = [ &
                cos(theta)*cos(phi), cos(theta)*sin(phi), sin(theta)]
            flux(vertex) = -dot_product(displacement, normal)/ &
                (4.0_dp*pi*radius**3)
        end do
        call evaluate_laplace_representation_torus_curved_3d( &
            parameters, triangles, major_radius, minor_radius, trace, flux, &
            target, 5, numerical, status, gradient)
        if (status /= 0) error stop "curved Laplace representation failed"
        errors(level) = abs(numerical - exact)/abs(exact)
        if (level == 2) then
            call record_condition(norm2(gradient - exact_gradient) < 3.0e-2_dp, &
                "curved Laplace representation recovers the exterior gradient")
        end if
        deallocate(vertices, triangles, parameters, trace, flux)
    end do

    write (*, '(a,2(es12.4,1x))') "relative field errors: ", errors
    call record_condition(errors(2) < 0.7_dp*errors(1), &
        "curved Laplace representation converges under refinement")
    call record_condition(errors(2) < 1.0e-1_dp, &
        "curved Laplace representation matches the fundamental solution")
    call check_summary("Curved torus Laplace representation gradient")
    if (.not. all_passed) error stop 1

contains

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_torus_curved_laplace_representation_gradient
