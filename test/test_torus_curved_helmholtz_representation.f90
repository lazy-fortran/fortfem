program test_torus_curved_helmholtz_representation
    use check, only: check_condition, check_summary
    use fortfem_boundary, only: evaluate_helmholtz_representation_torus_curved_3d
    use fortfem_core, only: generate_torus_surface_mesh
    use fortfem_kinds, only: dp
    implicit none

    real(dp), parameter :: major_radius = 2.0_dp, minor_radius = 0.6_dp
    real(dp), parameter :: wave_number = 0.7_dp
    integer, parameter :: counts(2) = [24, 36]
    integer, allocatable :: triangles(:, :)
    real(dp), allocatable :: parameters(:, :), vertices(:, :)
    complex(dp), allocatable :: flux(:), trace(:)
    complex(dp) :: exact, numerical, source_value
    real(dp) :: displacement(3), errors(2), normal(3), phi, radius, source(3)
    real(dp) :: target(3), theta
    integer :: level, status, vertex
    logical :: all_passed

    all_passed = .true.
    source = [major_radius, 0.0_dp, 0.0_dp]
    target = [0.0_dp, 0.0_dp, 0.0_dp]
    do level = 1, 2
        call generate_torus_surface_mesh( &
            major_radius, minor_radius, counts(level), counts(level) + 2, &
            vertices, triangles, parameters)
        allocate(trace(size(vertices, 2)), flux(size(vertices, 2)))
        do vertex = 1, size(vertices, 2)
            displacement = vertices(:, vertex) - source
            radius = norm2(displacement)
            trace(vertex) = &
                exp(cmplx(0.0_dp, wave_number*radius, dp))/radius
            theta = parameters(1, vertex)
            phi = parameters(2, vertex)
            normal = [ &
                cos(theta)*cos(phi), cos(theta)*sin(phi), sin(theta)]
            source_value = trace(vertex)
            flux(vertex) = source_value* &
                cmplx(-1.0_dp/radius, wave_number, dp)* &
                dot_product(displacement/radius, normal)
        end do
        call evaluate_helmholtz_representation_torus_curved_3d( &
            parameters, triangles, major_radius, minor_radius, trace, flux, &
            target, wave_number, 7, numerical, status)
        exact = exp(cmplx(0.0_dp, wave_number*major_radius, dp))/major_radius
        if (status /= 0) error stop "curved Helmholtz representation failed"
        errors(level) = abs(numerical - exact)/abs(exact)
        deallocate(vertices, triangles, parameters, trace, flux)
    end do

    write (*, '(a,2(es12.4,1x))') "relative field errors: ", errors
    call record_condition(errors(2) < 0.65_dp*errors(1), &
        "curved torus Helmholtz representation converges under refinement")
    call record_condition(errors(2) < 5.0e-2_dp, &
        "curved torus Helmholtz representation matches the outgoing field")
    call check_summary("Curved torus Helmholtz representation")
    if (.not. all_passed) error stop 1

contains

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_torus_curved_helmholtz_representation
