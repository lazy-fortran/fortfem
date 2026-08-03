program test_torus_curved_helmholtz_calderon_slow
    use check, only: check_condition, check_summary
    use fortfem_boundary, only: &
        assemble_helmholtz_torus_curved_calderon_3d, &
        solve_helmholtz_bem_dtn_torus_curved_3d
    use fortfem_core, only: generate_torus_surface_mesh
    use fortfem_kinds, only: dp
    implicit none

    real(dp), parameter :: major_radius = 2.0_dp, minor_radius = 0.6_dp
    real(dp), parameter :: wave_number = 0.7_dp
    integer, parameter :: counts(2) = [5, 7]
    integer, allocatable :: triangles(:, :)
    real(dp), allocatable :: mass(:, :), parameters(:, :), vertices(:, :)
    complex(dp), allocatable :: adjoint(:, :), double_layer(:, :)
    complex(dp), allocatable :: flux(:), hypersingular(:, :)
    complex(dp), allocatable :: single_layer(:, :), trace(:)
    complex(dp), allocatable :: first_term(:), second_term(:)
    real(dp) :: displacement(3), residuals(2), radius, source(3)
    integer :: level, status, vertex
    logical :: all_passed

    all_passed = .true.
    source = [major_radius, 0.0_dp, 0.0_dp]
    do level = 1, 2
        call generate_torus_surface_mesh( &
            major_radius, minor_radius, counts(level), counts(level) + 2, &
            vertices, triangles, parameters)
        allocate(trace(size(vertices, 2)))
        do vertex = 1, size(vertices, 2)
            displacement = vertices(:, vertex) - source
            radius = norm2(displacement)
            trace(vertex) = &
                exp(cmplx(0.0_dp, wave_number*radius, dp))/radius
        end do
        call solve_helmholtz_bem_dtn_torus_curved_3d( &
            parameters, triangles, major_radius, minor_radius, wave_number, &
            trace, 5, flux, status)
        if (status /= 0) error stop "Helmholtz Calderon flux solve failed"
        call assemble_helmholtz_torus_curved_calderon_3d( &
            parameters, triangles, major_radius, minor_radius, wave_number, 5, &
            single_layer, double_layer, adjoint, hypersingular, mass, status)
        if (status /= 0) error stop "curved Helmholtz Calderon assembly failed"
        allocate(first_term(size(trace)), second_term(size(trace)))
        first_term = matmul(hypersingular, trace)
        second_term = matmul( &
            adjoint + 0.5_dp*cmplx(transpose(mass), 0.0_dp, dp), flux)
        residuals(level) = complex_norm(first_term + second_term)/ &
            max(complex_norm(first_term) + complex_norm(second_term), &
            tiny(1.0_dp))
        if (level == 2) then
            call record_condition(maxval(abs(adjoint - &
                transpose(double_layer))) < 1.0e-13_dp, &
                "curved Helmholtz adjoint double layer is the transpose")
            call record_condition(maxval(abs(hypersingular - &
                transpose(hypersingular))) < 2.0e-12_dp, &
                "curved Helmholtz hypersingular operator is complex symmetric")
        end if
        deallocate( &
            vertices, triangles, parameters, trace, flux, single_layer, &
            double_layer, adjoint, hypersingular, mass, first_term, second_term)
    end do

    write (*, '(a,2(es12.4,1x))') "second Calderon residuals: ", residuals
    call record_condition(residuals(2) < residuals(1), &
        "curved Helmholtz second Calderon equation converges under refinement")
    call record_condition(residuals(2) < 3.0e-1_dp, &
        "curved Helmholtz second Calderon equation matches the outgoing trace")
    call check_summary("Curved torus Helmholtz Calderon")
    if (.not. all_passed) error stop 1

contains

    pure function complex_norm(values) result(value)
        complex(dp), intent(in) :: values(:)
        real(dp) :: value

        value = sqrt(sum(abs(values)**2))
    end function complex_norm

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_torus_curved_helmholtz_calderon_slow
