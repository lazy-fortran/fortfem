program test_torus_curved_helmholtz_costabel_slow
    use check, only: check_condition, check_summary
    use fortfem_boundary, only: &
        assemble_helmholtz_fem_bem_costabel_torus_curved_3d, &
        solve_helmholtz_fem_bem_costabel_torus_curved_3d
    use fortfem_core, only: generate_solid_torus_tetra_mesh
    use fortfem_kinds, only: dp
    implicit none

    real(dp), parameter :: major_radius = 2.0_dp, minor_radius = 0.6_dp
    real(dp), parameter :: wave_number = 0.7_dp
    integer, parameter :: counts(2) = [5, 7]
    integer, allocatable :: boundary_triangles(:, :), tetrahedra(:, :)
    real(dp), allocatable :: parameters(:, :), vertices(:, :)
    complex(dp), allocatable :: exact(:), flux(:), load(:), matrix(:, :)
    complex(dp), allocatable :: potential(:)
    logical, allocatable :: is_boundary(:)
    real(dp) :: displacement(3), errors(2), radius, source(3)
    integer :: level, local_node, node, panel, status
    logical :: all_passed

    all_passed = .true.
    do level = 1, 2
        call generate_solid_torus_tetra_mesh( &
            major_radius, minor_radius, 2, counts(level), counts(level) + 2, &
            vertices, tetrahedra, boundary_triangles, parameters)
        allocate(load(size(vertices, 2)))
        load = cmplx(0.0_dp, 0.0_dp, dp)
        load(1) = cmplx(1.0_dp, 0.0_dp, dp)
        source = vertices(:, 1)
        call solve_helmholtz_fem_bem_costabel_torus_curved_3d( &
            vertices, tetrahedra, parameters, boundary_triangles, &
            major_radius, minor_radius, wave_number, wave_number, load, 5, &
            potential, flux, status)
        if (status /= 0) error stop "curved torus Helmholtz Costabel failed"
        allocate(is_boundary(size(vertices, 2)), exact(size(vertices, 2)))
        is_boundary = .false.
        do panel = 1, size(boundary_triangles, 2)
            do local_node = 1, 3
                is_boundary(boundary_triangles(local_node, panel)) = .true.
            end do
        end do
        exact = cmplx(0.0_dp, 0.0_dp, dp)
        do node = 1, size(vertices, 2)
            if (is_boundary(node)) then
                displacement = vertices(:, node) - source
                radius = norm2(displacement)
                exact(node) = exp(cmplx(0.0_dp, wave_number*radius, dp))/ &
                    (4.0_dp*acos(-1.0_dp)*radius)
            end if
        end do
        errors(level) = complex_norm(pack(potential - exact, is_boundary))/ &
            complex_norm(pack(exact, is_boundary))
        if (level == 1) then
            call assemble_helmholtz_fem_bem_costabel_torus_curved_3d( &
                vertices, tetrahedra, parameters, boundary_triangles, &
                major_radius, minor_radius, wave_number, wave_number, 5, &
                matrix, status)
            call record_condition(status == 0 .and. maxval(abs( &
                matrix - transpose(matrix))) < 2.0e-11_dp, &
                "curved torus Helmholtz Costabel matrix is complex symmetric")
            deallocate(matrix)
        end if
        deallocate( &
            vertices, tetrahedra, boundary_triangles, parameters, load, &
            potential, flux, is_boundary, exact)
    end do

    write (*, '(a,2(es12.4,1x))') "boundary potential errors: ", errors
    call record_condition(errors(2) < 0.9_dp*errors(1), &
        "curved torus Helmholtz Costabel potential converges under refinement")
    call record_condition(errors(2) < 3.5e-1_dp, &
        "curved torus Helmholtz Costabel matches the outgoing point source")
    call check_summary("Curved torus Helmholtz Costabel coupling")
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

end program test_torus_curved_helmholtz_costabel_slow
