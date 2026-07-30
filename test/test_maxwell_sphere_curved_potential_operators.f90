program test_maxwell_sphere_curved_potential_operators
    use check, only: check_condition, check_summary
    use fortfem_api, only: &
        assemble_maxwell_sphere_curved_potential_operators_rwg_3d, &
        assemble_maxwell_sphere_curved_vector_potential_rwg_3d, &
        generate_sphere_surface_mesh
    use fortfem_kinds, only: dp
    implicit none

    integer, allocatable :: triangles(:, :)
    real(dp), allocatable :: scaled_vertices(:, :), vertices(:, :)
    complex(dp), allocatable :: scalar_potential(:, :), scaled_scalar(:, :)
    complex(dp), allocatable :: scaled_vector(:, :), vector_potential(:, :)
    complex(dp), allocatable :: vector_reference(:, :)
    real(dp) :: error, symmetry_error
    integer :: status
    logical :: all_passed

    all_passed = .true.
    call generate_sphere_surface_mesh(1.0_dp, 0, vertices, triangles)
    call assemble_maxwell_sphere_curved_potential_operators_rwg_3d( &
        vertices, triangles, 1.0_dp, 0.6_dp, 4, 1.0e-5_dp, 2, &
        vector_potential, scalar_potential, status)
    call record_condition(status == 0 .and. &
        maxval(abs(scalar_potential)) > 0.0_dp, &
        "curved RWG vector and scalar potentials assemble")
    symmetry_error = max( &
        maxval(abs(vector_potential - transpose(vector_potential))), &
        maxval(abs(scalar_potential - transpose(scalar_potential))))
    call record_condition(symmetry_error < 4.0e-12_dp, &
        "curved RWG potential blocks are complex symmetric")
    call assemble_maxwell_sphere_curved_vector_potential_rwg_3d( &
        vertices, triangles, 1.0_dp, 0.6_dp, 4, 1.0e-5_dp, 2, &
        vector_reference, status)
    call record_condition(maxval(abs(vector_reference - vector_potential)) < &
        3.0e-13_dp*max(1.0_dp, maxval(abs(vector_reference))), &
        "combined assembly reproduces the standalone vector block")

    call generate_sphere_surface_mesh(2.0_dp, 0, scaled_vertices, triangles)
    call assemble_maxwell_sphere_curved_potential_operators_rwg_3d( &
        scaled_vertices, triangles, 2.0_dp, 0.3_dp, 4, 1.0e-5_dp, 2, &
        scaled_vector, scaled_scalar, status)
    error = maxval(abs(scaled_vector - 8.0_dp*vector_potential))/ &
        maxval(abs(scaled_vector))
    call record_condition(status == 0 .and. error < 4.0e-11_dp, &
        "curved vector-potential block has analytical R-cubed scaling")
    error = maxval(abs(scaled_scalar - 2.0_dp*scalar_potential))/ &
        maxval(abs(scaled_scalar))
    call record_condition(error < 4.0e-11_dp, &
        "curved scalar-potential block has analytical linear-R scaling")

    call check_summary("Curved-sphere RWG potential operators")
    if (.not. all_passed) error stop 1

contains

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_maxwell_sphere_curved_potential_operators
