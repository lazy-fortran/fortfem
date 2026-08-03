program test_maxwell_sphere_curved_efie_bc_imaginary
    use check, only: check_condition, check_summary
    use fortfem_boundary, only: assemble_maxwell_sphere_curved_efie_bc_imaginary_3d
    use fortfem_core, only: generate_sphere_surface_mesh
    use fortfem_kinds, only: dp
    implicit none

    complex(dp), allocatable :: matrix(:, :), scaled(:, :)
    integer, allocatable :: triangles(:, :)
    real(dp), allocatable :: scaled_vertices(:, :), vertices(:, :)
    real(dp) :: error, symmetry_error
    integer :: status
    logical :: all_passed

    all_passed = .true.
    call generate_sphere_surface_mesh(1.0_dp, 0, vertices, triangles)
    call assemble_maxwell_sphere_curved_efie_bc_imaginary_3d( &
        vertices, triangles, 1.0_dp, 0.7_dp, 1.4_dp, 2, 2.0e-4_dp, 1, matrix, &
        status)
    symmetry_error = maxval(abs(matrix - transpose(matrix)))
    call record_condition(status == 0 .and. maxval(abs(matrix)) > 0.0_dp .and. &
        maxval(abs(aimag(matrix))) < 3.0e-13_dp .and. &
        symmetry_error < 3.0e-11_dp, &
        "curved BC imaginary-wave regularizer is real and symmetric")

    call generate_sphere_surface_mesh(2.0_dp, 0, scaled_vertices, triangles)
    call assemble_maxwell_sphere_curved_efie_bc_imaginary_3d( &
        scaled_vertices, triangles, 2.0_dp, 0.35_dp, 1.4_dp, 2, 2.0e-4_dp, 1, &
        scaled, status)
    error = maxval(abs(scaled - matrix))/maxval(abs(scaled))
    call record_condition(status == 0 .and. error < 2.0e-10_dp, &
        "curved BC regularizer is invariant under electromagnetic similarity")

    call check_summary("Curved-sphere BC imaginary-wave EFIE")
    if (.not. all_passed) error stop 1

contains

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_maxwell_sphere_curved_efie_bc_imaginary
