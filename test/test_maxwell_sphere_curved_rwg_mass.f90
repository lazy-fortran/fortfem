program test_maxwell_sphere_curved_rwg_mass
    use check, only: check_condition, check_summary
    use fortfem_boundary, only: &
        assemble_maxwell_sphere_curved_rwg_mass_matrix
    use fortfem_core, only: generate_sphere_surface_mesh
    use fortfem_kinds, only: dp
    implicit none

    integer, allocatable :: triangles(:, :)
    real(dp), allocatable :: matrix(:, :), scaled_matrix(:, :)
    real(dp), allocatable :: vertices(:, :), vector(:)
    real(dp), allocatable :: scaled_vertices(:, :)
    real(dp) :: energy, relative_scaling_error, symmetry_error
    integer :: basis, status
    logical :: all_passed

    all_passed = .true.
    call generate_sphere_surface_mesh(1.0_dp, 0, vertices, triangles)
    call assemble_maxwell_sphere_curved_rwg_mass_matrix( &
        vertices, triangles, 1.0_dp, 12, matrix, status)
    call record_condition(status == 0, "curved sphere RWG mass assembles")
    allocate(vector(size(matrix, 1)))
    do basis = 1, size(vector)
        vector(basis) = sin(real(3*basis + 1, dp))
    end do
    symmetry_error = maxval(abs(matrix - transpose(matrix)))
    energy = dot_product(vector, matmul(matrix, vector))
    call record_condition(symmetry_error < 3.0e-14_dp, &
        "curved RWG mass is symmetric")
    call record_condition(energy > 0.0_dp .and. minval(diagonal(matrix)) > &
        0.0_dp, "curved RWG mass has positive field energy")

    call generate_sphere_surface_mesh(2.0_dp, 0, scaled_vertices, triangles)
    call assemble_maxwell_sphere_curved_rwg_mass_matrix( &
        scaled_vertices, triangles, 2.0_dp, 12, scaled_matrix, status)
    relative_scaling_error = maxval(abs(scaled_matrix - 4.0_dp*matrix))/ &
        maxval(abs(scaled_matrix))
    call record_condition(status == 0 .and. relative_scaling_error < &
        3.0e-14_dp, "curved RWG mass obeys the analytical area scaling law")

    call check_summary("Curved-sphere RWG mass matrix")
    if (.not. all_passed) error stop 1

contains

    pure function diagonal(matrix) result(values)
        real(dp), intent(in) :: matrix(:, :)
        real(dp) :: values(min(size(matrix, 1), size(matrix, 2)))
        integer :: entry

        do entry = 1, size(values)
            values(entry) = matrix(entry, entry)
        end do
    end function diagonal

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_maxwell_sphere_curved_rwg_mass
