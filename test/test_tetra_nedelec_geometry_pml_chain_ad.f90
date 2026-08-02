program test_tetra_nedelec_geometry_pml_chain_ad
    use check, only: check_condition, check_summary
    use fortfem_api, only: &
        assemble_tetra_nedelec_geometry_pml_csc, &
        assemble_tetra_nedelec_geometry_pml_csc_jvp, &
        assemble_tetra_nedelec_geometry_pml_csc_vjp
    use fortfem_kinds, only: dp
    use fortsparse, only: csc_z_t, fortsparse_status_t
    implicit none

    integer, parameter :: order = 1, cell_count = 1
    integer, parameter :: tetrahedra(4, cell_count) = reshape([1, 2, 3, 4], &
        [4, cell_count])
    real(dp), parameter :: step = 2.0e-7_dp
    real(dp) :: vertices(3, 4), vertices_dot(3, 4), vertices_bar(3, 4)
    real(dp) :: origins(3, cell_count), origins_dot(3, cell_count)
    real(dp) :: origins_bar(3, cell_count)
    real(dp) :: normals(3, cell_count), normals_dot(3, cell_count)
    real(dp) :: normals_bar(3, cell_count)
    real(dp) :: width(cell_count), width_dot(cell_count), width_bar(cell_count)
    real(dp) :: wave_number, wave_number_dot, wave_number_bar
    real(dp) :: sigma_max, sigma_max_dot, sigma_max_bar
    complex(dp), allocatable :: matrix_values_bar(:)
    type(csc_z_t) :: matrix, matrix_dot, plus, minus
    type(fortsparse_status_t) :: status
    real(dp) :: lhs, rhs, relative_error
    integer :: entry
    logical :: all_passed

    all_passed = .true.
    vertices = reshape([ &
        0.0_dp, 0.0_dp, 0.0_dp, &
        1.0_dp, 0.0_dp, 0.0_dp, &
        0.0_dp, 1.0_dp, 0.0_dp, &
        0.0_dp, 0.0_dp, 1.0_dp], [3, 4])
    vertices_dot = reshape([ &
        0.01_dp, -0.02_dp, 0.015_dp, &
        -0.03_dp, 0.01_dp, 0.02_dp, &
        0.025_dp, -0.015_dp, 0.005_dp, &
        -0.01_dp, 0.03_dp, -0.02_dp], [3, 4])
    origins = reshape([0.0_dp, 0.0_dp, 0.0_dp], [3, cell_count])
    origins_dot = reshape([0.01_dp, -0.02_dp, 0.03_dp], [3, cell_count])
    normals = reshape([0.0_dp, 0.0_dp, 1.0_dp], [3, cell_count])
    normals_dot = reshape([0.08_dp, -0.04_dp, 0.0_dp], [3, cell_count])
    width = [1.2_dp]
    width_dot = [0.07_dp]
    wave_number = 1.7_dp
    wave_number_dot = 0.13_dp
    sigma_max = 0.8_dp
    sigma_max_dot = -0.06_dp

    call assemble_tetra_nedelec_geometry_pml_csc( &
        vertices, tetrahedra, order, origins, normals, width, wave_number, &
        sigma_max, 2, matrix, status)
    call record_condition(status%code == 0 .and. matrix%nnz > 0, &
        "geometry-generated curvilinear PML assembles")
    call assemble_tetra_nedelec_geometry_pml_csc_jvp( &
        vertices, tetrahedra, order, origins, normals, width, wave_number, &
        sigma_max, 2, vertices_dot, origins_dot, normals_dot, width_dot, &
        wave_number_dot, sigma_max_dot, matrix_dot, status)
    call record_condition(status%code == 0 .and. matrix_dot%nnz == matrix%nnz .and. &
        all(matrix_dot%col_ptr == matrix%col_ptr) .and. &
        all(matrix_dot%row_idx == matrix%row_idx), &
        "geometry-generated PML JVP preserves the CSC pattern")

    call assemble_tetra_nedelec_geometry_pml_csc( &
        vertices + step*vertices_dot, tetrahedra, order, &
        origins + step*origins_dot, normals + step*normals_dot, &
        width + step*width_dot, wave_number + step*wave_number_dot, &
        sigma_max + step*sigma_max_dot, 2, plus, status)
    call assemble_tetra_nedelec_geometry_pml_csc( &
        vertices - step*vertices_dot, tetrahedra, order, &
        origins - step*origins_dot, normals - step*normals_dot, &
        width - step*width_dot, wave_number - step*wave_number_dot, &
        sigma_max - step*sigma_max_dot, 2, minus, status)
    relative_error = maxval(abs(matrix_dot%val - &
        (plus%val - minus%val)/(2.0_dp*step)))/ &
        max(1.0_dp, maxval(abs(matrix_dot%val)))
    call record_condition(status%code == 0 .and. plus%nnz == matrix%nnz .and. &
        minus%nnz == matrix%nnz .and. relative_error < 8.0e-8_dp, &
        "geometry-generated PML JVP matches reassembly differences")

    allocate(matrix_values_bar(matrix%nnz))
    do entry = 1, matrix%nnz
        matrix_values_bar(entry) = cmplx( &
            0.003_dp*entry, -0.002_dp*entry, dp)
    end do
    call assemble_tetra_nedelec_geometry_pml_csc_vjp( &
        vertices, tetrahedra, order, origins, normals, width, wave_number, &
        sigma_max, 2, matrix_values_bar, vertices_bar, origins_bar, &
        normals_bar, width_bar, wave_number_bar, sigma_max_bar, status)
    lhs = real(sum(conjg(matrix_values_bar)*matrix_dot%val), dp)
    rhs = sum(vertices_bar*vertices_dot) + sum(origins_bar*origins_dot) + &
        sum(normals_bar*normals_dot) + dot_product(width_bar, width_dot) + &
        wave_number_bar*wave_number_dot + sigma_max_bar*sigma_max_dot
    call record_condition(status%code == 0, &
        "geometry-generated curvilinear PML VJP assembles")
    call record_condition(abs(lhs - rhs) < 2.0e-9_dp*max(1.0_dp, abs(lhs), abs(rhs)), &
        "geometry-generated PML VJP satisfies the complete real adjoint identity")

    call check_summary("geometry-generated curvilinear tetrahedral PML")
    if (.not. all_passed) error stop 1

contains

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_tetra_nedelec_geometry_pml_chain_ad
