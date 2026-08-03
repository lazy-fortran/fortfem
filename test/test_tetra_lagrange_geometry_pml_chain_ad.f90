program test_tetra_lagrange_geometry_pml_chain_ad
    use check, only: check_condition, check_summary
    use fortfem_feec, only: &
        assemble_tetra_lagrange_geometry_pml_csc, &
        assemble_tetra_lagrange_geometry_pml_csc_jvp, &
        assemble_tetra_lagrange_geometry_pml_csc_vjp
    use fortfem_kinds, only: dp
    use fortsparse, only: csc_z_t, fortsparse_status_t
    implicit none

    integer, parameter :: degree = 1, polynomial_degree = 2
    integer, parameter :: tetrahedra(4, 1) = reshape([1, 2, 3, 4], [4, 1])
    real(dp), parameter :: step = 2.0e-7_dp
    real(dp) :: vertices(3, 4), vertices_dot(3, 4), vertices_bar(3, 4)
    real(dp) :: origins(3, 1), origins_dot(3, 1), origins_bar(3, 1)
    real(dp) :: normals(3, 1), normals_dot(3, 1), normals_bar(3, 1)
    real(dp) :: width(1), width_dot(1), width_bar(1)
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
    origins = reshape([0.0_dp, 0.0_dp, 0.0_dp], [3, 1])
    origins_dot = reshape([0.01_dp, -0.02_dp, 0.03_dp], [3, 1])
    normals = reshape([0.0_dp, 0.0_dp, 1.0_dp], [3, 1])
    normals_dot = reshape([0.08_dp, -0.04_dp, 0.0_dp], [3, 1])
    width = [1.2_dp]
    width_dot = [0.07_dp]
    wave_number = 1.7_dp
    wave_number_dot = 0.13_dp
    sigma_max = 0.8_dp
    sigma_max_dot = -0.06_dp

    call assemble_tetra_lagrange_geometry_pml_csc( &
        vertices, tetrahedra, degree, origins, normals, width, wave_number, &
        sigma_max, polynomial_degree, matrix, status)
    call record_condition(status%code == 0 .and. matrix%nnz > 0, &
        "geometry-generated scalar PML assembles")
    call assemble_tetra_lagrange_geometry_pml_csc_jvp( &
        vertices, tetrahedra, degree, origins, normals, width, wave_number, &
        sigma_max, polynomial_degree, vertices_dot, origins_dot, normals_dot, &
        width_dot, wave_number_dot, sigma_max_dot, matrix_dot, status)
    call record_condition(status%code == 0 .and. matrix_dot%nnz == matrix%nnz .and. &
        all(matrix_dot%col_ptr == matrix%col_ptr) .and. &
        all(matrix_dot%row_idx == matrix%row_idx), &
        "geometry-generated scalar PML JVP preserves the CSC pattern")

    call assemble_tetra_lagrange_geometry_pml_csc( &
        vertices + step*vertices_dot, tetrahedra, degree, &
        origins + step*origins_dot, normals + step*normals_dot, &
        width + step*width_dot, wave_number + step*wave_number_dot, &
        sigma_max + step*sigma_max_dot, polynomial_degree, plus, status)
    call assemble_tetra_lagrange_geometry_pml_csc( &
        vertices - step*vertices_dot, tetrahedra, degree, &
        origins - step*origins_dot, normals - step*normals_dot, &
        width - step*width_dot, wave_number - step*wave_number_dot, &
        sigma_max - step*sigma_max_dot, polynomial_degree, minus, status)
    relative_error = maxval(abs(matrix_dot%val - &
        (plus%val - minus%val)/(2.0_dp*step)))/ &
        max(1.0_dp, maxval(abs(matrix_dot%val)))
    call record_condition(status%code == 0 .and. plus%nnz == matrix%nnz .and. &
        minus%nnz == matrix%nnz .and. relative_error < 8.0e-8_dp, &
        "geometry-generated scalar PML JVP matches reassembly differences")

    allocate(matrix_values_bar(matrix%nnz))
    do entry = 1, matrix%nnz
        matrix_values_bar(entry) = cmplx( &
            0.003_dp*entry, -0.002_dp*entry, dp)
    end do
    call assemble_tetra_lagrange_geometry_pml_csc_vjp( &
        vertices, tetrahedra, degree, origins, normals, width, wave_number, &
        sigma_max, polynomial_degree, matrix_values_bar, vertices_bar, origins_bar, &
        normals_bar, width_bar, wave_number_bar, sigma_max_bar, status)
    lhs = real(sum(conjg(matrix_values_bar)*matrix_dot%val), dp)
    rhs = sum(vertices_bar*vertices_dot) + sum(origins_bar*origins_dot) + &
        sum(normals_bar*normals_dot) + dot_product(width_bar, width_dot) + &
        wave_number_bar*wave_number_dot + sigma_max_bar*sigma_max_dot
    call record_condition(status%code == 0, &
        "geometry-generated scalar PML VJP assembles")
    call record_condition(abs(lhs - rhs) < 2.0e-9_dp*max(1.0_dp, abs(lhs), abs(rhs)), &
        "geometry-generated scalar PML VJP satisfies the complete adjoint identity")

    call check_summary("geometry-generated scalar tetrahedral PML")
    if (.not. all_passed) error stop 1

contains

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_tetra_lagrange_geometry_pml_chain_ad
