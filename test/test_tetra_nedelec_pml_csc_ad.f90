program test_tetra_nedelec_pml_csc_ad
    use check, only: check_condition, check_summary
    use fortfem_api, only: assemble_tetra_nedelec_pml_csc, &
        assemble_tetra_nedelec_pml_csc_jvp, &
        assemble_tetra_nedelec_pml_csc_vjp
    use fortfem_kinds, only: dp
    use fortsparse, only: csc_z_t, fortsparse_status_t
    implicit none

    integer, parameter :: order = 1
    integer, parameter :: tetrahedra(4, 2) = reshape([ &
        1, 2, 3, 4, 2, 3, 4, 5], [4, 2])
    real(dp), parameter :: step = 2.0e-7_dp
    real(dp) :: vertices(3, 5), vertices_dot(3, 5), vertices_bar(3, 5)
    complex(dp) :: stretch(3, 2), stretch_dot(3, 2), stretch_bar(3, 2)
    complex(dp), allocatable :: matrix_values_bar(:)
    type(csc_z_t) :: matrix, matrix_dot, plus, minus
    type(fortsparse_status_t) :: sparse_status
    real(dp) :: wave_number, wave_number_dot, wave_number_bar
    real(dp) :: lhs, rhs, relative_error
    integer :: entry

    vertices = reshape([ &
        0.0_dp, 0.0_dp, 0.0_dp, &
        1.0_dp, 0.0_dp, 0.0_dp, &
        0.0_dp, 1.0_dp, 0.0_dp, &
        0.0_dp, 0.0_dp, 1.0_dp, &
        1.0_dp, 1.0_dp, 1.0_dp], [3, 5])
    vertices_dot = reshape([ &
        0.01_dp, -0.02_dp, 0.015_dp, &
        -0.03_dp, 0.01_dp, 0.02_dp, &
        0.025_dp, -0.015_dp, 0.005_dp, &
        -0.01_dp, 0.03_dp, -0.02_dp, &
        0.02_dp, -0.01_dp, 0.025_dp], [3, 5])
    stretch(:, 1) = [ &
        cmplx(1.1_dp, 0.35_dp, dp), cmplx(0.95_dp, 0.12_dp, dp), &
        cmplx(1.05_dp, 0.27_dp, dp)]
    stretch(:, 2) = [ &
        cmplx(1.0_dp, 0.22_dp, dp), cmplx(1.08_dp, 0.18_dp, dp), &
        cmplx(0.97_dp, 0.31_dp, dp)]
    stretch_dot(:, 1) = [ &
        cmplx(0.04_dp, -0.03_dp, dp), cmplx(-0.02_dp, 0.01_dp, dp), &
        cmplx(0.03_dp, 0.02_dp, dp)]
    stretch_dot(:, 2) = [ &
        cmplx(-0.01_dp, 0.02_dp, dp), cmplx(0.025_dp, -0.015_dp, dp), &
        cmplx(0.02_dp, 0.01_dp, dp)]
    wave_number = 1.7_dp
    wave_number_dot = 0.13_dp

    call assemble_tetra_nedelec_pml_csc( &
        vertices, tetrahedra, order, stretch, wave_number, matrix, sparse_status)
    call assemble_tetra_nedelec_pml_csc_jvp( &
        vertices, tetrahedra, order, stretch, wave_number, vertices_dot, &
        stretch_dot, wave_number_dot, matrix_dot, sparse_status)
    call assemble_tetra_nedelec_pml_csc( &
        vertices + step*vertices_dot, tetrahedra, order, &
        stretch + step*stretch_dot, wave_number + step*wave_number_dot, plus, &
        sparse_status)
    call assemble_tetra_nedelec_pml_csc( &
        vertices - step*vertices_dot, tetrahedra, order, &
        stretch - step*stretch_dot, wave_number - step*wave_number_dot, minus, &
        sparse_status)
    call check_condition( &
        matrix_dot%nnz == matrix%nnz .and. plus%nnz == matrix%nnz .and. &
        minus%nnz == matrix%nnz .and. &
        all(matrix_dot%col_ptr == matrix%col_ptr) .and. &
        all(matrix_dot%row_idx == matrix%row_idx) .and. &
        all(plus%col_ptr == matrix%col_ptr) .and. &
        all(minus%col_ptr == matrix%col_ptr), &
        "Global PML JVP preserves the merged CSC pattern")
    relative_error = maxval(abs( &
        matrix_dot%val - (plus%val - minus%val)/(2.0_dp*step)))/ &
        max(1.0_dp, maxval(abs(matrix_dot%val)))
    call check_condition( &
        relative_error < 5.0e-8_dp, &
        "Global PML JVP matches complete merged-mesh reassembly")

    allocate(matrix_values_bar(matrix%nnz))
    do entry = 1, matrix%nnz
        matrix_values_bar(entry) = cmplx( &
            0.003_dp*entry, -0.002_dp*entry, dp)
    end do
    call assemble_tetra_nedelec_pml_csc_vjp( &
        vertices, tetrahedra, order, stretch, wave_number, matrix_values_bar, &
        vertices_bar, stretch_bar, wave_number_bar, sparse_status)
    lhs = real(sum(conjg(matrix_values_bar)*matrix_dot%val), dp)
    rhs = sum(vertices_bar*vertices_dot) + &
        real(sum(conjg(stretch_bar)*stretch_dot), dp) + &
        wave_number_bar*wave_number_dot
    call check_condition( &
        sparse_status%code == 0, "Global PML CSC VJP succeeds")
    call check_condition( &
        abs(lhs - rhs) < 3.0e-10_dp*max(1.0_dp, abs(lhs), abs(rhs)), &
        "Global PML products accumulate the shared-vertex adjoint")

    call check_summary("Global tetrahedral Nedelec PML derivatives")
end program test_tetra_nedelec_pml_csc_ad
