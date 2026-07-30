program test_tetra_nedelec_curl_mass_csc_ad
    use check, only: check_condition, check_summary
    use fortfem_api, only: assemble_tetra_nedelec_curl_mass_csc, &
        assemble_tetra_nedelec_curl_mass_csc_jvp, &
        assemble_tetra_nedelec_curl_mass_csc_vjp
    use fortfem_kinds, only: dp
    use fortsparse, only: csc_t, fortsparse_status_t
    implicit none

    integer, parameter :: order = 2
    integer, parameter :: tetrahedra(4, 2) = reshape([ &
        1, 2, 3, 4, 2, 3, 4, 5], [4, 2])
    real(dp), parameter :: step = 2.0e-7_dp
    real(dp) :: vertices(3, 5), vertices_dot(3, 5), vertices_bar(3, 5)
    real(dp), allocatable :: matrix_values_bar(:)
    type(csc_t) :: matrix, matrix_dot, plus, minus
    type(fortsparse_status_t) :: sparse_status
    real(dp) :: curl_coefficient, mass_coefficient
    real(dp) :: curl_coefficient_dot, mass_coefficient_dot
    real(dp) :: curl_coefficient_bar, mass_coefficient_bar
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
    curl_coefficient = 1.7_dp
    mass_coefficient = -0.8_dp
    curl_coefficient_dot = 0.13_dp
    mass_coefficient_dot = -0.09_dp

    call assemble_tetra_nedelec_curl_mass_csc( &
        vertices, tetrahedra, matrix, sparse_status, curl_coefficient, &
        mass_coefficient, order)
    call assemble_tetra_nedelec_curl_mass_csc_jvp( &
        vertices, tetrahedra, order, curl_coefficient, mass_coefficient, &
        vertices_dot, curl_coefficient_dot, mass_coefficient_dot, matrix_dot, &
        sparse_status)
    call assemble_tetra_nedelec_curl_mass_csc( &
        vertices + step*vertices_dot, tetrahedra, plus, sparse_status, &
        curl_coefficient + step*curl_coefficient_dot, &
        mass_coefficient + step*mass_coefficient_dot, order)
    call assemble_tetra_nedelec_curl_mass_csc( &
        vertices - step*vertices_dot, tetrahedra, minus, sparse_status, &
        curl_coefficient - step*curl_coefficient_dot, &
        mass_coefficient - step*mass_coefficient_dot, order)
    call check_condition( &
        matrix_dot%nnz == matrix%nnz .and. plus%nnz == matrix%nnz .and. &
        minus%nnz == matrix%nnz .and. &
        all(matrix_dot%col_ptr == matrix%col_ptr) .and. &
        all(matrix_dot%row_idx == matrix%row_idx), &
        "Global curl-mass JVP preserves the merged CSC pattern")
    relative_error = maxval(abs( &
        matrix_dot%val - (plus%val - minus%val)/(2.0_dp*step)))/ &
        max(1.0_dp, maxval(abs(matrix_dot%val)))
    call check_condition( &
        relative_error < 5.0e-8_dp, &
        "Global curl-mass JVP matches complete merged-mesh reassembly")

    allocate(matrix_values_bar(matrix%nnz))
    do entry = 1, matrix%nnz
        matrix_values_bar(entry) = 0.003_dp*entry - 0.00001_dp*entry**2
    end do
    call assemble_tetra_nedelec_curl_mass_csc_vjp( &
        vertices, tetrahedra, order, curl_coefficient, mass_coefficient, &
        matrix_values_bar, vertices_bar, curl_coefficient_bar, &
        mass_coefficient_bar, sparse_status)
    lhs = sum(matrix_values_bar*matrix_dot%val)
    rhs = sum(vertices_bar*vertices_dot) + &
        curl_coefficient_bar*curl_coefficient_dot + &
        mass_coefficient_bar*mass_coefficient_dot
    call check_condition( &
        sparse_status%code == 0, "Global curl-mass CSC VJP succeeds")
    call check_condition( &
        abs(lhs - rhs) < 3.0e-10_dp*max(1.0_dp, abs(lhs), abs(rhs)), &
        "Global curl-mass products accumulate the shared-vertex adjoint")

    call check_summary("Global tetrahedral Nedelec curl-mass derivatives")
end program test_tetra_nedelec_curl_mass_csc_ad
