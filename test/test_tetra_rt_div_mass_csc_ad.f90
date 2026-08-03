program test_tetra_rt_div_mass_csc_ad
    use check, only: check_condition, check_summary
    use fortfem_feec, only: assemble_tetra_rt_div_mass_csc, &
        assemble_tetra_rt_div_mass_csc_jvp, assemble_tetra_rt_div_mass_csc_vjp
    use fortfem_kinds, only: dp
    use fortsparse, only: csc_t, fortsparse_status_t
    implicit none

    integer, parameter :: degree = 1, quadrature_degree = 5
    integer, parameter :: tetrahedra(4, 2) = reshape([ &
        1, 2, 3, 4, 2, 3, 4, 5], [4, 2])
    real(dp), parameter :: step = 2.0e-7_dp
    real(dp) :: vertices(3, 5), vertices_dot(3, 5), vertices_bar(3, 5)
    real(dp), allocatable :: matrix_values_bar(:)
    type(csc_t) :: matrix, matrix_dot, plus, minus
    type(fortsparse_status_t) :: sparse_status
    real(dp) :: divergence_coefficient, mass_coefficient
    real(dp) :: divergence_coefficient_dot, mass_coefficient_dot
    real(dp) :: divergence_coefficient_bar, mass_coefficient_bar
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
    divergence_coefficient = 1.7_dp
    mass_coefficient = -0.8_dp
    divergence_coefficient_dot = 0.13_dp
    mass_coefficient_dot = -0.09_dp

    call assemble_tetra_rt_div_mass_csc( &
        vertices, tetrahedra, degree, quadrature_degree, matrix, &
        sparse_status, divergence_coefficient, mass_coefficient)
    call assemble_tetra_rt_div_mass_csc_jvp( &
        vertices, tetrahedra, degree, quadrature_degree, &
        divergence_coefficient, mass_coefficient, vertices_dot, &
        divergence_coefficient_dot, mass_coefficient_dot, matrix_dot, &
        sparse_status)
    call assemble_tetra_rt_div_mass_csc( &
        vertices + step*vertices_dot, tetrahedra, degree, quadrature_degree, &
        plus, sparse_status, &
        divergence_coefficient + step*divergence_coefficient_dot, &
        mass_coefficient + step*mass_coefficient_dot)
    call assemble_tetra_rt_div_mass_csc( &
        vertices - step*vertices_dot, tetrahedra, degree, quadrature_degree, &
        minus, sparse_status, &
        divergence_coefficient - step*divergence_coefficient_dot, &
        mass_coefficient - step*mass_coefficient_dot)
    call check_condition( &
        matrix_dot%nnz == matrix%nnz .and. plus%nnz == matrix%nnz .and. &
        minus%nnz == matrix%nnz .and. &
        all(matrix_dot%col_ptr == matrix%col_ptr) .and. &
        all(matrix_dot%row_idx == matrix%row_idx), &
        "Global RT JVP preserves the merged CSC pattern")
    relative_error = maxval(abs( &
        matrix_dot%val - (plus%val - minus%val)/(2.0_dp*step)))/ &
        max(1.0_dp, maxval(abs(matrix_dot%val)))
    call check_condition( &
        relative_error < 5.0e-8_dp, &
        "Global RT JVP matches independent complete mesh reassembly")

    allocate(matrix_values_bar(matrix%nnz))
    do entry = 1, matrix%nnz
        matrix_values_bar(entry) = 0.003_dp*entry - 0.00001_dp*entry**2
    end do
    call assemble_tetra_rt_div_mass_csc_vjp( &
        vertices, tetrahedra, degree, quadrature_degree, &
        divergence_coefficient, mass_coefficient, matrix_values_bar, &
        vertices_bar, divergence_coefficient_bar, mass_coefficient_bar, &
        sparse_status)
    lhs = sum(matrix_values_bar*matrix_dot%val)
    rhs = sum(vertices_bar*vertices_dot) + &
        divergence_coefficient_bar*divergence_coefficient_dot + &
        mass_coefficient_bar*mass_coefficient_dot
    call check_condition( &
        sparse_status%code == 0, "Global RT CSC VJP succeeds")
    call check_condition( &
        abs(lhs - rhs) < 3.0e-10_dp*max(1.0_dp, abs(lhs), abs(rhs)), &
        "Global RT products accumulate the shared-vertex adjoint")

    call check_summary("Global tetrahedral RT div-mass derivatives")
end program test_tetra_rt_div_mass_csc_ad
