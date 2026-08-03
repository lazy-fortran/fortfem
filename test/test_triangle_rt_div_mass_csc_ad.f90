program test_triangle_rt_div_mass_csc_ad
    use check, only: check_condition, check_summary
    use fortfem_feec, only: assemble_triangle_rt_div_mass_csc, &
        assemble_triangle_rt_div_mass_csc_jvp, &
        assemble_triangle_rt_div_mass_csc_vjp
    use fortfem_kinds, only: dp
    use fortfem_mesh_2d, only: mesh_2d_t
    use fortsparse, only: csc_t, fortsparse_status_t
    implicit none

    integer, parameter :: degree = 1, quadrature_degree = 5
    real(dp), parameter :: step = 2.0e-7_dp
    type(mesh_2d_t) :: mesh, minus_mesh, plus_mesh
    real(dp) :: vertices_dot(2, 4), vertices_bar(2, 4)
    real(dp), allocatable :: matrix_values_bar(:)
    type(csc_t) :: matrix, matrix_dot, minus, plus
    type(fortsparse_status_t) :: status
    real(dp) :: divergence_coefficient, divergence_coefficient_bar
    real(dp) :: divergence_coefficient_dot
    real(dp) :: mass_coefficient, mass_coefficient_bar, mass_coefficient_dot
    real(dp) :: lhs, relative_error, rhs
    integer :: entry

    call mesh%create_rectangular(2, 2, 0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp)
    call plus_mesh%create_rectangular(2, 2, 0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp)
    call minus_mesh%create_rectangular(2, 2, 0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp)
    vertices_dot = reshape([ &
        0.01_dp, -0.02_dp, -0.03_dp, 0.015_dp, &
        0.025_dp, -0.01_dp, 0.02_dp, 0.03_dp], [2, 4])
    plus_mesh%vertices = plus_mesh%vertices + step*vertices_dot
    minus_mesh%vertices = minus_mesh%vertices - step*vertices_dot
    divergence_coefficient = 1.7_dp
    mass_coefficient = -0.8_dp
    divergence_coefficient_dot = 0.13_dp
    mass_coefficient_dot = -0.09_dp

    call assemble_triangle_rt_div_mass_csc( &
        mesh, degree, quadrature_degree, matrix, status, &
        divergence_coefficient, mass_coefficient)
    call assemble_triangle_rt_div_mass_csc_jvp( &
        mesh, degree, quadrature_degree, divergence_coefficient, &
        mass_coefficient, vertices_dot, divergence_coefficient_dot, &
        mass_coefficient_dot, matrix_dot, status)
    call assemble_triangle_rt_div_mass_csc( &
        plus_mesh, degree, quadrature_degree, plus, status, &
        divergence_coefficient + step*divergence_coefficient_dot, &
        mass_coefficient + step*mass_coefficient_dot)
    call assemble_triangle_rt_div_mass_csc( &
        minus_mesh, degree, quadrature_degree, minus, status, &
        divergence_coefficient - step*divergence_coefficient_dot, &
        mass_coefficient - step*mass_coefficient_dot)
    call check_condition( &
        matrix_dot%nnz == matrix%nnz .and. plus%nnz == matrix%nnz .and. &
        minus%nnz == matrix%nnz .and. &
        all(matrix_dot%col_ptr == matrix%col_ptr) .and. &
        all(matrix_dot%row_idx == matrix%row_idx), &
        "global triangle RT JVP preserves the merged CSC pattern")
    relative_error = maxval(abs( &
        matrix_dot%val - (plus%val - minus%val)/(2.0_dp*step)))/ &
        max(1.0_dp, maxval(abs(matrix_dot%val)))
    call check_condition( &
        relative_error < 5.0e-8_dp, &
        "global triangle RT JVP matches independent mesh reassembly")

    allocate(matrix_values_bar(matrix%nnz))
    do entry = 1, matrix%nnz
        matrix_values_bar(entry) = 0.003_dp*entry - 0.00001_dp*entry**2
    end do
    call assemble_triangle_rt_div_mass_csc_vjp( &
        mesh, degree, quadrature_degree, divergence_coefficient, &
        mass_coefficient, matrix_values_bar, vertices_bar, &
        divergence_coefficient_bar, mass_coefficient_bar, status)
    lhs = sum(matrix_values_bar*matrix_dot%val)
    rhs = sum(vertices_bar*vertices_dot) + &
        divergence_coefficient_bar*divergence_coefficient_dot + &
        mass_coefficient_bar*mass_coefficient_dot
    call check_condition(status%code == 0, "global triangle RT VJP succeeds")
    call check_condition( &
        abs(lhs - rhs) < 3.0e-10_dp*max(1.0_dp, abs(lhs), abs(rhs)), &
        "global triangle RT products accumulate shared-vertex adjoints")

    call check_summary("Global triangle RT div-mass derivatives")
end program test_triangle_rt_div_mass_csc_ad
