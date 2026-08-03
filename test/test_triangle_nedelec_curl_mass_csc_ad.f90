program test_triangle_nedelec_curl_mass_csc_ad
    use check, only: check_condition, check_summary
    use fortfem_feec, only: assemble_triangle_nedelec_curl_mass_csc, &
        assemble_triangle_nedelec_curl_mass_csc_jvp, &
        assemble_triangle_nedelec_curl_mass_csc_vjp
    use fortfem_kinds, only: dp
    use fortfem_mesh_2d, only: mesh_2d_t
    use fortsparse, only: csc_t, fortsparse_status_t
    implicit none

    integer, parameter :: order = 2, quadrature_degree = 5
    real(dp), parameter :: step = 2.0e-7_dp
    type(mesh_2d_t) :: mesh, minus_mesh, plus_mesh
    real(dp) :: vertices_dot(2, 4), vertices_bar(2, 4)
    real(dp), allocatable :: matrix_values_bar(:)
    type(csc_t) :: matrix, matrix_dot, minus, plus
    type(fortsparse_status_t) :: status
    real(dp) :: curl_coefficient, curl_coefficient_bar, curl_coefficient_dot
    real(dp) :: mass_tensor(2, 2), mass_tensor_bar(2, 2)
    real(dp) :: mass_tensor_dot(2, 2)
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
    curl_coefficient = 1.7_dp
    curl_coefficient_dot = 0.13_dp
    mass_tensor = reshape([1.2_dp, 0.15_dp, -0.25_dp, 0.8_dp], [2, 2])
    mass_tensor_dot = reshape([ &
        -0.09_dp, 0.04_dp, 0.03_dp, 0.07_dp], [2, 2])

    call assemble_triangle_nedelec_curl_mass_csc( &
        mesh, order, quadrature_degree, matrix, status, curl_coefficient, &
        mass_tensor=mass_tensor)
    call assemble_triangle_nedelec_curl_mass_csc_jvp( &
        mesh, order, quadrature_degree, curl_coefficient, mass_tensor, &
        vertices_dot, curl_coefficient_dot, mass_tensor_dot, matrix_dot, &
        status)
    call assemble_triangle_nedelec_curl_mass_csc( &
        plus_mesh, order, quadrature_degree, plus, status, &
        curl_coefficient + step*curl_coefficient_dot, &
        mass_tensor=mass_tensor + step*mass_tensor_dot)
    call assemble_triangle_nedelec_curl_mass_csc( &
        minus_mesh, order, quadrature_degree, minus, status, &
        curl_coefficient - step*curl_coefficient_dot, &
        mass_tensor=mass_tensor - step*mass_tensor_dot)
    call check_condition( &
        matrix_dot%nnz == matrix%nnz .and. plus%nnz == matrix%nnz .and. &
        minus%nnz == matrix%nnz .and. &
        all(matrix_dot%col_ptr == matrix%col_ptr) .and. &
        all(matrix_dot%row_idx == matrix%row_idx), &
        "global triangle Nedelec JVP preserves the merged CSC pattern")
    relative_error = maxval(abs( &
        matrix_dot%val - (plus%val - minus%val)/(2.0_dp*step)))/ &
        max(1.0_dp, maxval(abs(matrix_dot%val)))
    call check_condition( &
        relative_error < 5.0e-8_dp, &
        "global triangle Nedelec JVP matches independent mesh reassembly")

    allocate(matrix_values_bar(matrix%nnz))
    do entry = 1, matrix%nnz
        matrix_values_bar(entry) = 0.003_dp*entry - 0.00001_dp*entry**2
    end do
    call assemble_triangle_nedelec_curl_mass_csc_vjp( &
        mesh, order, quadrature_degree, curl_coefficient, mass_tensor, &
        matrix_values_bar, vertices_bar, curl_coefficient_bar, &
        mass_tensor_bar, status)
    lhs = sum(matrix_values_bar*matrix_dot%val)
    rhs = sum(vertices_bar*vertices_dot) + &
        curl_coefficient_bar*curl_coefficient_dot + &
        sum(mass_tensor_bar*mass_tensor_dot)
    call check_condition( &
        status%code == 0, "global triangle Nedelec VJP succeeds")
    call check_condition( &
        abs(lhs - rhs) < 3.0e-10_dp*max(1.0_dp, abs(lhs), abs(rhs)), &
        "global triangle Nedelec products accumulate shared-vertex adjoints")

    call check_summary("Global triangle Nedelec curl-mass derivatives")
end program test_triangle_nedelec_curl_mass_csc_ad
