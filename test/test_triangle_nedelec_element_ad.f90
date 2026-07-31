program test_triangle_nedelec_element_ad
    use check, only: check_condition, check_summary
    use fortfem_api, only: assemble_triangle_nedelec_curl_mass_element, &
        assemble_triangle_nedelec_curl_mass_element_jvp, &
        assemble_triangle_nedelec_curl_mass_element_vjp
    use fortfem_kinds, only: dp
    implicit none

    integer, parameter :: order = 2, quadrature_degree = 6
    integer, parameter :: dof_count = order*(order + 2)
    real(dp), parameter :: step = 2.0e-7_dp
    real(dp) :: vertices(2, 3), vertices_dot(2, 3), vertices_bar(2, 3)
    real(dp), allocatable :: matrix(:, :), matrix_dot(:, :)
    real(dp), allocatable :: minus(:, :), plus(:, :)
    real(dp) :: matrix_bar(dof_count, dof_count)
    real(dp) :: curl_coefficient, curl_coefficient_bar, curl_coefficient_dot
    real(dp) :: mass_tensor(2, 2), mass_tensor_bar(2, 2)
    real(dp) :: mass_tensor_dot(2, 2)
    real(dp) :: lhs, relative_error, rhs
    integer :: column, row, status, status_minus, status_plus

    vertices = reshape([ &
        0.1_dp, -0.2_dp, 1.3_dp, -0.05_dp, 0.2_dp, 1.1_dp], [2, 3])
    vertices_dot = reshape([ &
        0.01_dp, -0.02_dp, -0.03_dp, 0.015_dp, 0.025_dp, -0.01_dp], [2, 3])
    curl_coefficient = 1.7_dp
    curl_coefficient_dot = 0.13_dp
    mass_tensor = reshape([1.2_dp, 0.15_dp, -0.25_dp, 0.8_dp], [2, 2])
    mass_tensor_dot = reshape([ &
        -0.09_dp, 0.04_dp, 0.03_dp, 0.07_dp], [2, 2])
    do column = 1, dof_count
        do row = 1, dof_count
            matrix_bar(row, column) = &
                0.002_dp*row - 0.0015_dp*column
        end do
    end do

    call assemble_triangle_nedelec_curl_mass_element( &
        vertices, order, quadrature_degree, matrix, status, &
        curl_coefficient, mass_tensor=mass_tensor)
    call assemble_triangle_nedelec_curl_mass_element_jvp( &
        vertices, order, quadrature_degree, curl_coefficient, mass_tensor, &
        vertices_dot, curl_coefficient_dot, mass_tensor_dot, matrix_dot, &
        status)
    call assemble_triangle_nedelec_curl_mass_element( &
        vertices + step*vertices_dot, order, quadrature_degree, plus, &
        status_plus, curl_coefficient + step*curl_coefficient_dot, &
        mass_tensor=mass_tensor + step*mass_tensor_dot)
    call assemble_triangle_nedelec_curl_mass_element( &
        vertices - step*vertices_dot, order, quadrature_degree, minus, &
        status_minus, curl_coefficient - step*curl_coefficient_dot, &
        mass_tensor=mass_tensor - step*mass_tensor_dot)
    relative_error = maxval(abs( &
        matrix_dot - (plus - minus)/(2.0_dp*step)))/ &
        max(1.0_dp, maxval(abs(matrix_dot)))
    call check_condition( &
        status == 0 .and. status_plus == 0 .and. status_minus == 0, &
        "triangle Nedelec element JVP accepts a valid direction")
    call check_condition( &
        relative_error < 5.0e-8_dp, &
        "triangle Nedelec element JVP matches independent full reassembly")

    call assemble_triangle_nedelec_curl_mass_element_vjp( &
        vertices, order, quadrature_degree, curl_coefficient, mass_tensor, &
        matrix_bar, vertices_bar, curl_coefficient_bar, mass_tensor_bar, &
        status)
    lhs = sum(matrix_bar*matrix_dot)
    rhs = sum(vertices_bar*vertices_dot) + &
        curl_coefficient_bar*curl_coefficient_dot + &
        sum(mass_tensor_bar*mass_tensor_dot)
    call check_condition(status == 0, "triangle Nedelec element VJP succeeds")
    call check_condition( &
        abs(lhs - rhs) < 2.0e-10_dp*max(1.0_dp, abs(lhs), abs(rhs)), &
        "triangle Nedelec element products obey the adjoint identity")

    call check_summary("Triangle Nedelec element derivatives")
end program test_triangle_nedelec_element_ad
