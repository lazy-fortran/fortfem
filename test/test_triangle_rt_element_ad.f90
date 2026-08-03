program test_triangle_rt_element_ad
    use check, only: check_condition, check_summary
    use fortfem_feec, only: assemble_triangle_rt_div_mass_element, &
        assemble_triangle_rt_div_mass_element_jvp, &
        assemble_triangle_rt_div_mass_element_vjp
    use fortfem_kinds, only: dp
    implicit none

    integer, parameter :: degree = 2, quadrature_degree = 6
    integer, parameter :: dof_count = (degree + 1)*(degree + 3)
    real(dp), parameter :: step = 2.0e-7_dp
    real(dp) :: vertices(2, 3), vertices_dot(2, 3), vertices_bar(2, 3)
    real(dp), allocatable :: matrix(:, :), matrix_dot(:, :)
    real(dp), allocatable :: minus(:, :), plus(:, :)
    real(dp) :: matrix_bar(dof_count, dof_count)
    real(dp) :: divergence_coefficient, divergence_coefficient_bar
    real(dp) :: divergence_coefficient_dot
    real(dp) :: mass_coefficient, mass_coefficient_bar, mass_coefficient_dot
    real(dp) :: lhs, relative_error, rhs
    integer :: column, row, status, status_minus, status_plus

    vertices = reshape([ &
        0.1_dp, -0.2_dp, 1.3_dp, -0.05_dp, 0.2_dp, 1.1_dp], [2, 3])
    vertices_dot = reshape([ &
        0.01_dp, -0.02_dp, -0.03_dp, 0.015_dp, 0.025_dp, -0.01_dp], [2, 3])
    divergence_coefficient = 1.7_dp
    mass_coefficient = -0.8_dp
    divergence_coefficient_dot = 0.13_dp
    mass_coefficient_dot = -0.09_dp
    do column = 1, dof_count
        do row = 1, dof_count
            matrix_bar(row, column) = &
                0.002_dp*row - 0.0015_dp*column
        end do
    end do

    call assemble_triangle_rt_div_mass_element( &
        vertices, degree, quadrature_degree, matrix, status, &
        divergence_coefficient, mass_coefficient)
    call assemble_triangle_rt_div_mass_element_jvp( &
        vertices, degree, quadrature_degree, divergence_coefficient, &
        mass_coefficient, vertices_dot, divergence_coefficient_dot, &
        mass_coefficient_dot, matrix_dot, status)
    call assemble_triangle_rt_div_mass_element( &
        vertices + step*vertices_dot, degree, quadrature_degree, plus, &
        status_plus, divergence_coefficient + &
        step*divergence_coefficient_dot, mass_coefficient + &
        step*mass_coefficient_dot)
    call assemble_triangle_rt_div_mass_element( &
        vertices - step*vertices_dot, degree, quadrature_degree, minus, &
        status_minus, divergence_coefficient - &
        step*divergence_coefficient_dot, mass_coefficient - &
        step*mass_coefficient_dot)
    relative_error = maxval(abs( &
        matrix_dot - (plus - minus)/(2.0_dp*step)))/ &
        max(1.0_dp, maxval(abs(matrix_dot)))
    call check_condition( &
        status == 0 .and. status_plus == 0 .and. status_minus == 0, &
        "triangle RT element JVP accepts a valid geometry direction")
    call check_condition( &
        relative_error < 5.0e-8_dp, &
        "triangle RT element JVP matches independent full reassembly")

    call assemble_triangle_rt_div_mass_element_vjp( &
        vertices, degree, quadrature_degree, divergence_coefficient, &
        mass_coefficient, matrix_bar, vertices_bar, &
        divergence_coefficient_bar, mass_coefficient_bar, status)
    lhs = sum(matrix_bar*matrix_dot)
    rhs = sum(vertices_bar*vertices_dot) + &
        divergence_coefficient_bar*divergence_coefficient_dot + &
        mass_coefficient_bar*mass_coefficient_dot
    call check_condition(status == 0, "triangle RT element VJP succeeds")
    call check_condition( &
        abs(lhs - rhs) < 2.0e-10_dp*max(1.0_dp, abs(lhs), abs(rhs)), &
        "triangle RT element products obey the adjoint identity")

    call check_summary("Triangle RT element derivatives")
end program test_triangle_rt_element_ad
