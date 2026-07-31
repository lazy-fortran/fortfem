program test_triangle_full_vector_element_ad
    use check, only: check_condition, check_summary
    use fortfem_api, only: &
        assemble_triangle_bdm_div_mass_element, &
        assemble_triangle_bdm_div_mass_element_jvp, &
        assemble_triangle_bdm_div_mass_element_vjp, &
        assemble_triangle_nedelec_second_curl_mass_element, &
        assemble_triangle_nedelec_second_curl_mass_element_jvp, &
        assemble_triangle_nedelec_second_curl_mass_element_vjp
    use fortfem_kinds, only: dp
    implicit none

    integer, parameter :: degree = 2, quadrature_degree = 6
    integer, parameter :: dof_count = (degree + 1)*(degree + 2)
    real(dp), parameter :: step = 2.0e-7_dp
    real(dp) :: vertices(2, 3), vertices_dot(2, 3), vertices_bar(2, 3)
    real(dp), allocatable :: matrix(:, :), matrix_dot(:, :)
    real(dp), allocatable :: minus(:, :), plus(:, :)
    real(dp) :: matrix_bar(dof_count, dof_count)
    real(dp) :: derivative_coefficient, derivative_coefficient_bar
    real(dp) :: derivative_coefficient_dot
    real(dp) :: mass_coefficient, mass_coefficient_bar, mass_coefficient_dot
    real(dp) :: lhs, relative_error, rhs
    integer :: column, row, status, status_minus, status_plus

    vertices = reshape([ &
        0.1_dp, -0.2_dp, 1.3_dp, -0.05_dp, 0.2_dp, 1.1_dp], [2, 3])
    vertices_dot = reshape([ &
        0.01_dp, -0.02_dp, -0.03_dp, 0.015_dp, 0.025_dp, -0.01_dp], [2, 3])
    derivative_coefficient = 1.7_dp
    mass_coefficient = -0.8_dp
    derivative_coefficient_dot = 0.13_dp
    mass_coefficient_dot = -0.09_dp
    do column = 1, dof_count
        do row = 1, dof_count
            matrix_bar(row, column) = &
                0.002_dp*row - 0.0015_dp*column
        end do
    end do

    call check_family(.true., "BDM")
    call check_family(.false., "second-kind Nedelec")
    call check_summary("Triangle full-vector element derivatives")

contains

    subroutine check_family(normal_family, label)
        logical, intent(in) :: normal_family
        character(*), intent(in) :: label

        if (normal_family) then
            call assemble_triangle_bdm_div_mass_element( &
                vertices, degree, quadrature_degree, matrix, status, &
                derivative_coefficient, mass_coefficient)
            call assemble_triangle_bdm_div_mass_element_jvp( &
                vertices, degree, quadrature_degree, derivative_coefficient, &
                mass_coefficient, vertices_dot, derivative_coefficient_dot, &
                mass_coefficient_dot, matrix_dot, status)
            call assemble_triangle_bdm_div_mass_element( &
                vertices + step*vertices_dot, degree, quadrature_degree, plus, &
                status_plus, derivative_coefficient + &
                step*derivative_coefficient_dot, mass_coefficient + &
                step*mass_coefficient_dot)
            call assemble_triangle_bdm_div_mass_element( &
                vertices - step*vertices_dot, degree, quadrature_degree, minus, &
                status_minus, derivative_coefficient - &
                step*derivative_coefficient_dot, mass_coefficient - &
                step*mass_coefficient_dot)
        else
            call assemble_triangle_nedelec_second_curl_mass_element( &
                vertices, degree, quadrature_degree, matrix, status, &
                derivative_coefficient, mass_coefficient)
            call assemble_triangle_nedelec_second_curl_mass_element_jvp( &
                vertices, degree, quadrature_degree, derivative_coefficient, &
                mass_coefficient, vertices_dot, derivative_coefficient_dot, &
                mass_coefficient_dot, matrix_dot, status)
            call assemble_triangle_nedelec_second_curl_mass_element( &
                vertices + step*vertices_dot, degree, quadrature_degree, plus, &
                status_plus, derivative_coefficient + &
                step*derivative_coefficient_dot, mass_coefficient + &
                step*mass_coefficient_dot)
            call assemble_triangle_nedelec_second_curl_mass_element( &
                vertices - step*vertices_dot, degree, quadrature_degree, minus, &
                status_minus, derivative_coefficient - &
                step*derivative_coefficient_dot, mass_coefficient - &
                step*mass_coefficient_dot)
        end if
        relative_error = maxval(abs( &
            matrix_dot - (plus - minus)/(2.0_dp*step)))/ &
            max(1.0_dp, maxval(abs(matrix_dot)))
        call check_condition( &
            status == 0 .and. status_plus == 0 .and. status_minus == 0, &
            label//" element JVP accepts a valid direction")
        call check_condition( &
            relative_error < 5.0e-8_dp, &
            label//" element JVP matches independent full reassembly")

        if (normal_family) then
            call assemble_triangle_bdm_div_mass_element_vjp( &
                vertices, degree, quadrature_degree, derivative_coefficient, &
                mass_coefficient, matrix_bar, vertices_bar, &
                derivative_coefficient_bar, mass_coefficient_bar, status)
        else
            call assemble_triangle_nedelec_second_curl_mass_element_vjp( &
                vertices, degree, quadrature_degree, derivative_coefficient, &
                mass_coefficient, matrix_bar, vertices_bar, &
                derivative_coefficient_bar, mass_coefficient_bar, status)
        end if
        lhs = sum(matrix_bar*matrix_dot)
        rhs = sum(vertices_bar*vertices_dot) + &
            derivative_coefficient_bar*derivative_coefficient_dot + &
            mass_coefficient_bar*mass_coefficient_dot
        call check_condition(status == 0, label//" element VJP succeeds")
        call check_condition( &
            abs(lhs - rhs) < 2.0e-10_dp*max(1.0_dp, abs(lhs), abs(rhs)), &
            label//" element products obey the adjoint identity")
    end subroutine check_family

end program test_triangle_full_vector_element_ad
