program test_triangle_full_vector_csc_ad
    use check, only: check_condition, check_summary
    use fortfem_api, only: &
        assemble_triangle_bdm_div_mass_csc, &
        assemble_triangle_bdm_div_mass_csc_jvp, &
        assemble_triangle_bdm_div_mass_csc_vjp, &
        assemble_triangle_nedelec_second_curl_mass_csc, &
        assemble_triangle_nedelec_second_curl_mass_csc_jvp, &
        assemble_triangle_nedelec_second_curl_mass_csc_vjp
    use fortfem_kinds, only: dp
    use fortfem_mesh_2d, only: mesh_2d_t
    use fortsparse, only: csc_t, fortsparse_status_t
    implicit none

    integer, parameter :: degree = 2, quadrature_degree = 5
    real(dp), parameter :: step = 2.0e-7_dp
    type(mesh_2d_t) :: mesh, minus_mesh, plus_mesh
    real(dp) :: vertices_dot(2, 4), vertices_bar(2, 4)
    real(dp), allocatable :: matrix_values_bar(:)
    type(csc_t) :: matrix, matrix_dot, minus, plus
    type(fortsparse_status_t) :: status
    real(dp) :: derivative_coefficient, derivative_coefficient_bar
    real(dp) :: derivative_coefficient_dot
    real(dp) :: mass_coefficient, mass_coefficient_bar, mass_coefficient_dot
    real(dp) :: lhs, relative_error, rhs
    integer :: entry

    vertices_dot = reshape([ &
        0.01_dp, -0.02_dp, -0.03_dp, 0.015_dp, &
        0.025_dp, -0.01_dp, 0.02_dp, 0.03_dp], [2, 4])
    derivative_coefficient = 1.7_dp
    mass_coefficient = -0.8_dp
    derivative_coefficient_dot = 0.13_dp
    mass_coefficient_dot = -0.09_dp

    call check_family(.true., "BDM")
    call check_family(.false., "second-kind Nedelec")
    call check_summary("Global triangle full-vector derivatives")

contains

    subroutine check_family(normal_family, label)
        logical, intent(in) :: normal_family
        character(*), intent(in) :: label

        call mesh%create_rectangular(2, 2, 0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp)
        call plus_mesh%create_rectangular( &
            2, 2, 0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp)
        call minus_mesh%create_rectangular( &
            2, 2, 0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp)
        plus_mesh%vertices = plus_mesh%vertices + step*vertices_dot
        minus_mesh%vertices = minus_mesh%vertices - step*vertices_dot
        if (normal_family) then
            call assemble_triangle_bdm_div_mass_csc( &
                mesh, degree, quadrature_degree, matrix, status, &
                derivative_coefficient, mass_coefficient)
            call assemble_triangle_bdm_div_mass_csc_jvp( &
                mesh, degree, quadrature_degree, derivative_coefficient, &
                mass_coefficient, vertices_dot, derivative_coefficient_dot, &
                mass_coefficient_dot, matrix_dot, status)
            call assemble_triangle_bdm_div_mass_csc( &
                plus_mesh, degree, quadrature_degree, plus, status, &
                derivative_coefficient + step*derivative_coefficient_dot, &
                mass_coefficient + step*mass_coefficient_dot)
            call assemble_triangle_bdm_div_mass_csc( &
                minus_mesh, degree, quadrature_degree, minus, status, &
                derivative_coefficient - step*derivative_coefficient_dot, &
                mass_coefficient - step*mass_coefficient_dot)
        else
            call assemble_triangle_nedelec_second_curl_mass_csc( &
                mesh, degree, quadrature_degree, matrix, status, &
                derivative_coefficient, mass_coefficient)
            call assemble_triangle_nedelec_second_curl_mass_csc_jvp( &
                mesh, degree, quadrature_degree, derivative_coefficient, &
                mass_coefficient, vertices_dot, derivative_coefficient_dot, &
                mass_coefficient_dot, matrix_dot, status)
            call assemble_triangle_nedelec_second_curl_mass_csc( &
                plus_mesh, degree, quadrature_degree, plus, status, &
                derivative_coefficient + step*derivative_coefficient_dot, &
                mass_coefficient + step*mass_coefficient_dot)
            call assemble_triangle_nedelec_second_curl_mass_csc( &
                minus_mesh, degree, quadrature_degree, minus, status, &
                derivative_coefficient - step*derivative_coefficient_dot, &
                mass_coefficient - step*mass_coefficient_dot)
        end if
        call check_condition( &
            matrix_dot%nnz == matrix%nnz .and. plus%nnz == matrix%nnz .and. &
            minus%nnz == matrix%nnz .and. &
            all(matrix_dot%col_ptr == matrix%col_ptr) .and. &
            all(matrix_dot%row_idx == matrix%row_idx), &
            label//" global JVP preserves the merged CSC pattern")
        relative_error = maxval(abs( &
            matrix_dot%val - (plus%val - minus%val)/(2.0_dp*step)))/ &
            max(1.0_dp, maxval(abs(matrix_dot%val)))
        call check_condition( &
            relative_error < 5.0e-8_dp, &
            label//" global JVP matches independent mesh reassembly")

        if (allocated(matrix_values_bar)) deallocate(matrix_values_bar)
        allocate(matrix_values_bar(matrix%nnz))
        do entry = 1, matrix%nnz
            matrix_values_bar(entry) = 0.003_dp*entry - 0.00001_dp*entry**2
        end do
        if (normal_family) then
            call assemble_triangle_bdm_div_mass_csc_vjp( &
                mesh, degree, quadrature_degree, derivative_coefficient, &
                mass_coefficient, matrix_values_bar, vertices_bar, &
                derivative_coefficient_bar, mass_coefficient_bar, status)
        else
            call assemble_triangle_nedelec_second_curl_mass_csc_vjp( &
                mesh, degree, quadrature_degree, derivative_coefficient, &
                mass_coefficient, matrix_values_bar, vertices_bar, &
                derivative_coefficient_bar, mass_coefficient_bar, status)
        end if
        lhs = sum(matrix_values_bar*matrix_dot%val)
        rhs = sum(vertices_bar*vertices_dot) + &
            derivative_coefficient_bar*derivative_coefficient_dot + &
            mass_coefficient_bar*mass_coefficient_dot
        call check_condition(status%code == 0, label//" global VJP succeeds")
        call check_condition( &
            abs(lhs - rhs) < 3.0e-10_dp*max(1.0_dp, abs(lhs), abs(rhs)), &
            label//" global products accumulate shared-vertex adjoints")
    end subroutine check_family

end program test_triangle_full_vector_csc_ad
