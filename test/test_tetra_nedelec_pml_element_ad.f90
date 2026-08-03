program test_tetra_nedelec_pml_element_ad
    use check, only: check_condition, check_summary
    use fortfem_feec, only: assemble_tetra_nedelec_pml_element, &
        assemble_tetra_nedelec_pml_element_jvp, &
        assemble_tetra_nedelec_pml_element_vjp
    use fortfem_kinds, only: dp
    implicit none

    integer, parameter :: order = 2, quadrature_degree = 6
    integer, parameter :: dof_count = order*(order + 2)*(order + 3)/2
    real(dp), parameter :: step = 2.0e-7_dp
    real(dp) :: vertices(3, 4), vertices_dot(3, 4), vertices_bar(3, 4)
    complex(dp) :: stretch(3), stretch_dot(3), stretch_bar(3)
    complex(dp), allocatable :: matrix(:, :), matrix_dot(:, :)
    complex(dp), allocatable :: plus(:, :), minus(:, :)
    complex(dp) :: matrix_bar(dof_count, dof_count)
    real(dp) :: wave_number, wave_number_dot, wave_number_bar
    real(dp) :: lhs, rhs, relative_error
    integer :: column, row, status, status_plus, status_minus

    vertices = reshape([ &
        0.1_dp, -0.2_dp, 0.3_dp, &
        1.3_dp, -0.1_dp, 0.25_dp, &
        0.25_dp, 0.9_dp, 0.35_dp, &
        0.0_dp, 0.1_dp, 1.4_dp], [3, 4])
    vertices_dot = reshape([ &
        0.01_dp, -0.02_dp, 0.015_dp, &
        -0.03_dp, 0.01_dp, 0.02_dp, &
        0.025_dp, -0.015_dp, 0.005_dp, &
        -0.01_dp, 0.03_dp, -0.02_dp], [3, 4])
    stretch = [ &
        cmplx(1.1_dp, 0.35_dp, dp), cmplx(0.95_dp, 0.12_dp, dp), &
        cmplx(1.05_dp, 0.27_dp, dp)]
    stretch_dot = [ &
        cmplx(0.04_dp, -0.03_dp, dp), cmplx(-0.02_dp, 0.01_dp, dp), &
        cmplx(0.03_dp, 0.02_dp, dp)]
    wave_number = 1.7_dp
    wave_number_dot = 0.13_dp
    do column = 1, dof_count
        do row = 1, dof_count
            matrix_bar(row, column) = cmplx( &
                0.002_dp*row - 0.0015_dp*column, &
                -0.001_dp*row + 0.0025_dp*column, dp)
        end do
    end do

    call assemble_tetra_nedelec_pml_element( &
        vertices, order, quadrature_degree, stretch, wave_number, matrix, status)
    call assemble_tetra_nedelec_pml_element_jvp( &
        vertices, order, quadrature_degree, stretch, wave_number, vertices_dot, &
        stretch_dot, wave_number_dot, matrix_dot, status)
    call assemble_tetra_nedelec_pml_element( &
        vertices + step*vertices_dot, order, quadrature_degree, &
        stretch + step*stretch_dot, wave_number + step*wave_number_dot, plus, &
        status_plus)
    call assemble_tetra_nedelec_pml_element( &
        vertices - step*vertices_dot, order, quadrature_degree, &
        stretch - step*stretch_dot, wave_number - step*wave_number_dot, minus, &
        status_minus)
    relative_error = maxval(abs( &
        matrix_dot - (plus - minus)/(2.0_dp*step)))/ &
        max(1.0_dp, maxval(abs(matrix_dot)))
    call check_condition( &
        status == 0 .and. status_plus == 0 .and. status_minus == 0, &
        "Nedelec PML element JVP accepts complex stretch and geometry")
    call check_condition( &
        relative_error < 6.0e-8_dp, &
        "Nedelec PML element JVP matches complete reassembly difference")

    call assemble_tetra_nedelec_pml_element_vjp( &
        vertices, order, quadrature_degree, stretch, wave_number, matrix_bar, &
        vertices_bar, stretch_bar, wave_number_bar, status)
    lhs = real(sum(conjg(matrix_bar)*matrix_dot), dp)
    rhs = sum(vertices_bar*vertices_dot) + &
        real(sum(conjg(stretch_bar)*stretch_dot), dp) + &
        wave_number_bar*wave_number_dot
    call check_condition(status == 0, "Nedelec PML element VJP succeeds")
    call check_condition( &
        abs(lhs - rhs) < 3.0e-10_dp*max(1.0_dp, abs(lhs), abs(rhs)), &
        "Nedelec PML element products satisfy the complex adjoint identity")

    call check_summary("Tetrahedral Nedelec PML element derivatives")
end program test_tetra_nedelec_pml_element_ad
