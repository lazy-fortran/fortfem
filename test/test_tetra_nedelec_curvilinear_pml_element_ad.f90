program test_tetra_nedelec_curvilinear_pml_element_ad
    use check, only: check_condition, check_summary
    use fortfem_api, only: &
        assemble_tetra_nedelec_curvilinear_pml_element, &
        assemble_tetra_nedelec_curvilinear_pml_element_jvp, &
        assemble_tetra_nedelec_curvilinear_pml_element_vjp, &
        assemble_tetra_nedelec_pml_element
    use fortfem_kinds, only: dp
    implicit none

    integer, parameter :: order = 1, quadrature_degree = 4
    integer, parameter :: dof_count = order*(order + 2)*(order + 3)/2
    real(dp), parameter :: step = 2.0e-7_dp
    real(dp) :: vertices(3, 4), vertices_dot(3, 4), vertices_bar(3, 4)
    complex(dp) :: stretch(3, 3), stretch_dot(3, 3), stretch_bar(3, 3)
    complex(dp) :: matrix_bar(dof_count, dof_count)
    complex(dp), allocatable :: matrix(:, :), matrix_dot(:, :)
    complex(dp), allocatable :: plus(:, :), minus(:, :)
    complex(dp) :: diagonal(3, 3), diagonal_vector(3)
    complex(dp), allocatable :: diagonal_matrix(:, :)
    real(dp) :: wave_number, wave_number_dot, wave_number_bar
    real(dp) :: lhs, rhs, relative_error
    integer :: column, row, status, status_plus, status_minus
    logical :: all_passed

    all_passed = .true.
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
    stretch = reshape([ &
        cmplx(1.1_dp, 0.2_dp, dp), cmplx(0.08_dp, -0.03_dp, dp), &
        cmplx(-0.02_dp, 0.01_dp, dp), cmplx(0.04_dp, 0.06_dp, dp), &
        cmplx(0.95_dp, 0.15_dp, dp), cmplx(0.05_dp, -0.02_dp, dp), &
        cmplx(0.03_dp, 0.01_dp, dp), cmplx(-0.04_dp, 0.02_dp, dp), &
        cmplx(1.05_dp, 0.18_dp, dp)], [3, 3])
    stretch_dot = reshape([ &
        cmplx(0.03_dp, -0.01_dp, dp), cmplx(-0.02_dp, 0.04_dp, dp), &
        cmplx(0.01_dp, 0.02_dp, dp), cmplx(0.05_dp, 0.01_dp, dp), &
        cmplx(-0.03_dp, 0.02_dp, dp), cmplx(0.02_dp, -0.03_dp, dp), &
        cmplx(-0.01_dp, 0.03_dp, dp), cmplx(0.04_dp, -0.02_dp, dp), &
        cmplx(0.01_dp, 0.01_dp, dp)], [3, 3])
    wave_number = 1.7_dp
    wave_number_dot = 0.13_dp
    do column = 1, dof_count
        do row = 1, dof_count
            matrix_bar(row, column) = cmplx( &
                0.002_dp*row - 0.0015_dp*column, &
                -0.001_dp*row + 0.0025_dp*column, dp)
        end do
    end do

    call assemble_tetra_nedelec_curvilinear_pml_element( &
        vertices, order, quadrature_degree, stretch, wave_number, matrix, status)
    call assemble_tetra_nedelec_curvilinear_pml_element_jvp( &
        vertices, order, quadrature_degree, stretch, wave_number, vertices_dot, &
        stretch_dot, wave_number_dot, matrix_dot, status)
    call assemble_tetra_nedelec_curvilinear_pml_element( &
        vertices + step*vertices_dot, order, quadrature_degree, &
        stretch + step*stretch_dot, wave_number + step*wave_number_dot, plus, &
        status_plus)
    call assemble_tetra_nedelec_curvilinear_pml_element( &
        vertices - step*vertices_dot, order, quadrature_degree, &
        stretch - step*stretch_dot, wave_number - step*wave_number_dot, minus, &
        status_minus)
    relative_error = maxval(abs( &
        matrix_dot - (plus - minus)/(2.0_dp*step)))/ &
        max(1.0_dp, maxval(abs(matrix_dot)))
    call record_condition( &
        status == 0 .and. status_plus == 0 .and. status_minus == 0, &
        "curvilinear tetrahedral PML element accepts a full stretch")
    call record_condition( &
        relative_error < 8.0e-8_dp, &
        "curvilinear tetrahedral PML JVP matches reassembly differences")

    call assemble_tetra_nedelec_curvilinear_pml_element_vjp( &
        vertices, order, quadrature_degree, stretch, wave_number, matrix_bar, &
        vertices_bar, stretch_bar, wave_number_bar, status)
    lhs = real(sum(conjg(matrix_bar)*matrix_dot), dp)
    rhs = sum(vertices_bar*vertices_dot) + &
        real(sum(conjg(stretch_bar)*stretch_dot), dp) + &
        wave_number_bar*wave_number_dot
    call record_condition(status == 0, "curvilinear tetrahedral PML VJP succeeds")
    call record_condition( &
        abs(lhs - rhs) < 4.0e-10_dp*max(1.0_dp, abs(lhs), abs(rhs)), &
        "curvilinear tetrahedral PML products satisfy the complex adjoint identity")

    diagonal = cmplx(0.0_dp, 0.0_dp, dp)
    diagonal(1, 1) = stretch(1, 1)
    diagonal(2, 2) = stretch(2, 2)
    diagonal(3, 3) = stretch(3, 3)
    diagonal_vector = [diagonal(1, 1), diagonal(2, 2), diagonal(3, 3)]
    call assemble_tetra_nedelec_curvilinear_pml_element( &
        vertices, order, quadrature_degree, diagonal, wave_number, &
        diagonal_matrix, status)
    call assemble_tetra_nedelec_pml_element( &
        vertices, order, quadrature_degree, diagonal_vector, wave_number, &
        matrix, status)
    call record_condition(status == 0, &
        "diagonal curvilinear tetrahedral PML assembly succeeds")
    call record_condition(maxval(abs(diagonal_matrix - matrix)) < 3.0e-12_dp, &
        "curvilinear tetrahedral PML reduces to Cartesian assembly")

    call check_summary("Curvilinear tetrahedral Nedelec PML derivatives")

    if (.not. all_passed) error stop 1

contains

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_tetra_nedelec_curvilinear_pml_element_ad
