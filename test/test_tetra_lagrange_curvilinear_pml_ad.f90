program test_tetra_lagrange_curvilinear_pml_ad
    use check, only: check_condition, check_summary
    use fortfem_api, only: &
        assemble_tetra_lagrange_curvilinear_pml_csc, &
        assemble_tetra_lagrange_curvilinear_pml_csc_jvp, &
        assemble_tetra_lagrange_curvilinear_pml_csc_vjp, &
        assemble_tetra_lagrange_curvilinear_pml_element, &
        assemble_tetra_lagrange_curvilinear_pml_element_jvp, &
        assemble_tetra_lagrange_curvilinear_pml_element_vjp, &
        curvilinear_scalar_helmholtz_pml_coefficients
    use fortfem_kinds, only: dp
    use fortsparse, only: csc_z_t, fortsparse_status_t
    implicit none

    real(dp), parameter :: difference_step = 1.0e-6_dp
    real(dp), parameter :: vertices(3, 4) = reshape([ &
        0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, &
        0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp], [3, 4])
    integer, parameter :: tetrahedra(4, 1) = reshape([1, 2, 3, 4], [4, 1])
    real(dp), parameter :: wave_number = 1.7_dp
    real(dp), parameter :: wave_number_dot = -0.13_dp
    integer, parameter :: degree = 1
    integer, parameter :: quadrature_degree = 2
    complex(dp) :: stretch(3, 3), stretch_dot(3, 3)
    complex(dp) :: stretch_minus(3, 3), stretch_plus(3, 3)
    complex(dp), allocatable :: matrix(:, :), matrix_dot(:, :)
    complex(dp), allocatable :: matrix_minus(:, :), matrix_plus(:, :)
    complex(dp), allocatable :: matrix_bar(:, :), stretch_bar(:, :, :)
    real(dp), allocatable :: vertices_bar(:, :)
    real(dp) :: vertices_dot(3, 4), wave_number_bar
    real(dp) :: forward_pairing, reverse_pairing
    real(dp) :: gradient(3, 4), reference_mass(4, 4)
    complex(dp) :: expected(4, 4)
    real(dp) :: volume
    complex(dp) :: gradient_coefficient(3, 3), mass_coefficient
    complex(dp), allocatable :: csc_matrix_bar(:), csc_stretch_bar(:, :, :)
    type(csc_z_t) :: csc_matrix, csc_matrix_dot, csc_minus, csc_plus
    type(fortsparse_status_t) :: sparse_status
    integer :: column, entry, row, status
    logical :: all_passed

    all_passed = .true.
    vertices_dot = reshape([ &
        0.03_dp, -0.02_dp, 0.04_dp, -0.01_dp, 0.05_dp, 0.02_dp, &
        -0.04_dp, 0.01_dp, 0.02_dp, 0.03_dp, -0.05_dp, 0.01_dp], [3, 4])
    stretch = cmplx(0.0_dp, 0.0_dp, dp)
    stretch(1, 1) = cmplx(1.0_dp, 0.08_dp, dp)
    stretch(2, 2) = cmplx(1.0_dp, 0.13_dp, dp)
    stretch(3, 3) = cmplx(1.0_dp, 0.05_dp, dp)
    stretch(1, 2) = cmplx(0.03_dp, -0.02_dp, dp)
    stretch(2, 1) = cmplx(-0.01_dp, 0.04_dp, dp)
    stretch_dot = reshape([ &
        cmplx(0.02_dp, -0.01_dp, dp), cmplx(-0.03_dp, 0.04_dp, dp), &
        cmplx(0.01_dp, 0.02_dp, dp), cmplx(0.03_dp, 0.01_dp, dp), &
        cmplx(-0.02_dp, 0.03_dp, dp), cmplx(0.01_dp, -0.04_dp, dp), &
        cmplx(0.02_dp, 0.01_dp, dp), cmplx(-0.01_dp, 0.02_dp, dp), &
        cmplx(0.03_dp, -0.02_dp, dp)], [3, 3])

    call assemble_tetra_lagrange_curvilinear_pml_element( &
        vertices, degree, quadrature_degree, stretch, wave_number, matrix, &
        status)
    call record_condition(status == 0, "curvilinear scalar PML element assembles")
    call curvilinear_scalar_helmholtz_pml_coefficients( &
        stretch, gradient_coefficient, mass_coefficient, status)
    volume = 1.0_dp/6.0_dp
    gradient(:, 1) = [-1.0_dp, -1.0_dp, -1.0_dp]
    gradient(:, 2) = [1.0_dp, 0.0_dp, 0.0_dp]
    gradient(:, 3) = [0.0_dp, 1.0_dp, 0.0_dp]
    gradient(:, 4) = [0.0_dp, 0.0_dp, 1.0_dp]
    reference_mass = 1.0_dp/20.0_dp
    do row = 1, 4
        reference_mass(row, row) = 1.0_dp/10.0_dp
    end do
    do column = 1, 4
        do row = 1, 4
            expected(row, column) = volume*dot_product(gradient(:, row), &
                matmul(gradient_coefficient, gradient(:, column))) - &
                volume*wave_number**2*mass_coefficient* &
                reference_mass(row, column)
        end do
    end do
    call record_condition(maxval(abs(matrix - expected)) < 2.0e-12_dp, &
        "curvilinear scalar PML element matches the P1 oracle")

    call assemble_tetra_lagrange_curvilinear_pml_element_jvp( &
        vertices, degree, quadrature_degree, stretch, wave_number, &
        vertices_dot, stretch_dot, wave_number_dot, matrix_dot, status)
    call assemble_tetra_lagrange_curvilinear_pml_element( &
        vertices - difference_step*vertices_dot, degree, quadrature_degree, &
        stretch - difference_step*stretch_dot, &
        wave_number - difference_step*wave_number_dot, matrix_minus, status)
    call assemble_tetra_lagrange_curvilinear_pml_element( &
        vertices + difference_step*vertices_dot, degree, quadrature_degree, &
        stretch + difference_step*stretch_dot, &
        wave_number + difference_step*wave_number_dot, matrix_plus, status)
    call record_condition(maxval(abs(matrix_dot - &
        (matrix_plus - matrix_minus)/(2.0_dp*difference_step))) < 4.0e-9_dp, &
        "curvilinear scalar PML element JVP matches central difference")

    allocate(matrix_bar(4, 4), vertices_bar(3, 4), stretch_bar(3, 3, 1))
    matrix_bar = reshape([ &
        (cmplx(0.03_dp*real(entry, dp), -0.02_dp*real(entry + 1, dp), dp), &
        entry=1, 16)], [4, 4])
    call assemble_tetra_lagrange_curvilinear_pml_element_vjp( &
        vertices, degree, quadrature_degree, stretch, wave_number, matrix_bar, &
        vertices_bar, stretch_bar(:, :, 1), wave_number_bar, status)
    forward_pairing = real(sum(conjg(matrix_bar)*matrix_dot), dp)
    reverse_pairing = real(sum(conjg(stretch_bar(:, :, 1))*stretch_dot), dp) + &
        sum(vertices_bar*vertices_dot) + wave_number_bar*wave_number_dot
    call record_condition(abs(forward_pairing - reverse_pairing) < 2.0e-11_dp, &
        "curvilinear scalar PML element VJP satisfies the complex dot identity")

    call assemble_tetra_lagrange_curvilinear_pml_csc( &
        vertices, tetrahedra, degree, stretch_reshape(stretch), wave_number, &
        csc_matrix, sparse_status)
    call record_condition(sparse_status%code == 0, &
        "curvilinear scalar PML CSC assembles")
    call assemble_tetra_lagrange_curvilinear_pml_csc_jvp( &
        vertices, tetrahedra, degree, stretch_reshape(stretch), wave_number, &
        vertices_dot, stretch_reshape(stretch_dot), wave_number_dot, &
        csc_matrix_dot, sparse_status)
    call assemble_tetra_lagrange_curvilinear_pml_csc( &
        vertices - difference_step*vertices_dot, tetrahedra, degree, &
        stretch_reshape(stretch_minus_value()), &
        wave_number - difference_step*wave_number_dot, csc_minus, sparse_status)
    call assemble_tetra_lagrange_curvilinear_pml_csc( &
        vertices + difference_step*vertices_dot, tetrahedra, degree, &
        stretch_reshape(stretch_plus_value()), &
        wave_number + difference_step*wave_number_dot, csc_plus, sparse_status)
    call record_condition(all(csc_matrix_dot%col_ptr == csc_matrix%col_ptr) .and. &
        maxval(abs(csc_matrix_dot%val - &
        (csc_plus%val - csc_minus%val)/(2.0_dp*difference_step))) < 4.0e-9_dp, &
        "curvilinear scalar PML CSC JVP preserves the pattern")

    allocate(csc_matrix_bar(csc_matrix%nnz), csc_stretch_bar(3, 3, 1))
    do column = 1, csc_matrix%ncol
        do entry = csc_matrix%col_ptr(column), csc_matrix%col_ptr(column + 1) - 1
            row = csc_matrix%row_idx(entry)
            csc_matrix_bar(entry) = matrix_bar(row, column)
        end do
    end do
    call assemble_tetra_lagrange_curvilinear_pml_csc_vjp( &
        vertices, tetrahedra, degree, stretch_reshape(stretch), wave_number, &
        csc_matrix_bar, vertices_bar, csc_stretch_bar, wave_number_bar, &
        sparse_status)
    forward_pairing = real(sum(conjg(csc_matrix_bar)*csc_matrix_dot%val), dp)
    reverse_pairing = real(sum(conjg(csc_stretch_bar)*stretch_reshape(stretch_dot)), &
        dp) + sum(vertices_bar*vertices_dot) + wave_number_bar*wave_number_dot
    call record_condition(abs(forward_pairing - reverse_pairing) < 2.0e-11_dp, &
        "curvilinear scalar PML CSC VJP satisfies the complex dot identity")

    call check_summary("Tetrahedral scalar curvilinear PML derivatives")
    if (.not. all_passed) error stop 1

contains

    function stretch_reshape(value) result(output)
        complex(dp), intent(in) :: value(3, 3)
        complex(dp) :: output(3, 3, 1)

        output(:, :, 1) = value
    end function stretch_reshape

    function stretch_minus_value() result(output)
        complex(dp) :: output(3, 3)

        output = stretch - difference_step*stretch_dot
    end function stretch_minus_value

    function stretch_plus_value() result(output)
        complex(dp) :: output(3, 3)

        output = stretch + difference_step*stretch_dot
    end function stretch_plus_value

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_tetra_lagrange_curvilinear_pml_ad
