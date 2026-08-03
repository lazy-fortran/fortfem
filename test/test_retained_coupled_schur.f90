program test_retained_coupled_schur
    use, intrinsic :: iso_fortran_env, only: real64
    use check, only: check_condition, check_summary
    use fortfem_feec, only: &
        assemble_retained_coupled_schur, &
        assemble_retained_coupled_schur_jvp, &
        assemble_retained_coupled_schur_vjp, &
        factor_retained_complex_field_split, factor_retained_field_split, &
        free_retained_complex_field_split, free_retained_field_split, &
        retained_complex_field_split_t, retained_field_split_t
    use fortsparse, only: csc_t, csc_z_t
    implicit none

    integer, parameter :: dp = real64
    type(csc_t) :: matrices(2), matrices_dot(2), matrices_bar(2)
    type(csc_z_t) :: complex_matrices(2), complex_matrices_dot(2)
    type(csc_z_t) :: complex_matrices_bar(2)
    type(retained_field_split_t) :: split
    type(retained_complex_field_split_t) :: complex_split
    real(dp) :: exterior_matrix(2, 2), exterior_to_fields(2, 2)
    real(dp) :: fields_to_exterior(2, 2), exterior_rhs(2), fields_rhs(2)
    real(dp) :: exterior_matrix_dot(2, 2), exterior_to_fields_dot(2, 2)
    real(dp) :: fields_to_exterior_dot(2, 2), exterior_rhs_dot(2), fields_rhs_dot(2)
    real(dp) :: schur_matrix(2, 2), schur_rhs(2), schur_matrix_dot(2, 2), schur_rhs_dot(2)
    real(dp) :: schur_matrix_ref(2, 2), schur_rhs_ref(2)
    real(dp) :: schur_plus(2, 2), schur_minus(2, 2), rhs_plus(2), rhs_minus(2)
    real(dp) :: schur_matrix_bar(2, 2), schur_rhs_bar(2)
    real(dp) :: exterior_matrix_bar(2, 2), exterior_to_fields_bar(2, 2)
    real(dp) :: fields_to_exterior_bar(2, 2), exterior_rhs_bar(2), fields_rhs_bar(2)
    real(dp) :: lhs, rhs_inner
    complex(dp) :: complex_exterior_matrix(2, 2), complex_exterior_to_fields(2, 2)
    complex(dp) :: complex_fields_to_exterior(2, 2), complex_exterior_rhs(2)
    complex(dp) :: complex_fields_rhs(2), complex_exterior_matrix_dot(2, 2)
    complex(dp) :: complex_exterior_to_fields_dot(2, 2)
    complex(dp) :: complex_fields_to_exterior_dot(2, 2)
    complex(dp) :: complex_exterior_rhs_dot(2), complex_fields_rhs_dot(2)
    complex(dp) :: complex_schur_matrix(2, 2), complex_schur_rhs(2)
    complex(dp) :: complex_schur_matrix_dot(2, 2), complex_schur_rhs_dot(2)
    complex(dp) :: complex_schur_matrix_ref(2, 2), complex_schur_rhs_ref(2)
    complex(dp) :: complex_schur_plus(2, 2), complex_schur_minus(2, 2)
    complex(dp) :: complex_rhs_plus(2), complex_rhs_minus(2)
    complex(dp) :: complex_schur_matrix_bar(2, 2), complex_schur_rhs_bar(2)
    complex(dp) :: complex_exterior_matrix_bar(2, 2)
    complex(dp) :: complex_exterior_to_fields_bar(2, 2)
    complex(dp) :: complex_fields_to_exterior_bar(2, 2)
    complex(dp) :: complex_exterior_rhs_bar(2), complex_fields_rhs_bar(2)
    real(dp) :: complex_lhs, complex_rhs_inner
    integer :: status
    logical :: all_passed
    real(dp), parameter :: step = 1.0e-6_dp

    all_passed = .true.
    call initialize_real_data( &
        matrices, matrices_dot, exterior_matrix, exterior_to_fields, &
        fields_to_exterior, exterior_rhs, fields_rhs, exterior_matrix_dot, &
        exterior_to_fields_dot, fields_to_exterior_dot, exterior_rhs_dot, &
        fields_rhs_dot)
    call factor_retained_field_split(matrices, split, status)
    call record_condition(status == 0, "real retained Schur factors are valid")
    call assemble_retained_coupled_schur( &
        split, exterior_matrix, exterior_to_fields, fields_to_exterior, &
        exterior_rhs, fields_rhs, schur_matrix, schur_rhs, status)
    call dense_real_oracle( &
        matrices(1)%val(1), matrices(2)%val(1), exterior_matrix, &
        exterior_to_fields, fields_to_exterior, exterior_rhs, fields_rhs, &
        schur_matrix_ref, schur_rhs_ref)
    call record_condition(status == 0 .and. maxval(abs(schur_matrix - schur_matrix_ref)) < &
        1.0e-13_dp .and. maxval(abs(schur_rhs - schur_rhs_ref)) < 1.0e-13_dp, &
        "real retained Schur value matches dense elimination oracle")
    call assemble_retained_coupled_schur_jvp( &
        split, exterior_matrix, exterior_to_fields, fields_to_exterior, &
        exterior_rhs, fields_rhs, matrices_dot, exterior_matrix_dot, &
        exterior_to_fields_dot, fields_to_exterior_dot, exterior_rhs_dot, &
        fields_rhs_dot, schur_matrix_dot, schur_rhs_dot, status)
    call dense_real_oracle( &
        matrices(1)%val(1) + step*matrices_dot(1)%val(1), &
        matrices(2)%val(1) + step*matrices_dot(2)%val(1), &
        exterior_matrix + step*exterior_matrix_dot, &
        exterior_to_fields + step*exterior_to_fields_dot, &
        fields_to_exterior + step*fields_to_exterior_dot, &
        exterior_rhs + step*exterior_rhs_dot, fields_rhs + step*fields_rhs_dot, &
        schur_plus, rhs_plus)
    call dense_real_oracle( &
        matrices(1)%val(1) - step*matrices_dot(1)%val(1), &
        matrices(2)%val(1) - step*matrices_dot(2)%val(1), &
        exterior_matrix - step*exterior_matrix_dot, &
        exterior_to_fields - step*exterior_to_fields_dot, &
        fields_to_exterior - step*fields_to_exterior_dot, &
        exterior_rhs - step*exterior_rhs_dot, fields_rhs - step*fields_rhs_dot, &
        schur_minus, rhs_minus)
    call record_condition(status == 0 .and. maxval(abs(schur_matrix_dot - &
        (schur_plus - schur_minus)/(2.0_dp*step))) < 1.0e-8_dp .and. &
        maxval(abs(schur_rhs_dot - (rhs_plus - rhs_minus)/(2.0_dp*step))) < 1.0e-8_dp, &
        "real retained Schur JVP matches independent reassembly")
    schur_matrix_bar = reshape([0.4_dp, -0.2_dp, 0.3_dp, -0.1_dp], [2, 2])
    schur_rhs_bar = [0.25_dp, -0.35_dp]
    call assemble_retained_coupled_schur_vjp( &
        split, exterior_matrix, exterior_to_fields, fields_to_exterior, &
        exterior_rhs, fields_rhs, schur_matrix_bar, schur_rhs_bar, &
        exterior_matrix_bar, exterior_to_fields_bar, fields_to_exterior_bar, &
        exterior_rhs_bar, fields_rhs_bar, matrices_bar, status)
    lhs = sum(schur_matrix_bar*schur_matrix_dot) + dot_product(schur_rhs_bar, schur_rhs_dot)
    rhs_inner = sum(exterior_matrix_bar*exterior_matrix_dot) + &
        sum(exterior_to_fields_bar*exterior_to_fields_dot) + &
        sum(fields_to_exterior_bar*fields_to_exterior_dot) + &
        dot_product(exterior_rhs_bar, exterior_rhs_dot) + &
        dot_product(fields_rhs_bar, fields_rhs_dot) + &
        dot_product(matrices_bar(1)%val, matrices_dot(1)%val) + &
        dot_product(matrices_bar(2)%val, matrices_dot(2)%val)
    call record_condition(status == 0 .and. abs(lhs - rhs_inner) < 1.0e-10_dp, &
        "real retained Schur VJP satisfies the adjoint identity")
    call free_retained_field_split(split)

    call initialize_complex_data( &
        complex_matrices, complex_matrices_dot, complex_exterior_matrix, &
        complex_exterior_to_fields, complex_fields_to_exterior, &
        complex_exterior_rhs, complex_fields_rhs, complex_exterior_matrix_dot, &
        complex_exterior_to_fields_dot, complex_fields_to_exterior_dot, &
        complex_exterior_rhs_dot, complex_fields_rhs_dot)
    call factor_retained_complex_field_split(complex_matrices, complex_split, status)
    call record_condition(status == 0, "complex retained Schur factors are valid")
    call assemble_retained_coupled_schur( &
        complex_split, complex_exterior_matrix, complex_exterior_to_fields, &
        complex_fields_to_exterior, complex_exterior_rhs, complex_fields_rhs, &
        complex_schur_matrix, complex_schur_rhs, status)
    call dense_complex_oracle( &
        complex_matrices(1)%val(1), complex_matrices(2)%val(1), &
        complex_exterior_matrix, complex_exterior_to_fields, &
        complex_fields_to_exterior, complex_exterior_rhs, complex_fields_rhs, &
        complex_schur_matrix_ref, complex_schur_rhs_ref)
    call record_condition(status == 0 .and. maxval(abs(complex_schur_matrix - &
        complex_schur_matrix_ref)) < 1.0e-13_dp .and. maxval(abs(complex_schur_rhs - &
        complex_schur_rhs_ref)) < 1.0e-13_dp, &
        "complex retained Schur value matches dense elimination oracle")
    call assemble_retained_coupled_schur_jvp( &
        complex_split, complex_exterior_matrix, complex_exterior_to_fields, &
        complex_fields_to_exterior, complex_exterior_rhs, complex_fields_rhs, &
        complex_matrices_dot, complex_exterior_matrix_dot, &
        complex_exterior_to_fields_dot, complex_fields_to_exterior_dot, &
        complex_exterior_rhs_dot, complex_fields_rhs_dot, complex_schur_matrix_dot, &
        complex_schur_rhs_dot, status)
    call dense_complex_oracle( &
        complex_matrices(1)%val(1) + step*complex_matrices_dot(1)%val(1), &
        complex_matrices(2)%val(1) + step*complex_matrices_dot(2)%val(1), &
        complex_exterior_matrix + step*complex_exterior_matrix_dot, &
        complex_exterior_to_fields + step*complex_exterior_to_fields_dot, &
        complex_fields_to_exterior + step*complex_fields_to_exterior_dot, &
        complex_exterior_rhs + step*complex_exterior_rhs_dot, &
        complex_fields_rhs + step*complex_fields_rhs_dot, &
        complex_schur_plus, complex_rhs_plus)
    call dense_complex_oracle( &
        complex_matrices(1)%val(1) - step*complex_matrices_dot(1)%val(1), &
        complex_matrices(2)%val(1) - step*complex_matrices_dot(2)%val(1), &
        complex_exterior_matrix - step*complex_exterior_matrix_dot, &
        complex_exterior_to_fields - step*complex_exterior_to_fields_dot, &
        complex_fields_to_exterior - step*complex_fields_to_exterior_dot, &
        complex_exterior_rhs - step*complex_exterior_rhs_dot, &
        complex_fields_rhs - step*complex_fields_rhs_dot, &
        complex_schur_minus, complex_rhs_minus)
    call record_condition(status == 0 .and. maxval(abs(complex_schur_matrix_dot - &
        (complex_schur_plus - complex_schur_minus)/(2.0_dp*step))) < 1.0e-8_dp .and. &
        maxval(abs(complex_schur_rhs_dot - (complex_rhs_plus - complex_rhs_minus)/ &
        (2.0_dp*step))) < 1.0e-8_dp, &
        "complex retained Schur JVP matches independent reassembly")
    complex_schur_matrix_bar = reshape([ &
        cmplx(0.4_dp, 0.1_dp, dp), cmplx(-0.2_dp, 0.3_dp, dp), &
        cmplx(0.3_dp, -0.2_dp, dp), cmplx(-0.1_dp, 0.05_dp, dp)], [2, 2])
    complex_schur_rhs_bar = [cmplx(0.25_dp, -0.1_dp, dp), cmplx(-0.35_dp, 0.2_dp, dp)]
    call assemble_retained_coupled_schur_vjp( &
        complex_split, complex_exterior_matrix, complex_exterior_to_fields, &
        complex_fields_to_exterior, complex_exterior_rhs, complex_fields_rhs, &
        complex_schur_matrix_bar, complex_schur_rhs_bar, &
        complex_exterior_matrix_bar, complex_exterior_to_fields_bar, &
        complex_fields_to_exterior_bar, complex_exterior_rhs_bar, &
        complex_fields_rhs_bar, complex_matrices_bar, status)
    complex_lhs = real(sum(conjg(complex_schur_matrix_bar)*complex_schur_matrix_dot) + &
        sum(conjg(complex_schur_rhs_bar)*complex_schur_rhs_dot), dp)
    complex_rhs_inner = real(sum(conjg(complex_exterior_matrix_bar)* &
        complex_exterior_matrix_dot) + sum(conjg(complex_exterior_to_fields_bar)* &
        complex_exterior_to_fields_dot) + sum(conjg(complex_fields_to_exterior_bar)* &
        complex_fields_to_exterior_dot) + sum(conjg(complex_exterior_rhs_bar)* &
        complex_exterior_rhs_dot) + sum(conjg(complex_fields_rhs_bar)* &
        complex_fields_rhs_dot) + sum(conjg(complex_matrices_bar(1)%val)* &
        complex_matrices_dot(1)%val) + sum(conjg(complex_matrices_bar(2)%val)* &
        complex_matrices_dot(2)%val), dp)
    call record_condition(status == 0 .and. abs(complex_lhs - complex_rhs_inner) < 1.0e-10_dp, &
        "complex retained Schur VJP satisfies the real adjoint identity")
    call free_retained_complex_field_split(complex_split)

    call check_summary("retained coupled Schur")
    if (.not. all_passed) error stop 1

contains

    subroutine initialize_real_data( &
            matrices, matrices_dot, exterior_matrix, exterior_to_fields, &
            fields_to_exterior, exterior_rhs, fields_rhs, exterior_matrix_dot, &
            exterior_to_fields_dot, fields_to_exterior_dot, exterior_rhs_dot, &
            fields_rhs_dot)
        type(csc_t), intent(out) :: matrices(:), matrices_dot(:)
        real(dp), intent(out) :: exterior_matrix(:, :), exterior_to_fields(:, :)
        real(dp), intent(out) :: fields_to_exterior(:, :), exterior_rhs(:), fields_rhs(:)
        real(dp), intent(out) :: exterior_matrix_dot(:, :), exterior_to_fields_dot(:, :)
        real(dp), intent(out) :: fields_to_exterior_dot(:, :)
        real(dp), intent(out) :: exterior_rhs_dot(:), fields_rhs_dot(:)

        call set_real_matrix(matrices(1), 2.0_dp)
        call set_real_matrix(matrices(2), 3.0_dp)
        matrices_dot = matrices
        matrices_dot(1)%val = [0.07_dp]
        matrices_dot(2)%val = [-0.04_dp]
        exterior_matrix = reshape([4.0_dp, 0.1_dp, 0.2_dp, 3.0_dp], [2, 2])
        exterior_to_fields = reshape([0.5_dp, 0.3_dp, -0.2_dp, 0.4_dp], [2, 2])
        fields_to_exterior = reshape([0.7_dp, -0.4_dp, 0.1_dp, 0.6_dp], [2, 2])
        exterior_rhs = [1.3_dp, -0.8_dp]
        fields_rhs = [0.9_dp, -1.2_dp]
        exterior_matrix_dot = reshape([0.03_dp, -0.02_dp, 0.01_dp, 0.04_dp], [2, 2])
        exterior_to_fields_dot = reshape([0.02_dp, 0.01_dp, -0.03_dp, 0.04_dp], [2, 2])
        fields_to_exterior_dot = reshape([-0.01_dp, 0.03_dp, 0.02_dp, -0.02_dp], [2, 2])
        exterior_rhs_dot = [0.05_dp, -0.04_dp]
        fields_rhs_dot = [-0.02_dp, 0.06_dp]
    end subroutine initialize_real_data

    subroutine initialize_complex_data( &
            matrices, matrices_dot, exterior_matrix, exterior_to_fields, &
            fields_to_exterior, exterior_rhs, fields_rhs, exterior_matrix_dot, &
            exterior_to_fields_dot, fields_to_exterior_dot, exterior_rhs_dot, &
            fields_rhs_dot)
        type(csc_z_t), intent(out) :: matrices(:), matrices_dot(:)
        complex(dp), intent(out) :: exterior_matrix(:, :), exterior_to_fields(:, :)
        complex(dp), intent(out) :: fields_to_exterior(:, :)
        complex(dp), intent(out) :: exterior_rhs(:), fields_rhs(:)
        complex(dp), intent(out) :: exterior_matrix_dot(:, :)
        complex(dp), intent(out) :: exterior_to_fields_dot(:, :)
        complex(dp), intent(out) :: fields_to_exterior_dot(:, :)
        complex(dp), intent(out) :: exterior_rhs_dot(:), fields_rhs_dot(:)

        call set_complex_matrix(matrices(1), cmplx(2.0_dp, 0.2_dp, dp))
        call set_complex_matrix(matrices(2), cmplx(3.0_dp, -0.3_dp, dp))
        matrices_dot = matrices
        matrices_dot(1)%val = [cmplx(0.07_dp, -0.01_dp, dp)]
        matrices_dot(2)%val = [cmplx(-0.04_dp, 0.02_dp, dp)]
        exterior_matrix = cmplx(reshape([4.0_dp, 0.1_dp, 0.2_dp, 3.0_dp], [2, 2]), &
            reshape([0.2_dp, -0.1_dp, 0.3_dp, 0.05_dp], [2, 2]), dp)
        exterior_to_fields = cmplx(reshape([0.5_dp, 0.3_dp, -0.2_dp, 0.4_dp], [2, 2]), &
            reshape([0.1_dp, -0.2_dp, 0.04_dp, 0.05_dp], [2, 2]), dp)
        fields_to_exterior = cmplx(reshape([0.7_dp, -0.4_dp, 0.1_dp, 0.6_dp], [2, 2]), &
            reshape([-0.03_dp, 0.02_dp, 0.01_dp, -0.04_dp], [2, 2]), dp)
        exterior_rhs = [cmplx(1.3_dp, 0.2_dp, dp), cmplx(-0.8_dp, 0.1_dp, dp)]
        fields_rhs = [cmplx(0.9_dp, -0.3_dp, dp), cmplx(-1.2_dp, 0.2_dp, dp)]
        exterior_matrix_dot = cmplx(reshape([0.03_dp, -0.02_dp, 0.01_dp, 0.04_dp], [2, 2]), &
            reshape([-0.04_dp, 0.02_dp, 0.03_dp, -0.01_dp], [2, 2]), dp)
        exterior_to_fields_dot = cmplx(reshape([0.02_dp, 0.01_dp, -0.03_dp, 0.04_dp], [2, 2]), &
            reshape([0.01_dp, -0.02_dp, 0.03_dp, 0.02_dp], [2, 2]), dp)
        fields_to_exterior_dot = cmplx(reshape([-0.01_dp, 0.03_dp, 0.02_dp, -0.02_dp], [2, 2]), &
            reshape([0.02_dp, 0.01_dp, -0.01_dp, 0.03_dp], [2, 2]), dp)
        exterior_rhs_dot = [cmplx(0.05_dp, -0.01_dp, dp), cmplx(-0.04_dp, 0.02_dp, dp)]
        fields_rhs_dot = [cmplx(-0.02_dp, 0.03_dp, dp), cmplx(0.06_dp, -0.01_dp, dp)]
    end subroutine initialize_complex_data

    subroutine set_real_matrix(matrix, value)
        type(csc_t), intent(out) :: matrix
        real(dp), intent(in) :: value
        matrix%nrow = 1
        matrix%ncol = 1
        matrix%nnz = 1
        matrix%col_ptr = [1, 2]
        matrix%row_idx = [1]
        matrix%val = [value]
    end subroutine set_real_matrix

    subroutine set_complex_matrix(matrix, value)
        type(csc_z_t), intent(out) :: matrix
        complex(dp), intent(in) :: value
        matrix%nrow = 1
        matrix%ncol = 1
        matrix%nnz = 1
        matrix%col_ptr = [1, 2]
        matrix%row_idx = [1]
        matrix%val = [value]
    end subroutine set_complex_matrix

    subroutine dense_real_oracle(d1, d2, e, c, f, rhs_e, rhs_f, schur, reduced_rhs)
        real(dp), intent(in) :: d1, d2, e(:, :), c(:, :), f(:, :)
        real(dp), intent(in) :: rhs_e(:), rhs_f(:)
        real(dp), intent(out) :: schur(:, :), reduced_rhs(:)
        real(dp) :: x(2), response(2, 2), d(2)

        d = [d1, d2]
        x = rhs_f/d
        response(1, :) = f(1, :)/d1
        response(2, :) = f(2, :)/d2
        schur = e - matmul(c, response)
        reduced_rhs = rhs_e - matmul(c, x)
    end subroutine dense_real_oracle

    subroutine dense_complex_oracle(d1, d2, e, c, f, rhs_e, rhs_f, schur, reduced_rhs)
        complex(dp), intent(in) :: d1, d2, e(:, :), c(:, :), f(:, :)
        complex(dp), intent(in) :: rhs_e(:), rhs_f(:)
        complex(dp), intent(out) :: schur(:, :), reduced_rhs(:)
        complex(dp) :: x(2), response(2, 2)

        x = rhs_f/[d1, d2]
        response(1, :) = f(1, :)/d1
        response(2, :) = f(2, :)/d2
        schur = e - matmul(c, response)
        reduced_rhs = rhs_e - matmul(c, x)
    end subroutine dense_complex_oracle

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_retained_coupled_schur
