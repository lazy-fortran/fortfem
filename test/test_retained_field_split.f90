program test_retained_field_split
    use, intrinsic :: iso_fortran_env, only: real64
    use check, only: check_condition, check_summary
    use fortfem_api, only: &
        apply_retained_complex_field_split, apply_retained_complex_field_split_jvp, &
        apply_retained_complex_field_split_vjp, apply_retained_field_split, &
        apply_retained_field_split_jvp, apply_retained_field_split_vjp, &
        factor_retained_complex_field_split, factor_retained_field_split, &
        free_retained_complex_field_split, free_retained_field_split, &
        retained_complex_field_split_t, retained_field_split_t
    use fortsparse, only: csc_t, csc_z_t
    implicit none

    integer, parameter :: dp = real64
    type(csc_t) :: matrices(2), matrices_dot(2), matrices_bar(2)
    type(csc_z_t) :: complex_matrices(2), complex_matrices_dot(2), complex_matrices_bar(2)
    type(retained_field_split_t) :: split
    type(retained_complex_field_split_t) :: complex_split
    real(dp) :: rhs(3), solution(3), rhs_dot(3), solution_dot(3), solution_bar(3)
    real(dp) :: rhs_bar(3), lhs, rhs_inner
    complex(dp) :: complex_rhs(2), complex_solution(2), complex_rhs_dot(2)
    complex(dp) :: complex_solution_dot(2), complex_solution_bar(2), complex_rhs_bar(2)
    real(dp) :: complex_lhs, complex_rhs_inner
    integer :: status
    logical :: all_passed

    all_passed = .true.
    call initialize_real_matrices(matrices, matrices_dot)
    rhs = [6.0_dp, 7.0_dp, 6.0_dp]
    rhs_dot = [0.3_dp, -0.4_dp, 0.2_dp]
    call factor_retained_field_split(matrices, split, status)
    call record_condition(status == 0, "real retained field split factors all blocks")
    call apply_retained_field_split(split, rhs, solution, status)
    call record_condition(status == 0 .and. maxval(abs(solution - [1.0_dp, 2.0_dp, 3.0_dp])) < &
        1.0e-13_dp, "real retained field split solves each independent field")
    call apply_retained_field_split_jvp( &
        split, matrices_dot, solution, rhs_dot, solution_dot, status)
    call record_condition(status == 0 .and. maxval(abs(solution_dot - &
        [1.95_dp/11.0_dp, -4.5_dp/11.0_dp, -0.125_dp])) < 1.0e-13_dp, &
        "real retained field split JVP matches the block solve oracle")
    solution_bar = [0.4_dp, -0.3_dp, 0.2_dp]
    call apply_retained_field_split_vjp( &
        split, solution, solution_bar, rhs_bar, matrices_bar, status)
    lhs = dot_product(solution_bar, solution_dot)
    rhs_inner = dot_product(rhs_bar, rhs_dot) + dot_product(matrices_bar(1)%val, matrices_dot(1)%val) + &
        dot_product(matrices_bar(2)%val, matrices_dot(2)%val)
    call record_condition(status == 0 .and. abs(lhs - rhs_inner) < 1.0e-12_dp, &
        "real retained field split VJP satisfies the adjoint identity")
    call free_retained_field_split(split)

    call initialize_complex_matrices(complex_matrices, complex_matrices_dot)
    complex_rhs = [ &
        complex_matrices(1)%val(1)*cmplx(1.0_dp, 0.1_dp, dp), &
        complex_matrices(2)%val(1)*cmplx(2.0_dp, -0.3_dp, dp)]
    complex_rhs_dot = [cmplx(0.03_dp, -0.01_dp, dp), cmplx(-0.02_dp, 0.02_dp, dp)]
    call factor_retained_complex_field_split(complex_matrices, complex_split, status)
    call record_condition(status == 0, "complex retained field split factors all blocks")
    call apply_retained_complex_field_split(complex_split, complex_rhs, complex_solution, status)
    call record_condition(status == 0 .and. maxval(abs(complex_solution - &
        [cmplx(1.0_dp, 0.1_dp, dp), cmplx(2.0_dp, -0.3_dp, dp)])) < 1.0e-13_dp, &
        "complex retained field split solves each independent field")
    call apply_retained_complex_field_split_jvp( &
        complex_split, complex_matrices_dot, complex_solution, complex_rhs_dot, &
        complex_solution_dot, status)
    complex_solution_bar = [cmplx(0.4_dp, 0.1_dp, dp), cmplx(-0.3_dp, 0.2_dp, dp)]
    call apply_retained_complex_field_split_vjp( &
        complex_split, complex_solution, complex_solution_bar, complex_rhs_bar, &
        complex_matrices_bar, status)
    complex_lhs = real(sum(conjg(complex_solution_bar)*complex_solution_dot), dp)
    complex_rhs_inner = real(sum(conjg(complex_rhs_bar)*complex_rhs_dot) + &
        sum(conjg(complex_matrices_bar(1)%val)*complex_matrices_dot(1)%val) + &
        sum(conjg(complex_matrices_bar(2)%val)*complex_matrices_dot(2)%val), dp)
    call record_condition(status == 0 .and. abs(complex_lhs - complex_rhs_inner) < 1.0e-12_dp, &
        "complex retained field split VJP satisfies the real adjoint identity")
    call free_retained_complex_field_split(complex_split)

    call check_summary("retained field split")
    if (.not. all_passed) error stop 1

contains

    subroutine initialize_real_matrices(matrices, matrices_dot)
        type(csc_t), intent(out) :: matrices(:), matrices_dot(:)

        matrices(1)%nrow = 2
        matrices(1)%ncol = 2
        matrices(1)%nnz = 4
        matrices(1)%col_ptr = [1, 3, 5]
        matrices(1)%row_idx = [1, 2, 1, 2]
        matrices(1)%val = [4.0_dp, 1.0_dp, 1.0_dp, 3.0_dp]
        matrices(2)%nrow = 1
        matrices(2)%ncol = 1
        matrices(2)%nnz = 1
        matrices(2)%col_ptr = [1, 2]
        matrices(2)%row_idx = [1]
        matrices(2)%val = [2.0_dp]
        matrices_dot = matrices
        matrices_dot(1)%val = [0.2_dp, 0.05_dp, -0.1_dp, 0.3_dp]
        matrices_dot(2)%val = [0.15_dp]
    end subroutine initialize_real_matrices

    subroutine initialize_complex_matrices(matrices, matrices_dot)
        type(csc_z_t), intent(out) :: matrices(:), matrices_dot(:)

        matrices(1)%nrow = 1
        matrices(1)%ncol = 1
        matrices(1)%nnz = 1
        matrices(1)%col_ptr = [1, 2]
        matrices(1)%row_idx = [1]
        matrices(1)%val = [cmplx(2.0_dp, 0.5_dp, dp)]
        matrices(2)%nrow = 1
        matrices(2)%ncol = 1
        matrices(2)%nnz = 1
        matrices(2)%col_ptr = [1, 2]
        matrices(2)%row_idx = [1]
        matrices(2)%val = [cmplx(3.0_dp, -0.2_dp, dp)]
        matrices_dot = matrices
        matrices_dot(1)%val = [cmplx(0.1_dp, -0.05_dp, dp)]
        matrices_dot(2)%val = [cmplx(0.2_dp, 0.03_dp, dp)]
    end subroutine initialize_complex_matrices

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_retained_field_split
