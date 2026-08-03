program test_fci_retained_cycle_composition
    use check, only: check_condition, check_summary
    use fortfem_feec, only: &
        evaluate_fci_retained_cycle_composition, &
        evaluate_fci_retained_cycle_composition_jvp, &
        evaluate_fci_retained_cycle_composition_vjp
    use fortfem_kinds, only: dp
    use fortfem_retained_field_split, only: &
        factor_retained_field_split, free_retained_field_split, &
        retained_field_split_t
    use fortsparse, only: csc_t, fortsparse_status_t
    implicit none

    type(csc_t) :: matrices(2), matrices_dot(2), matrices_bar(2)
    type(retained_field_split_t) :: split
    type(fortsparse_status_t) :: status
    real(dp), parameter :: residual(3) = [1.0_dp, 2.0_dp, 3.0_dp]
    real(dp), parameter :: parallel_action(3) = [0.6_dp, 0.3_dp, 1.2_dp]
    real(dp), parameter :: plane_action(3) = [0.2_dp, 0.4_dp, 0.1_dp]
    real(dp), parameter :: weights(3) = [0.5_dp, 0.25_dp, 0.75_dp]
    real(dp), parameter :: residual_dot(3) = [0.2_dp, -0.1_dp, 0.3_dp]
    real(dp), parameter :: parallel_action_dot(3) = [0.1_dp, 0.05_dp, -0.2_dp]
    real(dp), parameter :: plane_action_dot(3) = [-0.03_dp, 0.04_dp, 0.02_dp]
    real(dp), parameter :: weights_dot(3) = [0.02_dp, -0.01_dp, 0.03_dp]
    real(dp), parameter :: correction_bar(3) = [0.4_dp, -0.2_dp, 0.3_dp]
    real(dp), parameter :: work_bar(3) = [0.1_dp, -0.4_dp, 0.2_dp]
    real(dp), parameter :: total_bar = -0.15_dp
    real(dp) :: correction(3), work(3), total_work
    real(dp) :: correction_dot(3), work_dot(3), total_work_dot
    real(dp) :: expected_action(3), expected_action_dot(3)
    real(dp) :: expected_correction_dot(3), expected_work_dot(3)
    real(dp) :: residual_bar(3), parallel_action_bar(3), plane_action_bar(3)
    real(dp) :: weights_bar(3), lhs, rhs
    integer :: factor_status

    call initialize_matrices(matrices, matrices_dot)
    call factor_retained_field_split(matrices, split, factor_status)
    call check_condition(factor_status == 0, &
        "FCI retained cycle factors nonsymmetric caller blocks")

    call evaluate_fci_retained_cycle_composition( &
        split, residual, parallel_action, plane_action, weights, correction, &
        work, total_work, status)
    expected_action = solve_oracle(residual, matrices)
    call check_condition(status%code == 0 .and. &
        maxval(abs(expected_action - [0.1_dp, 0.6_dp, 1.5_dp])) < 1.0e-14_dp, &
        "independent retained-action oracle is the literal block solve")
    call check_condition( &
        maxval(abs(correction - [0.425_dp, 0.7_dp, 1.75_dp])) < 1.0e-14_dp .and. &
        maxval(abs(work - [2.4_dp, 0.325_dp, 4.35_dp])) < 1.0e-14_dp .and. &
        abs(total_work - 7.075_dp) < 1.0e-14_dp, &
        "FCI retained cycle matches the independent weighted-work oracle")

    call evaluate_fci_retained_cycle_composition_jvp( &
        split, matrices_dot, residual, parallel_action, plane_action, weights, &
        residual_dot, parallel_action_dot, plane_action_dot, weights_dot, &
        correction_dot, work_dot, total_work_dot, status)
    expected_action_dot = solve_oracle( &
        residual_dot - apply_oracle(matrices_dot, expected_action), matrices)
    expected_correction_dot = &
        weights_dot(1)*parallel_action + weights(1)*parallel_action_dot + &
        weights_dot(2)*plane_action + weights(2)*plane_action_dot + &
        weights_dot(3)*expected_action + weights(3)*expected_action_dot
    expected_work_dot(1) = weights_dot(1)*dot_product(residual, parallel_action) + &
        weights(1)*(dot_product(residual_dot, parallel_action) + &
        dot_product(residual, parallel_action_dot))
    expected_work_dot(2) = weights_dot(2)*dot_product(residual, plane_action) + &
        weights(2)*(dot_product(residual_dot, plane_action) + &
        dot_product(residual, plane_action_dot))
    expected_work_dot(3) = weights_dot(3)*dot_product(residual, expected_action) + &
        weights(3)*(dot_product(residual_dot, expected_action) + &
        dot_product(residual, expected_action_dot))
    call check_condition(status%code == 0 .and. &
        maxval(abs(correction_dot - expected_correction_dot)) < 2.0e-14_dp .and. &
        maxval(abs(work_dot - expected_work_dot)) < 2.0e-14_dp .and. &
        abs(total_work_dot - sum(expected_work_dot)) < 2.0e-14_dp, &
        "FCI retained-cycle JVP matches an independent inverse derivative")

    call evaluate_fci_retained_cycle_composition_vjp( &
        split, residual, parallel_action, plane_action, weights, correction_bar, &
        work_bar, total_bar, residual_bar, parallel_action_bar, plane_action_bar, &
        matrices_bar, weights_bar, status)
    lhs = dot_product(correction_bar, correction_dot) + &
        dot_product(work_bar, work_dot) + total_bar*total_work_dot
    rhs = dot_product(residual_bar, residual_dot) + &
        dot_product(parallel_action_bar, parallel_action_dot) + &
        dot_product(plane_action_bar, plane_action_dot) + &
        dot_product(weights_bar, weights_dot) + &
        dot_product(matrices_bar(1)%val, matrices_dot(1)%val) + &
        dot_product(matrices_bar(2)%val, matrices_dot(2)%val)
    call check_condition(status%code == 0 .and. abs(lhs - rhs) < 2.0e-13_dp, &
        "FCI retained-cycle VJP satisfies the nonsymmetric adjoint identity")

    call evaluate_fci_retained_cycle_composition( &
        split, residual, parallel_action, plane_action, [0.5_dp, 0.0_dp, 0.75_dp], &
        correction, work, total_work, status)
    call check_condition(status%code /= 0 .and. maxval(abs(correction)) == 0.0_dp, &
        "FCI retained cycle rejects a non-positive component weight")

    call free_retained_field_split(split)
    call check_summary("FCI retained cycle composition")

contains

    subroutine initialize_matrices(matrices, matrices_dot)
        type(csc_t), intent(out) :: matrices(:), matrices_dot(:)

        matrices(1)%nrow = 2
        matrices(1)%ncol = 2
        matrices(1)%nnz = 4
        matrices(1)%col_ptr = [1, 3, 5]
        matrices(1)%row_idx = [1, 2, 1, 2]
        matrices(1)%val = [4.0_dp, 2.0_dp, 1.0_dp, 3.0_dp]
        matrices(2)%nrow = 1
        matrices(2)%ncol = 1
        matrices(2)%nnz = 1
        matrices(2)%col_ptr = [1, 2]
        matrices(2)%row_idx = [1]
        matrices(2)%val = [2.0_dp]
        matrices_dot = matrices
        matrices_dot(1)%val = [0.1_dp, -0.05_dp, 0.02_dp, 0.08_dp]
        matrices_dot(2)%val = [0.04_dp]
    end subroutine initialize_matrices

    pure function solve_oracle(rhs, matrices) result(solution)
        real(dp), intent(in) :: rhs(:)
        type(csc_t), intent(in) :: matrices(:)
        real(dp) :: solution(3), determinant

        determinant = matrices(1)%val(1)*matrices(1)%val(4) - &
            matrices(1)%val(3)*matrices(1)%val(2)
        solution(1) = (matrices(1)%val(4)*rhs(1) - &
            matrices(1)%val(3)*rhs(2))/determinant
        solution(2) = (-matrices(1)%val(2)*rhs(1) + &
            matrices(1)%val(1)*rhs(2))/determinant
        solution(3) = rhs(3)/matrices(2)%val(1)
    end function solve_oracle

    pure function apply_oracle(matrices, field) result(product)
        type(csc_t), intent(in) :: matrices(:)
        real(dp), intent(in) :: field(:)
        real(dp) :: product(3)

        product(1) = matrices(1)%val(1)*field(1) + &
            matrices(1)%val(3)*field(2)
        product(2) = matrices(1)%val(2)*field(1) + &
            matrices(1)%val(4)*field(2)
        product(3) = matrices(2)%val(1)*field(3)
    end function apply_oracle

end program test_fci_retained_cycle_composition
