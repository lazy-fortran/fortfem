program test_enrichment_support_diagnostics
    use check, only: check_condition, check_summary
    use fortfem_api, only: evaluate_enrichment_support_gram, &
        evaluate_enrichment_support_gram_jvp, &
        evaluate_enrichment_support_gram_vjp, &
        evaluate_enrichment_support_rank_condition
    use fortfem_kinds, only: dp
    use fortsparse, only: fortsparse_status_t
    implicit none

    real(dp), parameter :: basis_values(4, 3) = reshape([ &
        1.0_dp, 0.0_dp, 1.0_dp, 2.0_dp, &
        0.0_dp, 1.0_dp, 1.0_dp, -1.0_dp, &
        1.0_dp, 1.0_dp, 2.0_dp, 1.0_dp], [4, 3])
    real(dp), parameter :: sample_weights(4) = [0.5_dp, 1.0_dp, 2.0_dp, 0.25_dp]
    logical, parameter :: active_mask(3) = [.true., .true., .false.]
    real(dp), parameter :: basis_dot(4, 3) = reshape([ &
        0.1_dp, -0.2_dp, 0.3_dp, 0.4_dp, &
        -0.3_dp, 0.2_dp, -0.1_dp, 0.5_dp, &
        0.6_dp, -0.4_dp, 0.2_dp, -0.2_dp], [4, 3])
    real(dp), parameter :: weights_dot(4) = [0.2_dp, -0.1_dp, 0.4_dp, 0.3_dp]
    real(dp), parameter :: gram_bar(3, 3) = reshape([ &
        0.2_dp, -0.4_dp, 0.1_dp, &
        0.3_dp, 0.5_dp, -0.2_dp, &
        -0.1_dp, 0.6_dp, 0.7_dp], [3, 3])
    real(dp) :: gram(3, 3), expected_gram(3, 3), gram_dot(3, 3)
    real(dp) :: basis_bar(4, 3), weights_bar(4)
    real(dp) :: gram_plus(3, 3), gram_minus(3, 3), finite_difference(3, 3)
    real(dp) :: lhs, rhs, epsilon_fd
    real(dp) :: eigenvalues(2), singular_eigenvalues(3)
    real(dp) :: expected_largest, expected_smallest
    real(dp) :: trace_value, discriminant, condition_estimate
    integer :: rank, i, j, point
    type(fortsparse_status_t) :: status
    logical :: all_passed

    all_passed = .true.
    epsilon_fd = 1.0e-7_dp

    expected_gram = 0.0_dp
    do i = 1, size(basis_values, 2)
        do j = 1, size(basis_values, 2)
            if (.not. active_mask(i) .or. .not. active_mask(j)) cycle
            do point = 1, size(basis_values, 1)
                expected_gram(i, j) = expected_gram(i, j) + &
                    sample_weights(point)*basis_values(point, i)* &
                    basis_values(point, j)
            end do
        end do
    end do

    call evaluate_enrichment_support_gram( &
        basis_values, sample_weights, active_mask, gram, status)
    call record(status%code == 0 .and. &
        maxval(abs(gram - expected_gram)) < 1.0e-14_dp, &
        "support Gram matrix matches the independent weighted oracle")

    call evaluate_enrichment_support_gram_jvp( &
        basis_values, sample_weights, active_mask, basis_dot, weights_dot, &
        gram_dot, status)
    call evaluate_enrichment_support_gram( &
        basis_values + epsilon_fd*basis_dot, &
        sample_weights + epsilon_fd*weights_dot, active_mask, gram_plus, status)
    call evaluate_enrichment_support_gram( &
        basis_values - epsilon_fd*basis_dot, &
        sample_weights - epsilon_fd*weights_dot, active_mask, gram_minus, status)
    finite_difference = (gram_plus - gram_minus)/(2.0_dp*epsilon_fd)
    call record(status%code == 0 .and. &
        maxval(abs(gram_dot - finite_difference)) < 2.0e-8_dp, &
        "support Gram JVP matches the fixed-activation finite difference")

    call evaluate_enrichment_support_gram_vjp( &
        basis_values, sample_weights, active_mask, gram_bar, basis_bar, &
        weights_bar, status)
    lhs = sum(gram_bar*gram_dot)
    rhs = sum(basis_bar*basis_dot) + sum(weights_bar*weights_dot)
    call record(status%code == 0 .and. abs(lhs - rhs) < 1.0e-13_dp, &
        "support Gram VJP satisfies the real dot-product identity")

    call evaluate_enrichment_support_rank_condition( &
        gram(1:2, 1:2), [.true., .true.], 1.0e-12_dp, rank, &
        condition_estimate, eigenvalues, status)
    trace_value = gram(1, 1) + gram(2, 2)
    discriminant = sqrt((gram(1, 1) - gram(2, 2))**2 + &
        4.0_dp*gram(1, 2)**2)
    expected_largest = 0.5_dp*(trace_value + discriminant)
    expected_smallest = 0.5_dp*(trace_value - discriminant)
    call record(status%code == 0 .and. rank == 2 .and. &
        abs(eigenvalues(1) - expected_largest) < 1.0e-12_dp .and. &
        abs(eigenvalues(2) - expected_smallest) < 1.0e-12_dp .and. &
        abs(condition_estimate - sqrt(expected_largest/expected_smallest)) < &
        1.0e-12_dp, "full-rank support condition estimate has an eigenvalue oracle")

    call evaluate_enrichment_support_rank_condition( &
        expected_gram, [ .true., .true., .true. ], 1.0e-12_dp, rank, &
        condition_estimate, singular_eigenvalues, status)
    call record(status%code == 0 .and. rank == 2 .and. &
        condition_estimate >= huge(1.0_dp), &
        "rank diagnostic marks the duplicated enrichment as singular")

    call evaluate_enrichment_support_gram( &
        basis_values, [-1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp], active_mask, gram, status)
    call record(status%code /= 0, "support Gram rejects negative quadrature weights")

    call check_summary("XFEM enrichment support diagnostics")
    if (.not. all_passed) error stop 1

contains

    subroutine record(condition, message)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: message

        call check_condition(condition, message)
        all_passed = all_passed .and. condition
    end subroutine record

end program test_enrichment_support_diagnostics
