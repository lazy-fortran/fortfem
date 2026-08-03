program test_enrichment_support_activation
    use check, only: check_condition, check_summary
    use fortfem_feec, only: &
        evaluate_enrichment_support_activation, &
        evaluate_enrichment_support_activation_jvp, &
        evaluate_enrichment_support_activation_vjp
    use fortfem_kinds, only: dp
    use fortsparse, only: fortsparse_status_t, FORTSPARSE_OK
    implicit none

    real(dp), parameter :: level_values(6) = [ &
        -2.0_dp, 4.0_dp, 3.0_dp, 2.0_dp, 3.0_dp, 4.0_dp]
    integer, parameter :: support_offsets(4) = [1, 4, 7, 9]
    integer, parameter :: support_nodes(8) = [1, 2, 3, 4, 5, 6, 1, 5]
    real(dp), parameter :: level_dot(6) = [ &
        0.3_dp, -0.2_dp, 0.1_dp, 0.4_dp, 0.2_dp, -0.5_dp]
    real(dp), parameter :: min_bar(3) = [0.2_dp, -0.4_dp, 0.3_dp]
    real(dp), parameter :: max_bar(3) = [-0.1_dp, 0.5_dp, -0.2_dp]
    real(dp), parameter :: margin_bar(3) = [0.7_dp, -0.6_dp, 0.4_dp]
    logical :: active_mask(3), expected_active(3)
    real(dp) :: support_min(3), support_max(3), margin(3)
    real(dp) :: expected_min(3), expected_max(3), expected_margin(3)
    real(dp) :: min_dot(3), max_dot(3), margin_dot(3)
    real(dp) :: min_plus(3), max_plus(3), margin_plus(3)
    real(dp) :: min_minus(3), max_minus(3), margin_minus(3)
    real(dp) :: level_bar(6)
    real(dp) :: tied_level_values(6)
    real(dp) :: lhs, rhs, epsilon_fd
    integer :: min_owner(3), max_owner(3), margin_branch(3)
    integer :: expected_min_owner(3), expected_max_owner(3)
    integer :: min_owner_plus(3), max_owner_plus(3), margin_branch_plus(3)
    integer :: min_owner_minus(3), max_owner_minus(3), margin_branch_minus(3)
    integer :: tied_min_owner(3), tied_max_owner(3), tied_margin_branch(3)
    integer :: invalid_min_owner(3), invalid_max_owner(3)
    integer :: basis, entry, first, last, owner
    type(fortsparse_status_t) :: status
    logical :: all_passed

    all_passed = .true.
    epsilon_fd = 1.0e-7_dp
    expected_active = [.true., .false., .true.]

    expected_min = 0.0_dp
    expected_max = 0.0_dp
    expected_margin = 0.0_dp
    expected_min_owner = 0
    expected_max_owner = 0
    do basis = 1, size(expected_active)
        first = support_offsets(basis)
        last = support_offsets(basis + 1) - 1
        expected_min(basis) = huge(1.0_dp)
        expected_max(basis) = -huge(1.0_dp)
        do entry = first, last
            owner = support_nodes(entry)
            if (level_values(owner) < expected_min(basis)) then
                expected_min(basis) = level_values(owner)
                expected_min_owner(basis) = owner
            end if
            if (level_values(owner) > expected_max(basis)) then
                expected_max(basis) = level_values(owner)
                expected_max_owner(basis) = owner
            end if
        end do
        expected_active(basis) = expected_min(basis) < 0.0_dp .and. &
            expected_max(basis) > 0.0_dp
        expected_margin(basis) = min(expected_max(basis), &
            -expected_min(basis))
    end do

    call evaluate_enrichment_support_activation( &
        level_values, support_offsets, support_nodes, active_mask, support_min, &
        support_max, margin, min_owner, max_owner, margin_branch, status)
    call record(status%code == FORTSPARSE_OK .and. &
        all(active_mask .eqv. expected_active) .and. &
        maxval(abs(support_min - expected_min)) < 1.0e-14_dp .and. &
        maxval(abs(support_max - expected_max)) < 1.0e-14_dp .and. &
        maxval(abs(margin - expected_margin)) < 1.0e-14_dp .and. &
        all(min_owner == expected_min_owner) .and. &
        all(max_owner == expected_max_owner), &
        "support activation matches the independent CSR sign oracle")

    call evaluate_enrichment_support_activation_jvp( &
        level_values, support_offsets, support_nodes, min_owner, max_owner, &
        margin_branch, level_dot, min_dot, max_dot, margin_dot, status)
    call evaluate_enrichment_support_activation( &
        level_values + epsilon_fd*level_dot, support_offsets, support_nodes, &
        active_mask, min_plus, max_plus, margin_plus, min_owner_plus, &
        max_owner_plus, margin_branch_plus, status)
    call evaluate_enrichment_support_activation( &
        level_values - epsilon_fd*level_dot, support_offsets, support_nodes, &
        active_mask, min_minus, max_minus, margin_minus, min_owner_minus, &
        max_owner_minus, margin_branch_minus, status)
    call record(status%code == FORTSPARSE_OK .and. &
        maxval(abs(min_dot - (min_plus - min_minus)/(2.0_dp*epsilon_fd))) < &
        1.0e-8_dp .and. &
        maxval(abs(max_dot - (max_plus - max_minus)/(2.0_dp*epsilon_fd))) < &
        1.0e-8_dp .and. &
        maxval(abs(margin_dot - (margin_plus - margin_minus)/ &
        (2.0_dp*epsilon_fd))) < 1.0e-8_dp, &
        "support activation JVP matches fixed-topology finite differences")

    call evaluate_enrichment_support_activation_vjp( &
        level_values, support_offsets, support_nodes, min_owner, max_owner, &
        margin_branch, min_bar, max_bar, margin_bar, level_bar, &
        status)
    lhs = sum(min_bar*min_dot) + sum(max_bar*max_dot) + &
        sum(margin_bar*margin_dot)
    rhs = sum(level_bar*level_dot)
    call record(status%code == FORTSPARSE_OK .and. abs(lhs - rhs) < 1.0e-13_dp, &
        "support activation VJP satisfies the real dot-product identity")

    call evaluate_enrichment_support_activation( &
        [0.0_dp, level_values(2:)], support_offsets, support_nodes, &
        active_mask, support_min, support_max, margin, min_owner, max_owner, &
        margin_branch, status)
    call record(status%code /= FORTSPARSE_OK, &
        "support activation rejects a level-set topology event")

    call evaluate_enrichment_support_activation( &
        level_values, [1, 4, 7, 10], support_nodes, active_mask, support_min, &
        support_max, margin, min_owner, max_owner, margin_branch, status)
    call record(status%code /= FORTSPARSE_OK, &
        "support activation rejects an out-of-range CSR offset")

    tied_level_values = level_values
    tied_level_values(5) = 2.0_dp
    call evaluate_enrichment_support_activation( &
        tied_level_values, support_offsets, support_nodes, active_mask, &
        support_min, support_max, margin, tied_min_owner, tied_max_owner, &
        tied_margin_branch, status)
    call evaluate_enrichment_support_activation_jvp( &
        tied_level_values, support_offsets, support_nodes, tied_min_owner, &
        tied_max_owner, tied_margin_branch, level_dot, min_dot, max_dot, &
        margin_dot, status)
    call record(status%code /= FORTSPARSE_OK, &
        "support activation JVP rejects a nondifferentiable margin tie")

    invalid_min_owner = min_owner
    invalid_min_owner(1) = 0
    call evaluate_enrichment_support_activation_jvp( &
        level_values, support_offsets, support_nodes, invalid_min_owner, &
        max_owner, margin_branch, level_dot, min_dot, max_dot, margin_dot, &
        status)
    call record(status%code /= FORTSPARSE_OK, &
        "support activation JVP rejects an out-of-range owner")

    invalid_max_owner = max_owner
    invalid_max_owner(2) = size(level_values) + 1
    call evaluate_enrichment_support_activation_vjp( &
        level_values, support_offsets, support_nodes, min_owner, &
        invalid_max_owner, margin_branch, min_bar, max_bar, margin_bar, &
        level_bar, status)
    call record(status%code /= FORTSPARSE_OK, &
        "support activation VJP rejects an out-of-range owner")

    call check_summary("XFEM enrichment support activation")
    if (.not. all_passed) error stop 1

contains

    subroutine record(condition, message)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: message

        call check_condition(condition, message)
        all_passed = all_passed .and. condition
    end subroutine record

end program test_enrichment_support_activation
