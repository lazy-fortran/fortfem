program test_physical_trace_reconciliation
    use check, only: check_condition, check_summary
    use fortfem_core, only: &
        initialize_physical_trace_ownership, &
        initialize_physical_trace_reconciliation, &
        physical_trace_reconciliation_maps, &
        physical_trace_ownership_t, &
        physical_trace_reconciliation_t, &
        reconcile_physical_trace_values, &
        reconcile_physical_trace_values_jvp, &
        reconcile_physical_trace_values_vjp
    use fortfem_kinds, only: dp
    use fortsparse, only: fortsparse_status_t
    implicit none

    real(dp), parameter :: reference_coordinates(2, 3) = reshape([ &
        0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, 1.0_dp], [2, 3])
    real(dp), parameter :: candidate_coordinates(2, 3) = reshape([ &
        0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp], [2, 3])
    real(dp), parameter :: bad_coordinates(2, 3) = reshape([ &
        0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, 1.0e-3_dp], [2, 3])
    integer, parameter :: reference_ids(3) = [11, 12, 13]
    integer, parameter :: candidate_ids(3) = [13, 11, 12]
    integer, parameter :: owners(3) = [0, 0, 1]
    logical, parameter :: reference_owned(3) = [.true., .true., .false.]
    logical, parameter :: candidate_owned(3) = [.false., .true., .true.]
    integer, parameter :: reference_orientation(3) = [1, -1, 1]
    integer, parameter :: candidate_orientation(3) = [1, -1, -1]
    real(dp), parameter :: candidate_values(3, 2) = reshape([ &
        10.0_dp, 1.0_dp, 20.0_dp, 2.0_dp, 30.0_dp, 3.0_dp], [3, 2])
    real(dp), parameter :: candidate_values_dot(3, 2) = reshape([ &
        0.3_dp, -0.1_dp, 0.5_dp, 0.2_dp, -0.4_dp, 0.6_dp], [3, 2])
    real(dp), parameter :: reference_bar(3, 2) = reshape([ &
        0.7_dp, -0.2_dp, 0.4_dp, 0.6_dp, -0.3_dp, 0.9_dp], [3, 2])
    real(dp), parameter :: step = 1.0e-5_dp
    real(dp) :: reference_values(3, 2), reference_values_dot(3, 2)
    real(dp) :: reference_values_plus(3, 2), reference_values_minus(3, 2)
    real(dp) :: candidate_values_bar(3, 2), lhs, rhs
    integer, allocatable :: candidate_index(:), orientation(:)
    real(dp) :: maximum_distance
    type(physical_trace_ownership_t) :: reference, candidate, bad
    type(physical_trace_reconciliation_t) :: reconciliation
    type(fortsparse_status_t) :: status
    logical :: all_passed

    all_passed = .true.
    call initialize_physical_trace_ownership( &
        reference, reference_coordinates, reference_ids, owners, reference_owned, &
        0, 1.0e-8_dp, status)
    call initialize_physical_trace_ownership( &
        candidate, candidate_coordinates, candidate_ids, [1, 0, 0], candidate_owned, &
        0, 1.0e-8_dp, status)
    call initialize_physical_trace_reconciliation( &
        reconciliation, reference, candidate, reference_orientation, &
        candidate_orientation, 1.0e-8_dp, status)
    call record_condition(status%code == 0, &
        "physical trace reconciliation accepts reordered coordinates and orientations")
    call physical_trace_reconciliation_maps( &
        reconciliation, candidate_index, orientation, maximum_distance, status)
    call record_condition(status%code == 0 .and. all(candidate_index == [2, 3, 1]) .and. &
        all(orientation == [-1, 1, 1]) .and. maximum_distance < 1.0e-14_dp, &
        "reconciliation exposes deterministic ID permutation, signs, and distance")

    call reconcile_physical_trace_values( &
        reconciliation, candidate_values, reference_values, status)
    call record_condition(status%code == 0 .and. all(reference_values == reshape([ &
        -1.0_dp, 20.0_dp, 10.0_dp, -30.0_dp, 3.0_dp, 2.0_dp], [3, 2])), &
        "reconciliation applies orientation signs while restoring reference order")
    call reconcile_physical_trace_values_jvp( &
        reconciliation, candidate_values_dot, reference_values_dot, status)
    call reconcile_physical_trace_values( &
        reconciliation, candidate_values + step*candidate_values_dot, reference_values_plus, status)
    call reconcile_physical_trace_values( &
        reconciliation, candidate_values - step*candidate_values_dot, reference_values_minus, status)
    call record_condition(status%code == 0 .and. maxval(abs(reference_values_dot - &
        (reference_values_plus - reference_values_minus)/(2.0_dp*step))) < 1.0e-8_dp, &
        "reconciliation JVP matches an independent central difference")

    call reconcile_physical_trace_values_vjp( &
        reconciliation, reference_bar, candidate_values_bar, status)
    lhs = sum(reference_bar*reference_values_dot)
    rhs = sum(candidate_values_bar*candidate_values_dot)
    call record_condition(status%code == 0 .and. abs(lhs - rhs) < 1.0e-14_dp, &
        "reconciliation VJP satisfies the real dot-product identity")

    call initialize_physical_trace_ownership( &
        bad, bad_coordinates, candidate_ids, [1, 0, 0], candidate_owned, &
        0, 1.0e-8_dp, status)
    call initialize_physical_trace_reconciliation( &
        reconciliation, reference, bad, reference_orientation, candidate_orientation, &
        1.0e-8_dp, status)
    call record_condition(status%code /= 0, &
        "reconciliation rejects coordinate mismatch beyond tolerance")
    call initialize_physical_trace_reconciliation( &
        reconciliation, reference, candidate, [1, 0, 1], candidate_orientation, &
        1.0e-8_dp, status)
    call record_condition(status%code /= 0, &
        "reconciliation rejects orientation entries other than plus or minus one")

    call check_summary("physical trace reconciliation")
    if (.not. all_passed) error stop 1

contains

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_physical_trace_reconciliation
