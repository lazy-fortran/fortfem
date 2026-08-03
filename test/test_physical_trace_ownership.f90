program test_physical_trace_ownership
    use check, only: check_condition, check_summary
    use fortfem_core, only: &
        compare_physical_trace_coordinates, &
        initialize_physical_trace_ownership, &
        physical_trace_ownership_dimension, &
        physical_trace_ownership_maps, &
        physical_trace_ownership_point_count, &
        physical_trace_ownership_rank, &
        physical_trace_ownership_t
    use fortfem_kinds, only: dp
    use fortsparse, only: fortsparse_status_t
    implicit none

    real(dp), parameter :: coordinates(2, 3) = reshape([ &
        0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, 1.0_dp], [2, 3])
    real(dp), parameter :: reordered_coordinates(2, 3) = reshape([ &
        0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp], [2, 3])
    real(dp), parameter :: perturbed_coordinates(2, 3) = reshape([ &
        0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, 1.0e-8_dp], [2, 3])
    integer, parameter :: global_ids(3) = [11, 12, 13]
    integer, parameter :: reordered_ids(3) = [13, 11, 12]
    integer, parameter :: owner_rank(3) = [0, 0, 1]
    logical, parameter :: owned(3) = [.true., .true., .false.]
    logical :: all_passed
    integer :: mismatch_count, copied_rank
    real(dp) :: maximum_distance
    real(dp), allocatable :: copied_coordinates(:, :)
    integer, allocatable :: copied_ids(:), copied_owners(:)
    logical, allocatable :: copied_owned(:)
    type(physical_trace_ownership_t) :: reference, candidate, bad
    type(fortsparse_status_t) :: status

    all_passed = .true.
    call initialize_physical_trace_ownership( &
        reference, coordinates, global_ids, owner_rank, owned, 0, 1.0e-7_dp, status)
    call record_condition(status%code == 0 .and. &
        physical_trace_ownership_dimension(reference) == 2 .and. &
        physical_trace_ownership_point_count(reference) == 3 .and. &
        physical_trace_ownership_rank(reference) == 0, &
        "physical trace ownership stores dimensions and rank")

    call physical_trace_ownership_maps( &
        reference, copied_coordinates, copied_ids, copied_owners, copied_owned, &
        copied_rank, status)
    call record_condition(status%code == 0 .and. copied_rank == 0 .and. &
        all(copied_coordinates == coordinates) .and. &
        all(copied_ids == global_ids) .and. &
        all(copied_owners == owner_rank) .and. &
        all(copied_owned .eqv. owned), &
        "physical trace ownership exports deterministic copies")

    call initialize_physical_trace_ownership( &
        candidate, reordered_coordinates, reordered_ids, [1, 0, 0], &
        [.false., .true., .true.], 0, 1.0e-7_dp, status)
    call record_condition(status%code == 0, &
        "independently ordered physical trace partition initializes")
    call compare_physical_trace_coordinates( &
        reference, candidate, 1.0e-7_dp, mismatch_count, maximum_distance, status)
    call record_condition(status%code == 0 .and. mismatch_count == 0 .and. &
        maximum_distance < 1.0e-14_dp, &
        "coordinate comparison matches rows independent of local ordering")

    call initialize_physical_trace_ownership( &
        candidate, perturbed_coordinates, reordered_ids, [1, 0, 0], &
        [.false., .true., .true.], 0, 1.0e-7_dp, status)
    call compare_physical_trace_coordinates( &
        reference, candidate, 1.0e-9_dp, mismatch_count, maximum_distance, status)
    call record_condition(status%code == 0 .and. mismatch_count == 1 .and. &
        maximum_distance > 1.0e-9_dp, &
        "coordinate comparison reports a physical mismatch")

    call initialize_physical_trace_ownership( &
        bad, coordinates, [11, 11, 13], owner_rank, owned, 0, 1.0e-7_dp, status)
    call record_condition(status%code /= 0, &
        "physical trace ownership rejects duplicate global IDs")
    call initialize_physical_trace_ownership( &
        bad, coordinates, global_ids, [1, 0, 1], owned, 0, 1.0e-7_dp, status)
    call record_condition(status%code /= 0, &
        "physical trace ownership rejects inconsistent local ownership")

    call check_summary("physical trace ownership")
    if (.not. all_passed) error stop 1

contains

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_physical_trace_ownership
