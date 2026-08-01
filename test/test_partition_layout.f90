program test_partition_layout
    use check, only: check_condition, check_summary
    use fortfem_api, only: assemble_partitioned_sum, assemble_partitioned_sum_jvp, &
        assemble_partitioned_sum_vjp, initialize_partition_layout, &
        partition_layout_ghost_count, partition_layout_global_count, &
        partition_layout_maps, partition_layout_owned_count, partition_layout_t
    use fortfem_kinds, only: dp
    use fortsparse, only: fortsparse_status_t
    implicit none

    integer, parameter :: global_count = 5, local_count = 3
    integer, parameter :: local_to_global(local_count) = [1, 3, 5]
    integer, parameter :: owner_rank(local_count) = [0, 1, 0]
    logical, parameter :: owned(local_count) = [.true., .false., .true.]
    real(dp), parameter :: epsilon_fd = 1.0e-7_dp
    real(dp), parameter :: local_values(local_count) = [2.0_dp, -4.0_dp, 3.0_dp]
    real(dp), parameter :: local_values_dot(local_count) = [-0.5_dp, 0.7_dp, 1.1_dp]
    real(dp), parameter :: expected_global(global_count) = [2.0_dp, 0.0_dp, -4.0_dp, 0.0_dp, 3.0_dp]
    real(dp) :: global_values(global_count), global_values_dot(global_count)
    real(dp) :: global_plus(global_count), global_minus(global_count)
    real(dp) :: local_plus(local_count), local_minus(local_count)
    real(dp) :: global_values_bar(global_count), local_values_bar(local_count)
    logical, allocatable :: owned_copy(:)
    integer, allocatable :: local_ids_copy(:), owners_copy(:)
    integer :: rank_copy
    integer, parameter :: bad_ids(2) = [1, 1]
    integer, parameter :: bad_owners(2) = [0, 0]
    logical, parameter :: bad_owned(2) = [.true., .false.]
    type(partition_layout_t) :: layout, bad_layout
    type(fortsparse_status_t) :: status
    real(dp) :: lhs, rhs
    logical :: all_passed

    all_passed = .true.
    call initialize_partition_layout( &
        layout, global_count, local_to_global, owner_rank, owned, 0, status)
    call record_condition(status%code == 0 .and. &
        partition_layout_global_count(layout) == global_count .and. &
        partition_layout_owned_count(layout) == 2 .and. &
        partition_layout_ghost_count(layout) == 1, &
        "partition layout records global, owned, and ghost counts")

    call partition_layout_maps( &
        layout, local_ids_copy, owners_copy, owned_copy, rank_copy, status)
    call record_condition(status%code == 0 .and. rank_copy == 0 .and. &
        all(local_ids_copy == local_to_global) .and. &
        all(owners_copy == owner_rank) .and. all(owned_copy .eqv. owned), &
        "partition layout exports deterministic owner and ghost metadata")

    call assemble_partitioned_sum(layout, local_values, global_values, status)
    call record_condition(status%code == 0 .and. &
        all(global_values == expected_global), &
        "serial partition reduction accumulates by global ID")
    call assemble_partitioned_sum_jvp(layout, local_values_dot, global_values_dot, status)
    local_plus = local_values + epsilon_fd*local_values_dot
    local_minus = local_values - epsilon_fd*local_values_dot
    call assemble_partitioned_sum(layout, local_plus, global_plus, status)
    call assemble_partitioned_sum(layout, local_minus, global_minus, status)
    call record_condition(status%code == 0 .and. &
        maxval(abs(global_values_dot - &
        (global_plus - global_minus)/(2.0_dp*epsilon_fd))) < 1.0e-8_dp, &
        "partition reduction JVP matches independent finite differences")

    global_values_bar = [0.3_dp, -0.2_dp, 0.5_dp, -0.7_dp, 0.9_dp]
    call assemble_partitioned_sum_vjp( &
        layout, global_values_bar, local_values_bar, status)
    lhs = sum(global_values_bar*global_values_dot)
    rhs = sum(local_values_bar*local_values_dot)
    call record_condition(status%code == 0 .and. abs(lhs - rhs) < 1.0e-14_dp, &
        "partition reduction VJP satisfies the real dot-product identity")

    call initialize_partition_layout( &
        bad_layout, global_count, bad_ids, bad_owners, bad_owned, 0, status)
    call record_condition(status%code /= 0, &
        "partition layout rejects duplicate local global IDs")
    if (allocated(local_ids_copy)) deallocate(local_ids_copy)
    if (allocated(owners_copy)) deallocate(owners_copy)
    if (allocated(owned_copy)) deallocate(owned_copy)
    call check_summary("partition layout")
    if (.not. all_passed) error stop 1

contains

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_partition_layout
