program test_distributed_trace_ownership
    use check, only: check_condition, check_summary
    use fortfem_feec, only: &
        assemble_distributed_trace_reduction, &
        assemble_distributed_trace_reduction_jvp, &
        assemble_distributed_trace_reduction_vjp, &
        distributed_trace_layout_component_count, &
        distributed_trace_layout_global_count, &
        distributed_trace_layout_local_count, &
        distributed_trace_layout_partition_count, &
        distributed_trace_layout_t, initialize_distributed_trace_layout, &
        initialize_partition_layout, partition_layout_t
    use fortfem_kinds, only: dp
    use fortsparse, only: fortsparse_status_t
    implicit none

    integer, parameter :: global_count = 4, component_count = 2
    integer, parameter :: partition_count = 2, local_count = 6
    integer, parameter :: local_to_global_0(3) = [1, 2, 4]
    integer, parameter :: owner_rank_0(3) = [0, 1, 0]
    logical, parameter :: owned_0(3) = [.true., .false., .true.]
    integer, parameter :: local_to_global_1(3) = [2, 3, 4]
    integer, parameter :: owner_rank_1(3) = [1, 1, 0]
    logical, parameter :: owned_1(3) = [.true., .true., .false.]
    real(dp), parameter :: epsilon_fd = 1.0e-5_dp
    real(dp), parameter :: local_values(local_count, component_count) = &
        reshape([ &
        1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp, 6.0_dp, &
        -1.0_dp, 2.0_dp, -3.0_dp, 4.0_dp, -5.0_dp, 6.0_dp], &
        [local_count, component_count])
    real(dp), parameter :: local_values_dot(local_count, component_count) = &
        reshape([ &
        -0.2_dp, 0.3_dp, -0.4_dp, 0.5_dp, -0.6_dp, 0.7_dp, &
        0.8_dp, -0.9_dp, 1.0_dp, -1.1_dp, 1.2_dp, -1.3_dp], &
        [local_count, component_count])
    type(partition_layout_t) :: partitions(partition_count)
    type(partition_layout_t) :: bad_partitions(partition_count)
    type(partition_layout_t) :: incomplete_partitions(partition_count)
    type(partition_layout_t) :: bad_partition, mismatched_partition
    type(distributed_trace_layout_t) :: layout, bad_layout, incomplete_layout
    type(fortsparse_status_t) :: status
    real(dp) :: global_values(global_count, component_count)
    real(dp) :: owner_values(global_count, component_count)
    real(dp) :: global_values_dot(global_count, component_count)
    real(dp) :: owner_values_dot(global_count, component_count)
    real(dp) :: global_plus(global_count, component_count)
    real(dp) :: owner_plus(global_count, component_count)
    real(dp) :: global_minus(global_count, component_count)
    real(dp) :: owner_minus(global_count, component_count)
    real(dp) :: global_values_bar(global_count, component_count)
    real(dp) :: owner_values_bar(global_count, component_count)
    real(dp) :: local_values_bar(local_count, component_count)
    real(dp) :: expected_global(global_count, component_count)
    real(dp) :: expected_owner(global_count, component_count)
    real(dp) :: lhs, rhs
    integer :: row, component, local, global_id
    logical :: all_passed

    all_passed = .true.
    call initialize_partition_layout( &
        partitions(1), global_count, local_to_global_0, owner_rank_0, &
        owned_0, 0, status)
    call record_condition(status%code == 0, &
        "rank-zero partition layout initializes")
    call initialize_partition_layout( &
        partitions(2), global_count, local_to_global_1, owner_rank_1, &
        owned_1, 1, status)
    call record_condition(status%code == 0, &
        "rank-one partition layout initializes")

    call initialize_distributed_trace_layout( &
        layout, partitions, component_count, status)
    call record_condition(status%code == 0 .and. &
        distributed_trace_layout_partition_count(layout) == partition_count .and. &
        distributed_trace_layout_global_count(layout) == global_count .and. &
        distributed_trace_layout_local_count(layout) == local_count .and. &
        distributed_trace_layout_component_count(layout) == component_count, &
        "distributed trace layout records packed ownership dimensions")

    expected_global = 0.0_dp
    expected_owner = 0.0_dp
    do row = 1, local_count
        if (row <= 3) then
            local = row
            global_id = local_to_global_0(local)
            do component = 1, component_count
                expected_global(global_id, component) = &
                    expected_global(global_id, component) + &
                    local_values(row, component)
                if (owned_0(local)) expected_owner(global_id, component) = &
                    expected_owner(global_id, component) + &
                    local_values(row, component)
            end do
        else
            local = row - 3
            global_id = local_to_global_1(local)
            do component = 1, component_count
                expected_global(global_id, component) = &
                    expected_global(global_id, component) + &
                    local_values(row, component)
                if (owned_1(local)) expected_owner(global_id, component) = &
                    expected_owner(global_id, component) + &
                    local_values(row, component)
            end do
        end if
    end do
    call assemble_distributed_trace_reduction( &
        layout, local_values, global_values, owner_values, status)
    call record_condition(status%code == 0 .and. &
        maxval(abs(global_values - expected_global)) < 1.0e-14_dp .and. &
        maxval(abs(owner_values - expected_owner)) < 1.0e-14_dp, &
        "trace reduction accumulates duplicate ghosts and owner rows")

    call assemble_distributed_trace_reduction_jvp( &
        layout, local_values_dot, global_values_dot, owner_values_dot, status)
    call assemble_distributed_trace_reduction( &
        layout, local_values + epsilon_fd*local_values_dot, &
        global_plus, owner_plus, status)
    call assemble_distributed_trace_reduction( &
        layout, local_values - epsilon_fd*local_values_dot, &
        global_minus, owner_minus, status)
    call record_condition(status%code == 0 .and. &
        maxval(abs(global_values_dot - &
        (global_plus - global_minus)/(2.0_dp*epsilon_fd))) < 1.0e-8_dp .and. &
        maxval(abs(owner_values_dot - &
        (owner_plus - owner_minus)/(2.0_dp*epsilon_fd))) < 1.0e-8_dp, &
        "trace reduction JVP matches an independent finite difference")

    global_values_bar = reshape([ &
        0.2_dp, -0.3_dp, 0.4_dp, -0.5_dp, 0.6_dp, -0.7_dp, 0.8_dp, -0.9_dp], &
        [global_count, component_count])
    owner_values_bar = reshape([ &
        -1.0_dp, 1.1_dp, -1.2_dp, 1.3_dp, -1.4_dp, 1.5_dp, -1.6_dp, 1.7_dp], &
        [global_count, component_count])
    call assemble_distributed_trace_reduction_vjp( &
        layout, global_values_bar, owner_values_bar, local_values_bar, status)
    lhs = sum(global_values_bar*global_values_dot) + &
        sum(owner_values_bar*owner_values_dot)
    rhs = sum(local_values_bar*local_values_dot)
    call record_condition(status%code == 0 .and. abs(lhs - rhs) < 1.0e-14_dp, &
        "trace reduction VJP satisfies the real transpose identity")

    call initialize_partition_layout( &
        bad_partition, global_count, [0, 2], [0, 1], [.true., .false.], 0, status)
    call record_condition(status%code /= 0, &
        "partition layouts reject invalid local-to-global IDs")
    call initialize_partition_layout( &
        bad_partition, global_count, [1, 2], [1, 1], [.true., .false.], 0, status)
    call record_condition(status%code /= 0, &
        "partition layouts reject an owned row with another owner")

    call initialize_partition_layout( &
        mismatched_partition, global_count, local_to_global_1, [1, 1, 2], &
        owned_1, 1, status)
    call record_condition(status%code == 0, &
        "a locally valid ghost owner is accepted before global validation")
    bad_partitions(1) = partitions(1)
    bad_partitions(2) = mismatched_partition
    call initialize_distributed_trace_layout( &
        bad_layout, bad_partitions, component_count, status)
    call record_condition(status%code /= 0, &
        "distributed trace layout rejects inconsistent ghost owner IDs")

    call initialize_partition_layout( &
        incomplete_partitions(1), global_count + 1, local_to_global_0, &
        owner_rank_0, owned_0, 0, status)
    call initialize_partition_layout( &
        incomplete_partitions(2), global_count + 1, local_to_global_1, &
        owner_rank_1, owned_1, 1, status)
    call initialize_distributed_trace_layout( &
        incomplete_layout, incomplete_partitions, component_count, status)
    call record_condition(status%code /= 0, &
        "distributed trace layout rejects an omitted global row")

    call check_summary("distributed trace ownership")
    if (.not. all_passed) error stop 1

contains

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_distributed_trace_ownership
