program test_mpi_trace_exchange
    use check, only: check_condition, check_summary
    use fortfem_mpi_trace_exchange, only: &
        initialize_mpi_trace_exchange_schedule, &
        pack_mpi_trace_exchange, pack_mpi_trace_exchange_jvp, &
        pack_mpi_trace_exchange_vjp, unpack_mpi_trace_exchange, &
        unpack_mpi_trace_exchange_jvp, unpack_mpi_trace_exchange_vjp, &
        mpi_trace_exchange_schedule_t
    use fortfem_kinds, only: dp
    use fortfem_partition_layout, only: &
        initialize_partition_layout, partition_layout_t
    use fortsparse, only: fortsparse_status_t
    implicit none

    integer, parameter :: global_count = 3, component_count = 2
    integer, parameter :: partition_count = 2, local_count = 5
    integer, parameter :: owner_count = 3, ghost_count = 2
    integer, parameter :: offsets(3) = [1, 4, 6]
    real(dp), parameter :: step = 1.0e-7_dp
    type(partition_layout_t) :: partitions(partition_count)
    type(partition_layout_t) :: bad_partitions(partition_count)
    type(mpi_trace_exchange_schedule_t) :: schedule
    integer, parameter :: local_to_global_0(3) = [1, 2, 3]
    integer, parameter :: owner_rank_0(3) = [0, 0, 1]
    logical, parameter :: owned_0(3) = [.true., .true., .false.]
    integer, parameter :: local_to_global_1(2) = [2, 3]
    integer, parameter :: owner_rank_1(2) = [0, 1]
    logical, parameter :: owned_1(2) = [.false., .true.]
    real(dp) :: local_values(local_count, component_count)
    real(dp) :: local_values_dot(local_count, component_count)
    real(dp) :: owner_values(owner_count, component_count), ghost_values(ghost_count, component_count)
    real(dp) :: owner_values_dot(owner_count, component_count), ghost_values_dot(ghost_count, component_count)
    real(dp) :: owner_values_plus(owner_count, component_count), ghost_values_plus(ghost_count, component_count)
    real(dp) :: owner_values_minus(owner_count, component_count), ghost_values_minus(ghost_count, component_count)
    real(dp) :: owner_bar(owner_count, component_count), ghost_bar(ghost_count, component_count)
    real(dp) :: local_bar(local_count, component_count)
    real(dp) :: received_ghost(ghost_count, component_count), local_unpacked(local_count, component_count)
    real(dp) :: received_ghost_dot(ghost_count, component_count), local_unpacked_dot(local_count, component_count)
    real(dp) :: local_unpacked_plus(local_count, component_count), local_unpacked_minus(local_count, component_count)
    real(dp) :: owner_unpacked_bar(owner_count, component_count), ghost_unpacked_bar(ghost_count, component_count)
    real(dp) :: lhs, rhs
    type(fortsparse_status_t) :: status
    logical :: all_passed

    all_passed = .true.
    call initialize_partition_layout( &
        partitions(1), global_count, local_to_global_0, owner_rank_0, owned_0, 0, status)
    call initialize_partition_layout( &
        partitions(2), global_count, local_to_global_1, owner_rank_1, owned_1, 1, status)
    call initialize_mpi_trace_exchange_schedule( &
        schedule, partitions, component_count, status)
    call record_condition(status%code == 0, &
        "MPI-ready trace schedule accepts complete owner/ghost partitions")

    local_values = reshape([1.0_dp, 2.0_dp, 3.0_dp, 0.5_dp, 0.6_dp, &
        0.7_dp, 0.8_dp, 0.9_dp, 1.1_dp, 1.2_dp], shape(local_values))
    local_values_dot = reshape([0.2_dp, -0.1_dp, 0.4_dp, 0.3_dp, -0.5_dp, &
        0.6_dp, 0.1_dp, -0.2_dp, 0.7_dp, -0.4_dp], shape(local_values_dot))
    call pack_mpi_trace_exchange(schedule, offsets, local_values, owner_values, ghost_values, status)
    call record_condition(status%code == 0 .and. &
        maxval(abs(owner_values - reshape([1.0_dp, 2.0_dp, 0.6_dp, 0.7_dp, 0.8_dp, 1.2_dp], [owner_count, component_count]))) < 1.0e-14_dp .and. &
        maxval(abs(ghost_values - reshape([3.0_dp, 0.5_dp, 0.9_dp, 1.1_dp], [ghost_count, component_count]))) < 1.0e-14_dp, &
        "trace packing follows deterministic owner/ghost row order")

    call pack_mpi_trace_exchange_jvp( &
        schedule, offsets, local_values_dot, owner_values_dot, ghost_values_dot, status)
    call pack_mpi_trace_exchange( &
        schedule, offsets, local_values + step*local_values_dot, owner_values_plus, ghost_values_plus, status)
    call pack_mpi_trace_exchange( &
        schedule, offsets, local_values - step*local_values_dot, owner_values_minus, ghost_values_minus, status)
    call record_condition(status%code == 0 .and. &
        maxval(abs(owner_values_dot - (owner_values_plus - owner_values_minus)/(2.0_dp*step))) < 2.0e-8_dp .and. &
        maxval(abs(ghost_values_dot - (ghost_values_plus - ghost_values_minus)/(2.0_dp*step))) < 2.0e-8_dp, &
        "trace packing JVP matches an independent central difference")
    owner_bar = reshape([0.4_dp, -0.2_dp, 0.7_dp, 0.3_dp, -0.5_dp, 0.1_dp], shape(owner_bar))
    ghost_bar = reshape([0.6_dp, -0.4_dp, 0.2_dp, 0.9_dp], shape(ghost_bar))
    call pack_mpi_trace_exchange_vjp( &
        schedule, offsets, owner_bar, ghost_bar, local_bar, status)
    lhs = sum(owner_bar*owner_values_dot) + sum(ghost_bar*ghost_values_dot)
    rhs = sum(local_bar*local_values_dot)
    call record_condition(status%code == 0 .and. abs(lhs - rhs) < 1.0e-13_dp, &
        "trace packing VJP satisfies the real dot-product identity")

    received_ghost = reshape([0.8_dp, 0.9_dp, 0.5_dp, 0.6_dp], shape(received_ghost))
    received_ghost_dot = reshape([0.7_dp, -0.2_dp, 0.1_dp, 0.4_dp], shape(received_ghost_dot))
    call unpack_mpi_trace_exchange( &
        schedule, offsets, owner_values, received_ghost, local_unpacked, status)
    call record_condition(status%code == 0 .and. &
        maxval(abs(local_unpacked - reshape([1.0_dp, 2.0_dp, 0.8_dp, 0.9_dp, 0.6_dp, 0.7_dp, &
        0.8_dp, 0.5_dp, 0.6_dp, 1.2_dp], [local_count, component_count]))) < 1.0e-14_dp, &
        "trace unpacking restores owner values and received ghost values")
    call unpack_mpi_trace_exchange_jvp( &
        schedule, offsets, owner_values_dot, received_ghost_dot, local_unpacked_dot, status)
    call unpack_mpi_trace_exchange( &
        schedule, offsets, owner_values + step*owner_values_dot, &
        received_ghost + step*received_ghost_dot, local_unpacked_plus, status)
    call unpack_mpi_trace_exchange( &
        schedule, offsets, owner_values - step*owner_values_dot, &
        received_ghost - step*received_ghost_dot, local_unpacked_minus, status)
    call record_condition(status%code == 0 .and. maxval(abs(local_unpacked_dot - &
        (local_unpacked_plus - local_unpacked_minus)/(2.0_dp*step))) < 2.0e-8_dp, &
        "trace unpacking JVP matches an independent central difference")
    call unpack_mpi_trace_exchange_vjp( &
        schedule, offsets, local_unpacked_dot, owner_unpacked_bar, ghost_unpacked_bar, status)
    lhs = sum(local_unpacked_dot*local_unpacked_dot)
    rhs = sum(owner_unpacked_bar*owner_values_dot) + sum(ghost_unpacked_bar*received_ghost_dot)
    call record_condition(status%code == 0 .and. abs(lhs - rhs) < 1.0e-13_dp, &
        "trace unpacking VJP satisfies the real dot-product identity")

    call pack_mpi_trace_exchange(schedule, [1, 4, 4], local_values, owner_values, ghost_values, status)
    call record_condition(status%code /= 0, "trace packing rejects inconsistent offsets")
    call initialize_partition_layout( &
        bad_partitions(1), global_count, local_to_global_0, owner_rank_0, owned_0, 0, status)
    call initialize_partition_layout( &
        bad_partitions(2), global_count, local_to_global_1, [1, 2], owned_1, 2, status)
    call initialize_mpi_trace_exchange_schedule( &
        schedule, bad_partitions, component_count, status)
    call record_condition(status%code /= 0, "trace schedule rejects inconsistent owner ranks")
    call check_summary("MPI-ready trace exchange")
    if (.not. all_passed) error stop 1

contains

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_mpi_trace_exchange
