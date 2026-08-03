program test_complex_mpi_trace_exchange
    use, intrinsic :: ieee_arithmetic, only: ieee_quiet_nan, ieee_value
    use check, only: check_condition, check_summary
    use fortfem_kinds, only: dp
    use fortfem_mpi_trace_exchange, only: &
        initialize_mpi_trace_exchange_schedule, mpi_trace_exchange_schedule_t, &
        pack_complex_mpi_trace_exchange, pack_complex_mpi_trace_exchange_jvp, &
        pack_complex_mpi_trace_exchange_vjp, unpack_complex_mpi_trace_exchange, &
        unpack_complex_mpi_trace_exchange_jvp, unpack_complex_mpi_trace_exchange_vjp
    use fortfem_partition_layout, only: initialize_partition_layout, partition_layout_t
    use fortsparse, only: fortsparse_status_t
    implicit none

    integer, parameter :: global_count = 3, component_count = 2
    integer, parameter :: partition_count = 2, local_count = 5
    integer, parameter :: owner_count = 3, ghost_count = 2
    integer, parameter :: offsets(3) = [1, 4, 6]
    real(dp), parameter :: step = 1.0e-7_dp
    type(partition_layout_t) :: partitions(partition_count)
    type(mpi_trace_exchange_schedule_t) :: schedule
    type(fortsparse_status_t) :: status
    integer, parameter :: local_to_global_0(3) = [1, 2, 3]
    integer, parameter :: owner_rank_0(3) = [0, 0, 1]
    logical, parameter :: owned_0(3) = [.true., .true., .false.]
    integer, parameter :: local_to_global_1(2) = [2, 3]
    integer, parameter :: owner_rank_1(2) = [0, 1]
    logical, parameter :: owned_1(2) = [.false., .true.]
    complex(dp) :: local_values(local_count, component_count)
    complex(dp) :: local_dot(local_count, component_count)
    complex(dp) :: owner_values(owner_count, component_count)
    complex(dp) :: ghost_values(ghost_count, component_count)
    complex(dp) :: owner_dot(owner_count, component_count)
    complex(dp) :: ghost_dot(ghost_count, component_count)
    complex(dp) :: owner_plus(owner_count, component_count)
    complex(dp) :: ghost_plus(ghost_count, component_count)
    complex(dp) :: owner_minus(owner_count, component_count)
    complex(dp) :: ghost_minus(ghost_count, component_count)
    complex(dp) :: local_plus(local_count, component_count)
    complex(dp) :: local_minus(local_count, component_count)
    complex(dp) :: owner_bar(owner_count, component_count)
    complex(dp) :: ghost_bar(ghost_count, component_count)
    complex(dp) :: local_bar(local_count, component_count)
    complex(dp) :: received_ghost(ghost_count, component_count)
    complex(dp) :: received_ghost_dot(ghost_count, component_count)
    complex(dp) :: unpacked(local_count, component_count)
    complex(dp) :: unpacked_dot(local_count, component_count)
    complex(dp) :: unpacked_plus(local_count, component_count)
    complex(dp) :: unpacked_minus(local_count, component_count)
    complex(dp) :: unpack_owner_plus(owner_count, component_count)
    complex(dp) :: unpack_owner_minus(owner_count, component_count)
    complex(dp) :: unpack_ghost_plus(ghost_count, component_count)
    complex(dp) :: unpack_ghost_minus(ghost_count, component_count)
    complex(dp) :: unpacked_bar(local_count, component_count)
    complex(dp) :: unpack_owner_bar(owner_count, component_count)
    complex(dp) :: unpack_ghost_bar(ghost_count, component_count)
    complex(dp) :: short_owner(owner_count - 1, component_count)
    real(dp) :: lhs, rhs
    logical :: all_passed

    all_passed = .true.
    call initialize_partition_layout( &
        partitions(1), global_count, local_to_global_0, owner_rank_0, &
        owned_0, 0, status)
    call initialize_partition_layout( &
        partitions(2), global_count, local_to_global_1, owner_rank_1, &
        owned_1, 1, status)
    call initialize_mpi_trace_exchange_schedule( &
        schedule, partitions, component_count, status)
    call record(status%code == 0, "complex trace test initializes a valid schedule")

    local_values(:, 1) = [ &
        cmplx(1.0_dp, -0.5_dp, dp), cmplx(-2.0_dp, 0.75_dp, dp), &
        cmplx(0.25_dp, 3.0_dp, dp), cmplx(4.0_dp, -1.5_dp, dp), &
        cmplx(-0.75_dp, -2.0_dp, dp)]
    local_values(:, 2) = [ &
        cmplx(-1.5_dp, 2.0_dp, dp), cmplx(0.5_dp, -3.0_dp, dp), &
        cmplx(2.25_dp, 0.125_dp, dp), cmplx(-4.5_dp, 1.0_dp, dp), &
        cmplx(3.25_dp, -0.625_dp, dp)]
    local_dot(:, 1) = [ &
        cmplx(0.3_dp, -0.2_dp, dp), cmplx(-0.4_dp, 0.7_dp, dp), &
        cmplx(0.8_dp, 0.1_dp, dp), cmplx(-0.6_dp, -0.5_dp, dp), &
        cmplx(0.2_dp, 0.9_dp, dp)]
    local_dot(:, 2) = [ &
        cmplx(-0.1_dp, 0.4_dp, dp), cmplx(0.9_dp, -0.3_dp, dp), &
        cmplx(-0.7_dp, 0.6_dp, dp), cmplx(0.5_dp, 0.2_dp, dp), &
        cmplx(-0.8_dp, -0.4_dp, dp)]

    call pack_complex_mpi_trace_exchange( &
        schedule, offsets, local_values, owner_values, ghost_values, status)
    call record(status%code == 0 .and. &
        maxval(abs(owner_values - local_values([1, 2, 5], :))) < 1.0e-14_dp .and. &
        maxval(abs(ghost_values - local_values([3, 4], :))) < 1.0e-14_dp, &
        "complex packing follows independently selected owner and ghost rows")

    call pack_complex_mpi_trace_exchange_jvp( &
        schedule, offsets, local_dot, owner_dot, ghost_dot, status)
    local_plus = local_values + step*local_dot
    local_minus = local_values - step*local_dot
    call pack_complex_mpi_trace_exchange( &
        schedule, offsets, local_plus, owner_plus, ghost_plus, status)
    call pack_complex_mpi_trace_exchange( &
        schedule, offsets, local_minus, owner_minus, ghost_minus, status)
    call record(status%code == 0 .and. &
        maxval(abs(owner_dot - (owner_plus - owner_minus)/(2.0_dp*step))) < &
        2.0e-8_dp .and. &
        maxval(abs(ghost_dot - (ghost_plus - ghost_minus)/(2.0_dp*step))) < 2.0e-8_dp, &
        "complex packing JVP matches central re-evaluation")

    owner_bar = cmplx(0.4_dp, -0.3_dp, dp)*owner_values + &
        cmplx(-0.2_dp, 0.1_dp, dp)
    ghost_bar = cmplx(-0.5_dp, 0.25_dp, dp)*ghost_values + &
        cmplx(0.1_dp, -0.4_dp, dp)
    call pack_complex_mpi_trace_exchange_vjp( &
        schedule, offsets, owner_bar, ghost_bar, local_bar, status)
    lhs = real(sum(conjg(owner_bar)*owner_dot) + sum(conjg(ghost_bar)*ghost_dot), dp)
    rhs = real(sum(conjg(local_bar)*local_dot), dp)
    call record(status%code == 0 .and. abs(lhs - rhs) < 1.0e-13_dp, &
        "complex packing VJP satisfies the real-part adjoint identity")

    received_ghost = ghost_values + cmplx(1.0_dp, -0.75_dp, dp)
    received_ghost_dot = ghost_dot + cmplx(-0.15_dp, 0.35_dp, dp)
    call unpack_complex_mpi_trace_exchange( &
        schedule, offsets, owner_values, received_ghost, unpacked, status)
    call record(status%code == 0 .and. &
        maxval(abs(unpacked([1, 2, 5], :) - owner_values)) < 1.0e-14_dp .and. &
        maxval(abs(unpacked([3, 4], :) - received_ghost)) < 1.0e-14_dp, &
        "complex unpacking restores independently selected transport rows")

    call unpack_complex_mpi_trace_exchange_jvp( &
        schedule, offsets, owner_dot, received_ghost_dot, unpacked_dot, status)
    unpack_owner_plus = owner_values + step*owner_dot
    unpack_owner_minus = owner_values - step*owner_dot
    unpack_ghost_plus = received_ghost + step*received_ghost_dot
    unpack_ghost_minus = received_ghost - step*received_ghost_dot
    call unpack_complex_mpi_trace_exchange( &
        schedule, offsets, unpack_owner_plus, unpack_ghost_plus, &
        unpacked_plus, status)
    call unpack_complex_mpi_trace_exchange( &
        schedule, offsets, unpack_owner_minus, unpack_ghost_minus, &
        unpacked_minus, status)
    call record(status%code == 0 .and. &
        maxval(abs(unpacked_dot - (unpacked_plus - unpacked_minus)/(2.0_dp*step))) < &
        2.0e-8_dp, "complex unpacking JVP matches central re-evaluation")

    unpacked_bar = cmplx(-0.3_dp, 0.2_dp, dp)*unpacked + &
        cmplx(0.6_dp, -0.1_dp, dp)
    call unpack_complex_mpi_trace_exchange_vjp( &
        schedule, offsets, unpacked_bar, unpack_owner_bar, unpack_ghost_bar, status)
    lhs = real(sum(conjg(unpacked_bar)*unpacked_dot), dp)
    rhs = real(sum(conjg(unpack_owner_bar)*owner_dot) + &
        sum(conjg(unpack_ghost_bar)*received_ghost_dot), dp)
    call record(status%code == 0 .and. abs(lhs - rhs) < 1.0e-13_dp, &
        "complex unpacking VJP satisfies the real-part adjoint identity")

    local_values(1, 1) = cmplx( &
        ieee_value(0.0_dp, ieee_quiet_nan), aimag(local_values(1, 1)), dp)
    call pack_complex_mpi_trace_exchange( &
        schedule, offsets, local_values, owner_values, ghost_values, status)
    call record(status%code /= 0, "complex packing rejects a non-finite real part")
    local_values(1, 1) = cmplx(1.0_dp, ieee_value(0.0_dp, ieee_quiet_nan), dp)
    call pack_complex_mpi_trace_exchange( &
        schedule, offsets, local_values, owner_values, ghost_values, status)
    call record(status%code /= 0, "complex packing rejects a non-finite imaginary part")
    call unpack_complex_mpi_trace_exchange( &
        schedule, offsets, short_owner, received_ghost, unpacked, status)
    call record(status%code /= 0, &
        "complex unpacking rejects an incompatible owner shape")

    call check_summary("complex MPI-ready trace exchange")
    if (.not. all_passed) error stop 1

contains

    subroutine record(condition, description)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record

end program test_complex_mpi_trace_exchange
