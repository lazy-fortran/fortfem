program test_complex_distributed_trace_reduction
    use, intrinsic :: ieee_arithmetic, only: ieee_quiet_nan, ieee_value
    use check, only: check_condition, check_summary
    use fortfem_feec, only: &
        assemble_complex_distributed_trace_reduction, &
        assemble_complex_distributed_trace_reduction_jvp, &
        assemble_complex_distributed_trace_reduction_vjp
    use fortfem_feec, only: &
        distributed_trace_layout_t, initialize_distributed_trace_layout
    use fortfem_kinds, only: dp
    use fortfem_feec, only: &
        initialize_partition_layout, partition_layout_t
    use fortsparse, only: FORTSPARSE_OK, fortsparse_status_t
    implicit none

    integer, parameter :: component_count = 2, global_count = 4
    integer, parameter :: local_count = 6, partition_count = 2
    real(dp), parameter :: difference_step = 1.0e-7_dp
    integer, parameter :: global_ids(local_count) = [1, 2, 4, 2, 3, 4]
    logical, parameter :: owner_mask(local_count) = &
        [.true., .false., .true., .true., .true., .false.]
    type(distributed_trace_layout_t) :: layout
    type(partition_layout_t) :: partitions(partition_count)
    type(fortsparse_status_t) :: status, status_minus, status_plus
    complex(dp) :: global_bar(global_count, component_count)
    complex(dp) :: global_dot(global_count, component_count)
    complex(dp) :: global_expected(global_count, component_count)
    complex(dp) :: global_minus(global_count, component_count)
    complex(dp) :: global_plus(global_count, component_count)
    complex(dp) :: global_values(global_count, component_count)
    complex(dp) :: local_bar(local_count, component_count)
    complex(dp) :: local_bar_expected(local_count, component_count)
    complex(dp) :: local_dot(local_count, component_count)
    complex(dp) :: local_dot_invalid(local_count, component_count)
    complex(dp) :: local_invalid(local_count, component_count)
    complex(dp) :: local_values(local_count, component_count)
    complex(dp) :: owner_bar(global_count, component_count)
    complex(dp) :: owner_dot(global_count, component_count)
    complex(dp) :: owner_expected(global_count, component_count)
    complex(dp) :: owner_minus(global_count, component_count)
    complex(dp) :: owner_plus(global_count, component_count)
    complex(dp) :: owner_values(global_count, component_count)
    complex(dp) :: short_global(global_count - 1, component_count)
    complex(dp) :: global_bar_invalid(global_count, component_count)
    real(dp) :: global_matrix(global_count, local_count)
    real(dp) :: owner_matrix(global_count, local_count)
    real(dp) :: forward_pairing, reverse_pairing
    integer :: component, row
    logical :: all_passed

    all_passed = .true.
    call initialize_partition_layout( &
        partitions(1), global_count, [1, 2, 4], [0, 1, 0], &
        [.true., .false., .true.], 0, status)
    call record_condition(status%code == FORTSPARSE_OK, &
        "first complex reduction partition initializes")
    call initialize_partition_layout( &
        partitions(2), global_count, [2, 3, 4], [1, 1, 0], &
        [.true., .true., .false.], 1, status)
    call record_condition(status%code == FORTSPARSE_OK, &
        "second complex reduction partition initializes")
    call initialize_distributed_trace_layout( &
        layout, partitions, component_count, status)
    call record_condition(status%code == FORTSPARSE_OK, &
        "complex reduction ownership ledger initializes")

    local_values = reshape([ &
        cmplx(0.4_dp, -0.2_dp, dp), cmplx(-0.7_dp, 0.3_dp, dp), &
        cmplx(0.1_dp, 0.8_dp, dp), cmplx(0.6_dp, -0.5_dp, dp), &
        cmplx(-0.3_dp, -0.4_dp, dp), cmplx(0.9_dp, 0.2_dp, dp), &
        cmplx(-0.2_dp, 0.5_dp, dp), cmplx(0.3_dp, -0.6_dp, dp), &
        cmplx(0.8_dp, 0.1_dp, dp), cmplx(-0.5_dp, 0.7_dp, dp), &
        cmplx(0.2_dp, -0.9_dp, dp), cmplx(0.7_dp, 0.4_dp, dp)], &
        shape(local_values))
    local_dot = reshape([ &
        cmplx(0.03_dp, -0.01_dp, dp), cmplx(-0.02_dp, 0.04_dp, dp), &
        cmplx(0.05_dp, 0.02_dp, dp), cmplx(-0.01_dp, -0.03_dp, dp), &
        cmplx(0.04_dp, -0.05_dp, dp), cmplx(-0.06_dp, 0.01_dp, dp), &
        cmplx(-0.02_dp, 0.03_dp, dp), cmplx(0.01_dp, -0.04_dp, dp), &
        cmplx(0.06_dp, 0.02_dp, dp), cmplx(-0.03_dp, 0.05_dp, dp), &
        cmplx(0.02_dp, -0.01_dp, dp), cmplx(0.04_dp, 0.03_dp, dp)], &
        shape(local_dot))

    global_matrix = 0.0_dp
    owner_matrix = 0.0_dp
    do row = 1, local_count
        global_matrix(global_ids(row), row) = 1.0_dp
        if (owner_mask(row)) owner_matrix(global_ids(row), row) = 1.0_dp
    end do
    do component = 1, component_count
        global_expected(:, component) = &
            matmul(global_matrix, local_values(:, component))
        owner_expected(:, component) = &
            matmul(owner_matrix, local_values(:, component))
    end do

    call assemble_complex_distributed_trace_reduction( &
        layout, local_values, global_values, owner_values, status)
    call record_condition(status%code == FORTSPARSE_OK .and. &
        maxval(abs(global_values - global_expected)) < 1.0e-14_dp .and. &
        maxval(abs(owner_values - owner_expected)) < 1.0e-14_dp, &
        "complex reduction matches independent aggregation matrices")

    call assemble_complex_distributed_trace_reduction_jvp( &
        layout, local_dot, global_dot, owner_dot, status)
    call assemble_complex_distributed_trace_reduction( &
        layout, local_values + difference_step*local_dot, global_plus, &
        owner_plus, status_plus)
    call assemble_complex_distributed_trace_reduction( &
        layout, local_values - difference_step*local_dot, global_minus, &
        owner_minus, status_minus)
    call record_condition(status%code == FORTSPARSE_OK .and. &
        status_plus%code == FORTSPARSE_OK .and. &
        status_minus%code == FORTSPARSE_OK .and. &
        maxval(abs(global_dot - &
        (global_plus - global_minus)/(2.0_dp*difference_step))) < 2.0e-9_dp .and. &
        maxval(abs(owner_dot - &
        (owner_plus - owner_minus)/(2.0_dp*difference_step))) < 2.0e-9_dp, &
        "complex reduction JVP matches complete central reassembly")

    global_bar = reshape([ &
        cmplx(0.2_dp, -0.1_dp, dp), cmplx(-0.3_dp, 0.5_dp, dp), &
        cmplx(0.4_dp, 0.2_dp, dp), cmplx(-0.6_dp, -0.4_dp, dp), &
        cmplx(0.7_dp, -0.2_dp, dp), cmplx(-0.5_dp, 0.3_dp, dp), &
        cmplx(0.1_dp, 0.6_dp, dp), cmplx(-0.4_dp, -0.7_dp, dp)], &
        shape(global_bar))
    owner_bar = reshape([ &
        cmplx(-0.8_dp, 0.2_dp, dp), cmplx(0.6_dp, -0.3_dp, dp), &
        cmplx(-0.1_dp, 0.4_dp, dp), cmplx(0.5_dp, 0.7_dp, dp), &
        cmplx(-0.2_dp, -0.6_dp, dp), cmplx(0.3_dp, 0.1_dp, dp), &
        cmplx(-0.7_dp, 0.5_dp, dp), cmplx(0.9_dp, -0.4_dp, dp)], &
        shape(owner_bar))
    call assemble_complex_distributed_trace_reduction_vjp( &
        layout, global_bar, owner_bar, local_bar, status)
    do component = 1, component_count
        local_bar_expected(:, component) = &
            matmul(transpose(global_matrix), global_bar(:, component)) + &
            matmul(transpose(owner_matrix), owner_bar(:, component))
    end do
    forward_pairing = real(sum(conjg(global_bar)*global_dot) + &
        sum(conjg(owner_bar)*owner_dot), dp)
    reverse_pairing = real(sum(conjg(local_bar)*local_dot), dp)
    call record_condition(status%code == FORTSPARSE_OK .and. &
        maxval(abs(local_bar - local_bar_expected)) < 1.0e-14_dp .and. &
        abs(forward_pairing - reverse_pairing) < 2.0e-14_dp, &
        "complex reduction VJP matches the matrix and real adjoint oracles")

    call assemble_complex_distributed_trace_reduction( &
        layout, local_values, short_global, owner_values, status)
    call record_condition(status%code /= FORTSPARSE_OK .and. &
        all(short_global == cmplx(0.0_dp, 0.0_dp, dp)), &
        "complex reduction rejects incompatible output shapes")

    local_invalid = local_values
    local_invalid(2, 1) = cmplx( &
        ieee_value(0.0_dp, ieee_quiet_nan), 0.0_dp, dp)
    call assemble_complex_distributed_trace_reduction( &
        layout, local_invalid, global_values, owner_values, status)
    call record_condition(status%code /= FORTSPARSE_OK .and. &
        all(global_values == cmplx(0.0_dp, 0.0_dp, dp)) .and. &
        all(owner_values == cmplx(0.0_dp, 0.0_dp, dp)), &
        "complex reduction rejects non-finite local values")

    local_dot_invalid = local_dot
    local_dot_invalid(4, 2) = cmplx( &
        0.0_dp, ieee_value(0.0_dp, ieee_quiet_nan), dp)
    call assemble_complex_distributed_trace_reduction_jvp( &
        layout, local_dot_invalid, global_dot, owner_dot, status)
    call record_condition(status%code /= FORTSPARSE_OK .and. &
        all(global_dot == cmplx(0.0_dp, 0.0_dp, dp)) .and. &
        all(owner_dot == cmplx(0.0_dp, 0.0_dp, dp)), &
        "complex reduction JVP rejects non-finite directions")

    global_bar_invalid = global_bar
    global_bar_invalid(1, 2) = cmplx( &
        ieee_value(0.0_dp, ieee_quiet_nan), 0.0_dp, dp)
    call assemble_complex_distributed_trace_reduction_vjp( &
        layout, global_bar_invalid, owner_bar, local_bar, status)
    call record_condition(status%code /= FORTSPARSE_OK .and. &
        all(local_bar == cmplx(0.0_dp, 0.0_dp, dp)), &
        "complex reduction VJP rejects non-finite adjoints")

    call check_summary("complex distributed trace reduction")
    if (.not. all_passed) error stop 1

contains

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_complex_distributed_trace_reduction
