program test_physical_trace_owner_selection
    use check, only: check_condition, check_summary
    use fortfem_core, only: &
        physical_trace_ownership_t, initialize_physical_trace_ownership, &
        physical_trace_owner_selection_t, initialize_physical_trace_owner_selection, &
        gather_physical_trace_values, gather_physical_trace_values_jvp, &
        gather_physical_trace_values_vjp
    use fortfem_kinds, only: dp
    use fortsparse, only: fortsparse_status_t
    implicit none

    type(physical_trace_ownership_t) :: ownership(2)
    type(physical_trace_owner_selection_t) :: selection
    real(dp) :: coordinates_1(2, 2), coordinates_2(2, 2)
    integer :: ids_1(2), ids_2(2), owner_1(2), owner_2(2)
    logical :: owned_1(2), owned_2(2)
    logical :: owned_bad(2)
    real(dp) :: local_values(4, 2), global_values(2, 2), global_dot(2, 2)
    real(dp) :: local_dot(4, 2), global_bar(2, 2), local_bar(4, 2)
    real(dp) :: lhs, rhs
    integer, parameter :: offsets(3) = [1, 3, 5]
    type(fortsparse_status_t) :: status

    coordinates_1 = reshape([0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp], [2, 2])
    coordinates_2 = reshape([1.0_dp + 1.0e-10_dp, 0.0_dp, 0.0_dp, 0.0_dp], [2, 2])
    ids_1 = [1, 2]
    ids_2 = [2, 1]
    owner_1 = [0, 1]
    owner_2 = [1, 0]
    owned_1 = [.true., .false.]
    owned_2 = [.true., .false.]
    call initialize_physical_trace_ownership(ownership(1), coordinates_1, ids_1, &
        owner_1, owned_1, 0, 1.0e-8_dp, status)
    call check_condition(status%code == 0, "first physical trace partition initializes")
    call initialize_physical_trace_ownership(ownership(2), coordinates_2, ids_2, &
        owner_2, owned_2, 1, 1.0e-8_dp, status)
    call check_condition(status%code == 0, "second reordered physical trace partition initializes")
    call initialize_physical_trace_owner_selection(selection, ownership, 1.0e-8_dp, status)
    call check_condition(status%code == 0, "owner selection finds one owner per physical trace ID")

    local_values = reshape([1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, &
        0.5_dp, 0.6_dp, 0.7_dp, 0.8_dp], [4, 2])
    call gather_physical_trace_values(selection, offsets, local_values, global_values, status)
    call check_condition(status%code == 0 .and. maxval(abs(global_values - reshape([ &
        1.0_dp, 3.0_dp, 0.5_dp, 0.7_dp], [2, 2]))) < 1.0e-14_dp, &
        "owner gather is deterministic in global-ID order")

    local_dot = reshape([0.2_dp, -0.1_dp, 0.4_dp, 0.3_dp, &
        -0.5_dp, 0.6_dp, 0.1_dp, -0.2_dp], [4, 2])
    call gather_physical_trace_values_jvp(selection, offsets, local_dot, global_dot, status)
    global_bar = reshape([0.4_dp, -0.2_dp, 0.7_dp, 0.3_dp], [2, 2])
    call gather_physical_trace_values_vjp(selection, offsets, global_bar, local_bar, status)
    lhs = sum(global_bar*global_dot)
    rhs = sum(local_bar*local_dot)
    call check_condition(status%code == 0 .and. abs(lhs - rhs) < 1.0e-14_dp, &
        "owner gather VJP satisfies the independent transpose oracle")
    call check_condition(maxval(abs(local_bar(2, :))) < 1.0e-14_dp .and. &
        maxval(abs(local_bar(4, :))) < 1.0e-14_dp, &
        "owner gather VJP leaves non-owner ghost rows zero")
    owned_bad = [.true., .true.]
    call initialize_physical_trace_ownership(ownership(2), coordinates_2, ids_2, &
        owner_2, owned_bad, 1, 1.0e-8_dp, status)
    call initialize_physical_trace_owner_selection(selection, ownership, 1.0e-8_dp, status)
    call check_condition(status%code /= 0, &
        "owner selection rejects duplicate owned copies")
    call check_summary("physical trace owner selection")
end program test_physical_trace_owner_selection
