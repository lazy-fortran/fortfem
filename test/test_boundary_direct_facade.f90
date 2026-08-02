program test_boundary_direct_facade
    use, intrinsic :: iso_fortran_env, only: real64
    use check, only: check_condition, check_summary
    use fortfem_boundary, only: &
        BOUNDARY_OPERATOR_BACKEND_DTN, &
        BOUNDARY_OPERATOR_TRACE_CHANNEL_SCALAR, &
        apply_planar_helmholtz_dtn, &
        boundary_operator_contract_t, &
        initialize_boundary_operator_contract, &
        initialize_boundary_operator_trace_metadata, &
        validate_boundary_operator_contract
    implicit none

    integer, parameter :: dp = real64
    integer, parameter :: sample_count = 8
    real(dp), parameter :: period = 2.0_dp*acos(-1.0_dp)
    real(dp), parameter :: wavenumber = 2.0_dp
    complex(dp) :: trace(sample_count), normal_derivative(sample_count)
    type(boundary_operator_contract_t) :: contract
    integer :: status

    trace = cmplx(1.0_dp, 0.0_dp, dp)
    call apply_planar_helmholtz_dtn( &
        trace, wavenumber, period, normal_derivative)
    call check_condition(maxval(abs(normal_derivative - &
        cmplx(0.0_dp, wavenumber, dp))) < 2.0e-13_dp, &
        "direct boundary facade returns the constant-mode DtN trace")

    call initialize_boundary_operator_contract( &
        contract, BOUNDARY_OPERATOR_BACKEND_DTN, "helmholtz", "H1-trace", &
        sample_count, sample_count, .true., .true., .true., .true., .true., &
        .true., "normal-length", "outgoing", &
        "planar Helmholtz DtN analytical smoke test", "periodic-line-1", &
        status)
    call check_condition(status == 0, &
        "direct boundary facade initializes DtN metadata")

    call initialize_boundary_operator_trace_metadata( &
        contract, BOUNDARY_OPERATOR_TRACE_CHANNEL_SCALAR, "complex-l2", status)
    call check_condition(status == 0 .and. &
        validate_boundary_operator_contract(contract, status) .and. &
        trim(contract%backend_name) == "DtN" .and. &
        trim(contract%topology_id) == "periodic-line-1", &
        "direct boundary facade validates trace metadata")

    call check_summary("direct boundary facade")
end program test_boundary_direct_facade
