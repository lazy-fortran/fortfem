program test_fci_field_line_tracer_jvp
    use check, only: check_condition, check_summary
    use fortfem_api, only: trace_fci_field_line_rk4_jvp
    use fortfem_kinds, only: dp
    use fortsparse, only: fortsparse_status_t
    implicit none

    real(dp), parameter :: initial_point(2) = [1.0_dp, -2.0_dp]
    real(dp), parameter :: initial_point_dot(2) = [0.5_dp, 0.25_dp]
    real(dp), parameter :: toroidal_span = 0.7_dp
    integer, parameter :: step_count = 32
    real(dp) :: endpoint(2), endpoint_dot(2)
    real(dp) :: expected_endpoint(2), expected_endpoint_dot(2)
    type(fortsparse_status_t) :: status

    call trace_fci_field_line_rk4_jvp( &
        initial_point, initial_point_dot, toroidal_span, step_count, &
        exponential_rhs_jvp, endpoint, endpoint_dot, status)
    expected_endpoint = exp(toroidal_span)*initial_point
    expected_endpoint_dot = exp(toroidal_span)*initial_point_dot
    call check_condition(status%code == 0, &
        "FCI tracer JVP accepts a tangent RHS callback")
    call check_condition(maxval(abs(endpoint - expected_endpoint)) < 1.0e-8_dp, &
        "FCI tracer JVP primal endpoint matches the exponential oracle")
    call check_condition(maxval(abs(endpoint_dot - expected_endpoint_dot)) < &
        1.0e-8_dp, "FCI tracer JVP endpoint tangent matches the oracle")

    call trace_fci_field_line_rk4_jvp( &
        initial_point, initial_point_dot, toroidal_span, 0, exponential_rhs_jvp, &
        endpoint, endpoint_dot, status)
    call check_condition(status%code /= 0, &
        "FCI tracer JVP rejects a non-positive step count")
    call check_summary("FCI field-line tracer JVP")

contains

    pure subroutine exponential_rhs_jvp( &
            phi, point, point_dot, derivative, derivative_dot)
        real(dp), intent(in) :: phi
        real(dp), intent(in) :: point(:), point_dot(:)
        real(dp), intent(out) :: derivative(:), derivative_dot(:)

        associate (unused_phi => phi)
            derivative = point
            derivative_dot = point_dot
        end associate
    end subroutine exponential_rhs_jvp

end program test_fci_field_line_tracer_jvp
