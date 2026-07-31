program test_fci_field_line_tracer
    use check, only: check_condition, check_summary
    use fortfem_api, only: trace_fci_field_line_rk4
    use fortfem_kinds, only: dp
    use fortsparse, only: fortsparse_status_t
    implicit none

    real(dp), parameter :: constant_initial(3) = [0.2_dp, 1.0_dp, -0.4_dp]
    real(dp), parameter :: constant_velocity(3) = [1.0_dp, 0.5_dp, -2.0_dp]
    real(dp), parameter :: exponential_initial(3) = [1.0_dp, 0.0_dp, 0.0_dp]
    real(dp) :: endpoint(3), endpoint_coarse(3), endpoint_fine(3)
    real(dp) :: coarse_error, fine_error
    type(fortsparse_status_t) :: status

    call trace_fci_field_line_rk4( &
        constant_initial, 2.0_dp, 4, constant_rhs, endpoint, status)
    call check_condition(status%code == 0, &
        "FCI field-line tracer accepts a positive step count")
    call check_condition(maxval(abs( &
        endpoint - (constant_initial + 2.0_dp*constant_velocity))) < &
        1.0e-14_dp, &
        "FCI field-line tracer is exact for a constant field-line velocity")

    call trace_fci_field_line_rk4( &
        exponential_initial, 1.0_dp, 4, exponential_rhs, endpoint_coarse, &
        status)
    call trace_fci_field_line_rk4( &
        exponential_initial, 1.0_dp, 8, exponential_rhs, endpoint_fine, &
        status)
    coarse_error = abs(endpoint_coarse(1) - exp(1.0_dp))
    fine_error = abs(endpoint_fine(1) - exp(1.0_dp))
    call check_condition(fine_error < coarse_error/10.0_dp, &
        "FCI RK4 field-line tracing converges against an exponential oracle")
    call check_condition(abs(endpoint_fine(2)) < 1.0e-14_dp .and. &
        abs(endpoint_fine(3)) < 1.0e-14_dp, &
        "FCI field-line tracing preserves inactive coordinates")

    call trace_fci_field_line_rk4( &
        constant_initial, 1.0_dp, 0, constant_rhs, endpoint, status)
    call check_condition(status%code /= 0, &
        "FCI field-line tracer rejects a nonpositive step count")
    call check_summary("PARALLAX-aligned FCI field-line tracer")

contains

    pure subroutine constant_rhs(phi, point, derivative)
        real(dp), intent(in) :: phi, point(:)
        real(dp), intent(out) :: derivative(:)

        associate (unused_phi => phi, unused_point => point)
            derivative = constant_velocity
        end associate
    end subroutine constant_rhs

    pure subroutine exponential_rhs(phi, point, derivative)
        real(dp), intent(in) :: phi, point(:)
        real(dp), intent(out) :: derivative(:)

        associate (unused_phi => phi)
            derivative = 0.0_dp
            derivative(1) = point(1)
        end associate
    end subroutine exponential_rhs

end program test_fci_field_line_tracer
