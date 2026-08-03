program test_deflated_residual
    use check, only: check_condition, check_summary
    use fortfem_feec, only: &
        assemble_deflated_residual, assemble_deflated_residual_jvp, &
        assemble_deflated_residual_vjp
    use fortfem_kinds, only: dp
    use fortsparse, only: fortsparse_status_t
    implicit none

    real(dp), parameter :: state(3) = [0.4_dp, -0.7_dp, 1.1_dp]
    real(dp), parameter :: state_dot(3) = [-0.2_dp, 0.3_dp, 0.5_dp]
    real(dp), parameter :: roots(3, 2) = reshape([ &
        0.1_dp, -0.2_dp, 0.8_dp, 1.0_dp, -0.5_dp, 1.4_dp], [3, 2])
    real(dp), parameter :: residual(3) = [0.6_dp, -0.4_dp, 0.9_dp]
    real(dp), parameter :: residual_dot(3) = [0.17_dp, -0.29_dp, 0.41_dp]
    real(dp), parameter :: scale = 0.8_dp, power = 2.5_dp, shift = 0.3_dp
    real(dp), parameter :: epsilon_fd = 1.0e-6_dp
    real(dp) :: output(3), output_plus(3), output_minus(3), output_dot(3)
    real(dp) :: output_bar(3), residual_bar(3), state_bar(3)
    real(dp) :: factor, factor_plus, factor_minus, factor_dot
    real(dp) :: lhs, rhs
    logical :: all_passed
    type(fortsparse_status_t) :: status

    all_passed = .true.
    call assemble_deflated_residual( &
        residual, state, roots, scale, power, shift, output, factor, status)
    call record_condition(status%code == 0 .and. factor > 1.0_dp .and. &
        maxval(abs(output - factor*residual)) < 1.0e-14_dp, &
        "deflated residual value matches the multiplicative oracle")

    call assemble_deflated_residual_jvp( &
        residual, state, roots, scale, power, shift, residual_dot, state_dot, &
        output_dot, factor_dot, status)
    call assemble_deflated_residual( &
        residual + epsilon_fd*residual_dot, state + epsilon_fd*state_dot, roots, &
        scale, power, shift, output_plus, factor_plus, status)
    call assemble_deflated_residual( &
        residual - epsilon_fd*residual_dot, state - epsilon_fd*state_dot, roots, &
        scale, power, shift, output_minus, factor_minus, status)
    call record_condition(status%code == 0 .and. &
        maxval(abs(output_dot - (output_plus - output_minus)/(2.0_dp*epsilon_fd))) &
        < 2.0e-8_dp .and. abs(factor_dot - &
        (factor_plus - factor_minus)/(2.0_dp*epsilon_fd)) < 2.0e-8_dp, &
        "deflated residual JVP matches an independent central difference")

    output_bar = [-0.3_dp, 0.6_dp, 1.2_dp]
    call assemble_deflated_residual_vjp( &
        residual, state, roots, scale, power, shift, output_bar, residual_bar, &
        state_bar, status)
    lhs = sum(output_bar*output_dot)
    rhs = sum(residual_bar*residual_dot) + sum(state_bar*state_dot)
    call record_condition(status%code == 0 .and. abs(lhs - rhs) < 2.0e-14_dp, &
        "deflated residual VJP satisfies the real adjoint identity")

    call assemble_deflated_residual( &
        residual, state, roots, scale, power, 0.0_dp, output, factor, status)
    call record_condition(status%code /= 0, "deflation rejects a zero shift")

    call check_summary("deflated residual")
    if (.not. all_passed) error stop 1

contains

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_deflated_residual
