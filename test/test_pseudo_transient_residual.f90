program test_pseudo_transient_residual
    use check, only: check_condition, check_summary
    use fortfem_api, only: &
        assemble_pseudo_transient_residual, &
        assemble_pseudo_transient_residual_jvp, &
        assemble_pseudo_transient_residual_vjp
    use fortfem_kinds, only: dp
    use fortsparse, only: FORTSPARSE_INVALID_MATRIX, FORTSPARSE_OK, &
        fortsparse_status_t
    implicit none

    real(dp), parameter :: epsilon_fd = 1.0e-7_dp
    logical :: all_passed

    all_passed = .true.
    call check_contract()
    call check_summary("Pseudo-transient residual")
    if (.not. all_passed) error stop 1

contains

    subroutine check_contract()
        integer, parameter :: state_size = 3
        real(dp) :: residual(state_size), mass(state_size, state_size)
        real(dp) :: state(state_size), previous_state(state_size), time_step
        real(dp) :: residual_dot(state_size), mass_dot(state_size, state_size)
        real(dp) :: state_dot(state_size), previous_state_dot(state_size), time_step_dot
        real(dp) :: augmented(state_size), augmented_dot(state_size)
        real(dp) :: augmented_plus(state_size), augmented_minus(state_size)
        real(dp) :: augmented_bar(state_size), residual_bar(state_size)
        real(dp) :: mass_bar(state_size, state_size), state_bar(state_size)
        real(dp) :: previous_state_bar(state_size), time_step_bar
        real(dp) :: expected(state_size), expected_dot(state_size)
        real(dp) :: lhs, rhs
        type(fortsparse_status_t) :: status

        residual = [0.4_dp, -0.2_dp, 0.7_dp]
        mass = reshape([1.2_dp, -0.3_dp, 0.4_dp, 0.5_dp, 1.1_dp, -0.2_dp, &
            -0.6_dp, 0.7_dp, 0.9_dp], [state_size, state_size])
        state = [1.1_dp, -0.6_dp, 0.3_dp]
        previous_state = [0.9_dp, -0.4_dp, 0.2_dp]
        time_step = 0.35_dp
        residual_dot = [-0.3_dp, 0.5_dp, -0.7_dp]
        mass_dot = reshape([0.2_dp, -0.1_dp, 0.3_dp, -0.4_dp, 0.5_dp, 0.6_dp, &
            0.7_dp, -0.8_dp, 0.1_dp], [state_size, state_size])
        state_dot = [0.2_dp, -0.8_dp, 0.4_dp]
        previous_state_dot = [0.1_dp, 0.2_dp, -0.5_dp]
        time_step_dot = -0.25_dp
        augmented_bar = [0.2_dp, -0.4_dp, 0.6_dp]

        expected = residual + matmul(mass, state - previous_state)/time_step
        call assemble_pseudo_transient_residual( &
            residual, mass, state, previous_state, time_step, augmented, status)
        call record(status%code == FORTSPARSE_OK .and. &
            maxval(abs(augmented - expected)) < 1.0e-14_dp, &
            "pseudo-transient residual value oracle")

        expected_dot = residual_dot + matmul(mass_dot, state - previous_state)/time_step + &
            matmul(mass, state_dot - previous_state_dot)/time_step - &
            matmul(mass, state - previous_state)*time_step_dot/time_step**2
        call assemble_pseudo_transient_residual_jvp( &
            residual, mass, state, previous_state, time_step, residual_dot, mass_dot, &
            state_dot, previous_state_dot, time_step_dot, augmented_dot, status)
        call assemble_pseudo_transient_residual( &
            residual + epsilon_fd*residual_dot, mass + epsilon_fd*mass_dot, &
            state + epsilon_fd*state_dot, previous_state + epsilon_fd*previous_state_dot, &
            time_step + epsilon_fd*time_step_dot, augmented_plus, status)
        call assemble_pseudo_transient_residual( &
            residual - epsilon_fd*residual_dot, mass - epsilon_fd*mass_dot, &
            state - epsilon_fd*state_dot, previous_state - epsilon_fd*previous_state_dot, &
            time_step - epsilon_fd*time_step_dot, augmented_minus, status)
        call record(status%code == FORTSPARSE_OK .and. &
            maxval(abs(augmented_dot - expected_dot)) < 1.0e-14_dp .and. &
            maxval(abs(augmented_dot - (augmented_plus - augmented_minus)/ &
            (2.0_dp*epsilon_fd))) < 2.0e-8_dp, &
            "pseudo-transient residual JVP oracle")

        call assemble_pseudo_transient_residual_vjp( &
            residual, mass, state, previous_state, time_step, augmented_bar, residual_bar, &
            mass_bar, state_bar, previous_state_bar, time_step_bar, status)
        lhs = dot_product(augmented_bar, augmented_dot)
        rhs = dot_product(residual_bar, residual_dot) + sum(mass_bar*mass_dot) + &
            dot_product(state_bar, state_dot) + dot_product(previous_state_bar, previous_state_dot) + &
            time_step_bar*time_step_dot
        call record(status%code == FORTSPARSE_OK .and. abs(lhs - rhs) < 1.0e-13_dp, &
            "pseudo-transient residual VJP dot-product oracle")

        call assemble_pseudo_transient_residual( &
            residual, mass, state, previous_state, 0.0_dp, augmented, status)
        call record(status%code == FORTSPARSE_INVALID_MATRIX, &
            "non-positive pseudo-transient step is rejected")
    end subroutine check_contract

    subroutine record(condition, message)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: message

        call check_condition(condition, message)
        all_passed = all_passed .and. condition
    end subroutine record

end program test_pseudo_transient_residual
