program test_pseudo_arclength_residual
    use check, only: check_condition, check_summary
    use fortfem_api, only: &
        assemble_pseudo_arclength_residual, &
        assemble_pseudo_arclength_residual_jvp, &
        assemble_pseudo_arclength_residual_vjp
    use fortfem_kinds, only: dp
    use fortsparse, only: FORTSPARSE_OK, fortsparse_status_t
    implicit none

    real(dp), parameter :: epsilon_fd = 1.0e-7_dp
    logical :: all_passed

    all_passed = .true.
    call check_contract()
    call check_summary("Pseudo-arclength continuation residual")
    if (.not. all_passed) error stop 1

contains

    subroutine check_contract()
        integer, parameter :: state_size = 3, parameter_size = 2
        real(dp) :: residual(state_size), state(state_size), parameter_value(parameter_size)
        real(dp) :: previous_state(state_size), previous_parameter(parameter_size)
        real(dp) :: tangent_state(state_size), tangent_parameter(parameter_size), step
        real(dp) :: residual_dot(state_size), state_dot(state_size), parameter_dot(parameter_size)
        real(dp) :: previous_state_dot(state_size), previous_parameter_dot(parameter_size)
        real(dp) :: tangent_state_dot(state_size), tangent_parameter_dot(parameter_size), step_dot
        real(dp) :: augmented(state_size + 1), augmented_dot(state_size + 1)
        real(dp) :: augmented_plus(state_size + 1), augmented_minus(state_size + 1)
        real(dp) :: augmented_bar(state_size + 1), residual_bar(state_size)
        real(dp) :: state_bar(state_size), parameter_bar(parameter_size)
        real(dp) :: previous_state_bar(state_size), previous_parameter_bar(parameter_size)
        real(dp) :: tangent_state_bar(state_size), tangent_parameter_bar(parameter_size), step_bar
        real(dp) :: expected(state_size + 1), expected_dot(state_size + 1)
        real(dp) :: lhs, rhs
        type(fortsparse_status_t) :: status

        residual = [0.4_dp, -0.2_dp, 0.7_dp]
        state = [1.1_dp, -0.6_dp, 0.3_dp]
        parameter_value = [0.8_dp, -0.5_dp]
        previous_state = [0.9_dp, -0.4_dp, 0.2_dp]
        previous_parameter = [0.7_dp, -0.3_dp]
        tangent_state = [0.6_dp, -0.1_dp, 0.8_dp]
        tangent_parameter = [0.2_dp, -0.4_dp]
        step = 0.35_dp
        residual_dot = [-0.3_dp, 0.5_dp, -0.7_dp]
        state_dot = [0.2_dp, -0.8_dp, 0.4_dp]
        parameter_dot = [-0.6_dp, 0.3_dp]
        previous_state_dot = [0.1_dp, 0.2_dp, -0.5_dp]
        previous_parameter_dot = [0.4_dp, -0.2_dp]
        tangent_state_dot = [-0.3_dp, 0.7_dp, 0.2_dp]
        tangent_parameter_dot = [0.5_dp, -0.1_dp]
        step_dot = -0.25_dp
        augmented_bar = [0.2_dp, -0.4_dp, 0.6_dp, -0.8_dp]

        expected(:state_size) = residual
        expected(state_size + 1) = dot_product(tangent_state, state - previous_state) + &
            dot_product(tangent_parameter, parameter_value - previous_parameter) - step
        call assemble_pseudo_arclength_residual( &
            residual, state, parameter_value, previous_state, previous_parameter, &
            tangent_state, tangent_parameter, step, augmented, status)
        call record(status%code == FORTSPARSE_OK .and. &
            maxval(abs(augmented - expected)) < 1.0e-14_dp, "continuation residual value oracle")

        expected_dot(:state_size) = residual_dot
        expected_dot(state_size + 1) = dot_product(tangent_state_dot, state - previous_state) + &
            dot_product(tangent_state, state_dot - previous_state_dot) + &
            dot_product(tangent_parameter_dot, parameter_value - previous_parameter) + &
            dot_product(tangent_parameter, parameter_dot - previous_parameter_dot) - step_dot
        call assemble_pseudo_arclength_residual_jvp( &
            residual, state, parameter_value, previous_state, previous_parameter, tangent_state, &
            tangent_parameter, step, residual_dot, state_dot, parameter_dot, previous_state_dot, &
            previous_parameter_dot, tangent_state_dot, tangent_parameter_dot, step_dot, augmented_dot, status)
        call assemble_pseudo_arclength_residual( &
            residual + epsilon_fd*residual_dot, state + epsilon_fd*state_dot, &
            parameter_value + epsilon_fd*parameter_dot, previous_state + epsilon_fd*previous_state_dot, &
            previous_parameter + epsilon_fd*previous_parameter_dot, tangent_state + epsilon_fd*tangent_state_dot, &
            tangent_parameter + epsilon_fd*tangent_parameter_dot, step + epsilon_fd*step_dot, &
            augmented_plus, status)
        call assemble_pseudo_arclength_residual( &
            residual - epsilon_fd*residual_dot, state - epsilon_fd*state_dot, &
            parameter_value - epsilon_fd*parameter_dot, previous_state - epsilon_fd*previous_state_dot, &
            previous_parameter - epsilon_fd*previous_parameter_dot, tangent_state - epsilon_fd*tangent_state_dot, &
            tangent_parameter - epsilon_fd*tangent_parameter_dot, step - epsilon_fd*step_dot, &
            augmented_minus, status)
        call record(status%code == FORTSPARSE_OK .and. &
            maxval(abs(augmented_dot - expected_dot)) < 1.0e-14_dp .and. &
            maxval(abs(augmented_dot - (augmented_plus - augmented_minus)/ &
            (2.0_dp*epsilon_fd))) < 2.0e-8_dp, "continuation residual JVP oracle")

        call assemble_pseudo_arclength_residual_vjp( &
            residual, state, parameter_value, previous_state, previous_parameter, tangent_state, &
            tangent_parameter, step, augmented_bar, residual_bar, state_bar, parameter_bar, &
            previous_state_bar, previous_parameter_bar, tangent_state_bar, tangent_parameter_bar, &
            step_bar, status)
        lhs = dot_product(augmented_bar, augmented_dot)
        rhs = dot_product(residual_bar, residual_dot) + dot_product(state_bar, state_dot) + &
            dot_product(parameter_bar, parameter_dot) + dot_product(previous_state_bar, previous_state_dot) + &
            dot_product(previous_parameter_bar, previous_parameter_dot) + dot_product(tangent_state_bar, tangent_state_dot) + &
            dot_product(tangent_parameter_bar, tangent_parameter_dot) + step_bar*step_dot
        call record(status%code == FORTSPARSE_OK .and. abs(lhs - rhs) < 1.0e-13_dp, &
            "continuation residual VJP dot-product oracle")
    end subroutine check_contract

    subroutine record(condition, message)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: message

        call check_condition(condition, message)
        all_passed = all_passed .and. condition
    end subroutine record

end program test_pseudo_arclength_residual
