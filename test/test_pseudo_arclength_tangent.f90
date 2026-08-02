program test_pseudo_arclength_tangent
    use check, only: check_condition, check_summary
    use fortfem_api, only: &
        normalize_pseudo_arclength_tangent, &
        normalize_pseudo_arclength_tangent_jvp, &
        normalize_pseudo_arclength_tangent_vjp
    use fortfem_kinds, only: dp
    use fortsparse, only: FORTSPARSE_INVALID_MATRIX, FORTSPARSE_OK, &
        fortsparse_status_t
    implicit none

    real(dp), parameter :: epsilon_fd = 1.0e-7_dp
    logical :: all_passed

    all_passed = .true.
    call check_contract()
    call check_summary("Pseudo-arclength tangent normalization")
    if (.not. all_passed) error stop 1

contains

    subroutine check_contract()
        integer, parameter :: state_size = 3, parameter_size = 2
        real(dp) :: tangent_state(state_size), tangent_parameter(parameter_size)
        real(dp) :: tangent_state_dot(state_size), tangent_parameter_dot(parameter_size)
        real(dp) :: normalized_state(state_size), normalized_parameter(parameter_size)
        real(dp) :: normalized_state_dot(state_size), normalized_parameter_dot(parameter_size)
        real(dp) :: normalized_state_plus(state_size), normalized_parameter_plus(parameter_size)
        real(dp) :: normalized_state_minus(state_size), normalized_parameter_minus(parameter_size)
        real(dp) :: tangent_state_bar(state_size), tangent_parameter_bar(parameter_size)
        real(dp) :: normalized_state_bar(state_size), normalized_parameter_bar(parameter_size)
        real(dp) :: norm, norm_dot, norm_plus, norm_minus, norm_bar
        real(dp) :: expected_state(state_size), expected_parameter(parameter_size)
        real(dp) :: expected_state_dot(state_size), expected_parameter_dot(parameter_size)
        real(dp) :: lhs, rhs
        type(fortsparse_status_t) :: status

        tangent_state = [0.6_dp, -0.8_dp, 0.3_dp]
        tangent_parameter = [0.4_dp, -0.2_dp]
        tangent_state_dot = [-0.3_dp, 0.5_dp, 0.7_dp]
        tangent_parameter_dot = [0.2_dp, -0.6_dp]
        norm = sqrt(sum(tangent_state**2) + sum(tangent_parameter**2))
        expected_state = tangent_state/norm
        expected_parameter = tangent_parameter/norm

        call normalize_pseudo_arclength_tangent( &
            tangent_state, tangent_parameter, normalized_state, normalized_parameter, &
            norm, status)
        call record(status%code == FORTSPARSE_OK .and. &
            maxval(abs(normalized_state - expected_state)) < 1.0e-14_dp .and. &
            maxval(abs(normalized_parameter - expected_parameter)) < 1.0e-14_dp, &
            "tangent normalization value oracle")

        norm_dot = (dot_product(tangent_state, tangent_state_dot) + &
            dot_product(tangent_parameter, tangent_parameter_dot))/norm
        expected_state_dot = (tangent_state_dot - expected_state*norm_dot)/norm
        expected_parameter_dot = (tangent_parameter_dot - expected_parameter*norm_dot)/norm
        call normalize_pseudo_arclength_tangent_jvp( &
            tangent_state, tangent_parameter, tangent_state_dot, tangent_parameter_dot, &
            normalized_state_dot, normalized_parameter_dot, norm_dot, status)
        call normalize_pseudo_arclength_tangent( &
            tangent_state + epsilon_fd*tangent_state_dot, &
            tangent_parameter + epsilon_fd*tangent_parameter_dot, &
            normalized_state_plus, normalized_parameter_plus, norm_plus, status)
        call normalize_pseudo_arclength_tangent( &
            tangent_state - epsilon_fd*tangent_state_dot, &
            tangent_parameter - epsilon_fd*tangent_parameter_dot, &
            normalized_state_minus, normalized_parameter_minus, norm_minus, status)
        call record(status%code == FORTSPARSE_OK .and. &
            maxval(abs(normalized_state_dot - expected_state_dot)) < 1.0e-14_dp .and. &
            maxval(abs(normalized_parameter_dot - expected_parameter_dot)) < 1.0e-14_dp .and. &
            maxval(abs(normalized_state_dot - (normalized_state_plus - normalized_state_minus)/ &
            (2.0_dp*epsilon_fd))) < 2.0e-8_dp .and. &
            maxval(abs(normalized_parameter_dot - (normalized_parameter_plus - normalized_parameter_minus)/ &
            (2.0_dp*epsilon_fd))) < 2.0e-8_dp .and. &
            abs(norm_dot - (norm_plus - norm_minus)/(2.0_dp*epsilon_fd)) < 2.0e-8_dp, &
            "tangent normalization JVP oracle")

        normalized_state_bar = [0.3_dp, -0.4_dp, 0.2_dp]
        normalized_parameter_bar = [-0.5_dp, 0.7_dp]
        norm_bar = -0.9_dp
        call normalize_pseudo_arclength_tangent_vjp( &
            tangent_state, tangent_parameter, normalized_state_bar, normalized_parameter_bar, &
            norm_bar, tangent_state_bar, tangent_parameter_bar, status)
        lhs = dot_product(normalized_state_bar, normalized_state_dot) + &
            dot_product(normalized_parameter_bar, normalized_parameter_dot) + norm_bar*norm_dot
        rhs = dot_product(tangent_state_bar, tangent_state_dot) + &
            dot_product(tangent_parameter_bar, tangent_parameter_dot)
        call record(status%code == FORTSPARSE_OK .and. abs(lhs - rhs) < 1.0e-13_dp, &
            "tangent normalization VJP dot-product oracle")

        call normalize_pseudo_arclength_tangent( &
            [0.0_dp, 0.0_dp, 0.0_dp], [0.0_dp, 0.0_dp], normalized_state, &
            normalized_parameter, norm, status)
        call record(status%code == FORTSPARSE_INVALID_MATRIX, &
            "zero tangent is rejected")
    end subroutine check_contract

    subroutine record(condition, message)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: message

        call check_condition(condition, message)
        all_passed = all_passed .and. condition
    end subroutine record

end program test_pseudo_arclength_tangent
