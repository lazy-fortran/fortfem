program test_toroidal_diagnostic_hooks
    use check, only: check_condition, check_summary
    use fortfem_kinds, only: dp
    use fortfem_toroidal_diagnostic_hooks, only: &
        near_axis_diagnostic_metadata_t, &
        evaluate_boozer_like_rotational_transform, &
        evaluate_boozer_like_rotational_transform_jvp, &
        evaluate_boozer_like_rotational_transform_vjp, &
        evaluate_near_axis_diagnostic_metadata, &
        evaluate_near_axis_diagnostic_metadata_jvp, &
        evaluate_near_axis_diagnostic_metadata_vjp
    use fortsparse, only: FORTSPARSE_OK, fortsparse_status_t
    implicit none

    integer, parameter :: sample_count = 4, coefficient_count = 5
    real(dp), parameter :: epsilon = 1.0e-7_dp
    real(dp), parameter :: tolerance = 3.0e-8_dp
    real(dp) :: poloidal_rate(sample_count), toroidal_rate(sample_count)
    real(dp) :: weights(sample_count), poloidal_rate_dot(sample_count)
    real(dp) :: toroidal_rate_dot(sample_count), weights_dot(sample_count)
    real(dp) :: iota, iota_dot, iota_plus, iota_minus
    real(dp) :: poloidal_integral, toroidal_integral
    real(dp) :: poloidal_integral_dot, toroidal_integral_dot
    real(dp) :: iota_bar, poloidal_rate_bar(sample_count)
    real(dp) :: toroidal_rate_bar(sample_count), weights_bar(sample_count)
    real(dp) :: lhs, rhs
    real(dp) :: weights_plus(sample_count), weights_minus(sample_count)
    real(dp) :: poloidal_plus(sample_count), poloidal_minus(sample_count)
    real(dp) :: toroidal_plus(sample_count), toroidal_minus(sample_count)
    integer, parameter :: radial_powers(coefficient_count) = [0, 2, 2, 3, 4]
    integer, parameter :: poloidal_modes(coefficient_count) = [0, 1, -1, 2, 0]
    complex(dp), parameter :: coefficients(coefficient_count) = [ &
        cmplx(1.0_dp, 0.2_dp, dp), cmplx(-0.4_dp, 0.1_dp, dp), &
        cmplx(0.3_dp, -0.7_dp, dp), cmplx(0.2_dp, 0.5_dp, dp), &
        cmplx(-0.1_dp, 0.4_dp, dp)]
    complex(dp), parameter :: coefficients_dot(coefficient_count) = [ &
        cmplx(-0.2_dp, 0.3_dp, dp), cmplx(0.1_dp, 0.4_dp, dp), &
        cmplx(-0.3_dp, -0.2_dp, dp), cmplx(0.6_dp, -0.1_dp, dp), &
        cmplx(0.2_dp, -0.5_dp, dp)]
    complex(dp) :: coefficients_plus(coefficient_count), coefficients_minus(coefficient_count)
    complex(dp) :: axis_value_dot, axis_value_plus, axis_value_minus
    complex(dp) :: axis_value_bar, coefficients_bar(coefficient_count)
    real(dp) :: coefficient_norm_dot, coefficient_norm_plus, coefficient_norm_minus
    real(dp) :: axis_norm_dot, axis_norm_plus, axis_norm_minus
    real(dp) :: coefficient_norm_bar, axis_norm_bar
    type(near_axis_diagnostic_metadata_t) :: metadata, metadata_plus, metadata_minus
    type(fortsparse_status_t) :: status
    logical :: all_passed

    all_passed = .true.
    poloidal_rate = [1.0_dp, 2.0_dp, -0.5_dp, 3.0_dp]
    toroidal_rate = [2.0_dp, 4.0_dp, 1.0_dp, 2.5_dp]
    weights = [1.0_dp, 2.0_dp, 1.5_dp, 0.5_dp]
    poloidal_rate_dot = [0.2_dp, -0.3_dp, 0.4_dp, 0.1_dp]
    toroidal_rate_dot = [-0.1_dp, 0.4_dp, -0.2_dp, 0.3_dp]
    weights_dot = [0.1_dp, -0.2_dp, 0.05_dp, 0.15_dp]

    call evaluate_boozer_like_rotational_transform( &
        poloidal_rate, toroidal_rate, weights, iota, poloidal_integral, &
        toroidal_integral, status)
    call record(status%code == FORTSPARSE_OK .and. &
        abs(iota - sum(weights*poloidal_rate)/sum(weights*toroidal_rate)) < tolerance .and. &
        abs(poloidal_integral - sum(weights*poloidal_rate)) < tolerance, &
        "rotational-transform reduction matches independent ratio oracle")

    call evaluate_boozer_like_rotational_transform_jvp( &
        poloidal_rate, toroidal_rate, weights, poloidal_rate_dot, &
        toroidal_rate_dot, weights_dot, iota_dot, poloidal_integral_dot, &
        toroidal_integral_dot, status)
    poloidal_plus = poloidal_rate + epsilon*poloidal_rate_dot
    poloidal_minus = poloidal_rate - epsilon*poloidal_rate_dot
    toroidal_plus = toroidal_rate + epsilon*toroidal_rate_dot
    toroidal_minus = toroidal_rate - epsilon*toroidal_rate_dot
    weights_plus = weights + epsilon*weights_dot
    weights_minus = weights - epsilon*weights_dot
    call evaluate_boozer_like_rotational_transform( &
        poloidal_plus, toroidal_plus, weights_plus, iota_plus, &
        poloidal_integral, toroidal_integral, status)
    call evaluate_boozer_like_rotational_transform( &
        poloidal_minus, toroidal_minus, weights_minus, iota_minus, &
        poloidal_integral, toroidal_integral, status)
    call record(abs(iota_dot - (iota_plus - iota_minus)/(2.0_dp*epsilon)) < tolerance, &
        "rotational-transform JVP matches central differences")

    iota_bar = -0.8_dp
    call evaluate_boozer_like_rotational_transform_vjp( &
        poloidal_rate, toroidal_rate, weights, iota_bar, poloidal_rate_bar, &
        toroidal_rate_bar, weights_bar, status)
    lhs = iota_bar*iota_dot
    rhs = dot_product(poloidal_rate_bar, poloidal_rate_dot) + &
        dot_product(toroidal_rate_bar, toroidal_rate_dot) + &
        dot_product(weights_bar, weights_dot)
    call record(status%code == FORTSPARSE_OK .and. abs(lhs - rhs) < tolerance, &
        "rotational-transform VJP satisfies the real dot-product oracle")

    call evaluate_near_axis_diagnostic_metadata( &
        radial_powers, poloidal_modes, coefficients, metadata, status)
    call record(status%code == FORTSPARSE_OK .and. metadata%coefficient_count == 5 .and. &
        metadata%regular_count == 2 .and. metadata%irregular_count == 3 .and. &
        metadata%axis_mode_count == 1 .and. .not. metadata%axis_regular .and. &
        abs(metadata%coefficient_l2_norm - sqrt(sum(abs(coefficients)**2))) < tolerance .and. &
        abs(metadata%axis_value - coefficients(1)) < tolerance, &
        "near-axis metadata reports supplied mode regularity and norms")

    coefficients_plus = coefficients + epsilon*coefficients_dot
    coefficients_minus = coefficients - epsilon*coefficients_dot
    call evaluate_near_axis_diagnostic_metadata_jvp( &
        radial_powers, poloidal_modes, coefficients, coefficients_dot, &
        coefficient_norm_dot, axis_value_dot, axis_norm_dot, status)
    call evaluate_near_axis_diagnostic_metadata( &
        radial_powers, poloidal_modes, coefficients_plus, metadata_plus, status)
    call evaluate_near_axis_diagnostic_metadata( &
        radial_powers, poloidal_modes, coefficients_minus, metadata_minus, status)
    coefficient_norm_plus = metadata_plus%coefficient_l2_norm
    coefficient_norm_minus = metadata_minus%coefficient_l2_norm
    axis_value_plus = metadata_plus%axis_value
    axis_value_minus = metadata_minus%axis_value
    axis_norm_plus = metadata_plus%axis_value_norm
    axis_norm_minus = metadata_minus%axis_value_norm
    call record(abs(coefficient_norm_dot - (coefficient_norm_plus - coefficient_norm_minus)/ &
        (2.0_dp*epsilon)) < tolerance .and. &
        abs(axis_value_dot - (axis_value_plus - axis_value_minus)/ &
        (2.0_dp*epsilon)) < tolerance .and. &
        abs(axis_norm_dot - (axis_norm_plus - axis_norm_minus)/(2.0_dp*epsilon)) < tolerance, &
        "near-axis smooth metadata JVP matches central differences")

    coefficient_norm_bar = 0.7_dp
    axis_value_bar = cmplx(-0.5_dp, 0.4_dp, dp)
    axis_norm_bar = -0.3_dp
    call evaluate_near_axis_diagnostic_metadata_vjp( &
        radial_powers, poloidal_modes, coefficients, coefficient_norm_bar, &
        axis_value_bar, axis_norm_bar, coefficients_bar, status)
    lhs = coefficient_norm_bar*coefficient_norm_dot + &
        real(conjg(axis_value_bar)*axis_value_dot, dp) + axis_norm_bar*axis_norm_dot
    rhs = real(sum(conjg(coefficients_bar)*coefficients_dot), dp)
    call record(status%code == FORTSPARSE_OK .and. abs(lhs - rhs) < tolerance, &
        "near-axis coefficient VJP satisfies the complex real-part oracle")

    call evaluate_boozer_like_rotational_transform( &
        poloidal_rate, [-2.0_dp, -4.0_dp, -1.0_dp, -2.5_dp], weights, iota, &
        poloidal_integral, toroidal_integral, status)
    call record(status%code == FORTSPARSE_OK .and. iota < 0.0_dp .and. &
        abs(iota - sum(weights*poloidal_rate)/toroidal_integral) < tolerance, &
        "rotational-transform orientation follows supplied rate signs")

    call evaluate_boozer_like_rotational_transform( &
        poloidal_rate, [0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp], weights, iota, &
        poloidal_integral, toroidal_integral, status)
    call record(status%code /= FORTSPARSE_OK, &
        "singular toroidal rate is rejected")

    call evaluate_near_axis_diagnostic_metadata( &
        [-1, 2, 2, 3, 4], poloidal_modes, coefficients, metadata, status)
    call record(status%code /= FORTSPARSE_OK, &
        "negative radial powers are rejected")

    call check_summary("Toroidal diagnostic hooks")
    if (.not. all_passed) error stop 1

contains

    subroutine record(condition, description)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record

end program test_toroidal_diagnostic_hooks
