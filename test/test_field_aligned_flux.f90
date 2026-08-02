program test_field_aligned_flux
    use check, only: check_condition, check_summary
    use fortfem_feec, only: &
        evaluate_field_aligned_flux, evaluate_field_aligned_flux_jvp, &
        evaluate_field_aligned_flux_vjp
    use fortfem_kinds, only: dp
    use fortsparse, only: fortsparse_status_t
    implicit none

    real(dp), parameter :: parallel_coefficient = 100.0_dp
    real(dp), parameter :: perpendicular_coefficient = 1.0_dp
    real(dp), parameter :: direction(3) = [0.6_dp, 0.8_dp, 0.0_dp]
    real(dp), parameter :: gradient(3) = [2.0_dp, -1.0_dp, 3.0_dp]
    real(dp), parameter :: expected_flux(3) = [25.76_dp, 30.68_dp, 3.0_dp]
    real(dp), parameter :: coefficient_dot = -0.7_dp
    real(dp), parameter :: perpendicular_dot = 0.2_dp
    real(dp), parameter :: direction_dot(3) = [0.01_dp, -0.0075_dp, 0.0_dp]
    real(dp), parameter :: gradient_dot(3) = [-0.3_dp, 0.4_dp, 0.1_dp]
    real(dp) :: flux(3), flux_dot(3), flux_plus(3), flux_minus(3)
    real(dp) :: parallel_bar, perpendicular_bar, direction_bar(3)
    real(dp) :: gradient_bar(3), left, right, epsilon
    real(dp), parameter :: flux_bar(3) = [0.4_dp, -0.6_dp, 0.8_dp]
    real(dp), parameter :: bad_direction(3) = [0.6_dp, 0.8_dp, 0.1_dp]
    type(fortsparse_status_t) :: status

    call evaluate_field_aligned_flux( &
        parallel_coefficient, perpendicular_coefficient, direction, gradient, &
        flux, status)
    call check_condition(status%code == 0, &
        "field-aligned flux accepts a unit direction")
    call check_condition(maxval(abs(flux - expected_flux)) < 1.0e-13_dp, &
        "field-aligned flux matches the anisotropic tensor oracle")

    call evaluate_field_aligned_flux_jvp( &
        parallel_coefficient, perpendicular_coefficient, direction, gradient, &
        coefficient_dot, perpendicular_dot, direction_dot, gradient_dot, &
        flux_dot, status)
    epsilon = 1.0e-6_dp
    call evaluate_field_aligned_flux( &
        parallel_coefficient + epsilon*coefficient_dot, &
        perpendicular_coefficient + epsilon*perpendicular_dot, &
        direction + epsilon*direction_dot, gradient + epsilon*gradient_dot, &
        flux_plus, status)
    call evaluate_field_aligned_flux( &
        parallel_coefficient - epsilon*coefficient_dot, &
        perpendicular_coefficient - epsilon*perpendicular_dot, &
        direction - epsilon*direction_dot, gradient - epsilon*gradient_dot, &
        flux_minus, status)
    call check_condition(maxval(abs(flux_dot - &
        (flux_plus - flux_minus)/(2.0_dp*epsilon))) < 1.0e-7_dp, &
        "field-aligned flux JVP matches central differences")

    call evaluate_field_aligned_flux_vjp( &
        parallel_coefficient, perpendicular_coefficient, direction, gradient, &
        flux_bar, parallel_bar, perpendicular_bar, direction_bar, gradient_bar, &
        status)
    left = dot_product(flux_bar, flux_dot)
    right = parallel_bar*coefficient_dot + perpendicular_bar*perpendicular_dot + &
        dot_product(direction_bar, direction_dot) + &
        dot_product(gradient_bar, gradient_dot)
    call check_condition(status%code == 0 .and. abs(left - right) < 1.0e-12_dp, &
        "field-aligned flux VJP satisfies the real dot-product identity")

    call evaluate_field_aligned_flux( &
        parallel_coefficient, perpendicular_coefficient, bad_direction, gradient, &
        flux, status)
    call check_condition(status%code /= 0, &
        "field-aligned flux rejects a non-unit direction")
    call check_summary("field-aligned flux")
end program test_field_aligned_flux
