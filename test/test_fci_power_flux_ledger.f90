program test_fci_power_flux_ledger
    use check, only: check_condition, check_summary
    use fortfem_feec, only: evaluate_fci_power_flux_ledger, &
        evaluate_fci_power_flux_ledger_jvp, evaluate_fci_power_flux_ledger_vjp
    use fortfem_kinds, only: dp
    use fortsparse, only: fortsparse_status_t
    implicit none

    real(dp), parameter :: gradient(2) = [1.0_dp, -2.0_dp]
    real(dp), parameter :: parallel_coefficient(2) = [2.0_dp, 3.0_dp]
    real(dp), parameter :: staggered_volumes(2) = [4.0_dp, 5.0_dp]
    real(dp), parameter :: field(3) = [1.5_dp, -0.5_dp, 2.0_dp]
    real(dp), parameter :: perpendicular_action(3) = [-2.0_dp, 1.0_dp, -3.0_dp]
    real(dp), parameter :: canonical_volumes(3) = [2.0_dp, 3.0_dp, 4.0_dp]
    real(dp), parameter :: gradient_dot(2) = [0.2_dp, -0.4_dp]
    real(dp), parameter :: parallel_coefficient_dot(2) = [-0.1_dp, 0.3_dp]
    real(dp), parameter :: staggered_volumes_dot(2) = [0.5_dp, -0.2_dp]
    real(dp), parameter :: field_dot(3) = [0.1_dp, -0.2_dp, 0.3_dp]
    real(dp), parameter :: perpendicular_action_dot(3) = [0.4_dp, 0.2_dp, -0.1_dp]
    real(dp), parameter :: canonical_volumes_dot(3) = [-0.1_dp, 0.2_dp, 0.3_dp]
    real(dp), parameter :: parallel_flux_bar(2) = [0.4_dp, -0.8_dp]
    real(dp), parameter :: parallel_power_bar = 0.7_dp
    real(dp), parameter :: perpendicular_power_bar = -0.3_dp
    real(dp), parameter :: total_power_bar = 0.2_dp
    real(dp), parameter :: step = 1.0e-7_dp
    real(dp) :: parallel_flux(2), parallel_power, perpendicular_power, total_power
    real(dp) :: parallel_flux_dot(2), parallel_power_dot, perpendicular_power_dot
    real(dp) :: total_power_dot
    real(dp) :: parallel_flux_plus(2), parallel_flux_minus(2)
    real(dp) :: parallel_power_plus, perpendicular_power_plus, total_power_plus
    real(dp) :: parallel_power_minus, perpendicular_power_minus, total_power_minus
    real(dp) :: gradient_bar(2), parallel_coefficient_bar(2)
    real(dp) :: staggered_volumes_bar(2), field_bar(3)
    real(dp) :: perpendicular_action_bar(3), canonical_volumes_bar(3)
    real(dp) :: left, right
    real(dp) :: expected_parallel_power, expected_perpendicular_power
    real(dp) :: bad_coefficient(2)
    type(fortsparse_status_t) :: status

    expected_parallel_power = -sum(staggered_volumes*parallel_coefficient*gradient**2)
    expected_perpendicular_power = sum( &
        canonical_volumes*field*perpendicular_action)
    call evaluate_fci_power_flux_ledger( &
        gradient, parallel_coefficient, staggered_volumes, field, &
        perpendicular_action, canonical_volumes, parallel_flux, parallel_power, &
        perpendicular_power, total_power, status)
    call check_condition(status%code == 0, &
        "FCI power/flux ledger accepts positive FCI data")
    call check_condition(maxval(abs(parallel_flux - &
        [-2.0_dp, 6.0_dp])) < 1.0e-14_dp .and. &
        abs(parallel_power - expected_parallel_power) < 1.0e-14_dp .and. &
        abs(perpendicular_power - expected_perpendicular_power) < 1.0e-14_dp .and. &
        abs(total_power - expected_parallel_power - expected_perpendicular_power) < &
        1.0e-14_dp, "parallel flux and split powers match an independent oracle")

    call evaluate_fci_power_flux_ledger_jvp( &
        gradient, parallel_coefficient, staggered_volumes, field, &
        perpendicular_action, canonical_volumes, gradient_dot, &
        parallel_coefficient_dot, staggered_volumes_dot, field_dot, &
        perpendicular_action_dot, canonical_volumes_dot, parallel_flux_dot, &
        parallel_power_dot, perpendicular_power_dot, total_power_dot, status)
    call evaluate_fci_power_flux_ledger( &
        gradient + step*gradient_dot, parallel_coefficient + step*parallel_coefficient_dot, &
        staggered_volumes + step*staggered_volumes_dot, field + step*field_dot, &
        perpendicular_action + step*perpendicular_action_dot, &
        canonical_volumes + step*canonical_volumes_dot, parallel_flux_plus, &
        parallel_power_plus, perpendicular_power_plus, total_power_plus, status)
    call evaluate_fci_power_flux_ledger( &
        gradient - step*gradient_dot, parallel_coefficient - step*parallel_coefficient_dot, &
        staggered_volumes - step*staggered_volumes_dot, field - step*field_dot, &
        perpendicular_action - step*perpendicular_action_dot, &
        canonical_volumes - step*canonical_volumes_dot, parallel_flux_minus, &
        parallel_power_minus, perpendicular_power_minus, total_power_minus, status)
    call check_condition(maxval(abs(parallel_flux_dot - &
        (parallel_flux_plus - parallel_flux_minus)/(2.0_dp*step))) < 2.0e-8_dp .and. &
        maxval(abs([parallel_power_dot, perpendicular_power_dot, total_power_dot] - &
        [(parallel_power_plus - parallel_power_minus)/(2.0_dp*step), &
        (perpendicular_power_plus - perpendicular_power_minus)/(2.0_dp*step), &
        (total_power_plus - total_power_minus)/(2.0_dp*step)])) < 5.0e-8_dp, &
        "FCI power/flux ledger JVP matches independent central differences")

    call evaluate_fci_power_flux_ledger_vjp( &
        gradient, parallel_coefficient, staggered_volumes, field, &
        perpendicular_action, canonical_volumes, parallel_flux_bar, &
        parallel_power_bar, perpendicular_power_bar, total_power_bar, gradient_bar, &
        parallel_coefficient_bar, staggered_volumes_bar, field_bar, &
        perpendicular_action_bar, canonical_volumes_bar, status)
    left = dot_product(parallel_flux_bar, parallel_flux_dot) + &
        parallel_power_bar*parallel_power_dot + &
        perpendicular_power_bar*perpendicular_power_dot + total_power_bar*total_power_dot
    right = dot_product(gradient_bar, gradient_dot) + &
        dot_product(parallel_coefficient_bar, parallel_coefficient_dot) + &
        dot_product(staggered_volumes_bar, staggered_volumes_dot) + &
        dot_product(field_bar, field_dot) + &
        dot_product(perpendicular_action_bar, perpendicular_action_dot) + &
        dot_product(canonical_volumes_bar, canonical_volumes_dot)
    call check_condition(status%code == 0 .and. abs(left - right) < 1.0e-12_dp, &
        "FCI power/flux ledger VJP satisfies the real transpose oracle")

    bad_coefficient = parallel_coefficient
    bad_coefficient(1) = -1.0_dp
    call evaluate_fci_power_flux_ledger( &
        gradient, bad_coefficient, staggered_volumes, field, perpendicular_action, &
        canonical_volumes, parallel_flux, parallel_power, perpendicular_power, &
        total_power, status)
    call check_condition(status%code /= 0, &
        "FCI power/flux ledger rejects a negative parallel coefficient")
    call check_summary("FCI power/flux ledger")
end program test_fci_power_flux_ledger
