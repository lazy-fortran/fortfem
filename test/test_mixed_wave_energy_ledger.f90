program test_mixed_wave_energy_ledger
    use check, only: check_condition, check_summary
    use fortfem_kinds, only: dp
    use fortfem_mixed_wave_energy_ledger, only: &
        evaluate_mixed_wave_energy_helicity_port_ledger, &
        evaluate_mixed_wave_energy_helicity_port_ledger_jvp, &
        evaluate_mixed_wave_energy_helicity_port_ledger_vjp
    use fortsparse, only: fortsparse_status_t
    implicit none

    real(dp), parameter :: weights(3) = [0.5_dp, 1.0_dp, 1.5_dp]
    real(dp), parameter :: pressure(3) = [1.0_dp, -0.5_dp, 0.8_dp]
    real(dp), parameter :: velocity(3) = [0.2_dp, 0.6_dp, -0.4_dp]
    real(dp), parameter :: displacement(3) = [-0.3_dp, 0.7_dp, 0.1_dp]
    real(dp), parameter :: momentum(3) = [0.9_dp, -0.2_dp, 0.5_dp]
    real(dp), parameter :: stress(3) = [0.4_dp, -0.8_dp, 0.3_dp]
    real(dp), parameter :: boundary_power(3) = [-2.0_dp, 0.5_dp, -0.25_dp]
    real(dp), parameter :: weights_dot(3) = [0.1_dp, -0.2_dp, 0.05_dp]
    real(dp), parameter :: pressure_dot(3) = [0.3_dp, -0.1_dp, 0.2_dp]
    real(dp), parameter :: velocity_dot(3) = [-0.2_dp, 0.4_dp, 0.1_dp]
    real(dp), parameter :: displacement_dot(3) = [0.15_dp, 0.05_dp, -0.1_dp]
    real(dp), parameter :: momentum_dot(3) = [0.1_dp, -0.3_dp, 0.2_dp]
    real(dp), parameter :: stress_dot(3) = [-0.05_dp, 0.2_dp, 0.1_dp]
    real(dp), parameter :: boundary_power_dot(3) = [0.2_dp, -0.1_dp, 0.3_dp]
    real(dp), parameter :: finite_difference_step = 1.0e-7_dp
    real(dp) :: energy, helicity, port_power, total
    real(dp) :: energy_dot, helicity_dot, port_power_dot, total_dot
    real(dp) :: energy_plus, helicity_plus, port_plus, total_plus
    real(dp) :: energy_minus, helicity_minus, port_minus, total_minus
    real(dp) :: weights_bar(3), pressure_bar(3), velocity_bar(3)
    real(dp) :: displacement_bar(3), momentum_bar(3), stress_bar(3), boundary_power_bar(3)
    real(dp), parameter :: output_bar(4) = [0.7_dp, -0.4_dp, 0.3_dp, 0.6_dp]
    real(dp) :: lhs, rhs, expected_energy, expected_helicity, expected_port
    real(dp), parameter :: invalid_weights(3) = [0.5_dp, 0.0_dp, 1.5_dp]
    type(fortsparse_status_t) :: status
    integer :: sample

    expected_energy = 0.0_dp
    expected_helicity = 0.0_dp
    expected_port = 0.0_dp
    do sample = 1, 3
        expected_energy = expected_energy + 0.5_dp*weights(sample)*(pressure(sample)**2 + &
            velocity(sample)**2 + displacement(sample)**2 + momentum(sample)**2 + stress(sample)**2)
        expected_helicity = expected_helicity + weights(sample)*(pressure(sample)*displacement(sample) + &
            velocity(sample)*momentum(sample))
        expected_port = expected_port + weights(sample)*boundary_power(sample)
    end do
    call evaluate_mixed_wave_energy_helicity_port_ledger(weights, pressure, velocity, displacement, &
        momentum, stress, boundary_power, energy, helicity, port_power, total, status)
    call check_condition(status%code == 0 .and. abs(energy - expected_energy) < 1.0e-14_dp .and. &
        abs(helicity - expected_helicity) < 1.0e-14_dp .and. abs(port_power - expected_port) < 1.0e-14_dp .and. &
        abs(total - expected_energy - expected_helicity - expected_port) < 1.0e-14_dp, &
        "mixed-wave ledger matches an independent nested-loop oracle")
    call check_condition(port_power < 0.0_dp, "signed boundary power is preserved")

    call evaluate_mixed_wave_energy_helicity_port_ledger_jvp(weights, pressure, velocity, displacement, momentum, &
        stress, boundary_power, weights_dot, pressure_dot, velocity_dot, displacement_dot, momentum_dot, stress_dot, &
        boundary_power_dot, energy_dot, helicity_dot, port_power_dot, total_dot, status)
    call evaluate_mixed_wave_energy_helicity_port_ledger(weights + finite_difference_step*weights_dot, &
        pressure + finite_difference_step*pressure_dot, velocity + finite_difference_step*velocity_dot, &
        displacement + finite_difference_step*displacement_dot, momentum + finite_difference_step*momentum_dot, &
        stress + finite_difference_step*stress_dot, boundary_power + finite_difference_step*boundary_power_dot, &
        energy_plus, helicity_plus, port_plus, total_plus, status)
    call evaluate_mixed_wave_energy_helicity_port_ledger(weights - finite_difference_step*weights_dot, &
        pressure - finite_difference_step*pressure_dot, velocity - finite_difference_step*velocity_dot, &
        displacement - finite_difference_step*displacement_dot, momentum - finite_difference_step*momentum_dot, &
        stress - finite_difference_step*stress_dot, boundary_power - finite_difference_step*boundary_power_dot, &
        energy_minus, helicity_minus, port_minus, total_minus, status)
    call check_condition(maxval(abs([energy_dot, helicity_dot, port_power_dot, total_dot] - &
        [(energy_plus - energy_minus)/(2.0_dp*finite_difference_step), &
         (helicity_plus - helicity_minus)/(2.0_dp*finite_difference_step), &
         (port_plus - port_minus)/(2.0_dp*finite_difference_step), &
         (total_plus - total_minus)/(2.0_dp*finite_difference_step)])) < 1.0e-8_dp, &
        "mixed-wave ledger JVP matches central differences")

    call evaluate_mixed_wave_energy_helicity_port_ledger_vjp(weights, pressure, velocity, displacement, momentum, &
        stress, boundary_power, output_bar(1), output_bar(2), output_bar(3), output_bar(4), weights_bar, &
        pressure_bar, velocity_bar, displacement_bar, momentum_bar, stress_bar, boundary_power_bar, status)
    lhs = output_bar(1)*energy_dot + output_bar(2)*helicity_dot + output_bar(3)*port_power_dot + &
        output_bar(4)*total_dot
    rhs = dot_product(weights_bar, weights_dot) + dot_product(pressure_bar, pressure_dot) + &
        dot_product(velocity_bar, velocity_dot) + dot_product(displacement_bar, displacement_dot) + &
        dot_product(momentum_bar, momentum_dot) + dot_product(stress_bar, stress_dot) + &
        dot_product(boundary_power_bar, boundary_power_dot)
    call check_condition(status%code == 0 .and. abs(lhs - rhs) < 1.0e-13_dp, &
        "mixed-wave ledger VJP satisfies the real transpose oracle")
    call evaluate_mixed_wave_energy_helicity_port_ledger(invalid_weights, pressure, velocity, displacement, &
        momentum, stress, boundary_power, energy, helicity, port_power, total, status)
    call check_condition(status%code /= 0, "mixed-wave ledger rejects nonpositive weights")
    call check_summary("mixed-wave energy helicity port ledger")
end program test_mixed_wave_energy_ledger
