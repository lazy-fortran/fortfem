module fortfem_mixed_wave_energy_ledger
    !! Neutral mixed-wave energy, helicity, and signed-port ledger.
    !!
    !! Samples are caller-owned scalar contractions.  The ledger intentionally
    !! imposes no constitutive closure: positive quadrature weights combine the
    !! quadratic energy, canonical pressure/displacement and velocity/momentum
    !! pairings, and the supplied signed boundary power.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortfem_kinds, only: dp
    use fortsparse, only: fortsparse_status_t, status_set, &
        FORTSPARSE_INVALID_MATRIX, FORTSPARSE_OK
    implicit none
    private

    public :: evaluate_mixed_wave_energy_helicity_port_ledger
    public :: evaluate_mixed_wave_energy_helicity_port_ledger_jvp
    public :: evaluate_mixed_wave_energy_helicity_port_ledger_vjp

contains

    subroutine evaluate_mixed_wave_energy_helicity_port_ledger( &
            weights, pressure, velocity, displacement, momentum, stress, &
            boundary_power, energy, helicity, port_power, total, status)
        real(dp), intent(in) :: weights(:), pressure(:), velocity(:), displacement(:)
        real(dp), intent(in) :: momentum(:), stress(:), boundary_power(:)
        real(dp), intent(out) :: energy, helicity, port_power, total
        type(fortsparse_status_t), intent(out) :: status
        integer :: sample

        energy = 0.0_dp
        helicity = 0.0_dp
        port_power = 0.0_dp
        total = 0.0_dp
        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            "mixed-wave ledger received incompatible samples")
        if (.not. valid_values(weights, pressure, velocity, displacement, momentum, &
            stress, boundary_power)) return
        do sample = 1, size(weights)
            energy = energy + 0.5_dp*weights(sample)*(pressure(sample)**2 + &
                velocity(sample)**2 + displacement(sample)**2 + momentum(sample)**2 + &
                stress(sample)**2)
            helicity = helicity + weights(sample)*(pressure(sample)*displacement(sample) + &
                velocity(sample)*momentum(sample))
            port_power = port_power + weights(sample)*boundary_power(sample)
        end do
        total = energy + helicity + port_power
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine evaluate_mixed_wave_energy_helicity_port_ledger

    subroutine evaluate_mixed_wave_energy_helicity_port_ledger_jvp( &
            weights, pressure, velocity, displacement, momentum, stress, boundary_power, &
            weights_dot, pressure_dot, velocity_dot, displacement_dot, momentum_dot, &
            stress_dot, boundary_power_dot, energy_dot, helicity_dot, port_power_dot, &
            total_dot, status)
        real(dp), intent(in) :: weights(:), pressure(:), velocity(:), displacement(:)
        real(dp), intent(in) :: momentum(:), stress(:), boundary_power(:)
        real(dp), intent(in) :: weights_dot(:), pressure_dot(:), velocity_dot(:)
        real(dp), intent(in) :: displacement_dot(:), momentum_dot(:), stress_dot(:)
        real(dp), intent(in) :: boundary_power_dot(:)
        real(dp), intent(out) :: energy_dot, helicity_dot, port_power_dot, total_dot
        type(fortsparse_status_t), intent(out) :: status
        integer :: sample
        real(dp) :: energy_sample, helicity_sample

        energy_dot = 0.0_dp
        helicity_dot = 0.0_dp
        port_power_dot = 0.0_dp
        total_dot = 0.0_dp
        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            "mixed-wave ledger JVP received incompatible samples")
        if (.not. valid_values(weights, pressure, velocity, displacement, momentum, &
            stress, boundary_power) .or. .not. valid_directions(weights_dot, pressure_dot, &
            velocity_dot, displacement_dot, momentum_dot, stress_dot, boundary_power_dot, &
            size(weights))) return
        do sample = 1, size(weights)
            energy_sample = 0.5_dp*(pressure(sample)**2 + velocity(sample)**2 + &
                displacement(sample)**2 + momentum(sample)**2 + stress(sample)**2)
            helicity_sample = pressure(sample)*displacement(sample) + &
                velocity(sample)*momentum(sample)
            energy_dot = energy_dot + weights_dot(sample)*energy_sample + weights(sample)*( &
                pressure(sample)*pressure_dot(sample) + velocity(sample)*velocity_dot(sample) + &
                displacement(sample)*displacement_dot(sample) + momentum(sample)*momentum_dot(sample) + &
                stress(sample)*stress_dot(sample))
            helicity_dot = helicity_dot + weights_dot(sample)*helicity_sample + weights(sample)*( &
                pressure_dot(sample)*displacement(sample) + pressure(sample)*displacement_dot(sample) + &
                velocity_dot(sample)*momentum(sample) + velocity(sample)*momentum_dot(sample))
            port_power_dot = port_power_dot + weights_dot(sample)*boundary_power(sample) + &
                weights(sample)*boundary_power_dot(sample)
        end do
        total_dot = energy_dot + helicity_dot + port_power_dot
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine evaluate_mixed_wave_energy_helicity_port_ledger_jvp

    subroutine evaluate_mixed_wave_energy_helicity_port_ledger_vjp( &
            weights, pressure, velocity, displacement, momentum, stress, boundary_power, &
            energy_bar, helicity_bar, port_power_bar, total_bar, weights_bar, pressure_bar, &
            velocity_bar, displacement_bar, momentum_bar, stress_bar, boundary_power_bar, status)
        real(dp), intent(in) :: weights(:), pressure(:), velocity(:), displacement(:)
        real(dp), intent(in) :: momentum(:), stress(:), boundary_power(:)
        real(dp), intent(in) :: energy_bar, helicity_bar, port_power_bar, total_bar
        real(dp), intent(out) :: weights_bar(:), pressure_bar(:), velocity_bar(:)
        real(dp), intent(out) :: displacement_bar(:), momentum_bar(:), stress_bar(:)
        real(dp), intent(out) :: boundary_power_bar(:)
        type(fortsparse_status_t), intent(out) :: status
        integer :: sample
        real(dp) :: effective_energy_bar, effective_helicity_bar, effective_port_bar
        real(dp) :: energy_sample, helicity_sample

        weights_bar = 0.0_dp
        pressure_bar = 0.0_dp
        velocity_bar = 0.0_dp
        displacement_bar = 0.0_dp
        momentum_bar = 0.0_dp
        stress_bar = 0.0_dp
        boundary_power_bar = 0.0_dp
        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            "mixed-wave ledger VJP received incompatible samples")
        if (.not. valid_values(weights, pressure, velocity, displacement, momentum, &
            stress, boundary_power) .or. size(weights_bar) /= size(weights) .or. &
            size(pressure_bar) /= size(weights) .or. size(velocity_bar) /= size(weights) .or. &
            size(displacement_bar) /= size(weights) .or. size(momentum_bar) /= size(weights) .or. &
            size(stress_bar) /= size(weights) .or. size(boundary_power_bar) /= size(weights) .or. &
            .not. ieee_is_finite(energy_bar) .or. .not. ieee_is_finite(helicity_bar) .or. &
            .not. ieee_is_finite(port_power_bar) .or. .not. ieee_is_finite(total_bar)) return
        effective_energy_bar = energy_bar + total_bar
        effective_helicity_bar = helicity_bar + total_bar
        effective_port_bar = port_power_bar + total_bar
        do sample = 1, size(weights)
            energy_sample = 0.5_dp*(pressure(sample)**2 + velocity(sample)**2 + &
                displacement(sample)**2 + momentum(sample)**2 + stress(sample)**2)
            helicity_sample = pressure(sample)*displacement(sample) + &
                velocity(sample)*momentum(sample)
            weights_bar(sample) = effective_energy_bar*energy_sample + &
                effective_helicity_bar*helicity_sample + effective_port_bar*boundary_power(sample)
            pressure_bar(sample) = weights(sample)*(effective_energy_bar*pressure(sample) + &
                effective_helicity_bar*displacement(sample))
            displacement_bar(sample) = weights(sample)*(effective_energy_bar*displacement(sample) + &
                effective_helicity_bar*pressure(sample))
            velocity_bar(sample) = weights(sample)*(effective_energy_bar*velocity(sample) + &
                effective_helicity_bar*momentum(sample))
            momentum_bar(sample) = weights(sample)*(effective_energy_bar*momentum(sample) + &
                effective_helicity_bar*velocity(sample))
            stress_bar(sample) = weights(sample)*effective_energy_bar*stress(sample)
            boundary_power_bar(sample) = weights(sample)*effective_port_bar
        end do
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine evaluate_mixed_wave_energy_helicity_port_ledger_vjp

    logical function valid_values(weights, pressure, velocity, displacement, momentum, &
            stress, boundary_power) result(valid)
        real(dp), intent(in) :: weights(:), pressure(:), velocity(:), displacement(:)
        real(dp), intent(in) :: momentum(:), stress(:), boundary_power(:)

        valid = size(weights) > 0 .and. size(pressure) == size(weights) .and. &
            size(velocity) == size(weights) .and. size(displacement) == size(weights) .and. &
            size(momentum) == size(weights) .and. size(stress) == size(weights) .and. &
            size(boundary_power) == size(weights) .and. all(ieee_is_finite(weights)) .and. &
            all(ieee_is_finite(pressure)) .and. all(ieee_is_finite(velocity)) .and. &
            all(ieee_is_finite(displacement)) .and. all(ieee_is_finite(momentum)) .and. &
            all(ieee_is_finite(stress)) .and. all(ieee_is_finite(boundary_power)) .and. &
            all(weights > 0.0_dp)
    end function valid_values

    logical function valid_directions(weights_dot, pressure_dot, velocity_dot, &
            displacement_dot, momentum_dot, stress_dot, boundary_power_dot, count) result(valid)
        real(dp), intent(in) :: weights_dot(:), pressure_dot(:), velocity_dot(:)
        real(dp), intent(in) :: displacement_dot(:), momentum_dot(:), stress_dot(:)
        real(dp), intent(in) :: boundary_power_dot(:)
        integer, intent(in) :: count

        valid = size(weights_dot) == count .and. size(pressure_dot) == count .and. &
            size(velocity_dot) == count .and. size(displacement_dot) == count .and. &
            size(momentum_dot) == count .and. size(stress_dot) == count .and. &
            size(boundary_power_dot) == count .and. all(ieee_is_finite(weights_dot)) .and. &
            all(ieee_is_finite(pressure_dot)) .and. all(ieee_is_finite(velocity_dot)) .and. &
            all(ieee_is_finite(displacement_dot)) .and. all(ieee_is_finite(momentum_dot)) .and. &
            all(ieee_is_finite(stress_dot)) .and. all(ieee_is_finite(boundary_power_dot))
    end function valid_directions

end module fortfem_mixed_wave_energy_ledger
