program test_mixed_wave_long_time
    !! Long-time invariant campaign for the ideal mixed-wave split.
    !!
    !! The oracle is the caller-owned quadratic Hamiltonian, not an internal
    !! quantity returned by the time-stepper.  A second pass with the signed
    !! step checks that accumulated forward evolution is reversible.
    use check, only: check_condition, check_summary
    use fortfem_api, only: advance_mixed_wave_strang
    use fortfem_kinds, only: dp
    use fortsparse, only: fortsparse_status_t
    implicit none

    real(dp), parameter :: mass_q(2, 2) = reshape([ &
        2.0_dp, 0.3_dp, 0.3_dp, 1.5_dp], [2, 2])
    real(dp), parameter :: mass_v(2, 2) = reshape([ &
        1.4_dp, -0.2_dp, -0.2_dp, 1.8_dp], [2, 2])
    real(dp), parameter :: coupling_a(2, 2) = reshape([ &
        0.7_dp, -0.3_dp, 0.4_dp, 1.1_dp], [2, 2])
    real(dp), parameter :: coupling_b(2, 2) = reshape([ &
        -0.2_dp, 0.8_dp, -0.6_dp, 0.5_dp], [2, 2])
    real(dp), parameter :: initial_q(2) = [0.7_dp, -0.4_dp]
    real(dp), parameter :: initial_v(2) = [-0.2_dp, 0.9_dp]
    real(dp), parameter :: time_step = 0.17_dp
    integer, parameter :: step_count = 2000
    real(dp), parameter :: energy_tolerance = 5.0e-10_dp
    real(dp), parameter :: reversibility_tolerance = 5.0e-10_dp
    real(dp) :: q(2), v(2), q_initial(2), v_initial(2)
    real(dp) :: energy_initial, energy_error, maximum_energy_error
    integer :: step
    type(fortsparse_status_t) :: status

    q_initial = initial_q
    v_initial = initial_v
    q = q_initial
    v = v_initial
    energy_initial = wave_energy(q, v)
    maximum_energy_error = 0.0_dp

    do step = 1, step_count
        call advance_mixed_wave_strang( &
            mass_q, mass_v, coupling_a, coupling_b, time_step, q, v, status)
        call check_condition(status%code == 0, &
            "long-time Strang evolution accepts every forward step")
        if (status%code /= 0) error stop 1
        energy_error = abs(wave_energy(q, v) - energy_initial)
        maximum_energy_error = max(maximum_energy_error, energy_error)
    end do

    call check_condition(maximum_energy_error < energy_tolerance, &
        "long-time Strang evolution preserves the independent quadratic energy")

    do step = 1, step_count
        call advance_mixed_wave_strang( &
            mass_q, mass_v, coupling_a, coupling_b, -time_step, q, v, status)
        call check_condition(status%code == 0, &
            "long-time signed-step reversal accepts every backward step")
        if (status%code /= 0) error stop 1
    end do

    call check_condition( &
        maxval(abs(q - q_initial)) < reversibility_tolerance .and. &
        maxval(abs(v - v_initial)) < reversibility_tolerance, &
        "long-time Strang evolution reverses to the initial state")
    call check_summary("Long-time mixed-wave structure preservation")

contains

    pure function wave_energy(q_state, v_state) result(energy)
        real(dp), intent(in) :: q_state(:), v_state(:)
        real(dp) :: energy

        energy = 0.5_dp*(dot_product(q_state, matmul(mass_q, q_state)) + &
            dot_product(v_state, matmul(mass_v, v_state)))
    end function wave_energy

end program test_mixed_wave_long_time
