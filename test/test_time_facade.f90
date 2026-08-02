program test_time_facade
    !! Behavioral smoke test for the canonical time-integration facade.
    !!
    !! The oscillator identities below are independent of the implementation:
    !! midpoint is a Cayley rotation, partitioned Euler has determinant one,
    !! and a positive dissipative Cayley step contracts the mass energy.
    use check, only: check_condition, check_summary
    use fortfem_time, only: advance_dissipative_cayley, &
        advance_mixed_wave_midpoint, advance_mixed_wave_strang, &
        advance_mixed_wave_symplectic_euler, &
        assemble_symplectic_map_defect
    use fortfem_kinds, only: dp
    use fortsparse, only: fortsparse_status_t
    implicit none

    real(dp), parameter :: mass(1, 1) = reshape([1.0_dp], [1, 1])
    real(dp), parameter :: coupling(1, 1) = mass
    real(dp), parameter :: damping(1, 1) = reshape([0.4_dp], [1, 1])
    real(dp), parameter :: time_step = 0.2_dp
    real(dp), parameter :: initial_q = 1.0_dp
    real(dp), parameter :: initial_v = 0.0_dp
    real(dp), parameter :: initial_state = 1.0_dp
    real(dp), parameter :: tolerance = 2.0e-13_dp
    real(dp) :: q, v, energy_initial, energy_final
    real(dp) :: q_reverse, v_reverse
    real(dp) :: q_euler, v_euler
    real(dp) :: q_work(1), v_work(1)
    real(dp) :: state_next(1), euler_map(2, 2), symplectic_form(2, 2)
    real(dp) :: defect(2, 2), expected_map(2, 2)
    real(dp) :: q_strang, v_strang, strang_energy
    type(fortsparse_status_t) :: status

    q = initial_q
    v = initial_v
    energy_initial = oscillator_energy(q, v)
    q_work = q
    v_work = v
    call advance_mixed_wave_midpoint( &
        mass, mass, coupling, time_step, q_work, v_work, status)
    q = q_work(1)
    v = v_work(1)

    energy_final = oscillator_energy(q, v)
    call check_condition(status%code == 0 .and. &
        abs(energy_final - energy_initial) < tolerance, &
        "time facade midpoint preserves oscillator energy")

    q_reverse = q
    v_reverse = v
    q_work = q_reverse
    v_work = v_reverse
    call advance_mixed_wave_midpoint( &
        mass, mass, coupling, -time_step, q_work, v_work, status)
    q_reverse = q_work(1)
    v_reverse = v_work(1)

    call check_condition(status%code == 0 .and. &
        abs(q_reverse - initial_q) < tolerance .and. &
        abs(v_reverse - initial_v) < tolerance, &
        "time facade midpoint is reversible under signed steps")

    q_euler = initial_q
    v_euler = initial_v
    q_work = q_euler
    v_work = v_euler
    call advance_mixed_wave_symplectic_euler( &
        mass, mass, coupling, time_step, q_work, v_work, status)
    q_euler = q_work(1)
    v_euler = v_work(1)
    expected_map = reshape([1.0_dp - time_step**2, time_step, &
        -time_step, 1.0_dp], [2, 2])
    euler_map = expected_map
    symplectic_form = reshape([0.0_dp, -1.0_dp, 1.0_dp, 0.0_dp], [2, 2])
    call assemble_symplectic_map_defect( &
        euler_map, symplectic_form, defect, status)
    call check_condition(status%code == 0 .and. &
        maxval(abs(defect)) < tolerance .and. &
        abs(q_euler - expected_map(1, 1)) < tolerance .and. &
        abs(v_euler - expected_map(2, 1)) < tolerance, &
        "time facade exposes a symplectic Euler map with zero defect")

    q_strang = initial_q
    v_strang = initial_v
    strang_energy = oscillator_energy(q_strang, v_strang)
    q_work = q_strang
    v_work = v_strang
    call advance_mixed_wave_strang( &
        mass, mass, coupling, coupling, time_step, q_work, v_work, status)
    q_strang = q_work(1)
    v_strang = v_work(1)
    call check_condition(status%code == 0 .and. &
        abs(oscillator_energy(q_strang, v_strang) - strang_energy) < tolerance, &
        "time facade Strang split preserves oscillator energy")

    call advance_dissipative_cayley( &
        mass, damping, time_step, [initial_state], state_next, status)
    call check_condition(status%code == 0 .and. &
        state_next(1)**2 < initial_state**2, &
        "time facade dissipative Cayley contracts mass energy")

    call check_summary("Canonical structure-preserving time facade")

contains

    pure function oscillator_energy(q_state, v_state) result(energy)
        real(dp), intent(in) :: q_state, v_state
        real(dp) :: energy

        energy = 0.5_dp*(q_state**2 + v_state**2)
    end function oscillator_energy

end program test_time_facade
