program test_mixed_wave_midpoint
    use check, only: check_condition, check_summary
    use fortfem_api, only: advance_mixed_wave_midpoint
    use fortfem_kinds, only: dp
    use fortsparse, only: fortsparse_status_t
    implicit none

    real(dp), parameter :: mass_q(2, 2) = reshape([ &
        1.0_dp, 0.0_dp, 0.0_dp, 1.0_dp], [2, 2])
    real(dp), parameter :: mass_v(2, 2) = mass_q
    real(dp), parameter :: coupling(2, 2) = mass_q
    real(dp), parameter :: time_step = 0.2_dp
    real(dp), parameter :: initial_q(2) = [1.0_dp, 0.0_dp]
    real(dp), parameter :: initial_v(2) = [0.0_dp, 1.0_dp]
    real(dp), parameter :: cayley_cosine = 0.9801980198019802_dp
    real(dp), parameter :: cayley_sine = 0.1980198019801980_dp
    real(dp) :: q(2), v(2), energy_initial, energy_final
    type(fortsparse_status_t) :: status

    q = initial_q
    v = initial_v
    energy_initial = 0.5_dp*(dot_product(q, matmul(mass_q, q)) + &
        dot_product(v, matmul(mass_v, v)))
    call advance_mixed_wave_midpoint( &
        mass_q, mass_v, coupling, time_step, q, v, status)
    energy_final = 0.5_dp*(dot_product(q, matmul(mass_q, q)) + &
        dot_product(v, matmul(mass_v, v)))
    call check_condition(status%code == 0, &
        "mixed first-order midpoint accepts compatible wave blocks")
    call check_condition(maxval(abs(q - [cayley_cosine, -cayley_sine])) < &
        2.0e-14_dp .and. maxval(abs(v - [cayley_sine, cayley_cosine])) < &
        2.0e-14_dp, &
        "mixed midpoint matches the independent Cayley oscillator oracle")
    call check_condition(abs(energy_final - energy_initial) < 2.0e-14_dp, &
        "mixed midpoint preserves the discrete wave energy")

    call advance_mixed_wave_midpoint( &
        mass_q, mass_v, coupling, -time_step, q, v, status)
    call check_condition(maxval(abs(q - initial_q)) < 2.0e-14_dp .and. &
        maxval(abs(v - initial_v)) < 2.0e-14_dp, &
        "mixed midpoint is exactly time reversible to roundoff")

    call check_summary("Structure-preserving mixed first-order wave")
end program test_mixed_wave_midpoint
