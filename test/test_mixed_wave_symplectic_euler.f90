program test_mixed_wave_symplectic_euler
    use check, only: check_condition, check_summary
    use fortfem_api, only: advance_mixed_wave_symplectic_euler
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
    real(dp), parameter :: expected_q(2) = [0.96_dp, -0.2_dp]
    real(dp), parameter :: expected_v(2) = [0.2_dp, 1.0_dp]
    real(dp) :: q(2), v(2), q_other(2), v_other(2), q_bad(1)
    real(dp) :: symplectic_initial, symplectic_final
    type(fortsparse_status_t) :: status

    q = initial_q
    v = initial_v
    call advance_mixed_wave_symplectic_euler( &
        mass_q, mass_v, coupling, time_step, q, v, status)
    call check_condition(status%code == 0, &
        "mixed symplectic Euler accepts compatible wave blocks")
    call check_condition(maxval(abs(q - expected_q)) < 1.0e-14_dp .and. &
        maxval(abs(v - expected_v)) < 1.0e-14_dp, &
        "mixed symplectic Euler matches the independent partitioned oracle")

    q_other = [0.3_dp, -0.7_dp]
    v_other = [-0.4_dp, 0.9_dp]
    symplectic_initial = dot_product(initial_q, v_other) - &
        dot_product(initial_v, q_other)
    call advance_mixed_wave_symplectic_euler( &
        mass_q, mass_v, coupling, time_step, q_other, v_other, status)
    symplectic_final = dot_product(q, v_other) - dot_product(v, q_other)
    call check_condition(abs(symplectic_final - symplectic_initial) < &
        1.0e-14_dp, "mixed symplectic Euler preserves the canonical two-state form")

    q_bad = 0.0_dp
    call advance_mixed_wave_symplectic_euler( &
        mass_q, mass_v, coupling, time_step, q_bad, v, status)
    call check_condition(status%code /= 0, &
        "mixed symplectic Euler rejects incompatible state dimensions")
    call check_summary("Structure-preserving mixed symplectic Euler")
end program test_mixed_wave_symplectic_euler
