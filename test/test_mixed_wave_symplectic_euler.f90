program test_mixed_wave_symplectic_euler
    use check, only: check_condition, check_summary
    use fortfem_api, only: advance_mixed_wave_symplectic_euler, &
        advance_mixed_wave_symplectic_euler_jvp, &
        advance_mixed_wave_symplectic_euler_vjp
    use fortfem_kinds, only: dp
    use fortsparse, only: fortsparse_status_t
    implicit none

    real(dp), parameter :: mass_q(2, 2) = reshape([ &
        1.0_dp, 0.0_dp, 0.0_dp, 1.0_dp], [2, 2])
    real(dp), parameter :: mass_v(2, 2) = mass_q
    real(dp), parameter :: coupling(2, 2) = mass_q
    real(dp), parameter :: mass_q_deriv(2, 2) = reshape([ &
        2.0_dp, 0.0_dp, 0.0_dp, 3.0_dp], [2, 2])
    real(dp), parameter :: mass_v_deriv(2, 2) = reshape([ &
        4.0_dp, 0.0_dp, 0.0_dp, 5.0_dp], [2, 2])
    real(dp), parameter :: coupling_deriv(2, 2) = reshape([ &
        1.0_dp, -0.5_dp, 2.0_dp, 1.5_dp], [2, 2])
    real(dp), parameter :: time_step = 0.2_dp
    real(dp), parameter :: initial_q(2) = [1.0_dp, 0.0_dp]
    real(dp), parameter :: initial_v(2) = [0.0_dp, 1.0_dp]
    real(dp), parameter :: expected_q(2) = [0.96_dp, -0.2_dp]
    real(dp), parameter :: expected_v(2) = [0.2_dp, 1.0_dp]
    real(dp) :: q(2), v(2), q_other(2), v_other(2), q_bad(1)
    real(dp) :: q_dot(2), v_dot(2), q_next_dot(2), v_next_dot(2)
    real(dp) :: q_plus(2), v_plus(2), q_minus(2), v_minus(2)
    real(dp) :: q_bar(2), v_bar(2), time_step_bar
    real(dp) :: q_next_bar(2), v_next_bar(2), vjp_left, vjp_right
    real(dp), parameter :: time_step_dot = -0.07_dp
    real(dp), parameter :: finite_difference_step = 1.0e-7_dp
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

    q = initial_q
    v = initial_v
    q_dot = [0.25_dp, -0.4_dp]
    v_dot = [-0.3_dp, 0.15_dp]
    call advance_mixed_wave_symplectic_euler_jvp( &
        mass_q_deriv, mass_v_deriv, coupling_deriv, time_step, q, v, &
        time_step_dot, q_dot, v_dot, q_next_dot, v_next_dot, status)
    q_plus = initial_q + finite_difference_step*q_dot
    v_plus = initial_v + finite_difference_step*v_dot
    call advance_mixed_wave_symplectic_euler( &
        mass_q_deriv, mass_v_deriv, coupling_deriv, &
        time_step + finite_difference_step*time_step_dot, q_plus, v_plus, status)
    q_minus = initial_q - finite_difference_step*q_dot
    v_minus = initial_v - finite_difference_step*v_dot
    call advance_mixed_wave_symplectic_euler( &
        mass_q_deriv, mass_v_deriv, coupling_deriv, &
        time_step - finite_difference_step*time_step_dot, q_minus, v_minus, status)
    call check_condition(maxval(abs((q_plus - q_minus)/(2.0_dp* &
        finite_difference_step) - q_next_dot)) < 2.0e-8_dp .and. &
        maxval(abs((v_plus - v_minus)/(2.0_dp*finite_difference_step) - &
        v_next_dot)) < 2.0e-8_dp, &
        "mixed symplectic Euler JVP matches an independent central difference")

    q_next_bar = [0.6_dp, -0.2_dp]
    v_next_bar = [-0.5_dp, 0.8_dp]
    call advance_mixed_wave_symplectic_euler_vjp( &
        mass_q_deriv, mass_v_deriv, coupling_deriv, time_step, initial_q, initial_v, &
        q_next_bar, v_next_bar, q_bar, v_bar, time_step_bar, status)
    vjp_left = dot_product(q_next_bar, q_next_dot) + &
        dot_product(v_next_bar, v_next_dot)
    vjp_right = dot_product(q_bar, q_dot) + dot_product(v_bar, v_dot) + &
        time_step_bar*time_step_dot
    call check_condition(abs(vjp_left - vjp_right) < 2.0e-13_dp, &
        "mixed symplectic Euler VJP satisfies the real adjoint identity")
    call check_summary("Structure-preserving mixed symplectic Euler")
end program test_mixed_wave_symplectic_euler
