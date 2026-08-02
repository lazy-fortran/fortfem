program test_mixed_wave_wall_midpoint
    use, intrinsic :: iso_fortran_env, only: real64
    use check, only: check_condition, check_summary
    use fortfem_time, only: &
        advance_mixed_wave_wall_midpoint, &
        advance_mixed_wave_wall_midpoint_jvp, &
        advance_mixed_wave_wall_midpoint_vjp, &
        evaluate_mixed_wave_wall_energy_balance
    use fortfem_kinds, only: dp
    use fortsparse, only: fortsparse_status_t
    implicit none

    real(dp), parameter :: h = 0.15_dp, eps = 1.0e-7_dp
    real(dp) :: mq(2, 2), mv(2, 2), inductance(1, 1), resistance(1, 1)
    real(dp) :: coupling(2, 2), port(1, 2)
    real(dp) :: mq_dot(2, 2), mv_dot(2, 2), l_dot(1, 1), r_dot(1, 1)
    real(dp) :: coupling_dot(2, 2), port_dot(1, 2), h_dot
    real(dp) :: mq_bar(2, 2), mv_bar(2, 2), l_bar(1, 1), r_bar(1, 1)
    real(dp) :: coupling_bar(2, 2), port_bar(1, 2)
    real(dp) :: q(2), v(2), current(1), q_dot(2), v_dot(2), current_dot(1)
    real(dp) :: q_next_dot(2), v_next_dot(2), current_next_dot(1)
    real(dp) :: q_plus(2), v_plus(2), current_plus(1), q_minus(2), v_minus(2), current_minus(1)
    real(dp) :: q_next_bar(2), v_next_bar(2), current_next_bar(1)
    real(dp) :: q_bar(2), v_bar(2), current_bar(1), h_bar
    real(dp) :: energy_n, energy_next, dissipation, balance
    real(dp) :: lhs, rhs
    type(fortsparse_status_t) :: status

    mq = reshape([2.0_dp, 0.1_dp, 0.1_dp, 1.5_dp], shape(mq))
    mv = reshape([1.4_dp, 0.0_dp, 0.0_dp, 1.1_dp], shape(mv))
    inductance(1, 1) = 0.9_dp
    resistance(1, 1) = 0.35_dp
    coupling = reshape([1.0_dp, -0.2_dp, 0.3_dp, 0.8_dp], shape(coupling))
    port(1, :) = [0.6_dp, -0.4_dp]
    q = [0.7_dp, -0.3_dp]
    v = [-0.2_dp, 0.5_dp]
    current(1) = 0.4_dp
    call advance_mixed_wave_wall_midpoint( &
        mq, mv, inductance, resistance, coupling, port, h, q, v, current, status)
    call check_condition(status%code == 0, "coupled wave-wall midpoint accepts compatible blocks")
    call check_condition(abs(current(1)) > 1.0e-8_dp, &
        "coupled midpoint produces a nontrivial wall current")
    call evaluate_mixed_wave_wall_energy_balance( &
        mq, mv, inductance, resistance, coupling, port, h, &
        [0.7_dp, -0.3_dp], [-0.2_dp, 0.5_dp], [0.4_dp], q, v, current, &
        energy_n, energy_next, dissipation, balance, status)
    call check_condition(status%code == 0 .and. dissipation >= 0.0_dp .and. &
        abs(balance) < 2.0e-12_dp, &
        "coupled midpoint satisfies the independent energy-dissipation ledger")

    mq_dot = reshape([0.03_dp, -0.01_dp, -0.01_dp, 0.02_dp], shape(mq_dot))
    mv_dot = reshape([0.01_dp, 0.02_dp, 0.02_dp, -0.03_dp], shape(mv_dot))
    l_dot(1, 1) = 0.04_dp
    r_dot(1, 1) = -0.02_dp
    coupling_dot = reshape([0.02_dp, -0.03_dp, 0.01_dp, 0.04_dp], shape(coupling_dot))
    port_dot(1, :) = [0.03_dp, 0.01_dp]
    h_dot = -0.06_dp
    q_dot = [0.1_dp, -0.2_dp]
    v_dot = [-0.15_dp, 0.08_dp]
    current_dot(1) = 0.05_dp
    call advance_mixed_wave_wall_midpoint_jvp( &
        mq, mv, inductance, resistance, coupling, port, h, q, v, current, &
        mq_dot, mv_dot, l_dot, r_dot, coupling_dot, port_dot, h_dot, &
        q_dot, v_dot, current_dot, q_next_dot, v_next_dot, current_next_dot, status)
    q_plus = q + eps*q_dot
    v_plus = v + eps*v_dot
    current_plus = current + eps*current_dot
    call advance_mixed_wave_wall_midpoint( &
        mq + eps*mq_dot, mv + eps*mv_dot, inductance + eps*l_dot, &
        resistance + eps*r_dot, coupling + eps*coupling_dot, port + eps*port_dot, &
        h + eps*h_dot, q_plus, v_plus, current_plus, status)
    q_minus = q - eps*q_dot
    v_minus = v - eps*v_dot
    current_minus = current - eps*current_dot
    call advance_mixed_wave_wall_midpoint( &
        mq - eps*mq_dot, mv - eps*mv_dot, inductance - eps*l_dot, &
        resistance - eps*r_dot, coupling - eps*coupling_dot, port - eps*port_dot, &
        h - eps*h_dot, q_minus, v_minus, current_minus, status)
    call check_condition(maxval(abs(q_next_dot - (q_plus - q_minus)/(2.0_dp*eps))) < 3.0e-8_dp .and. &
        maxval(abs(v_next_dot - (v_plus - v_minus)/(2.0_dp*eps))) < 3.0e-8_dp .and. &
        maxval(abs(current_next_dot - (current_plus - current_minus)/(2.0_dp*eps))) < 3.0e-8_dp, &
        "coupled wave-wall JVP matches central differences")

    q_next_bar = [0.4_dp, -0.1_dp]
    v_next_bar = [-0.3_dp, 0.2_dp]
    current_next_bar(1) = 0.5_dp
    call advance_mixed_wave_wall_midpoint_vjp( &
        mq, mv, inductance, resistance, coupling, port, h, q, v, current, &
        q_next_bar, v_next_bar, current_next_bar, mq_bar, mv_bar, l_bar, r_bar, &
        coupling_bar, port_bar, q_bar, v_bar, current_bar, h_bar, status)
    lhs = dot_product(q_next_bar, q_next_dot) + dot_product(v_next_bar, v_next_dot) + &
        dot_product(current_next_bar, current_next_dot)
    rhs = sum(mq_bar*mq_dot) + sum(mv_bar*mv_dot) + sum(l_bar*l_dot) + &
        sum(r_bar*r_dot) + sum(coupling_bar*coupling_dot) + sum(port_bar*port_dot) + &
        dot_product(q_bar, q_dot) + &
        dot_product(v_bar, v_dot) + dot_product(current_bar, current_dot) + h_bar*h_dot
    call check_condition(status%code == 0 .and. abs(lhs - rhs) < 3.0e-11_dp, &
        "coupled wave-wall VJP satisfies the real adjoint identity")
    call check_summary("Structure-preserving coupled wave-wall midpoint")
end program test_mixed_wave_wall_midpoint
