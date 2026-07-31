program test_dissipative_cayley
    use check, only: check_condition, check_summary
    use fortfem_api, only: &
        advance_dissipative_cayley, advance_dissipative_cayley_jvp, &
        advance_dissipative_cayley_vjp
    use fortfem_kinds, only: dp
    use fortsparse, only: fortsparse_status_t
    implicit none

    real(dp), parameter :: mass(2, 2) = reshape([ &
        2.0_dp, 0.2_dp, 0.2_dp, 1.5_dp], [2, 2])
    real(dp), parameter :: damping(2, 2) = reshape([ &
        0.5_dp, 0.1_dp, 0.1_dp, 0.3_dp], [2, 2])
    real(dp), parameter :: mass_dot(2, 2) = reshape([ &
        0.03_dp, -0.02_dp, -0.02_dp, 0.04_dp], [2, 2])
    real(dp), parameter :: damping_dot(2, 2) = reshape([ &
        -0.01_dp, 0.02_dp, 0.02_dp, -0.03_dp], [2, 2])
    real(dp), parameter :: state(2) = [1.0_dp, -0.4_dp]
    real(dp), parameter :: state_dot(2) = [0.2_dp, -0.1_dp]
    real(dp), parameter :: state_bar(2) = [0.7_dp, -0.5_dp]
    real(dp), parameter :: time_step = 0.6_dp
    real(dp), parameter :: time_step_dot = -0.08_dp
    real(dp), parameter :: eps = 1.0e-7_dp
    real(dp) :: state_next(2), state_next_dot(2), state_plus(2), state_minus(2)
    real(dp) :: oracle(2), lhs, rhs
    real(dp) :: mass_bar(2, 2), damping_bar(2, 2), state_input_bar(2)
    real(dp) :: time_step_bar
    real(dp) :: matrix_a(2, 2), matrix_b(2, 2), rhs_vector(2)
    real(dp) :: energy_initial, energy_final
    real(dp) :: bad_state(1)
    type(fortsparse_status_t) :: status

    matrix_a = mass + 0.5_dp*time_step*damping
    matrix_b = mass - 0.5_dp*time_step*damping
    rhs_vector = matmul(matrix_b, state)
    oracle = solve_two_by_two(matrix_a, rhs_vector)
    call advance_dissipative_cayley( &
        mass, damping, time_step, state, state_next, status)
    energy_initial = dot_product(state, matmul(mass, state))
    energy_final = dot_product(state_next, matmul(mass, state_next))
    call check_condition(status%code == 0 .and. &
        maxval(abs(state_next - oracle)) < 1.0e-14_dp .and. &
        energy_final <= energy_initial + 1.0e-13_dp, &
        "dissipative Cayley step matches oracle and decreases M-energy")

    call advance_dissipative_cayley_jvp( &
        mass, damping, time_step, state, mass_dot, damping_dot, time_step_dot, &
        state_dot, state_next_dot, status)
    call advance_dissipative_cayley( &
        mass + eps*mass_dot, damping + eps*damping_dot, &
        time_step + eps*time_step_dot, state + eps*state_dot, state_plus, status)
    call advance_dissipative_cayley( &
        mass - eps*mass_dot, damping - eps*damping_dot, &
        time_step - eps*time_step_dot, state - eps*state_dot, state_minus, status)
    call check_condition(status%code == 0 .and. &
        maxval(abs((state_plus - state_minus)/(2.0_dp*eps) - state_next_dot)) &
        < 1.0e-8_dp, "dissipative Cayley JVP matches central differences")

    call advance_dissipative_cayley_vjp( &
        mass, damping, time_step, state, state_next, state_bar, mass_bar, &
        damping_bar, time_step_bar, state_input_bar, status)
    lhs = dot_product(state_bar, state_next_dot)
    rhs = sum(mass_bar*mass_dot) + sum(damping_bar*damping_dot) + &
        time_step_bar*time_step_dot + dot_product(state_input_bar, state_dot)
    call check_condition(status%code == 0 .and. abs(lhs - rhs) < 1.0e-12_dp, &
        "dissipative Cayley VJP satisfies the real dot-product identity")

    call advance_dissipative_cayley( &
        mass, damping, time_step, state, bad_state, status)
    call check_condition(status%code /= 0, &
        "dissipative Cayley rejects an incompatible state shape")
    call advance_dissipative_cayley( &
        mass, damping, -time_step, state, state_next, status)
    call check_condition(status%code /= 0, &
        "dissipative Cayley rejects a negative time step")
    call check_summary("Dissipative Cayley")

contains

    pure function solve_two_by_two(matrix, vector) result(solution)
        real(dp), intent(in) :: matrix(2, 2), vector(2)
        real(dp) :: solution(2), determinant

        determinant = matrix(1, 1)*matrix(2, 2) - &
            matrix(1, 2)*matrix(2, 1)
        solution = [ &
            (matrix(2, 2)*vector(1) - matrix(1, 2)*vector(2))/determinant, &
            (-matrix(2, 1)*vector(1) + matrix(1, 1)*vector(2))/determinant]
    end function solve_two_by_two

end program test_dissipative_cayley
