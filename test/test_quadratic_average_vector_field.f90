program test_quadratic_average_vector_field
    use check, only: check_condition, check_summary
    use fortfem_api, only: advance_quadratic_avf, &
        advance_quadratic_avf_jvp, advance_quadratic_avf_vjp
    use fortfem_kinds, only: dp
    use fortsparse, only: fortsparse_status_t
    implicit none

    real(dp), parameter :: hamiltonian(3, 3) = reshape([ &
        2.0_dp, 0.15_dp, -0.10_dp, 0.15_dp, 1.5_dp, 0.20_dp, &
        -0.10_dp, 0.20_dp, 1.2_dp], [3, 3])
    real(dp), parameter :: interconnection(3, 3) = reshape([ &
        0.0_dp, 0.7_dp, -0.2_dp, -0.7_dp, 0.0_dp, 0.4_dp, &
        0.2_dp, -0.4_dp, 0.0_dp], [3, 3])
    real(dp), parameter :: linear_term(3) = [0.3_dp, -0.4_dp, 0.2_dp]
    real(dp), parameter :: hamiltonian_dot(3, 3) = reshape([ &
        0.04_dp, -0.01_dp, 0.02_dp, -0.01_dp, 0.03_dp, 0.01_dp, &
        0.02_dp, 0.01_dp, -0.02_dp], [3, 3])
    real(dp), parameter :: interconnection_dot(3, 3) = reshape([ &
        0.0_dp, -0.03_dp, 0.02_dp, 0.03_dp, 0.0_dp, -0.01_dp, &
        -0.02_dp, 0.01_dp, 0.0_dp], [3, 3])
    real(dp), parameter :: linear_term_dot(3) = [-0.02_dp, 0.05_dp, -0.04_dp]
    real(dp), parameter :: state(3) = [0.8_dp, -0.5_dp, 0.35_dp]
    real(dp), parameter :: state_dot(3) = [0.2_dp, -0.15_dp, 0.1_dp]
    real(dp), parameter :: state_next_bar(3) = [0.6_dp, -0.2_dp, 0.9_dp]
    real(dp), parameter :: time_step = 0.35_dp
    real(dp), parameter :: time_step_dot = -0.07_dp
    real(dp), parameter :: finite_difference_step = 1.0e-6_dp
    real(dp) :: state_next(3), state_plus(3), state_minus(3)
    real(dp) :: state_next_dot(3), state_next_oracle(3)
    real(dp) :: hamiltonian_bar(3, 3), interconnection_bar(3, 3)
    real(dp) :: linear_term_bar(3), state_bar(3), time_step_bar
    real(dp) :: identity_matrix(3, 3), matrix_a(3, 3), matrix_b(3, 3)
    real(dp) :: rhs(3), initial_energy, final_energy, lhs, rhs_adjoint
    real(dp) :: state_reversed(3), state_reverse_next(3), oracle_error
    type(fortsparse_status_t) :: status

    identity_matrix = 0.0_dp
    identity_matrix(1, 1) = 1.0_dp
    identity_matrix(2, 2) = 1.0_dp
    identity_matrix(3, 3) = 1.0_dp
    matrix_a = identity_matrix - 0.5_dp*time_step*matmul( &
        interconnection, hamiltonian)
    matrix_b = identity_matrix + 0.5_dp*time_step*matmul( &
        interconnection, hamiltonian)
    rhs = matmul(matrix_b, state) + time_step*matmul(interconnection, linear_term)
    state_next_oracle = solve_three_by_three(matrix_a, rhs)

    call advance_quadratic_avf(hamiltonian, interconnection, linear_term, &
        time_step, state, state_next, status)
    oracle_error = maxval(abs(state_next - state_next_oracle))
    call check_condition(status%code == 0 .and. oracle_error < 2.0e-14_dp, &
        "quadratic AVF matches the independent linear-solve oracle")

    initial_energy = 0.5_dp*dot_product(state, matmul(hamiltonian, state)) + &
        dot_product(linear_term, state)
    final_energy = 0.5_dp*dot_product(state_next, &
        matmul(hamiltonian, state_next)) + dot_product(linear_term, state_next)
    call check_condition(abs(final_energy - initial_energy) < 2.0e-13_dp, &
        "quadratic AVF preserves the quadratic Hamiltonian")

    state_reversed = state_next
    call advance_quadratic_avf(hamiltonian, interconnection, linear_term, &
        -time_step, state_reversed, state_reverse_next, status)
    state_reversed = state_reverse_next
    call check_condition(status%code == 0 .and. &
        maxval(abs(state_reversed - state)) < 2.0e-13_dp, &
        "quadratic AVF is reversible under signed step reversal")

    call advance_quadratic_avf_jvp(hamiltonian, interconnection, linear_term, &
        time_step, state, hamiltonian_dot, interconnection_dot, &
        linear_term_dot, time_step_dot, state_dot, state_next_dot, status)
    call advance_quadratic_avf( &
        hamiltonian + finite_difference_step*hamiltonian_dot, &
        interconnection + finite_difference_step*interconnection_dot, &
        linear_term + finite_difference_step*linear_term_dot, &
        time_step + finite_difference_step*time_step_dot, &
        state + finite_difference_step*state_dot, state_plus, status)
    call advance_quadratic_avf( &
        hamiltonian - finite_difference_step*hamiltonian_dot, &
        interconnection - finite_difference_step*interconnection_dot, &
        linear_term - finite_difference_step*linear_term_dot, &
        time_step - finite_difference_step*time_step_dot, &
        state - finite_difference_step*state_dot, state_minus, status)
    call check_condition(maxval(abs((state_plus - state_minus) / &
        (2.0_dp*finite_difference_step) - state_next_dot)) < 2.0e-8_dp, &
        "quadratic AVF JVP matches an independent central difference")

    call advance_quadratic_avf_vjp(hamiltonian, interconnection, linear_term, &
        time_step, state, state_next, state_next_bar, hamiltonian_bar, &
        interconnection_bar, linear_term_bar, time_step_bar, state_bar, status)
    lhs = dot_product(state_next_bar, state_next_dot)
    rhs_adjoint = sum(hamiltonian_bar*hamiltonian_dot) + &
        sum(interconnection_bar*interconnection_dot) + &
        dot_product(linear_term_bar, linear_term_dot) + &
        time_step_bar*time_step_dot + dot_product(state_bar, state_dot)
    call check_condition(abs(lhs - rhs_adjoint) < 2.0e-12_dp, &
        "quadratic AVF VJP satisfies the real adjoint identity")

    call check_summary("Quadratic average-vector-field step")

contains

    pure function solve_three_by_three(matrix, vector) result(solution)
        real(dp), intent(in) :: matrix(3, 3), vector(3)
        real(dp) :: solution(3), work(3, 4), factor
        integer :: pivot, row

        work(:, :3) = matrix
        work(:, 4) = vector
        do pivot = 1, 2
            do row = pivot + 1, 3
                factor = work(row, pivot)/work(pivot, pivot)
                work(row, pivot:4) = work(row, pivot:4) - &
                    factor*work(pivot, pivot:4)
            end do
        end do
        solution(3) = work(3, 4)/work(3, 3)
        solution(2) = (work(2, 4) - work(2, 3)*solution(3))/work(2, 2)
        solution(1) = (work(1, 4) - work(1, 2)*solution(2) - &
            work(1, 3)*solution(3))/work(1, 1)
    end function solve_three_by_three

end program test_quadratic_average_vector_field
