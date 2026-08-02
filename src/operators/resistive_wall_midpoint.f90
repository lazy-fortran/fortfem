module fortfem_resistive_wall_midpoint
    !! Structure-preserving implicit-midpoint update for a neutral wall RL block.
    !!
    !! The caller supplies inductance, resistance, and voltage data. The
    !! update is an operator-level building block for conducting-wall and
    !! response-matrix couplings; it contains no wall geometry or application
    !! normalization.
    use fortnum_linalg, only: dense_solve
    use fortfem_kinds, only: dp
    implicit none
    private

    public :: advance_resistive_wall_midpoint
    public :: advance_resistive_wall_midpoint_jvp
    public :: advance_resistive_wall_midpoint_vjp
    public :: evaluate_resistive_wall_energy_balance

contains

    subroutine advance_resistive_wall_midpoint( &
            inductance, resistance, step_size, current, voltage_n, voltage_next, &
            current_next, status)
        real(dp), intent(in) :: inductance(:, :), resistance(:, :), step_size
        real(dp), intent(in) :: current(:), voltage_n(:), voltage_next(:)
        real(dp), intent(out) :: current_next(:)
        integer, intent(out) :: status

        real(dp), allocatable :: matrix_a(:, :), matrix_b(:, :), right_hand_side(:)
        real(dp) :: half_step
        integer :: info

        current_next = 0.0_dp
        status = 1
        if (.not. valid_inputs( &
            inductance, resistance, step_size, current, voltage_n, voltage_next, &
            current_next)) return
        half_step = 0.5_dp*step_size
        allocate( &
            matrix_a(size(current), size(current)), &
            matrix_b(size(current), size(current)), right_hand_side(size(current)))
        matrix_a = inductance + half_step*resistance
        matrix_b = inductance - half_step*resistance
        right_hand_side = matmul(matrix_b, current) + &
            step_size*0.5_dp*(voltage_n + voltage_next)
        call dense_solve(matrix_a, right_hand_side, current_next, info)
        if (info /= 0) then
            status = 2
            return
        end if
        status = 0
    end subroutine advance_resistive_wall_midpoint

    subroutine advance_resistive_wall_midpoint_jvp( &
            inductance, resistance, step_size, current, voltage_n, voltage_next, &
            inductance_dot, resistance_dot, step_size_dot, current_dot, &
            voltage_n_dot, voltage_next_dot, current_next_dot, status)
        real(dp), intent(in) :: inductance(:, :), resistance(:, :), step_size
        real(dp), intent(in) :: current(:), voltage_n(:), voltage_next(:)
        real(dp), intent(in) :: inductance_dot(:, :), resistance_dot(:, :)
        real(dp), intent(in) :: step_size_dot, current_dot(:)
        real(dp), intent(in) :: voltage_n_dot(:), voltage_next_dot(:)
        real(dp), intent(out) :: current_next_dot(:)
        integer, intent(out) :: status

        real(dp), allocatable :: matrix_a(:, :), matrix_b(:, :)
        real(dp), allocatable :: matrix_a_dot(:, :), matrix_b_dot(:, :)
        real(dp), allocatable :: right_hand_side(:), right_hand_side_dot(:)
        real(dp), allocatable :: current_next(:)
        real(dp) :: half_step, half_step_dot
        integer :: info

        current_next_dot = 0.0_dp
        status = 1
        if (.not. valid_inputs( &
            inductance, resistance, step_size, current, voltage_n, voltage_next, &
            current_next_dot) .or. &
            any(shape(inductance_dot) /= shape(inductance)) .or. &
            any(shape(resistance_dot) /= shape(resistance)) .or. &
            size(current_dot) /= size(current) .or. &
            size(voltage_n_dot) /= size(voltage_n) .or. &
            size(voltage_next_dot) /= size(voltage_next)) return
        allocate(current_next(size(current)))
        call advance_resistive_wall_midpoint( &
            inductance, resistance, step_size, current, voltage_n, voltage_next, &
            current_next, status)
        if (status /= 0) return
        half_step = 0.5_dp*step_size
        half_step_dot = 0.5_dp*step_size_dot
        allocate( &
            matrix_a(size(current), size(current)), &
            matrix_b(size(current), size(current)), &
            matrix_a_dot(size(current), size(current)), &
            matrix_b_dot(size(current), size(current)), &
            right_hand_side(size(current)), right_hand_side_dot(size(current)))
        matrix_a = inductance + half_step*resistance
        matrix_b = inductance - half_step*resistance
        matrix_a_dot = inductance_dot + half_step_dot*resistance + &
            half_step*resistance_dot
        matrix_b_dot = inductance_dot - half_step_dot*resistance - &
            half_step*resistance_dot
        right_hand_side = matmul(matrix_b, current) + &
            step_size*0.5_dp*(voltage_n + voltage_next)
        right_hand_side_dot = matmul(matrix_b_dot, current) + &
            matmul(matrix_b, current_dot) + &
            step_size_dot*0.5_dp*(voltage_n + voltage_next) + &
            step_size*0.5_dp*(voltage_n_dot + voltage_next_dot)
        call dense_solve(matrix_a, right_hand_side_dot - &
            matmul(matrix_a_dot, current_next), current_next_dot, info)
        if (info /= 0) then
            status = 2
            return
        end if
        status = 0
    end subroutine advance_resistive_wall_midpoint_jvp

    subroutine advance_resistive_wall_midpoint_vjp( &
            inductance, resistance, step_size, current, voltage_n, voltage_next, &
            current_next_bar, inductance_bar, resistance_bar, step_size_bar, &
            current_bar, voltage_n_bar, voltage_next_bar, status)
        real(dp), intent(in) :: inductance(:, :), resistance(:, :), step_size
        real(dp), intent(in) :: current(:), voltage_n(:), voltage_next(:)
        real(dp), intent(in) :: current_next_bar(:)
        real(dp), intent(out) :: inductance_bar(:, :), resistance_bar(:, :)
        real(dp), intent(out) :: step_size_bar, current_bar(:), voltage_n_bar(:)
        real(dp), intent(out) :: voltage_next_bar(:)
        integer, intent(out) :: status

        real(dp), allocatable :: matrix_a(:, :), matrix_b(:, :), current_next(:)
        real(dp), allocatable :: right_hand_side(:), adjoint(:)
        real(dp), allocatable :: matrix_a_bar(:, :), matrix_b_bar(:, :)
        real(dp) :: half_step
        integer :: info

        inductance_bar = 0.0_dp
        resistance_bar = 0.0_dp
        step_size_bar = 0.0_dp
        current_bar = 0.0_dp
        voltage_n_bar = 0.0_dp
        voltage_next_bar = 0.0_dp
        status = 1
        if (.not. valid_inputs( &
            inductance, resistance, step_size, current, voltage_n, voltage_next, &
            current_next_bar) .or. size(current_next_bar) /= size(current) .or. &
            any(shape(inductance_bar) /= shape(inductance)) .or. &
            any(shape(resistance_bar) /= shape(resistance)) .or. &
            size(current_bar) /= size(current) .or. &
            size(voltage_n_bar) /= size(voltage_n) .or. &
            size(voltage_next_bar) /= size(voltage_next)) return
        allocate( &
            matrix_a(size(current), size(current)), &
            matrix_b(size(current), size(current)), current_next(size(current)), &
            right_hand_side(size(current)), adjoint(size(current)), &
            matrix_a_bar(size(current), size(current)), &
            matrix_b_bar(size(current), size(current)))
        half_step = 0.5_dp*step_size
        matrix_a = inductance + half_step*resistance
        matrix_b = inductance - half_step*resistance
        right_hand_side = matmul(matrix_b, current) + &
            step_size*0.5_dp*(voltage_n + voltage_next)
        call dense_solve(matrix_a, right_hand_side, current_next, info)
        if (info /= 0) then
            status = 2
            return
        end if
        call dense_solve(transpose(matrix_a), current_next_bar, adjoint, info)
        if (info /= 0) then
            status = 2
            return
        end if
        matrix_a_bar = -outer_product(adjoint, current_next)
        matrix_b_bar = outer_product(adjoint, current)
        inductance_bar = matrix_a_bar + matrix_b_bar
        resistance_bar = half_step*(matrix_a_bar - matrix_b_bar)
        step_size_bar = 0.5_dp*sum(resistance*(matrix_a_bar - matrix_b_bar)) + &
            0.5_dp*dot_product(adjoint, voltage_n + voltage_next)
        current_bar = matmul(transpose(matrix_b), adjoint)
        voltage_n_bar = 0.5_dp*step_size*adjoint
        voltage_next_bar = voltage_n_bar
        status = 0
    end subroutine advance_resistive_wall_midpoint_vjp

    subroutine evaluate_resistive_wall_energy_balance( &
            inductance, resistance, step_size, current, current_next, voltage_n, &
            voltage_next, energy_n, energy_next, input_work, dissipation, balance, &
            status)
        real(dp), intent(in) :: inductance(:, :), resistance(:, :), step_size
        real(dp), intent(in) :: current(:), current_next(:)
        real(dp), intent(in) :: voltage_n(:), voltage_next(:)
        real(dp), intent(out) :: energy_n, energy_next, input_work, dissipation
        real(dp), intent(out) :: balance
        integer, intent(out) :: status

        real(dp) :: current_mid(size(current)), voltage_mid(size(current))

        energy_n = 0.0_dp
        energy_next = 0.0_dp
        input_work = 0.0_dp
        dissipation = 0.0_dp
        balance = huge(1.0_dp)
        status = 1
        if (.not. valid_inputs( &
            inductance, resistance, step_size, current, voltage_n, voltage_next, &
            current_next) .or. size(current_next) /= size(current)) return
        current_mid = 0.5_dp*(current + current_next)
        voltage_mid = 0.5_dp*(voltage_n + voltage_next)
        energy_n = 0.5_dp*dot_product(current, matmul(inductance, current))
        energy_next = 0.5_dp*dot_product( &
            current_next, matmul(inductance, current_next))
        input_work = step_size*dot_product(voltage_mid, current_mid)
        dissipation = step_size*dot_product( &
            current_mid, matmul(resistance, current_mid))
        balance = energy_next - energy_n - input_work + dissipation
        status = 0
    end subroutine evaluate_resistive_wall_energy_balance

    logical function valid_inputs( &
            inductance, resistance, step_size, current, voltage_n, voltage_next, &
            output) result(valid)
        real(dp), intent(in) :: inductance(:, :), resistance(:, :), step_size
        real(dp), intent(in) :: current(:), voltage_n(:), voltage_next(:), output(:)

        valid = .false.
        if (step_size <= 0.0_dp) return
        if (size(inductance, 1) < 1 .or. &
            size(inductance, 1) /= size(inductance, 2) .or. &
            any(shape(resistance) /= shape(inductance))) return
        if (size(current) /= size(inductance, 1) .or. &
            size(voltage_n) /= size(current) .or. &
            size(voltage_next) /= size(current) .or. &
            size(output) /= size(current)) return
        valid = .true.
    end function valid_inputs

    pure function outer_product(first, second) result(product)
        real(dp), intent(in) :: first(:), second(:)
        real(dp) :: product(size(first), size(second))

        product = spread(first, 2, size(second))* &
            spread(second, 1, size(first))
    end function outer_product

end module fortfem_resistive_wall_midpoint
