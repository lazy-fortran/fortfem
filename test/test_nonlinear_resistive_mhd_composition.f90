program test_nonlinear_resistive_mhd_composition
    use check, only: check_condition, check_summary
    use fortfem_api, only: &
        RESISTIVE_MHD_AMPERE, &
        RESISTIVE_MHD_ANISOTROPIC_TRANSPORT, &
        RESISTIVE_MHD_BLOCK_COUNT, &
        RESISTIVE_MHD_FARADAY, &
        RESISTIVE_MHD_FREE_BOUNDARY, &
        RESISTIVE_MHD_MOMENTUM, &
        RESISTIVE_MHD_PRESSURE, &
        RESISTIVE_MHD_TENSOR, &
        RESISTIVE_MHD_WALL, &
        assemble_nonlinear_resistive_mhd_residual, &
        assemble_nonlinear_resistive_mhd_residual_jvp, &
        assemble_nonlinear_resistive_mhd_residual_vjp, &
        nonlinear_resistive_mhd_energy_ledger_t
    use fortfem_kinds, only: dp
    use fortsparse, only: FORTSPARSE_OK, fortsparse_status_t
    implicit none

    integer, parameter :: state_size = 3
    real(dp), parameter :: epsilon = 1.0e-6_dp
    real(dp) :: state(state_size), state_dot(state_size)
    real(dp) :: residual(state_size), residual_dot(state_size)
    real(dp) :: residual_plus(state_size), residual_minus(state_size)
    real(dp) :: state_bar(state_size), residual_bar(state_size)
    real(dp) :: expected(state_size), expected_dot(state_size)
    real(dp) :: lhs, rhs
    type(nonlinear_resistive_mhd_energy_ledger_t) :: ledger, ledger_dot
    type(nonlinear_resistive_mhd_energy_ledger_t) :: ledger_plus, ledger_minus
    type(nonlinear_resistive_mhd_energy_ledger_t) :: ledger_bar
    type(fortsparse_status_t) :: status
    logical :: callback_seen(RESISTIVE_MHD_BLOCK_COUNT)
    integer :: block, component

    state = [0.7_dp, -0.35_dp, 0.9_dp]
    state_dot = [-0.2_dp, 0.4_dp, 0.1_dp]
    residual_bar = [0.3_dp, -0.6_dp, 0.2_dp]
    ledger_bar%stored_energy = 0.4_dp
    ledger_bar%input_power = -0.7_dp
    ledger_bar%dissipation = 0.5_dp
    ledger_bar%balance = 0.9_dp
    callback_seen = .false.

    call assemble_nonlinear_resistive_mhd_residual( &
        state, nonlinear_value, residual, ledger, status)
    call independent_value_oracle(state, expected, ledger)
    call check_condition(status%code == FORTSPARSE_OK, &
        "nonlinear resistive-MHD block composition accepts caller callbacks")
    call check_condition(all(callback_seen), &
        "all Faraday/Ampere, force, transport, wall, and boundary blocks are visited")
    call check_condition(maxval(abs(residual - expected)) < 1.0e-14_dp, &
        "nonlinear block residual matches independent oracle")

    call assemble_nonlinear_resistive_mhd_residual_jvp( &
        state, state_dot, nonlinear_value, nonlinear_jvp, residual_dot, ledger_dot, status)
    call assemble_nonlinear_resistive_mhd_residual( &
        state + epsilon*state_dot, nonlinear_value, residual_plus, ledger_plus, status)
    call assemble_nonlinear_resistive_mhd_residual( &
        state - epsilon*state_dot, nonlinear_value, residual_minus, ledger_minus, status)
    call independent_jvp_oracle(state, state_dot, expected_dot, ledger_dot)
    call check_condition(status%code == FORTSPARSE_OK .and. &
        maxval(abs(residual_dot - expected_dot)) < 1.0e-13_dp .and. &
        maxval(abs(residual_dot - (residual_plus - residual_minus) / &
            (2.0_dp*epsilon))) < 2.0e-8_dp .and. &
        abs(ledger_dot%stored_energy - (ledger_plus%stored_energy - &
            ledger_minus%stored_energy)/(2.0_dp*epsilon)) < 2.0e-8_dp .and. &
        abs(ledger_dot%input_power - (ledger_plus%input_power - &
            ledger_minus%input_power)/(2.0_dp*epsilon)) < 2.0e-8_dp .and. &
        abs(ledger_dot%dissipation - (ledger_plus%dissipation - &
            ledger_minus%dissipation)/(2.0_dp*epsilon)) < 2.0e-8_dp .and. &
        abs(ledger_dot%balance - (ledger_plus%balance - &
            ledger_minus%balance)/(2.0_dp*epsilon)) < 2.0e-8_dp, &
        "nonlinear block JVP matches independent and central-difference oracles")

    call assemble_nonlinear_resistive_mhd_residual_vjp( &
        state, nonlinear_value, nonlinear_vjp, residual_bar, ledger_bar, state_bar, status)
    lhs = dot_product(residual_bar, residual_dot) + &
        ledger_bar%stored_energy*ledger_dot%stored_energy + &
        ledger_bar%input_power*ledger_dot%input_power + &
        ledger_bar%dissipation*ledger_dot%dissipation + &
        ledger_bar%balance*ledger_dot%balance
    rhs = dot_product(state_bar, state_dot)
    call check_condition(status%code == FORTSPARSE_OK .and. abs(lhs - rhs) < 1.0e-12_dp, &
        "nonlinear block VJP satisfies the full energy-ledger dot-product identity")

    call assemble_nonlinear_resistive_mhd_residual( &
        state, rejecting_value, residual, ledger, status)
    call check_condition(status%code /= FORTSPARSE_OK .and. maxval(abs(residual)) == 0.0_dp, &
        "callback failure is rejected without a partial residual")

    call check_summary("nonlinear resistive-MHD composition")

contains

    subroutine nonlinear_value(state_in, block_id, block_residual, stored_energy, &
            input_power, dissipation, local_status)
        real(dp), intent(in) :: state_in(:)
        integer, intent(in) :: block_id
        real(dp), intent(out) :: block_residual(:)
        real(dp), intent(out) :: stored_energy, input_power, dissipation
        integer, intent(out) :: local_status
        real(dp) :: coefficient

        block_residual = 0.0_dp
        stored_energy = 0.0_dp
        input_power = 0.0_dp
        dissipation = 0.0_dp
        local_status = 1
        if (block_id < 1 .or. block_id > RESISTIVE_MHD_BLOCK_COUNT) return
        if (size(state_in) /= state_size .or. size(block_residual) /= state_size) return
        coefficient = real(block_id, dp)
        block_residual = coefficient*state_in**2 + &
            0.1_dp*coefficient*state_in
        stored_energy = 0.5_dp*coefficient*sum(state_in**2)
        input_power = coefficient*sum(state_in)
        dissipation = 0.25_dp*coefficient*sum(state_in**2)
        callback_seen(block_id) = .true.
        local_status = 0
    end subroutine nonlinear_value

    subroutine nonlinear_jvp(state_in, state_direction, block_id, block_residual_dot, &
            stored_energy_dot, input_power_dot, dissipation_dot, local_status)
        real(dp), intent(in) :: state_in(:), state_direction(:)
        integer, intent(in) :: block_id
        real(dp), intent(out) :: block_residual_dot(:)
        real(dp), intent(out) :: stored_energy_dot, input_power_dot, dissipation_dot
        integer, intent(out) :: local_status
        real(dp) :: coefficient

        block_residual_dot = 0.0_dp
        stored_energy_dot = 0.0_dp
        input_power_dot = 0.0_dp
        dissipation_dot = 0.0_dp
        local_status = 1
        if (block_id < 1 .or. block_id > RESISTIVE_MHD_BLOCK_COUNT) return
        if (size(state_in) /= state_size .or. size(state_direction) /= state_size .or. &
            size(block_residual_dot) /= state_size) return
        coefficient = real(block_id, dp)
        block_residual_dot = (2.0_dp*coefficient*state_in + &
            0.1_dp*coefficient)*state_direction
        stored_energy_dot = coefficient*dot_product(state_in, state_direction)
        input_power_dot = coefficient*sum(state_direction)
        dissipation_dot = 0.5_dp*coefficient*dot_product(state_in, state_direction)
        local_status = 0
    end subroutine nonlinear_jvp

    subroutine nonlinear_vjp(state_in, block_id, block_residual_bar, stored_energy_bar, &
            input_power_bar, dissipation_bar, balance_bar, state_bar_out, local_status)
        real(dp), intent(in) :: state_in(:), block_residual_bar(:)
        integer, intent(in) :: block_id
        real(dp), intent(in) :: stored_energy_bar, input_power_bar, dissipation_bar
        real(dp), intent(in) :: balance_bar
        real(dp), intent(out) :: state_bar_out(:)
        integer, intent(out) :: local_status
        real(dp) :: coefficient

        state_bar_out = 0.0_dp
        local_status = 1
        if (block_id < 1 .or. block_id > RESISTIVE_MHD_BLOCK_COUNT) return
        if (size(state_in) /= state_size .or. size(block_residual_bar) /= state_size .or. &
            size(state_bar_out) /= state_size) return
        coefficient = real(block_id, dp)
        state_bar_out = block_residual_bar*(2.0_dp*coefficient*state_in + &
            0.1_dp*coefficient) + (stored_energy_bar + 0.5_dp*dissipation_bar)* &
            coefficient*state_in + (input_power_bar + balance_bar)*coefficient - &
            balance_bar*0.5_dp*coefficient*state_in
        local_status = 0
    end subroutine nonlinear_vjp

    subroutine rejecting_value(state_in, block_id, block_residual, stored_energy, &
            input_power, dissipation, local_status)
        real(dp), intent(in) :: state_in(:)
        integer, intent(in) :: block_id
        real(dp), intent(out) :: block_residual(:)
        real(dp), intent(out) :: stored_energy, input_power, dissipation
        integer, intent(out) :: local_status

        block_residual = 0.0_dp
        stored_energy = 0.0_dp
        input_power = 0.0_dp
        dissipation = 0.0_dp
        local_status = merge(7, 0, block_id == RESISTIVE_MHD_AMPERE)
    end subroutine rejecting_value

    subroutine independent_value_oracle(state_in, residual_out, ledger_out)
        real(dp), intent(in) :: state_in(:)
        real(dp), intent(out) :: residual_out(:)
        type(nonlinear_resistive_mhd_energy_ledger_t), intent(out) :: ledger_out
        real(dp) :: coefficient

        residual_out = 0.0_dp
        ledger_out = nonlinear_resistive_mhd_energy_ledger_t()
        do block = 1, RESISTIVE_MHD_BLOCK_COUNT
            coefficient = real(block, dp)
            residual_out = residual_out + coefficient*state_in**2 + &
                0.1_dp*coefficient*state_in
            ledger_out%stored_energy = ledger_out%stored_energy + &
                0.5_dp*coefficient*sum(state_in**2)
            ledger_out%input_power = ledger_out%input_power + coefficient*sum(state_in)
            ledger_out%dissipation = ledger_out%dissipation + &
                0.25_dp*coefficient*sum(state_in**2)
        end do
        ledger_out%balance = ledger_out%input_power - ledger_out%dissipation
    end subroutine independent_value_oracle

    subroutine independent_jvp_oracle(state_in, direction, residual_out, ledger_out)
        real(dp), intent(in) :: state_in(:), direction(:)
        real(dp), intent(out) :: residual_out(:)
        type(nonlinear_resistive_mhd_energy_ledger_t), intent(out) :: ledger_out
        real(dp) :: coefficient

        residual_out = 0.0_dp
        ledger_out = nonlinear_resistive_mhd_energy_ledger_t()
        do block = 1, RESISTIVE_MHD_BLOCK_COUNT
            coefficient = real(block, dp)
            residual_out = residual_out + (2.0_dp*coefficient*state_in + &
                0.1_dp*coefficient)*direction
            ledger_out%stored_energy = ledger_out%stored_energy + &
                coefficient*dot_product(state_in, direction)
            ledger_out%input_power = ledger_out%input_power + coefficient*sum(direction)
            ledger_out%dissipation = ledger_out%dissipation + &
                0.5_dp*coefficient*dot_product(state_in, direction)
        end do
        ledger_out%balance = ledger_out%input_power - ledger_out%dissipation
    end subroutine independent_jvp_oracle

end program test_nonlinear_resistive_mhd_composition
