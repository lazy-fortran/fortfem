program test_nonlinear_resistive_mhd_properties
    !! Seeded property campaign for the closure-neutral nonlinear MHD composer.
    !!
    !! The callbacks below own small caller-side block residuals and their
    !! derivatives.  Independent summed-loop oracles keep this test separate
    !! from the implementation's block accumulation and ledger bookkeeping.
    use, intrinsic :: iso_fortran_env, only: int32
    use check, only: check_condition, check_property, check_summary, &
        property_random_integer, property_random_unit, property_rng_initialize, &
        property_rng_t
    use fortfem_interop, only: &
        RESISTIVE_MHD_BLOCK_COUNT, &
        assemble_nonlinear_resistive_mhd_residual, &
        assemble_nonlinear_resistive_mhd_residual_jvp, &
        assemble_nonlinear_resistive_mhd_residual_vjp, &
        nonlinear_resistive_mhd_energy_ledger_t
    use fortfem_kinds, only: dp
    use fortsparse, only: FORTSPARSE_OK, fortsparse_status_t
    implicit none

    logical :: all_passed
    integer(int32) :: first_failed_seed, shrunk_seed
    integer :: active_component
    real(dp), allocatable :: active_coefficients(:)

    call check_property( &
        "random nonlinear resistive-MHD composition preserves ledger and adjoint", &
        20260802_int32, 20, mhd_case, all_passed, first_failed_seed, shrunk_seed)
    call check_condition(all_passed .and. first_failed_seed == 0_int32 .and. &
        shrunk_seed == 0_int32, &
        "random nonlinear resistive-MHD property reports no failure seed")
    call check_summary("random nonlinear resistive-MHD composition properties")
    if (.not. all_passed) error stop 1

contains

    logical function mhd_case(case_seed)
        integer(int32), intent(in) :: case_seed
        type(property_rng_t) :: rng
        type(fortsparse_status_t) :: status
        type(nonlinear_resistive_mhd_energy_ledger_t) :: ledger
        type(nonlinear_resistive_mhd_energy_ledger_t) :: ledger_expected
        type(nonlinear_resistive_mhd_energy_ledger_t) :: ledger_dot
        type(nonlinear_resistive_mhd_energy_ledger_t) :: ledger_plus
        type(nonlinear_resistive_mhd_energy_ledger_t) :: ledger_minus
        type(nonlinear_resistive_mhd_energy_ledger_t) :: ledger_bar
        real(dp), allocatable :: state(:), state_dot(:)
        real(dp), allocatable :: residual(:), residual_expected(:)
        real(dp), allocatable :: residual_dot(:), residual_dot_expected(:)
        real(dp), allocatable :: residual_plus(:), residual_minus(:)
        real(dp), allocatable :: residual_bar(:), state_bar(:)
        real(dp) :: epsilon, lhs, rhs
        integer :: block

        mhd_case = .false.
        call property_rng_initialize(rng, case_seed)
        active_component = property_random_integer(rng, 1, 4)
        if (allocated(active_coefficients)) deallocate(active_coefficients)
        allocate(active_coefficients(RESISTIVE_MHD_BLOCK_COUNT))
        allocate(state(active_component), state_dot(active_component), &
            residual(active_component), residual_expected(active_component), &
            residual_dot(active_component), residual_dot_expected(active_component), &
            residual_plus(active_component), residual_minus(active_component), &
            residual_bar(active_component), state_bar(active_component))
        do block = 1, active_component
            state(block) = 1.6_dp*property_random_unit(rng) - 0.8_dp
            state_dot(block) = 1.2_dp*property_random_unit(rng) - 0.6_dp
            residual_bar(block) = 1.4_dp*property_random_unit(rng) - 0.7_dp
        end do
        do block = 1, RESISTIVE_MHD_BLOCK_COUNT
            active_coefficients(block) = &
                0.2_dp + 1.8_dp*property_random_unit(rng)
        end do
        ledger_bar%stored_energy = 1.4_dp*property_random_unit(rng) - 0.7_dp
        ledger_bar%input_power = 1.4_dp*property_random_unit(rng) - 0.7_dp
        ledger_bar%dissipation = 1.4_dp*property_random_unit(rng) - 0.7_dp
        ledger_bar%balance = 1.4_dp*property_random_unit(rng) - 0.7_dp
        epsilon = 1.0e-6_dp

        call assemble_nonlinear_resistive_mhd_residual( &
            state, value_law, residual, ledger, status)
        call independent_value_oracle( &
            state, active_coefficients, residual_expected, ledger_expected)
        if (status%code /= FORTSPARSE_OK .or. &
            maxval(abs(residual - residual_expected)) > 2.0e-13_dp .or. &
            .not. ledger_close(ledger, ledger_expected, 2.0e-13_dp) .or. &
            ledger%dissipation < 0.0_dp .or. &
            abs(ledger%balance - (ledger%input_power - ledger%dissipation)) > &
                2.0e-13_dp) return

        call assemble_nonlinear_resistive_mhd_residual_jvp( &
            state, state_dot, value_law, jvp_law, residual_dot, ledger_dot, status)
        call independent_jvp_oracle( &
            state, state_dot, active_coefficients, residual_dot_expected, ledger_expected)
        if (status%code /= FORTSPARSE_OK .or. &
            maxval(abs(residual_dot - residual_dot_expected)) > 3.0e-13_dp .or. &
            .not. ledger_close(ledger_dot, ledger_expected, 3.0e-13_dp) .or. &
            abs(ledger_dot%balance - (ledger_dot%input_power - &
                ledger_dot%dissipation)) > 3.0e-13_dp) return

        call assemble_nonlinear_resistive_mhd_residual( &
            state + epsilon*state_dot, value_law, residual_plus, ledger_plus, status)
        if (status%code /= FORTSPARSE_OK) return
        call assemble_nonlinear_resistive_mhd_residual( &
            state - epsilon*state_dot, value_law, residual_minus, ledger_minus, status)
        if (status%code /= FORTSPARSE_OK) return
        if (maxval(abs(residual_dot - (residual_plus - residual_minus) / &
                (2.0_dp*epsilon))) > 3.0e-8_dp .or. &
            abs(ledger_dot%stored_energy - (ledger_plus%stored_energy - &
                ledger_minus%stored_energy)/(2.0_dp*epsilon)) > 3.0e-8_dp .or. &
            abs(ledger_dot%input_power - (ledger_plus%input_power - &
                ledger_minus%input_power)/(2.0_dp*epsilon)) > 3.0e-8_dp .or. &
            abs(ledger_dot%dissipation - (ledger_plus%dissipation - &
                ledger_minus%dissipation)/(2.0_dp*epsilon)) > 3.0e-8_dp .or. &
            abs(ledger_dot%balance - (ledger_plus%balance - &
                ledger_minus%balance)/(2.0_dp*epsilon)) > 3.0e-8_dp) return

        call assemble_nonlinear_resistive_mhd_residual_vjp( &
            state, value_law, vjp_law, residual_bar, ledger_bar, state_bar, status)
        lhs = dot_product(residual_bar, residual_dot) + &
            ledger_bar%stored_energy*ledger_dot%stored_energy + &
            ledger_bar%input_power*ledger_dot%input_power + &
            ledger_bar%dissipation*ledger_dot%dissipation + &
            ledger_bar%balance*ledger_dot%balance
        rhs = dot_product(state_bar, state_dot)
        if (status%code /= FORTSPARSE_OK .or. abs(lhs - rhs) > 2.0e-12_dp) return
        mhd_case = .true.

    end function mhd_case

    subroutine value_law(state_in, block_id, block_residual, stored_energy, &
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
        if (block_id < 1 .or. block_id > RESISTIVE_MHD_BLOCK_COUNT .or. &
            .not. allocated(active_coefficients) .or. &
            size(state_in) /= active_component .or. &
            size(block_residual) /= active_component) return
        coefficient = active_coefficients(block_id)
        block_residual = coefficient*state_in**2 + &
            0.1_dp*coefficient*state_in
        stored_energy = 0.5_dp*coefficient*sum(state_in**2)
        input_power = coefficient*sum(state_in)
        dissipation = 0.25_dp*coefficient*sum(state_in**2)
        local_status = 0
    end subroutine value_law

    subroutine jvp_law(state_in, state_direction, block_id, block_residual_dot, &
            stored_energy_dot, input_power_dot, dissipation_dot, local_status)
        real(dp), intent(in) :: state_in(:), state_direction(:)
        integer, intent(in) :: block_id
        real(dp), intent(out) :: block_residual_dot(:)
        real(dp), intent(out) :: stored_energy_dot, input_power_dot
        real(dp), intent(out) :: dissipation_dot
        integer, intent(out) :: local_status
        real(dp) :: coefficient

        block_residual_dot = 0.0_dp
        stored_energy_dot = 0.0_dp
        input_power_dot = 0.0_dp
        dissipation_dot = 0.0_dp
        local_status = 1
        if (block_id < 1 .or. block_id > RESISTIVE_MHD_BLOCK_COUNT .or. &
            .not. allocated(active_coefficients) .or. &
            size(state_in) /= active_component .or. &
            size(state_direction) /= active_component .or. &
            size(block_residual_dot) /= active_component) return
        coefficient = active_coefficients(block_id)
        block_residual_dot = (2.0_dp*coefficient*state_in + &
            0.1_dp*coefficient)*state_direction
        stored_energy_dot = coefficient*dot_product(state_in, state_direction)
        input_power_dot = coefficient*sum(state_direction)
        dissipation_dot = 0.5_dp*coefficient*dot_product(state_in, state_direction)
        local_status = 0
    end subroutine jvp_law

    subroutine vjp_law(state_in, block_id, block_residual_bar, stored_energy_bar, &
            input_power_bar, dissipation_bar, balance_bar, state_bar_out, &
            local_status)
        real(dp), intent(in) :: state_in(:), block_residual_bar(:)
        integer, intent(in) :: block_id
        real(dp), intent(in) :: stored_energy_bar, input_power_bar
        real(dp), intent(in) :: dissipation_bar, balance_bar
        real(dp), intent(out) :: state_bar_out(:)
        integer, intent(out) :: local_status
        real(dp) :: coefficient

        state_bar_out = 0.0_dp
        local_status = 1
        if (block_id < 1 .or. block_id > RESISTIVE_MHD_BLOCK_COUNT .or. &
            .not. allocated(active_coefficients) .or. &
            size(state_in) /= active_component .or. &
            size(block_residual_bar) /= active_component .or. &
            size(state_bar_out) /= active_component) return
        coefficient = active_coefficients(block_id)
        state_bar_out = block_residual_bar*(2.0_dp*coefficient*state_in + &
            0.1_dp*coefficient) + (stored_energy_bar + 0.5_dp*dissipation_bar)* &
            coefficient*state_in + (input_power_bar + balance_bar)*coefficient - &
            balance_bar*0.5_dp*coefficient*state_in
        local_status = 0
    end subroutine vjp_law

    subroutine independent_value_oracle(state, coefficients, residual, ledger)
        real(dp), intent(in) :: state(:), coefficients(:)
        real(dp), intent(out) :: residual(:)
        type(nonlinear_resistive_mhd_energy_ledger_t), intent(out) :: ledger
        integer :: block

        residual = 0.0_dp
        ledger = nonlinear_resistive_mhd_energy_ledger_t()
        do block = 1, RESISTIVE_MHD_BLOCK_COUNT
            residual = residual + coefficients(block)*state**2 + &
                0.1_dp*coefficients(block)*state
            ledger%stored_energy = ledger%stored_energy + &
                0.5_dp*coefficients(block)*sum(state**2)
            ledger%input_power = ledger%input_power + &
                coefficients(block)*sum(state)
            ledger%dissipation = ledger%dissipation + &
                0.25_dp*coefficients(block)*sum(state**2)
        end do
        ledger%balance = ledger%input_power - ledger%dissipation
    end subroutine independent_value_oracle

    subroutine independent_jvp_oracle(state, direction, coefficients, residual, ledger)
        real(dp), intent(in) :: state(:), direction(:), coefficients(:)
        real(dp), intent(out) :: residual(:)
        type(nonlinear_resistive_mhd_energy_ledger_t), intent(out) :: ledger
        integer :: block

        residual = 0.0_dp
        ledger = nonlinear_resistive_mhd_energy_ledger_t()
        do block = 1, RESISTIVE_MHD_BLOCK_COUNT
            residual = residual + (2.0_dp*coefficients(block)*state + &
                0.1_dp*coefficients(block))*direction
            ledger%stored_energy = ledger%stored_energy + &
                coefficients(block)*dot_product(state, direction)
            ledger%input_power = ledger%input_power + &
                coefficients(block)*sum(direction)
            ledger%dissipation = ledger%dissipation + &
                0.5_dp*coefficients(block)*dot_product(state, direction)
        end do
        ledger%balance = ledger%input_power - ledger%dissipation
    end subroutine independent_jvp_oracle

    pure logical function ledger_close(left, right, tolerance) result(close)
        type(nonlinear_resistive_mhd_energy_ledger_t), intent(in) :: left, right
        real(dp), intent(in) :: tolerance

        close = abs(left%stored_energy - right%stored_energy) <= tolerance .and. &
            abs(left%input_power - right%input_power) <= tolerance .and. &
            abs(left%dissipation - right%dissipation) <= tolerance .and. &
            abs(left%balance - right%balance) <= tolerance
    end function ledger_close

end program test_nonlinear_resistive_mhd_properties
