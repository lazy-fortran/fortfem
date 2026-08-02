module fortfem_nonlinear_resistive_mhd_composition
    !! Closure-neutral nonlinear resistive-MHD block composition.
    !!
    !! A client supplies one nonlinear callback for each named block.  The
    !! callbacks own constitutive laws, units, geometry, state selection, and
    !! any Faraday/Ampere-compatible FEEC, tensor, FEM/BEM, DtN, PML, wall, or
    !! free-boundary discretization.  FortFEM only composes their common state
    !! residuals and an explicit energy/input/dissipation ledger.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortfem_kinds, only: dp
    use fortsparse, only: FORTSPARSE_INVALID_MATRIX, FORTSPARSE_OK, &
        fortsparse_status_t, status_set
    implicit none
    private

    integer, parameter, public :: RESISTIVE_MHD_FARADAY = 1
    integer, parameter, public :: RESISTIVE_MHD_AMPERE = 2
    integer, parameter, public :: RESISTIVE_MHD_MOMENTUM = 3
    integer, parameter, public :: RESISTIVE_MHD_PRESSURE = 4
    integer, parameter, public :: RESISTIVE_MHD_TENSOR = 5
    integer, parameter, public :: RESISTIVE_MHD_ANISOTROPIC_TRANSPORT = 6
    integer, parameter, public :: RESISTIVE_MHD_WALL = 7
    integer, parameter, public :: RESISTIVE_MHD_FREE_BOUNDARY = 8
    integer, parameter, public :: RESISTIVE_MHD_BLOCK_COUNT = 8

    type, public :: nonlinear_resistive_mhd_energy_ledger_t
        !! Additive diagnostics returned by all nonlinear blocks.
        !!
        !! `balance = input_power - dissipation`; a caller may compare this
        !! instantaneous power balance with its own time derivative of
        !! `stored_energy`.  FortFEM does not infer a time integrator or a
        !! plasma closure from these values.
        real(dp) :: stored_energy = 0.0_dp
        real(dp) :: input_power = 0.0_dp
        real(dp) :: dissipation = 0.0_dp
        real(dp) :: balance = 0.0_dp
    end type nonlinear_resistive_mhd_energy_ledger_t

    abstract interface
        subroutine nonlinear_resistive_mhd_value_law( &
                state, block_id, block_residual, stored_energy, input_power, &
                dissipation, status)
            import dp
            real(dp), intent(in) :: state(:)
            integer, intent(in) :: block_id
            real(dp), intent(out) :: block_residual(:)
            real(dp), intent(out) :: stored_energy, input_power, dissipation
            integer, intent(out) :: status
        end subroutine nonlinear_resistive_mhd_value_law

        subroutine nonlinear_resistive_mhd_jvp_law( &
                state, state_dot, block_id, block_residual_dot, stored_energy_dot, &
                input_power_dot, dissipation_dot, status)
            import dp
            real(dp), intent(in) :: state(:), state_dot(:)
            integer, intent(in) :: block_id
            real(dp), intent(out) :: block_residual_dot(:)
            real(dp), intent(out) :: stored_energy_dot, input_power_dot
            real(dp), intent(out) :: dissipation_dot
            integer, intent(out) :: status
        end subroutine nonlinear_resistive_mhd_jvp_law

        subroutine nonlinear_resistive_mhd_vjp_law( &
                state, block_id, block_residual_bar, stored_energy_bar, &
                input_power_bar, dissipation_bar, balance_bar, state_bar, status)
            import dp
            real(dp), intent(in) :: state(:), block_residual_bar(:)
            integer, intent(in) :: block_id
            real(dp), intent(in) :: stored_energy_bar, input_power_bar
            real(dp), intent(in) :: dissipation_bar, balance_bar
            real(dp), intent(out) :: state_bar(:)
            integer, intent(out) :: status
        end subroutine nonlinear_resistive_mhd_vjp_law
    end interface

    public :: assemble_nonlinear_resistive_mhd_residual
    public :: assemble_nonlinear_resistive_mhd_residual_jvp
    public :: assemble_nonlinear_resistive_mhd_residual_vjp

contains

    subroutine assemble_nonlinear_resistive_mhd_residual( &
            state, value_law, residual, ledger, status)
        !! Compose all named nonlinear block values and diagnostics.
        real(dp), intent(in) :: state(:)
        procedure(nonlinear_resistive_mhd_value_law) :: value_law
        real(dp), intent(out) :: residual(:)
        type(nonlinear_resistive_mhd_energy_ledger_t), intent(out) :: ledger
        type(fortsparse_status_t), intent(out) :: status

        real(dp) :: block_residual(size(residual))
        real(dp) :: block_energy, block_input, block_dissipation
        integer :: block, local_status

        residual = 0.0_dp
        ledger = nonlinear_resistive_mhd_energy_ledger_t()
        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            "nonlinear resistive-MHD composition received incompatible state")
        if (.not. valid_value_inputs(state, residual)) return

        do block = 1, RESISTIVE_MHD_BLOCK_COUNT
            block_residual = 0.0_dp
            block_energy = 0.0_dp
            block_input = 0.0_dp
            block_dissipation = 0.0_dp
            call value_law(state, block, block_residual, block_energy, block_input, &
                block_dissipation, local_status)
            if (local_status /= 0 .or. any(.not. ieee_is_finite(block_residual)) .or. &
                .not. ieee_is_finite(block_energy) .or. &
                .not. ieee_is_finite(block_input) .or. &
                .not. ieee_is_finite(block_dissipation) .or. block_dissipation < 0.0_dp) then
                residual = 0.0_dp
                ledger = nonlinear_resistive_mhd_energy_ledger_t()
                return
            end if
            residual = residual + block_residual
            ledger%stored_energy = ledger%stored_energy + block_energy
            ledger%input_power = ledger%input_power + block_input
            ledger%dissipation = ledger%dissipation + block_dissipation
        end do
        ledger%balance = ledger%input_power - ledger%dissipation
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine assemble_nonlinear_resistive_mhd_residual

    subroutine assemble_nonlinear_resistive_mhd_residual_jvp( &
            state, state_dot, value_law, jvp_law, residual_dot, ledger_dot, status)
        !! Compose exact directional derivatives of all nonlinear blocks.
        real(dp), intent(in) :: state(:), state_dot(:)
        procedure(nonlinear_resistive_mhd_value_law) :: value_law
        procedure(nonlinear_resistive_mhd_jvp_law) :: jvp_law
        real(dp), intent(out) :: residual_dot(:)
        type(nonlinear_resistive_mhd_energy_ledger_t), intent(out) :: ledger_dot
        type(fortsparse_status_t), intent(out) :: status

        real(dp) :: block_residual(size(residual_dot))
        real(dp) :: block_residual_dot(size(residual_dot))
        real(dp) :: block_energy, block_input, block_dissipation
        real(dp) :: block_energy_dot, block_input_dot, block_dissipation_dot
        integer :: block, local_status

        residual_dot = 0.0_dp
        ledger_dot = nonlinear_resistive_mhd_energy_ledger_t()
        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            "nonlinear resistive-MHD JVP received incompatible state")
        if (.not. valid_value_inputs(state, residual_dot) .or. &
            size(state_dot) /= size(state) .or. any(.not. ieee_is_finite(state_dot))) return

        do block = 1, RESISTIVE_MHD_BLOCK_COUNT
            block_residual = 0.0_dp
            block_energy = 0.0_dp
            block_input = 0.0_dp
            block_dissipation = 0.0_dp
            call value_law(state, block, block_residual, block_energy, block_input, &
                block_dissipation, local_status)
            if (local_status /= 0 .or. any(.not. ieee_is_finite(block_residual)) .or. &
                .not. ieee_is_finite(block_energy) .or. &
                .not. ieee_is_finite(block_input) .or. &
                .not. ieee_is_finite(block_dissipation) .or. block_dissipation < 0.0_dp) then
                residual_dot = 0.0_dp
                ledger_dot = nonlinear_resistive_mhd_energy_ledger_t()
                return
            end if
            block_residual_dot = 0.0_dp
            block_energy_dot = 0.0_dp
            block_input_dot = 0.0_dp
            block_dissipation_dot = 0.0_dp
            call jvp_law(state, state_dot, block, block_residual_dot, block_energy_dot, &
                block_input_dot, block_dissipation_dot, local_status)
            if (local_status /= 0 .or. any(.not. ieee_is_finite(block_residual_dot)) .or. &
                .not. ieee_is_finite(block_energy_dot) .or. &
                .not. ieee_is_finite(block_input_dot) .or. &
                .not. ieee_is_finite(block_dissipation_dot)) then
                residual_dot = 0.0_dp
                ledger_dot = nonlinear_resistive_mhd_energy_ledger_t()
                return
            end if
            residual_dot = residual_dot + block_residual_dot
            ledger_dot%stored_energy = ledger_dot%stored_energy + block_energy_dot
            ledger_dot%input_power = ledger_dot%input_power + block_input_dot
            ledger_dot%dissipation = ledger_dot%dissipation + block_dissipation_dot
        end do
        ledger_dot%balance = ledger_dot%input_power - ledger_dot%dissipation
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine assemble_nonlinear_resistive_mhd_residual_jvp

    subroutine assemble_nonlinear_resistive_mhd_residual_vjp( &
            state, value_law, vjp_law, residual_bar, ledger_bar, state_bar, status)
        !! Compose the real adjoint of residual and energy-ledger outputs.
        real(dp), intent(in) :: state(:)
        procedure(nonlinear_resistive_mhd_value_law) :: value_law
        procedure(nonlinear_resistive_mhd_vjp_law) :: vjp_law
        real(dp), intent(in) :: residual_bar(:)
        type(nonlinear_resistive_mhd_energy_ledger_t), intent(in) :: ledger_bar
        real(dp), intent(out) :: state_bar(:)
        type(fortsparse_status_t), intent(out) :: status

        real(dp) :: block_residual(size(residual_bar)), block_state_bar(size(state))
        real(dp) :: block_energy, block_input, block_dissipation
        integer :: block, local_status

        state_bar = 0.0_dp
        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            "nonlinear resistive-MHD VJP received incompatible state")
        if (.not. valid_value_inputs(state, residual_bar) .or. &
            size(state_bar) /= size(state) .or. &
            .not. finite_ledger(ledger_bar)) return

        do block = 1, RESISTIVE_MHD_BLOCK_COUNT
            block_residual = 0.0_dp
            block_energy = 0.0_dp
            block_input = 0.0_dp
            block_dissipation = 0.0_dp
            call value_law(state, block, block_residual, block_energy, block_input, &
                block_dissipation, local_status)
            if (local_status /= 0 .or. any(.not. ieee_is_finite(block_residual)) .or. &
                .not. ieee_is_finite(block_energy) .or. &
                .not. ieee_is_finite(block_input) .or. &
                .not. ieee_is_finite(block_dissipation) .or. block_dissipation < 0.0_dp) then
                state_bar = 0.0_dp
                return
            end if
            block_state_bar = 0.0_dp
            call vjp_law(state, block, residual_bar, ledger_bar%stored_energy, &
                ledger_bar%input_power, ledger_bar%dissipation, ledger_bar%balance, &
                block_state_bar, local_status)
            if (local_status /= 0 .or. any(.not. ieee_is_finite(block_state_bar))) then
                state_bar = 0.0_dp
                return
            end if
            state_bar = state_bar + block_state_bar
        end do
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine assemble_nonlinear_resistive_mhd_residual_vjp

    pure logical function valid_value_inputs(state, output) result(valid)
        real(dp), intent(in) :: state(:), output(:)

        valid = size(state) > 0 .and. size(output) == size(state) .and. &
            all(ieee_is_finite(state)) .and. all(ieee_is_finite(output))
    end function valid_value_inputs

    pure logical function finite_ledger(ledger) result(valid)
        type(nonlinear_resistive_mhd_energy_ledger_t), intent(in) :: ledger

        valid = ieee_is_finite(ledger%stored_energy) .and. &
            ieee_is_finite(ledger%input_power) .and. &
            ieee_is_finite(ledger%dissipation) .and. ieee_is_finite(ledger%balance)
    end function finite_ledger

end module fortfem_nonlinear_resistive_mhd_composition
