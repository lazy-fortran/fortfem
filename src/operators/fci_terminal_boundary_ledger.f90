module fortfem_fci_terminal_boundary_ledger
    !! Vector-valued conservative terminal ledger for FCI support operators.
    !!
    !! Each terminal event belongs to one canonical cell.  For every conserved
    !! component, its signed terminal integral is divided by the canonical
    !! volume for the returned cell contribution.  The separate global ledger
    !! retains the integrated event sum, so clients can compare it directly
    !! with volume and material balances without introducing a model-specific
    !! species, energy, charge, or closure convention.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortfem_kinds, only: dp
    use fortsparse, only: fortsparse_status_t, status_set, &
        FORTSPARSE_INVALID_MATRIX, FORTSPARSE_OK
    implicit none
    private

    public :: assemble_fci_terminal_boundary_ledger
    public :: assemble_fci_terminal_boundary_ledger_jvp
    public :: assemble_fci_terminal_boundary_ledger_vjp

contains

    subroutine assemble_fci_terminal_boundary_ledger( &
            terminal_owners, terminal_weights, terminal_flux, canonical_volumes, &
            contribution, global_ledger, status)
        integer, intent(in) :: terminal_owners(:)
        real(dp), intent(in) :: terminal_weights(:), terminal_flux(:, :)
        real(dp), intent(in) :: canonical_volumes(:)
        real(dp), intent(out) :: contribution(:, :), global_ledger(:)
        type(fortsparse_status_t), intent(out) :: status
        integer :: event, component, owner

        contribution = 0.0_dp
        global_ledger = 0.0_dp
        if (.not. valid_value_inputs(terminal_owners, terminal_weights, &
            terminal_flux, canonical_volumes, contribution, global_ledger, &
            status)) return

        do event = 1, size(terminal_owners)
            owner = terminal_owners(event)
            do component = 1, size(terminal_flux, 2)
                contribution(owner, component) = contribution(owner, component) + &
                    terminal_weights(event)*terminal_flux(event, component)/ &
                    canonical_volumes(owner)
                global_ledger(component) = global_ledger(component) + &
                    terminal_weights(event)*terminal_flux(event, component)
            end do
        end do
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine assemble_fci_terminal_boundary_ledger

    subroutine assemble_fci_terminal_boundary_ledger_jvp( &
            terminal_owners, terminal_weights, terminal_flux, canonical_volumes, &
            terminal_weights_dot, terminal_flux_dot, canonical_volumes_dot, &
            contribution_dot, global_ledger_dot, status)
        integer, intent(in) :: terminal_owners(:)
        real(dp), intent(in) :: terminal_weights(:), terminal_flux(:, :)
        real(dp), intent(in) :: canonical_volumes(:)
        real(dp), intent(in) :: terminal_weights_dot(:), terminal_flux_dot(:, :)
        real(dp), intent(in) :: canonical_volumes_dot(:)
        real(dp), intent(out) :: contribution_dot(:, :), global_ledger_dot(:)
        type(fortsparse_status_t), intent(out) :: status
        integer :: event, component, owner
        real(dp) :: weighted_flux, weighted_flux_dot

        contribution_dot = 0.0_dp
        global_ledger_dot = 0.0_dp
        if (.not. valid_value_inputs(terminal_owners, terminal_weights, &
            terminal_flux, canonical_volumes, contribution_dot, &
            global_ledger_dot, status)) return
        if (size(terminal_weights_dot) /= size(terminal_owners) .or. &
            size(terminal_flux_dot, 1) /= size(terminal_flux, 1) .or. &
            size(terminal_flux_dot, 2) /= size(terminal_flux, 2) .or. &
            size(canonical_volumes_dot) /= size(canonical_volumes) .or. &
            any(.not. ieee_is_finite(terminal_weights_dot)) .or. &
            any(.not. ieee_is_finite(terminal_flux_dot)) .or. &
            any(.not. ieee_is_finite(canonical_volumes_dot))) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "vector FCI terminal ledger JVP has incompatible tangents")
            return
        end if

        do event = 1, size(terminal_owners)
            owner = terminal_owners(event)
            do component = 1, size(terminal_flux, 2)
                weighted_flux = terminal_weights(event)*terminal_flux(event, component)
                weighted_flux_dot = terminal_weights_dot(event)* &
                    terminal_flux(event, component) + terminal_weights(event)* &
                    terminal_flux_dot(event, component)
                contribution_dot(owner, component) = &
                    contribution_dot(owner, component) + weighted_flux_dot/ &
                    canonical_volumes(owner) - weighted_flux* &
                    canonical_volumes_dot(owner)/canonical_volumes(owner)**2
                global_ledger_dot(component) = global_ledger_dot(component) + &
                    weighted_flux_dot
            end do
        end do
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine assemble_fci_terminal_boundary_ledger_jvp

    subroutine assemble_fci_terminal_boundary_ledger_vjp( &
            terminal_owners, terminal_weights, terminal_flux, canonical_volumes, &
            contribution_bar, global_ledger_bar, terminal_weights_bar, &
            terminal_flux_bar, canonical_volumes_bar, status)
        integer, intent(in) :: terminal_owners(:)
        real(dp), intent(in) :: terminal_weights(:), terminal_flux(:, :)
        real(dp), intent(in) :: canonical_volumes(:)
        real(dp), intent(in) :: contribution_bar(:, :), global_ledger_bar(:)
        real(dp), intent(out) :: terminal_weights_bar(:), terminal_flux_bar(:, :)
        real(dp), intent(out) :: canonical_volumes_bar(:)
        type(fortsparse_status_t), intent(out) :: status
        integer :: event, component, owner
        real(dp) :: cotangent, inverse_volume

        terminal_weights_bar = 0.0_dp
        terminal_flux_bar = 0.0_dp
        canonical_volumes_bar = 0.0_dp
        if (.not. valid_value_inputs(terminal_owners, terminal_weights, &
            terminal_flux, canonical_volumes, contribution_bar, &
            global_ledger_bar, status)) return
        if (size(terminal_weights_bar) /= size(terminal_owners) .or. &
            size(terminal_flux_bar, 1) /= size(terminal_flux, 1) .or. &
            size(terminal_flux_bar, 2) /= size(terminal_flux, 2) .or. &
            size(canonical_volumes_bar) /= size(canonical_volumes) .or. &
            any(.not. ieee_is_finite(contribution_bar)) .or. &
            any(.not. ieee_is_finite(global_ledger_bar))) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "vector FCI terminal ledger VJP has incompatible cotangents")
            return
        end if

        do event = 1, size(terminal_owners)
            owner = terminal_owners(event)
            inverse_volume = 1.0_dp/canonical_volumes(owner)
            do component = 1, size(terminal_flux, 2)
                cotangent = contribution_bar(owner, component)
                terminal_weights_bar(event) = terminal_weights_bar(event) + &
                    cotangent*terminal_flux(event, component)*inverse_volume + &
                    global_ledger_bar(component)*terminal_flux(event, component)
                terminal_flux_bar(event, component) = &
                    terminal_flux_bar(event, component) + &
                    cotangent*terminal_weights(event)*inverse_volume + &
                    global_ledger_bar(component)*terminal_weights(event)
                canonical_volumes_bar(owner) = canonical_volumes_bar(owner) - &
                    cotangent*terminal_weights(event)*terminal_flux(event, component)* &
                    inverse_volume**2
            end do
        end do
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine assemble_fci_terminal_boundary_ledger_vjp

    logical function valid_value_inputs( &
            terminal_owners, terminal_weights, terminal_flux, canonical_volumes, &
            contribution, global_ledger, status)
        integer, intent(in) :: terminal_owners(:)
        real(dp), intent(in) :: terminal_weights(:), terminal_flux(:, :)
        real(dp), intent(in) :: canonical_volumes(:)
        real(dp), intent(in) :: contribution(:, :), global_ledger(:)
        type(fortsparse_status_t), intent(out) :: status

        valid_value_inputs = .false.
        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            "vector FCI terminal ledger has incompatible arrays")
        if (size(terminal_owners) < 1 .or. size(canonical_volumes) < 1 .or. &
            size(terminal_weights) /= size(terminal_owners) .or. &
            size(terminal_flux, 1) /= size(terminal_owners) .or. &
            size(terminal_flux, 2) < 1 .or. &
            size(contribution, 1) /= size(canonical_volumes) .or. &
            size(contribution, 2) /= size(terminal_flux, 2) .or. &
            size(global_ledger) /= size(terminal_flux, 2)) return
        if (any(terminal_owners < 1) .or. &
            any(terminal_owners > size(canonical_volumes))) return
        if (any(.not. ieee_is_finite(terminal_weights)) .or. &
            any(.not. ieee_is_finite(terminal_flux)) .or. &
            any(.not. ieee_is_finite(canonical_volumes)) .or. &
            any(.not. ieee_is_finite(contribution)) .or. &
            any(.not. ieee_is_finite(global_ledger))) return
        if (any(terminal_weights <= 0.0_dp) .or. &
            any(canonical_volumes <= 0.0_dp)) return
        valid_value_inputs = .true.
        call status_set(status, FORTSPARSE_OK, "")
    end function valid_value_inputs

end module fortfem_fci_terminal_boundary_ledger
