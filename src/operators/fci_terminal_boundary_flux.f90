module fortfem_fci_terminal_boundary_flux
    !! Conservative terminal boundary contribution for FCI support operators.
    !!
    !! Each terminal event belongs to one canonical cell.  Its signed outward
    !! flux, multiplied by a positive terminal measure, is divided by that
    !! cell's canonical volume.  Therefore the canonical-volume weighted sum of
    !! the returned contribution equals the integrated terminal flux.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortfem_kinds, only: dp
    use fortsparse, only: fortsparse_status_t, status_set, &
        FORTSPARSE_INVALID_MATRIX, FORTSPARSE_OK
    implicit none

    private

    public :: assemble_fci_terminal_boundary_flux
    public :: assemble_fci_terminal_boundary_flux_jvp
    public :: assemble_fci_terminal_boundary_flux_vjp

contains

    subroutine assemble_fci_terminal_boundary_flux( &
            terminal_owners, terminal_weights, terminal_flux, canonical_volumes, &
            contribution, status)
        !! Assemble `W_c^{-1}` times signed terminal flux integrals.
        integer, intent(in) :: terminal_owners(:)
        real(dp), intent(in) :: terminal_weights(:), terminal_flux(:)
        real(dp), intent(in) :: canonical_volumes(:)
        real(dp), intent(out) :: contribution(:)
        type(fortsparse_status_t), intent(out) :: status

        integer :: event, owner

        contribution = 0.0_dp
        if (.not. valid_terminal_inputs( &
            terminal_owners, terminal_weights, terminal_flux, canonical_volumes, &
            contribution, status)) return
        do event = 1, size(terminal_owners)
            owner = terminal_owners(event)
            contribution(owner) = contribution(owner) + &
                terminal_weights(event)*terminal_flux(event)/canonical_volumes(owner)
        end do
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine assemble_fci_terminal_boundary_flux

    subroutine assemble_fci_terminal_boundary_flux_jvp( &
            terminal_owners, terminal_weights, terminal_flux, canonical_volumes, &
            terminal_weights_dot, terminal_flux_dot, canonical_volumes_dot, &
            contribution_dot, status)
        !! Apply the fixed-owner JVP of the terminal flux contribution.
        integer, intent(in) :: terminal_owners(:)
        real(dp), intent(in) :: terminal_weights(:), terminal_flux(:)
        real(dp), intent(in) :: canonical_volumes(:)
        real(dp), intent(in) :: terminal_weights_dot(:), terminal_flux_dot(:)
        real(dp), intent(in) :: canonical_volumes_dot(:)
        real(dp), intent(out) :: contribution_dot(:)
        type(fortsparse_status_t), intent(out) :: status

        integer :: event, owner
        real(dp) :: weighted_flux, weighted_flux_dot

        contribution_dot = 0.0_dp
        if (.not. valid_terminal_inputs( &
            terminal_owners, terminal_weights, terminal_flux, canonical_volumes, &
            contribution_dot, status)) return
        if (size(terminal_weights_dot) /= size(terminal_owners) .or. &
            size(terminal_flux_dot) /= size(terminal_owners) .or. &
            size(canonical_volumes_dot) /= size(canonical_volumes)) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "FCI terminal flux JVP received incompatible tangents")
            return
        end if
        if (any(.not. ieee_is_finite(terminal_weights_dot)) .or. &
            any(.not. ieee_is_finite(terminal_flux_dot)) .or. &
            any(.not. ieee_is_finite(canonical_volumes_dot))) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "FCI terminal flux JVP received non-finite tangents")
            return
        end if
        do event = 1, size(terminal_owners)
            owner = terminal_owners(event)
            weighted_flux = terminal_weights(event)*terminal_flux(event)
            weighted_flux_dot = terminal_weights_dot(event)*terminal_flux(event) + &
                terminal_weights(event)*terminal_flux_dot(event)
            contribution_dot(owner) = contribution_dot(owner) + &
                weighted_flux_dot/canonical_volumes(owner) - &
                weighted_flux*canonical_volumes_dot(owner)/canonical_volumes(owner)**2
        end do
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine assemble_fci_terminal_boundary_flux_jvp

    subroutine assemble_fci_terminal_boundary_flux_vjp( &
            terminal_owners, terminal_weights, terminal_flux, canonical_volumes, &
            contribution_bar, terminal_weights_bar, terminal_flux_bar, &
            canonical_volumes_bar, status)
        !! Apply the real VJP of the terminal flux contribution.
        integer, intent(in) :: terminal_owners(:)
        real(dp), intent(in) :: terminal_weights(:), terminal_flux(:)
        real(dp), intent(in) :: canonical_volumes(:), contribution_bar(:)
        real(dp), intent(out) :: terminal_weights_bar(:), terminal_flux_bar(:)
        real(dp), intent(out) :: canonical_volumes_bar(:)
        type(fortsparse_status_t), intent(out) :: status

        integer :: event, owner
        real(dp) :: cotangent, inverse_volume

        terminal_weights_bar = 0.0_dp
        terminal_flux_bar = 0.0_dp
        canonical_volumes_bar = 0.0_dp
        if (.not. valid_terminal_inputs( &
            terminal_owners, terminal_weights, terminal_flux, canonical_volumes, &
            contribution_bar, status)) return
        if (size(terminal_weights_bar) /= size(terminal_owners) .or. &
            size(terminal_flux_bar) /= size(terminal_owners) .or. &
            size(canonical_volumes_bar) /= size(canonical_volumes)) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "FCI terminal flux VJP received incompatible cotangents")
            return
        end if
        if (any(.not. ieee_is_finite(contribution_bar))) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "FCI terminal flux VJP received non-finite cotangents")
            return
        end if
        do event = 1, size(terminal_owners)
            owner = terminal_owners(event)
            cotangent = contribution_bar(owner)
            inverse_volume = 1.0_dp/canonical_volumes(owner)
            terminal_weights_bar(event) = terminal_weights_bar(event) + &
                cotangent*terminal_flux(event)*inverse_volume
            terminal_flux_bar(event) = terminal_flux_bar(event) + &
                cotangent*terminal_weights(event)*inverse_volume
            canonical_volumes_bar(owner) = canonical_volumes_bar(owner) - &
                cotangent*terminal_weights(event)*terminal_flux(event)* &
                inverse_volume**2
        end do
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine assemble_fci_terminal_boundary_flux_vjp

    logical function valid_terminal_inputs( &
            terminal_owners, terminal_weights, terminal_flux, canonical_volumes, &
            output, status)
        integer, intent(in) :: terminal_owners(:)
        real(dp), intent(in) :: terminal_weights(:), terminal_flux(:)
        real(dp), intent(in) :: canonical_volumes(:), output(:)
        type(fortsparse_status_t), intent(out) :: status

        valid_terminal_inputs = .false.
        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            "FCI terminal flux received incompatible arrays")
        if (size(canonical_volumes) < 1 .or. size(output) /= size(canonical_volumes)) return
        if (size(terminal_weights) /= size(terminal_owners) .or. &
            size(terminal_flux) /= size(terminal_owners)) return
        if (any(terminal_owners < 1) .or. &
            any(terminal_owners > size(canonical_volumes))) return
        if (any(.not. ieee_is_finite(terminal_weights)) .or. &
            any(.not. ieee_is_finite(terminal_flux)) .or. &
            any(.not. ieee_is_finite(canonical_volumes))) return
        if (any(.not. ieee_is_finite(output))) return
        if (any(terminal_weights <= 0.0_dp) .or. &
            any(canonical_volumes <= 0.0_dp)) return
        valid_terminal_inputs = .true.
        call status_set(status, FORTSPARSE_OK, "")
    end function valid_terminal_inputs

end module fortfem_fci_terminal_boundary_flux
