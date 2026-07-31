module fortfem_surface_current_constraints
    !! Fixed-topology linear constraints for integrated sheet currents.
    !!
    !! A loop row is an oriented integer combination of manifold currents.
    !! The residual is B_loop K - K_target.  Open-edge balance is supplied by
    !! the junction ledger; this contract supplies closed-loop or prescribed
    !! cycle constraints without assigning physical units or closure laws.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortfem_kinds, only: dp
    use fortsparse, only: fortsparse_status_t, status_set, &
        FORTSPARSE_INVALID_MATRIX, FORTSPARSE_OK
    implicit none
    private

    public :: assemble_surface_current_loop_constraints
    public :: assemble_surface_current_loop_constraints_jvp
    public :: assemble_surface_current_loop_constraints_vjp

contains

    subroutine assemble_surface_current_loop_constraints( &
            loop_basis, manifold_current, target_loop_current, residual, status)
        integer, intent(in) :: loop_basis(:, :)
        real(dp), intent(in) :: manifold_current(:, :), target_loop_current(:, :)
        real(dp), intent(out) :: residual(:, :)
        type(fortsparse_status_t), intent(out) :: status
        integer :: loop, manifold

        residual = 0.0_dp
        call validate_loop_inputs( &
            loop_basis, manifold_current, target_loop_current, residual, status)
        if (status%code /= FORTSPARSE_OK) return

        do loop = 1, size(loop_basis, 1)
            residual(:, loop) = -target_loop_current(:, loop)
            do manifold = 1, size(loop_basis, 2)
                residual(:, loop) = residual(:, loop) + &
                    real(loop_basis(loop, manifold), dp) * &
                    manifold_current(:, manifold)
            end do
        end do
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine assemble_surface_current_loop_constraints

    subroutine assemble_surface_current_loop_constraints_jvp( &
            loop_basis, manifold_current, target_loop_current, &
            manifold_current_dot, target_loop_current_dot, residual_dot, status)
        integer, intent(in) :: loop_basis(:, :)
        real(dp), intent(in) :: manifold_current(:, :), target_loop_current(:, :)
        real(dp), intent(in) :: manifold_current_dot(:, :), target_loop_current_dot(:, :)
        real(dp), intent(out) :: residual_dot(:, :)
        type(fortsparse_status_t), intent(out) :: status
        integer :: loop, manifold

        residual_dot = 0.0_dp
        call validate_loop_inputs( &
            loop_basis, manifold_current, target_loop_current, residual_dot, status)
        if (status%code /= FORTSPARSE_OK) return
        if (size(manifold_current_dot, 1) /= 3 .or. &
            size(manifold_current_dot, 2) /= size(loop_basis, 2) .or. &
            size(target_loop_current_dot, 1) /= 3 .or. &
            size(target_loop_current_dot, 2) /= size(loop_basis, 1) .or. &
            any(.not. ieee_is_finite(manifold_current_dot)) .or. &
            any(.not. ieee_is_finite(target_loop_current_dot))) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "surface-current loop JVP has incompatible increments")
            return
        end if

        do loop = 1, size(loop_basis, 1)
            residual_dot(:, loop) = -target_loop_current_dot(:, loop)
            do manifold = 1, size(loop_basis, 2)
                residual_dot(:, loop) = residual_dot(:, loop) + &
                    real(loop_basis(loop, manifold), dp) * &
                    manifold_current_dot(:, manifold)
            end do
        end do
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine assemble_surface_current_loop_constraints_jvp

    subroutine assemble_surface_current_loop_constraints_vjp( &
            loop_basis, manifold_current, target_loop_current, residual_bar, &
            manifold_current_bar, target_loop_current_bar, status)
        integer, intent(in) :: loop_basis(:, :)
        real(dp), intent(in) :: manifold_current(:, :), target_loop_current(:, :)
        real(dp), intent(in) :: residual_bar(:, :)
        real(dp), intent(out) :: manifold_current_bar(:, :)
        real(dp), intent(out) :: target_loop_current_bar(:, :)
        type(fortsparse_status_t), intent(out) :: status
        integer :: loop, manifold

        manifold_current_bar = 0.0_dp
        target_loop_current_bar = 0.0_dp
        call validate_loop_inputs( &
            loop_basis, manifold_current, target_loop_current, residual_bar, status)
        if (status%code /= FORTSPARSE_OK) return
        if (size(manifold_current_bar, 1) /= 3 .or. &
            size(manifold_current_bar, 2) /= size(loop_basis, 2) .or. &
            size(target_loop_current_bar, 1) /= 3 .or. &
            size(target_loop_current_bar, 2) /= size(loop_basis, 1) .or. &
            any(.not. ieee_is_finite(residual_bar))) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "surface-current loop VJP has incompatible cotangents")
            return
        end if

        do loop = 1, size(loop_basis, 1)
            target_loop_current_bar(:, loop) = -residual_bar(:, loop)
            do manifold = 1, size(loop_basis, 2)
                manifold_current_bar(:, manifold) = &
                    manifold_current_bar(:, manifold) + &
                    real(loop_basis(loop, manifold), dp) * residual_bar(:, loop)
            end do
        end do
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine assemble_surface_current_loop_constraints_vjp

    subroutine validate_loop_inputs( &
            loop_basis, manifold_current, target_loop_current, residual, status)
        integer, intent(in) :: loop_basis(:, :)
        real(dp), intent(in) :: manifold_current(:, :), target_loop_current(:, :)
        real(dp), intent(in) :: residual(:, :)
        type(fortsparse_status_t), intent(out) :: status

        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            "surface-current loop constraint has incompatible arrays")
        if (size(manifold_current, 1) /= 3 .or. &
            size(target_loop_current, 1) /= 3 .or. &
            size(manifold_current, 2) /= size(loop_basis, 2) .or. &
            size(target_loop_current, 2) /= size(loop_basis, 1) .or. &
            size(residual, 1) /= 3 .or. &
            size(residual, 2) /= size(loop_basis, 1)) return
        if (any(.not. ieee_is_finite(manifold_current)) .or. &
            any(.not. ieee_is_finite(target_loop_current)) .or. &
            any(.not. ieee_is_finite(residual))) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "surface-current loop constraint received non-finite data")
            return
        end if
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine validate_loop_inputs

end module fortfem_surface_current_constraints
