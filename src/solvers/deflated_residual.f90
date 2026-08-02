module fortfem_deflated_residual
    !! Smooth root-deflation wrapper for caller-owned nonlinear residuals.
    !!
    !! For fixed reference states z_j, this module returns
    !!
    !!   F_d(x) = D(x) F(x),
    !!   D(x) = 1 + scale sum_j (||x-z_j||^2 + shift^2)^(-power/2).
    !!
    !! The base residual and its directional derivative remain caller-owned.
    !! A positive shift makes the wrapper smooth at every state, including a
    !! previously converged root.  Root selection, continuation policy, and
    !! nonlinear solver acceptance are deliberately outside this contract.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortfem_kinds, only: dp
    use fortsparse, only: fortsparse_status_t, status_set, &
        FORTSPARSE_INVALID_MATRIX, FORTSPARSE_OK
    implicit none
    private

    public :: assemble_deflated_residual
    public :: assemble_deflated_residual_jvp
    public :: assemble_deflated_residual_vjp

contains

    subroutine assemble_deflated_residual( &
            residual, state, reference_states, scale, power, shift, &
            deflated_residual, deflation_factor, status)
        real(dp), intent(in) :: residual(:), state(:), reference_states(:, :)
        real(dp), intent(in) :: scale, power, shift
        real(dp), intent(out) :: deflated_residual(:), deflation_factor
        type(fortsparse_status_t), intent(out) :: status
        real(dp) :: gradient(size(state))

        deflated_residual = 0.0_dp
        deflation_factor = 0.0_dp
        if (.not. valid_inputs( &
            residual, state, reference_states, scale, power, shift, &
            deflated_residual)) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "deflated residual has incompatible inputs")
            return
        end if
        call evaluate_factor_and_gradient( &
            state, reference_states, scale, power, shift, &
            deflation_factor, gradient)
        deflated_residual = deflation_factor*residual
        if (.not. ieee_is_finite(deflation_factor) .or. &
            .not. all(ieee_is_finite(deflated_residual))) then
            deflated_residual = 0.0_dp
            deflation_factor = 0.0_dp
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "deflated residual is non-finite")
            return
        end if
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine assemble_deflated_residual

    subroutine assemble_deflated_residual_jvp( &
            residual, state, reference_states, scale, power, shift, residual_dot, &
            state_dot, deflated_residual_dot, deflation_factor_dot, status)
        real(dp), intent(in) :: residual(:), state(:), reference_states(:, :)
        real(dp), intent(in) :: scale, power, shift, residual_dot(:), state_dot(:)
        real(dp), intent(out) :: deflated_residual_dot(:), deflation_factor_dot
        type(fortsparse_status_t), intent(out) :: status
        real(dp) :: deflation_factor, gradient(size(state))

        deflated_residual_dot = 0.0_dp
        deflation_factor_dot = 0.0_dp
        if (.not. valid_inputs( &
            residual, state, reference_states, scale, power, shift, &
            deflated_residual_dot) .or. size(residual_dot) /= size(residual) .or. &
            size(state_dot) /= size(state) .or. .not. all(ieee_is_finite(residual_dot)) &
            .or. .not. all(ieee_is_finite(state_dot))) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "deflated residual JVP has incompatible inputs")
            return
        end if
        call evaluate_factor_and_gradient( &
            state, reference_states, scale, power, shift, deflation_factor, gradient)
        deflation_factor_dot = dot_product(gradient, state_dot)
        deflated_residual_dot = deflation_factor*residual_dot + &
            deflation_factor_dot*residual
        if (.not. ieee_is_finite(deflation_factor_dot) .or. &
            .not. all(ieee_is_finite(deflated_residual_dot))) then
            deflated_residual_dot = 0.0_dp
            deflation_factor_dot = 0.0_dp
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "deflated residual JVP is non-finite")
            return
        end if
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine assemble_deflated_residual_jvp

    subroutine assemble_deflated_residual_vjp( &
            residual, state, reference_states, scale, power, shift, &
            deflated_residual_bar, residual_bar, state_bar, status)
        real(dp), intent(in) :: residual(:), state(:), reference_states(:, :)
        real(dp), intent(in) :: scale, power, shift, deflated_residual_bar(:)
        real(dp), intent(out) :: residual_bar(:), state_bar(:)
        type(fortsparse_status_t), intent(out) :: status
        real(dp) :: deflation_factor, gradient(size(state)), scalar_bar

        residual_bar = 0.0_dp
        state_bar = 0.0_dp
        if (.not. valid_inputs( &
            residual, state, reference_states, scale, power, shift, residual_bar) .or. &
            size(deflated_residual_bar) /= size(residual) .or. &
            size(state_bar) /= size(state) .or. &
            .not. all(ieee_is_finite(deflated_residual_bar))) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "deflated residual VJP has incompatible inputs")
            return
        end if
        call evaluate_factor_and_gradient( &
            state, reference_states, scale, power, shift, deflation_factor, gradient)
        residual_bar = deflation_factor*deflated_residual_bar
        scalar_bar = dot_product(deflated_residual_bar, residual)
        state_bar = scalar_bar*gradient
        if (.not. all(ieee_is_finite(residual_bar)) .or. &
            .not. all(ieee_is_finite(state_bar))) then
            residual_bar = 0.0_dp
            state_bar = 0.0_dp
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "deflated residual VJP is non-finite")
            return
        end if
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine assemble_deflated_residual_vjp

    subroutine evaluate_factor_and_gradient( &
            state, reference_states, scale, power, shift, factor, gradient)
        real(dp), intent(in) :: state(:), reference_states(:, :)
        real(dp), intent(in) :: scale, power, shift
        real(dp), intent(out) :: factor, gradient(:)
        real(dp) :: delta(size(state)), squared_distance, denominator
        integer :: root

        factor = 1.0_dp
        gradient = 0.0_dp
        do root = 1, size(reference_states, 2)
            delta = state - reference_states(:, root)
            squared_distance = dot_product(delta, delta)
            denominator = squared_distance + shift**2
            factor = factor + scale*denominator**(-0.5_dp*power)
            gradient = gradient - scale*power* &
                denominator**(-0.5_dp*power - 1.0_dp)*delta
        end do
    end subroutine evaluate_factor_and_gradient

    logical function valid_inputs( &
            residual, state, reference_states, scale, power, shift, target) result(valid)
        real(dp), intent(in) :: residual(:), state(:), reference_states(:, :)
        real(dp), intent(in) :: scale, power, shift, target(:)

        valid = size(residual) > 0 .and. size(state) == size(residual) .and. &
            size(reference_states, 1) == size(state) .and. &
            size(reference_states, 2) > 0 .and. scale >= 0.0_dp .and. &
            power > 0.0_dp .and. shift > 0.0_dp .and. size(target) == size(residual)
        if (.not. valid) return
        valid = ieee_is_finite(scale) .and. ieee_is_finite(power) .and. &
            ieee_is_finite(shift) .and. all(ieee_is_finite(residual)) .and. &
            all(ieee_is_finite(state)) .and. all(ieee_is_finite(reference_states))
    end function valid_inputs

end module fortfem_deflated_residual
