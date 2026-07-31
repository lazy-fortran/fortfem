module fortfem_surface_current_trace_residual
    !! Work-conjugate residual for a caller-owned surface-current trace space.
    !!
    !! For independent test and trial vector traces, coefficients c, and a
    !! target current K_target, the residual is
    !!
    !!   r_i = sum_q w_q dot(T_i(q), sum_j S_j(q)c_j - K_target(q)).
    !!
    !! This is a neutral surface L2 pairing.  It gives fitted, cut, DG, and
    !! IGA callers an explicit ownership point for a tangential sheet-current
    !! unknown without selecting a constitutive law or a physical closure.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortfem_kinds, only: dp
    use fortsparse, only: fortsparse_status_t, status_set, &
        FORTSPARSE_INVALID_MATRIX, FORTSPARSE_OK
    implicit none
    private

    public :: assemble_surface_current_trace_residual
    public :: assemble_surface_current_trace_residual_jvp
    public :: assemble_surface_current_trace_residual_vjp

contains

    subroutine assemble_surface_current_trace_residual( &
            test_basis, trial_basis, surface_weights, coefficients, &
            target_current, residual, status)
        !! Assemble the weighted trace-space current residual.
        real(dp), intent(in) :: test_basis(:, :, :), trial_basis(:, :, :)
        real(dp), intent(in) :: surface_weights(:), coefficients(:)
        real(dp), intent(in) :: target_current(:, :)
        real(dp), intent(out) :: residual(:)
        type(fortsparse_status_t), intent(out) :: status

        integer :: quadrature_count, test_count, trial_count
        integer :: quadrature, test_dof, trial_dof
        real(dp) :: current(3), state(3)

        residual = 0.0_dp
        call validate_trace_residual_inputs( &
            test_basis, trial_basis, surface_weights, coefficients, &
            target_current, quadrature_count, test_count, trial_count, status)
        if (status%code /= FORTSPARSE_OK) return
        if (size(residual) /= test_count) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "surface-current trace residual output has incompatible size")
            return
        end if

        do quadrature = 1, quadrature_count
            current = 0.0_dp
            do trial_dof = 1, trial_count
                current = current + trial_basis(quadrature, trial_dof, :)* &
                    coefficients(trial_dof)
            end do
            state = current - target_current(quadrature, :)
            do test_dof = 1, test_count
                residual(test_dof) = residual(test_dof) + &
                    surface_weights(quadrature)*dot_product( &
                    test_basis(quadrature, test_dof, :), state)
            end do
        end do
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine assemble_surface_current_trace_residual

    subroutine assemble_surface_current_trace_residual_jvp( &
            test_basis, trial_basis, surface_weights, coefficients, &
            target_current, test_basis_dot, trial_basis_dot, &
            surface_weights_dot, coefficients_dot, target_current_dot, &
            residual_dot, status)
        !! Apply the full product-rule JVP of the trace residual.
        real(dp), intent(in) :: test_basis(:, :, :), trial_basis(:, :, :)
        real(dp), intent(in) :: surface_weights(:), coefficients(:)
        real(dp), intent(in) :: target_current(:, :)
        real(dp), intent(in) :: test_basis_dot(:, :, :), trial_basis_dot(:, :, :)
        real(dp), intent(in) :: surface_weights_dot(:), coefficients_dot(:)
        real(dp), intent(in) :: target_current_dot(:, :)
        real(dp), intent(out) :: residual_dot(:)
        type(fortsparse_status_t), intent(out) :: status

        integer :: quadrature_count, test_count, trial_count
        integer :: quadrature, test_dof, trial_dof
        real(dp) :: current(3), state(3), current_dot(3), state_dot(3)

        residual_dot = 0.0_dp
        call validate_trace_residual_inputs( &
            test_basis, trial_basis, surface_weights, coefficients, &
            target_current, quadrature_count, test_count, trial_count, status)
        if (status%code /= FORTSPARSE_OK) return
        if (size(residual_dot) /= test_count .or. &
            size(test_basis_dot, 1) /= quadrature_count .or. &
            size(test_basis_dot, 2) /= test_count .or. &
            size(test_basis_dot, 3) /= 3 .or. &
            size(trial_basis_dot, 1) /= quadrature_count .or. &
            size(trial_basis_dot, 2) /= trial_count .or. &
            size(trial_basis_dot, 3) /= 3 .or. &
            size(surface_weights_dot) /= quadrature_count .or. &
            size(coefficients_dot) /= trial_count .or. &
            size(target_current_dot, 1) /= quadrature_count .or. &
            size(target_current_dot, 2) /= 3 .or. &
            any(.not. ieee_is_finite(test_basis_dot)) .or. &
            any(.not. ieee_is_finite(trial_basis_dot)) .or. &
            any(.not. ieee_is_finite(surface_weights_dot)) .or. &
            any(.not. ieee_is_finite(coefficients_dot)) .or. &
            any(.not. ieee_is_finite(target_current_dot))) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "surface-current trace JVP received incompatible increments")
            return
        end if

        do quadrature = 1, quadrature_count
            current = 0.0_dp
            current_dot = 0.0_dp
            do trial_dof = 1, trial_count
                current = current + trial_basis(quadrature, trial_dof, :)* &
                    coefficients(trial_dof)
                current_dot = current_dot + &
                    trial_basis_dot(quadrature, trial_dof, :)* &
                    coefficients(trial_dof) + &
                    trial_basis(quadrature, trial_dof, :)* &
                    coefficients_dot(trial_dof)
            end do
            state = current - target_current(quadrature, :)
            state_dot = current_dot - target_current_dot(quadrature, :)
            do test_dof = 1, test_count
                residual_dot(test_dof) = residual_dot(test_dof) + &
                    surface_weights_dot(quadrature)*dot_product( &
                    test_basis(quadrature, test_dof, :), state) + &
                    surface_weights(quadrature)*dot_product( &
                    test_basis_dot(quadrature, test_dof, :), state) + &
                    surface_weights(quadrature)*dot_product( &
                    test_basis(quadrature, test_dof, :), state_dot)
            end do
        end do
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine assemble_surface_current_trace_residual_jvp

    subroutine assemble_surface_current_trace_residual_vjp( &
            test_basis, trial_basis, surface_weights, coefficients, &
            target_current, residual_bar, test_basis_bar, trial_basis_bar, &
            surface_weights_bar, coefficients_bar, target_current_bar, status)
        !! Apply the real transpose of the trace residual.
        real(dp), intent(in) :: test_basis(:, :, :), trial_basis(:, :, :)
        real(dp), intent(in) :: surface_weights(:), coefficients(:)
        real(dp), intent(in) :: target_current(:, :), residual_bar(:)
        real(dp), intent(out) :: test_basis_bar(:, :, :)
        real(dp), intent(out) :: trial_basis_bar(:, :, :)
        real(dp), intent(out) :: surface_weights_bar(:), coefficients_bar(:)
        real(dp), intent(out) :: target_current_bar(:, :)
        type(fortsparse_status_t), intent(out) :: status

        integer :: quadrature_count, test_count, trial_count
        integer :: quadrature, test_dof, trial_dof
        real(dp) :: current(3), state(3), state_bar(3)

        test_basis_bar = 0.0_dp
        trial_basis_bar = 0.0_dp
        surface_weights_bar = 0.0_dp
        coefficients_bar = 0.0_dp
        target_current_bar = 0.0_dp
        call validate_trace_residual_inputs( &
            test_basis, trial_basis, surface_weights, coefficients, &
            target_current, quadrature_count, test_count, trial_count, status)
        if (status%code /= FORTSPARSE_OK) return
        if (size(residual_bar) /= test_count .or. &
            size(test_basis_bar, 1) /= quadrature_count .or. &
            size(test_basis_bar, 2) /= test_count .or. &
            size(test_basis_bar, 3) /= 3 .or. &
            size(trial_basis_bar, 1) /= quadrature_count .or. &
            size(trial_basis_bar, 2) /= trial_count .or. &
            size(trial_basis_bar, 3) /= 3 .or. &
            size(surface_weights_bar) /= quadrature_count .or. &
            size(coefficients_bar) /= trial_count .or. &
            size(target_current_bar, 1) /= quadrature_count .or. &
            size(target_current_bar, 2) /= 3 .or. &
            any(.not. ieee_is_finite(residual_bar))) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "surface-current trace VJP received incompatible cotangents")
            return
        end if

        do quadrature = 1, quadrature_count
            current = 0.0_dp
            do trial_dof = 1, trial_count
                current = current + trial_basis(quadrature, trial_dof, :)* &
                    coefficients(trial_dof)
            end do
            state = current - target_current(quadrature, :)
            state_bar = 0.0_dp
            do test_dof = 1, test_count
                test_basis_bar(quadrature, test_dof, :) = &
                    test_basis_bar(quadrature, test_dof, :) + &
                    surface_weights(quadrature)*residual_bar(test_dof)*state
                surface_weights_bar(quadrature) = &
                    surface_weights_bar(quadrature) + residual_bar(test_dof)* &
                    dot_product(test_basis(quadrature, test_dof, :), state)
                state_bar = state_bar + residual_bar(test_dof)* &
                    test_basis(quadrature, test_dof, :)
            end do
            target_current_bar(quadrature, :) = &
                target_current_bar(quadrature, :) - &
                surface_weights(quadrature)*state_bar
            do trial_dof = 1, trial_count
                trial_basis_bar(quadrature, trial_dof, :) = &
                    trial_basis_bar(quadrature, trial_dof, :) + &
                    surface_weights(quadrature)*coefficients(trial_dof)*state_bar
                coefficients_bar(trial_dof) = coefficients_bar(trial_dof) + &
                    surface_weights(quadrature)*dot_product( &
                    trial_basis(quadrature, trial_dof, :), state_bar)
            end do
        end do
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine assemble_surface_current_trace_residual_vjp

    subroutine validate_trace_residual_inputs( &
            test_basis, trial_basis, surface_weights, coefficients, &
            target_current, quadrature_count, test_count, trial_count, status)
        real(dp), intent(in) :: test_basis(:, :, :), trial_basis(:, :, :)
        real(dp), intent(in) :: surface_weights(:), coefficients(:)
        real(dp), intent(in) :: target_current(:, :)
        integer, intent(out) :: quadrature_count, test_count, trial_count
        type(fortsparse_status_t), intent(out) :: status

        quadrature_count = size(test_basis, 1)
        test_count = size(test_basis, 2)
        trial_count = size(trial_basis, 2)
        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            "surface-current trace residual received incompatible arrays")
        if (quadrature_count < 1 .or. test_count < 1 .or. trial_count < 1) return
        if (size(test_basis, 3) /= 3 .or. &
            size(trial_basis, 1) /= quadrature_count .or. &
            size(trial_basis, 3) /= 3 .or. &
            size(surface_weights) /= quadrature_count .or. &
            size(coefficients) /= trial_count .or. &
            size(target_current, 1) /= quadrature_count .or. &
            size(target_current, 2) /= 3) return
        if (any(.not. ieee_is_finite(test_basis)) .or. &
            any(.not. ieee_is_finite(trial_basis)) .or. &
            any(.not. ieee_is_finite(surface_weights)) .or. &
            any(.not. ieee_is_finite(coefficients)) .or. &
            any(.not. ieee_is_finite(target_current)) .or. &
            any(surface_weights <= 0.0_dp)) return
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine validate_trace_residual_inputs

end module fortfem_surface_current_trace_residual
