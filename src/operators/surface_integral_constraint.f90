module fortfem_surface_integral_constraint
    !! Neutral weighted fixed-topology scalar surface-integral constraint.
    !!
    !! For caller-owned samples s(q), positive surface weights w(q), and a
    !! target t, the constraint value is
    !!
    !!   c = sum_q w(q) s(q) - t.
    !!
    !! Sample meaning, units, quadrature, topology, and target selection remain
    !! external.  The same contract can therefore serve area, flux, volume,
    !! loop, or shape ledgers without selecting a physical model.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortfem_kinds, only: dp
    use fortsparse, only: FORTSPARSE_INVALID_MATRIX, FORTSPARSE_OK, &
        fortsparse_status_t, status_set
    implicit none
    private

    public :: evaluate_surface_integral_constraint
    public :: evaluate_surface_integral_constraint_jvp
    public :: evaluate_surface_integral_constraint_vjp

contains

    subroutine evaluate_surface_integral_constraint( &
            samples, weights, target, constraint, status)
        real(dp), intent(in) :: samples(:), weights(:), target
        real(dp), intent(out) :: constraint
        type(fortsparse_status_t), intent(out) :: status

        constraint = 0.0_dp
        if (.not. validate_inputs(samples, weights, target, status)) return
        constraint = dot_product(samples, weights) - target
        if (.not. ieee_is_finite(constraint)) then
            constraint = 0.0_dp
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "surface-integral constraint is non-finite")
            return
        end if
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine evaluate_surface_integral_constraint

    subroutine evaluate_surface_integral_constraint_jvp( &
            samples, weights, target, samples_dot, weights_dot, target_dot, &
            constraint_dot, status)
        real(dp), intent(in) :: samples(:), weights(:), target
        real(dp), intent(in) :: samples_dot(:), weights_dot(:), target_dot
        real(dp), intent(out) :: constraint_dot
        type(fortsparse_status_t), intent(out) :: status

        constraint_dot = 0.0_dp
        if (.not. validate_inputs(samples, weights, target, status)) return
        if (size(samples_dot) /= size(samples)) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "surface-integral JVP has an incompatible sample tangent")
            return
        end if
        if (size(weights_dot) /= size(weights)) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "surface-integral JVP has an incompatible weight tangent")
            return
        end if
        if (.not. ieee_is_finite(target_dot)) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "surface-integral JVP target tangent is non-finite")
            return
        end if
        if (.not. all(ieee_is_finite(samples_dot))) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "surface-integral JVP sample tangent is non-finite")
            return
        end if
        if (.not. all(ieee_is_finite(weights_dot))) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "surface-integral JVP weight tangent is non-finite")
            return
        end if
        constraint_dot = dot_product(samples_dot, weights) + &
            dot_product(samples, weights_dot) - target_dot
        if (.not. ieee_is_finite(constraint_dot)) then
            constraint_dot = 0.0_dp
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "surface-integral constraint JVP is non-finite")
            return
        end if
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine evaluate_surface_integral_constraint_jvp

    subroutine evaluate_surface_integral_constraint_vjp( &
            samples, weights, target, constraint_bar, samples_bar, weights_bar, &
            target_bar, status)
        real(dp), intent(in) :: samples(:), weights(:), target, constraint_bar
        real(dp), intent(out) :: samples_bar(:), weights_bar(:), target_bar
        type(fortsparse_status_t), intent(out) :: status

        samples_bar = 0.0_dp
        weights_bar = 0.0_dp
        target_bar = 0.0_dp
        if (.not. validate_inputs(samples, weights, target, status)) return
        if (size(samples_bar) /= size(samples)) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "surface-integral VJP has an incompatible sample cotangent")
            return
        end if
        if (size(weights_bar) /= size(weights)) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "surface-integral VJP has an incompatible weight cotangent")
            return
        end if
        if (.not. ieee_is_finite(constraint_bar)) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "surface-integral VJP constraint cotangent is non-finite")
            return
        end if
        samples_bar = constraint_bar*weights
        weights_bar = constraint_bar*samples
        target_bar = -constraint_bar
        if (.not. all(ieee_is_finite(samples_bar))) then
            samples_bar = 0.0_dp
            weights_bar = 0.0_dp
            target_bar = 0.0_dp
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "surface-integral sample VJP is non-finite")
            return
        end if
        if (.not. all(ieee_is_finite(weights_bar))) then
            samples_bar = 0.0_dp
            weights_bar = 0.0_dp
            target_bar = 0.0_dp
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "surface-integral weight VJP is non-finite")
            return
        end if
        if (.not. ieee_is_finite(target_bar)) then
            samples_bar = 0.0_dp
            weights_bar = 0.0_dp
            target_bar = 0.0_dp
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "surface-integral target VJP is non-finite")
            return
        end if
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine evaluate_surface_integral_constraint_vjp

    logical function validate_inputs(samples, weights, target, status) result(valid)
        real(dp), intent(in) :: samples(:), weights(:), target
        type(fortsparse_status_t), intent(out) :: status

        valid = .false.
        if (size(samples) < 1) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "surface-integral constraint requires samples")
            return
        end if
        if (size(weights) /= size(samples)) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "surface-integral constraint has incompatible topology")
            return
        end if
        if (.not. ieee_is_finite(target)) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "surface-integral constraint target is non-finite")
            return
        end if
        if (.not. all(ieee_is_finite(samples))) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "surface-integral samples are non-finite")
            return
        end if
        if (.not. all(ieee_is_finite(weights))) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "surface-integral weights are non-finite")
            return
        end if
        if (.not. all(weights > 0.0_dp)) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "surface-integral weights must be positive")
            return
        end if
        valid = .true.
        call status_set(status, FORTSPARSE_OK, "")
    end function validate_inputs

end module fortfem_surface_integral_constraint
