module fortfem_mortar_trace_coupling
    !! Cross-mass coupling for independently discretized interface traces.
    !!
    !! M_ij = integral T_test_i T_trial_j dS.  The two trace spaces may have
    !! different numbers of degrees of freedom; only their quadrature rows
    !! must agree.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortfem_kinds, only: dp
    use fortsparse, only: fortsparse_status_t, status_set, &
        FORTSPARSE_INVALID_MATRIX, FORTSPARSE_OK
    implicit none
    private

    public :: assemble_mortar_trace_coupling
    public :: assemble_mortar_trace_coupling_jvp
    public :: assemble_mortar_trace_coupling_vjp

contains

    subroutine assemble_mortar_trace_coupling( &
            test_trace, trial_trace, surface_weights, matrix, status)
        !! Assemble the weighted test/trial trace cross-mass matrix.
        real(dp), intent(in) :: test_trace(:, :), trial_trace(:, :)
        real(dp), intent(in) :: surface_weights(:)
        real(dp), intent(out) :: matrix(:, :)
        type(fortsparse_status_t), intent(out) :: status

        integer :: quadrature_count, test_count, trial_count
        integer :: quadrature, test_dof, trial_dof
        real(dp) :: scale

        matrix = 0.0_dp
        call validate_mortar_inputs( &
            test_trace, trial_trace, surface_weights, matrix, status)
        if (status%code /= FORTSPARSE_OK) return
        quadrature_count = size(test_trace, 1)
        test_count = size(test_trace, 2)
        trial_count = size(trial_trace, 2)

        do quadrature = 1, quadrature_count
            scale = surface_weights(quadrature)
            do test_dof = 1, test_count
                do trial_dof = 1, trial_count
                    matrix(test_dof, trial_dof) = &
                        matrix(test_dof, trial_dof) + scale* &
                        test_trace(quadrature, test_dof)* &
                        trial_trace(quadrature, trial_dof)
                end do
            end do
        end do
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine assemble_mortar_trace_coupling

    subroutine assemble_mortar_trace_coupling_jvp( &
            test_trace, trial_trace, surface_weights, test_trace_dot, &
            trial_trace_dot, surface_weights_dot, matrix_dot, status)
        !! Apply the product-rule JVP of the weighted cross-mass matrix.
        real(dp), intent(in) :: test_trace(:, :), trial_trace(:, :)
        real(dp), intent(in) :: surface_weights(:)
        real(dp), intent(in) :: test_trace_dot(:, :), trial_trace_dot(:, :)
        real(dp), intent(in) :: surface_weights_dot(:)
        real(dp), intent(out) :: matrix_dot(:, :)
        type(fortsparse_status_t), intent(out) :: status

        integer :: quadrature_count, test_count, trial_count
        integer :: quadrature, test_dof, trial_dof

        matrix_dot = 0.0_dp
        call validate_mortar_inputs( &
            test_trace, trial_trace, surface_weights, matrix_dot, status)
        if (status%code /= FORTSPARSE_OK) return
        quadrature_count = size(test_trace, 1)
        test_count = size(test_trace, 2)
        trial_count = size(trial_trace, 2)
        if (size(test_trace_dot, 1) /= quadrature_count .or. &
            size(test_trace_dot, 2) /= test_count .or. &
            size(trial_trace_dot, 1) /= quadrature_count .or. &
            size(trial_trace_dot, 2) /= trial_count .or. &
            size(surface_weights_dot) /= quadrature_count .or. &
            any(.not. ieee_is_finite(test_trace_dot)) .or. &
            any(.not. ieee_is_finite(trial_trace_dot)) .or. &
            any(.not. ieee_is_finite(surface_weights_dot))) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "mortar trace coupling JVP has incompatible increments")
            return
        end if
        do quadrature = 1, quadrature_count
            do test_dof = 1, test_count
                do trial_dof = 1, trial_count
                    matrix_dot(test_dof, trial_dof) = &
                        matrix_dot(test_dof, trial_dof) + &
                        surface_weights_dot(quadrature)* &
                        test_trace(quadrature, test_dof)* &
                        trial_trace(quadrature, trial_dof) + &
                        surface_weights(quadrature)*test_trace_dot( &
                        quadrature, test_dof)*trial_trace(quadrature, trial_dof) + &
                        surface_weights(quadrature)*test_trace(quadrature, test_dof)* &
                        trial_trace_dot(quadrature, trial_dof)
                end do
            end do
        end do
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine assemble_mortar_trace_coupling_jvp

    subroutine assemble_mortar_trace_coupling_vjp( &
            test_trace, trial_trace, surface_weights, matrix_bar, test_trace_bar, &
            trial_trace_bar, surface_weights_bar, status)
        !! Apply the reverse product of the weighted cross-mass matrix.
        real(dp), intent(in) :: test_trace(:, :), trial_trace(:, :)
        real(dp), intent(in) :: surface_weights(:), matrix_bar(:, :)
        real(dp), intent(out) :: test_trace_bar(:, :), trial_trace_bar(:, :)
        real(dp), intent(out) :: surface_weights_bar(:)
        type(fortsparse_status_t), intent(out) :: status

        integer :: quadrature_count, test_count, trial_count
        integer :: quadrature, test_dof, trial_dof

        test_trace_bar = 0.0_dp
        trial_trace_bar = 0.0_dp
        surface_weights_bar = 0.0_dp
        call validate_mortar_inputs( &
            test_trace, trial_trace, surface_weights, matrix_bar, status)
        if (status%code /= FORTSPARSE_OK) return
        quadrature_count = size(test_trace, 1)
        test_count = size(test_trace, 2)
        trial_count = size(trial_trace, 2)
        if (size(test_trace_bar, 1) /= quadrature_count .or. &
            size(test_trace_bar, 2) /= test_count .or. &
            size(trial_trace_bar, 1) /= quadrature_count .or. &
            size(trial_trace_bar, 2) /= trial_count .or. &
            size(surface_weights_bar) /= quadrature_count) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "mortar trace coupling VJP has incompatible cotangents")
            return
        end if
        do quadrature = 1, quadrature_count
            do test_dof = 1, test_count
                do trial_dof = 1, trial_count
                    test_trace_bar(quadrature, test_dof) = &
                        test_trace_bar(quadrature, test_dof) + &
                        surface_weights(quadrature)*matrix_bar(test_dof, trial_dof)* &
                        trial_trace(quadrature, trial_dof)
                    trial_trace_bar(quadrature, trial_dof) = &
                        trial_trace_bar(quadrature, trial_dof) + &
                        surface_weights(quadrature)*matrix_bar(test_dof, trial_dof)* &
                        test_trace(quadrature, test_dof)
                    surface_weights_bar(quadrature) = &
                        surface_weights_bar(quadrature) + &
                        matrix_bar(test_dof, trial_dof)* &
                        test_trace(quadrature, test_dof)* &
                        trial_trace(quadrature, trial_dof)
                end do
            end do
        end do
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine assemble_mortar_trace_coupling_vjp

    subroutine validate_mortar_inputs( &
            test_trace, trial_trace, surface_weights, matrix, status)
        real(dp), intent(in) :: test_trace(:, :), trial_trace(:, :)
        real(dp), intent(in) :: surface_weights(:), matrix(:, :)
        type(fortsparse_status_t), intent(out) :: status
        integer :: quadrature_count, test_count, trial_count

        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            "mortar trace coupling received incompatible arrays")
        quadrature_count = size(test_trace, 1)
        test_count = size(test_trace, 2)
        trial_count = size(trial_trace, 2)
        if (quadrature_count < 1 .or. test_count < 1 .or. trial_count < 1) return
        if (size(trial_trace, 1) /= quadrature_count .or. &
            size(surface_weights) /= quadrature_count .or. &
            size(matrix, 1) /= test_count .or. size(matrix, 2) /= trial_count) return
        if (any(.not. ieee_is_finite(test_trace)) .or. &
            any(.not. ieee_is_finite(trial_trace)) .or. &
            any(.not. ieee_is_finite(surface_weights)) .or. &
            any(.not. ieee_is_finite(matrix))) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "mortar trace coupling received non-finite data")
            return
        end if
        if (any(surface_weights <= 0.0_dp)) return
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine validate_mortar_inputs

end module fortfem_mortar_trace_coupling
