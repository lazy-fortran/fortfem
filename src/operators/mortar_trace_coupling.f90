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

        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            "mortar trace coupling received incompatible arrays")
        matrix = 0.0_dp
        quadrature_count = size(test_trace, 1)
        test_count = size(test_trace, 2)
        trial_count = size(trial_trace, 2)
        if (quadrature_count < 1 .or. test_count < 1 .or. trial_count < 1) return
        if (size(trial_trace, 1) /= quadrature_count .or. &
            size(surface_weights) /= quadrature_count .or. &
            size(matrix, 1) /= test_count .or. &
            size(matrix, 2) /= trial_count) return
        if (any(.not. ieee_is_finite(test_trace)) .or. &
            any(.not. ieee_is_finite(trial_trace)) .or. &
            any(.not. ieee_is_finite(surface_weights))) return
        if (any(surface_weights <= 0.0_dp)) return

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

end module fortfem_mortar_trace_coupling
