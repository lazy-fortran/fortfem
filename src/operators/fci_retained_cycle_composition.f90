module fortfem_fci_retained_cycle_composition
    !! Compose a retained coarse action with the neutral FCI split ledger.
    !!
    !! The retained block matrices and their transpose factors are owned by a
    !! `retained_field_split_t`.  This layer applies those factors without
    !! assembling a global preconditioner, then delegates all additive action
    !! and residual-work bookkeeping to the coupled FCI ledger.
    use fortfem_fci_coupled_field_split_ledger, only: &
        evaluate_fci_coupled_field_split_ledger, &
        evaluate_fci_coupled_field_split_ledger_jvp, &
        evaluate_fci_coupled_field_split_ledger_vjp
    use fortfem_kinds, only: dp
    use fortfem_retained_field_split, only: &
        apply_retained_field_split, apply_retained_field_split_jvp, &
        apply_retained_field_split_vjp, retained_field_split_t
    use fortsparse, only: csc_t, fortsparse_status_t, status_set, &
        FORTSPARSE_INVALID_MATRIX, FORTSPARSE_OK
    implicit none
    private

    public :: evaluate_fci_retained_cycle_composition
    public :: evaluate_fci_retained_cycle_composition_jvp
    public :: evaluate_fci_retained_cycle_composition_vjp

contains

    subroutine evaluate_fci_retained_cycle_composition( &
            split, residual, parallel_action, plane_action, weights, correction, &
            work_ledger, total_work, status)
        type(retained_field_split_t), intent(inout) :: split
        real(dp), intent(in) :: residual(:), parallel_action(:), plane_action(:)
        real(dp), intent(in) :: weights(:)
        real(dp), intent(out) :: correction(:), work_ledger(:), total_work
        type(fortsparse_status_t), intent(out) :: status
        real(dp), allocatable :: retained_actions(:, :)
        integer :: solve_status

        correction = 0.0_dp
        work_ledger = 0.0_dp
        total_work = 0.0_dp
        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            "FCI retained cycle received incompatible inputs")
        allocate(retained_actions(size(residual), 1))
        call apply_retained_field_split( &
            split, residual, retained_actions(:, 1), solve_status)
        if (solve_status /= FORTSPARSE_OK) return
        call evaluate_fci_coupled_field_split_ledger( &
            residual, parallel_action, plane_action, retained_actions, weights, &
            correction, work_ledger, total_work, status)
    end subroutine evaluate_fci_retained_cycle_composition

    subroutine evaluate_fci_retained_cycle_composition_jvp( &
            split, matrices_dot, residual, parallel_action, plane_action, weights, &
            residual_dot, parallel_action_dot, plane_action_dot, weights_dot, &
            correction_dot, work_ledger_dot, total_work_dot, status)
        type(retained_field_split_t), intent(inout) :: split
        type(csc_t), intent(in) :: matrices_dot(:)
        real(dp), intent(in) :: residual(:), parallel_action(:), plane_action(:)
        real(dp), intent(in) :: weights(:), residual_dot(:)
        real(dp), intent(in) :: parallel_action_dot(:), plane_action_dot(:)
        real(dp), intent(in) :: weights_dot(:)
        real(dp), intent(out) :: correction_dot(:), work_ledger_dot(:)
        real(dp), intent(out) :: total_work_dot
        type(fortsparse_status_t), intent(out) :: status
        real(dp), allocatable :: retained_actions(:, :)
        real(dp), allocatable :: retained_actions_dot(:, :)
        integer :: solve_status

        correction_dot = 0.0_dp
        work_ledger_dot = 0.0_dp
        total_work_dot = 0.0_dp
        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            "FCI retained-cycle JVP received incompatible inputs")
        allocate(retained_actions(size(residual), 1))
        allocate(retained_actions_dot(size(residual), 1))
        call apply_retained_field_split( &
            split, residual, retained_actions(:, 1), solve_status)
        if (solve_status /= FORTSPARSE_OK) return
        call apply_retained_field_split_jvp( &
            split, matrices_dot, retained_actions(:, 1), residual_dot, &
            retained_actions_dot(:, 1), solve_status)
        if (solve_status /= FORTSPARSE_OK) return
        call evaluate_fci_coupled_field_split_ledger_jvp( &
            residual, parallel_action, plane_action, retained_actions, weights, &
            residual_dot, parallel_action_dot, plane_action_dot, &
            retained_actions_dot, weights_dot, correction_dot, work_ledger_dot, &
            total_work_dot, status)
    end subroutine evaluate_fci_retained_cycle_composition_jvp

    subroutine evaluate_fci_retained_cycle_composition_vjp( &
            split, residual, parallel_action, plane_action, weights, &
            correction_bar, work_ledger_bar, total_work_bar, residual_bar, &
            parallel_action_bar, plane_action_bar, matrices_bar, weights_bar, &
            status)
        type(retained_field_split_t), intent(inout) :: split
        real(dp), intent(in) :: residual(:), parallel_action(:), plane_action(:)
        real(dp), intent(in) :: weights(:), correction_bar(:), work_ledger_bar(:)
        real(dp), intent(in) :: total_work_bar
        real(dp), intent(out) :: residual_bar(:), parallel_action_bar(:)
        real(dp), intent(out) :: plane_action_bar(:), weights_bar(:)
        type(csc_t), intent(out) :: matrices_bar(:)
        type(fortsparse_status_t), intent(out) :: status
        real(dp), allocatable :: direct_residual_bar(:), retained_actions(:, :)
        real(dp), allocatable :: retained_actions_bar(:, :), retained_rhs_bar(:)
        integer :: solve_status

        residual_bar = 0.0_dp
        parallel_action_bar = 0.0_dp
        plane_action_bar = 0.0_dp
        weights_bar = 0.0_dp
        call zero_matrix_bars(split, matrices_bar)
        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            "FCI retained-cycle VJP received incompatible inputs")
        allocate(direct_residual_bar(size(residual)))
        allocate(retained_actions(size(residual), 1))
        allocate(retained_actions_bar(size(residual), 1))
        allocate(retained_rhs_bar(size(residual)))
        call apply_retained_field_split( &
            split, residual, retained_actions(:, 1), solve_status)
        if (solve_status /= FORTSPARSE_OK) return
        call evaluate_fci_coupled_field_split_ledger_vjp( &
            residual, parallel_action, plane_action, retained_actions, weights, &
            correction_bar, work_ledger_bar, total_work_bar, direct_residual_bar, &
            parallel_action_bar, plane_action_bar, retained_actions_bar, &
            weights_bar, status)
        if (status%code /= FORTSPARSE_OK) return
        call apply_retained_field_split_vjp( &
            split, retained_actions(:, 1), retained_actions_bar(:, 1), &
            retained_rhs_bar, matrices_bar, solve_status)
        if (solve_status /= FORTSPARSE_OK) then
            residual_bar = 0.0_dp
            parallel_action_bar = 0.0_dp
            plane_action_bar = 0.0_dp
            weights_bar = 0.0_dp
            call zero_matrix_bars(split, matrices_bar)
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "FCI retained-cycle transpose solve failed")
            return
        end if
        residual_bar = direct_residual_bar + retained_rhs_bar
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine evaluate_fci_retained_cycle_composition_vjp

    subroutine zero_matrix_bars(split, matrices_bar)
        type(retained_field_split_t), intent(in) :: split
        type(csc_t), intent(out) :: matrices_bar(:)
        integer :: field

        if (.not. allocated(split%matrices)) return
        if (size(matrices_bar) /= size(split%matrices)) return
        matrices_bar = split%matrices
        do field = 1, size(matrices_bar)
            matrices_bar(field)%val = 0.0_dp
        end do
    end subroutine zero_matrix_bars

end module fortfem_fci_retained_cycle_composition
