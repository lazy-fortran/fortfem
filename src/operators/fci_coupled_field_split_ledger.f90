module fortfem_fci_coupled_field_split_ledger
    !! Neutral additive composition of caller-owned FCI split actions.
    !!
    !! The first two actions conventionally come from a parallel Jacobi solve
    !! and a ragged perpendicular-plane cycle.  Additional columns may come
    !! from retained coarse, plane, or coupled-field factors.  This module does
    !! not own those solvers: it combines their actions, records the weighted
    !! residual work of every contribution, and exposes exact derivatives.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortfem_kinds, only: dp
    use fortsparse, only: fortsparse_status_t, status_set, &
        FORTSPARSE_INVALID_MATRIX, FORTSPARSE_OK
    implicit none
    private

    public :: evaluate_fci_coupled_field_split_ledger
    public :: evaluate_fci_coupled_field_split_ledger_jvp
    public :: evaluate_fci_coupled_field_split_ledger_vjp

contains

    subroutine evaluate_fci_coupled_field_split_ledger( &
            residual, parallel_action, plane_action, retained_actions, weights, &
            correction, work_ledger, total_work, status)
        real(dp), intent(in) :: residual(:), parallel_action(:), plane_action(:)
        real(dp), intent(in) :: retained_actions(:, :), weights(:)
        real(dp), intent(out) :: correction(:), work_ledger(:), total_work
        type(fortsparse_status_t), intent(out) :: status
        integer :: block
        real(dp), allocatable :: unweighted_work(:)

        correction = 0.0_dp
        work_ledger = 0.0_dp
        total_work = 0.0_dp
        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            "FCI coupled field split received incompatible actions")
        if (.not. valid_values( &
            residual, parallel_action, plane_action, retained_actions, weights, &
            correction, work_ledger)) return

        allocate(unweighted_work(size(weights)))
        unweighted_work(1) = dot_product(residual, parallel_action)
        unweighted_work(2) = dot_product(residual, plane_action)
        do block = 1, size(retained_actions, 2)
            unweighted_work(block + 2) = &
                dot_product(residual, retained_actions(:, block))
        end do
        if (any(unweighted_work < 0.0_dp)) return

        correction = weights(1)*parallel_action + weights(2)*plane_action
        do block = 1, size(retained_actions, 2)
            correction = correction + weights(block + 2)*retained_actions(:, block)
        end do
        work_ledger = weights*unweighted_work
        total_work = sum(work_ledger)
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine evaluate_fci_coupled_field_split_ledger

    subroutine evaluate_fci_coupled_field_split_ledger_jvp( &
            residual, parallel_action, plane_action, retained_actions, weights, &
            residual_dot, parallel_action_dot, plane_action_dot, &
            retained_actions_dot, weights_dot, correction_dot, work_ledger_dot, &
            total_work_dot, status)
        real(dp), intent(in) :: residual(:), parallel_action(:), plane_action(:)
        real(dp), intent(in) :: retained_actions(:, :), weights(:)
        real(dp), intent(in) :: residual_dot(:), parallel_action_dot(:)
        real(dp), intent(in) :: plane_action_dot(:), retained_actions_dot(:, :)
        real(dp), intent(in) :: weights_dot(:)
        real(dp), intent(out) :: correction_dot(:), work_ledger_dot(:)
        real(dp), intent(out) :: total_work_dot
        type(fortsparse_status_t), intent(out) :: status
        integer :: block
        real(dp) :: work, work_dot

        correction_dot = 0.0_dp
        work_ledger_dot = 0.0_dp
        total_work_dot = 0.0_dp
        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            "FCI coupled field-split JVP received incompatible actions")
        if (.not. valid_values( &
            residual, parallel_action, plane_action, retained_actions, weights, &
            correction_dot, work_ledger_dot)) return
        if (.not. valid_directions( &
            residual_dot, parallel_action_dot, plane_action_dot, &
            retained_actions_dot, weights_dot, size(residual), &
            size(retained_actions, 2))) return
        if (.not. nonnegative_component_work( &
            residual, parallel_action, plane_action, retained_actions)) return

        correction_dot = weights_dot(1)*parallel_action + &
            weights(1)*parallel_action_dot + weights_dot(2)*plane_action + &
            weights(2)*plane_action_dot
        work = dot_product(residual, parallel_action)
        work_dot = dot_product(residual_dot, parallel_action) + &
            dot_product(residual, parallel_action_dot)
        work_ledger_dot(1) = weights_dot(1)*work + weights(1)*work_dot
        work = dot_product(residual, plane_action)
        work_dot = dot_product(residual_dot, plane_action) + &
            dot_product(residual, plane_action_dot)
        work_ledger_dot(2) = weights_dot(2)*work + weights(2)*work_dot
        do block = 1, size(retained_actions, 2)
            correction_dot = correction_dot + &
                weights_dot(block + 2)*retained_actions(:, block) + &
                weights(block + 2)*retained_actions_dot(:, block)
            work = dot_product(residual, retained_actions(:, block))
            work_dot = dot_product(residual_dot, retained_actions(:, block)) + &
                dot_product(residual, retained_actions_dot(:, block))
            work_ledger_dot(block + 2) = weights_dot(block + 2)*work + &
                weights(block + 2)*work_dot
        end do
        total_work_dot = sum(work_ledger_dot)
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine evaluate_fci_coupled_field_split_ledger_jvp

    subroutine evaluate_fci_coupled_field_split_ledger_vjp( &
            residual, parallel_action, plane_action, retained_actions, weights, &
            correction_bar, work_ledger_bar, total_work_bar, residual_bar, &
            parallel_action_bar, plane_action_bar, retained_actions_bar, &
            weights_bar, status)
        real(dp), intent(in) :: residual(:), parallel_action(:), plane_action(:)
        real(dp), intent(in) :: retained_actions(:, :), weights(:)
        real(dp), intent(in) :: correction_bar(:), work_ledger_bar(:)
        real(dp), intent(in) :: total_work_bar
        real(dp), intent(out) :: residual_bar(:), parallel_action_bar(:)
        real(dp), intent(out) :: plane_action_bar(:), retained_actions_bar(:, :)
        real(dp), intent(out) :: weights_bar(:)
        type(fortsparse_status_t), intent(out) :: status
        integer :: block
        real(dp) :: effective_work_bar, work

        residual_bar = 0.0_dp
        parallel_action_bar = 0.0_dp
        plane_action_bar = 0.0_dp
        retained_actions_bar = 0.0_dp
        weights_bar = 0.0_dp
        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            "FCI coupled field-split VJP received incompatible actions")
        if (.not. valid_values( &
            residual, parallel_action, plane_action, retained_actions, weights, &
            residual_bar, weights_bar)) return
        if (.not. valid_cotangents( &
            correction_bar, work_ledger_bar, total_work_bar, residual_bar, &
            parallel_action_bar, plane_action_bar, retained_actions_bar, &
            weights_bar, size(residual), size(retained_actions, 2))) return
        if (.not. nonnegative_component_work( &
            residual, parallel_action, plane_action, retained_actions)) return

        effective_work_bar = work_ledger_bar(1) + total_work_bar
        work = dot_product(residual, parallel_action)
        residual_bar = effective_work_bar*weights(1)*parallel_action
        parallel_action_bar = weights(1)*correction_bar + &
            effective_work_bar*weights(1)*residual
        weights_bar(1) = dot_product(correction_bar, parallel_action) + &
            effective_work_bar*work

        effective_work_bar = work_ledger_bar(2) + total_work_bar
        work = dot_product(residual, plane_action)
        residual_bar = residual_bar + effective_work_bar*weights(2)*plane_action
        plane_action_bar = weights(2)*correction_bar + &
            effective_work_bar*weights(2)*residual
        weights_bar(2) = dot_product(correction_bar, plane_action) + &
            effective_work_bar*work

        do block = 1, size(retained_actions, 2)
            effective_work_bar = work_ledger_bar(block + 2) + total_work_bar
            work = dot_product(residual, retained_actions(:, block))
            residual_bar = residual_bar + effective_work_bar*weights(block + 2)* &
                retained_actions(:, block)
            retained_actions_bar(:, block) = weights(block + 2)*correction_bar + &
                effective_work_bar*weights(block + 2)*residual
            weights_bar(block + 2) = &
                dot_product(correction_bar, retained_actions(:, block)) + &
                effective_work_bar*work
        end do
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine evaluate_fci_coupled_field_split_ledger_vjp

    logical function valid_values( &
            residual, parallel_action, plane_action, retained_actions, weights, &
            correction, work_ledger) result(valid)
        real(dp), intent(in) :: residual(:), parallel_action(:), plane_action(:)
        real(dp), intent(in) :: retained_actions(:, :), weights(:)
        real(dp), intent(in) :: correction(:), work_ledger(:)
        integer :: vector_size, retained_count

        valid = .false.
        vector_size = size(residual)
        retained_count = size(retained_actions, 2)
        if (vector_size < 1 .or. retained_count < 1) return
        if (size(parallel_action) /= vector_size .or. &
            size(plane_action) /= vector_size) return
        if (size(retained_actions, 1) /= vector_size) return
        if (size(weights) /= retained_count + 2) return
        if (size(correction) /= vector_size .or. &
            size(work_ledger) /= retained_count + 2) return
        if (any(.not. ieee_is_finite(residual)) .or. &
            any(.not. ieee_is_finite(parallel_action)) .or. &
            any(.not. ieee_is_finite(plane_action)) .or. &
            any(.not. ieee_is_finite(retained_actions)) .or. &
            any(.not. ieee_is_finite(weights))) return
        if (any(weights <= 0.0_dp)) return
        valid = .true.
    end function valid_values

    logical function nonnegative_component_work( &
            residual, parallel_action, plane_action, retained_actions) result(valid)
        real(dp), intent(in) :: residual(:), parallel_action(:), plane_action(:)
        real(dp), intent(in) :: retained_actions(:, :)
        integer :: block

        valid = dot_product(residual, parallel_action) >= 0.0_dp
        if (.not. valid) return
        valid = dot_product(residual, plane_action) >= 0.0_dp
        if (.not. valid) return
        do block = 1, size(retained_actions, 2)
            valid = dot_product(residual, retained_actions(:, block)) >= 0.0_dp
            if (.not. valid) return
        end do
    end function nonnegative_component_work

    logical function valid_directions( &
            residual_dot, parallel_action_dot, plane_action_dot, &
            retained_actions_dot, weights_dot, vector_size, retained_count) &
            result(valid)
        real(dp), intent(in) :: residual_dot(:), parallel_action_dot(:)
        real(dp), intent(in) :: plane_action_dot(:), retained_actions_dot(:, :)
        real(dp), intent(in) :: weights_dot(:)
        integer, intent(in) :: vector_size, retained_count

        valid = .false.
        if (size(residual_dot) /= vector_size .or. &
            size(parallel_action_dot) /= vector_size .or. &
            size(plane_action_dot) /= vector_size) return
        if (size(retained_actions_dot, 1) /= vector_size .or. &
            size(retained_actions_dot, 2) /= retained_count) return
        if (size(weights_dot) /= retained_count + 2) return
        if (any(.not. ieee_is_finite(residual_dot)) .or. &
            any(.not. ieee_is_finite(parallel_action_dot)) .or. &
            any(.not. ieee_is_finite(plane_action_dot)) .or. &
            any(.not. ieee_is_finite(retained_actions_dot)) .or. &
            any(.not. ieee_is_finite(weights_dot))) return
        valid = .true.
    end function valid_directions

    logical function valid_cotangents( &
            correction_bar, work_ledger_bar, total_work_bar, residual_bar, &
            parallel_action_bar, plane_action_bar, retained_actions_bar, &
            weights_bar, vector_size, retained_count) result(valid)
        real(dp), intent(in) :: correction_bar(:), work_ledger_bar(:)
        real(dp), intent(in) :: total_work_bar
        real(dp), intent(in) :: residual_bar(:), parallel_action_bar(:)
        real(dp), intent(in) :: plane_action_bar(:), retained_actions_bar(:, :)
        real(dp), intent(in) :: weights_bar(:)
        integer, intent(in) :: vector_size, retained_count

        valid = .false.
        if (size(correction_bar) /= vector_size .or. &
            size(work_ledger_bar) /= retained_count + 2) return
        if (size(residual_bar) /= vector_size .or. &
            size(parallel_action_bar) /= vector_size .or. &
            size(plane_action_bar) /= vector_size) return
        if (size(retained_actions_bar, 1) /= vector_size .or. &
            size(retained_actions_bar, 2) /= retained_count) return
        if (size(weights_bar) /= retained_count + 2) return
        if (any(.not. ieee_is_finite(correction_bar)) .or. &
            any(.not. ieee_is_finite(work_ledger_bar)) .or. &
            .not. ieee_is_finite(total_work_bar)) return
        valid = .true.
    end function valid_cotangents

end module fortfem_fci_coupled_field_split_ledger
