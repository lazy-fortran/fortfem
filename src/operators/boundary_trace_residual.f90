module fortfem_boundary_trace_residual
    !! Real normal/tangential boundary-port residual.
    !!
    !! For caller-owned trace maps N and T, positive work weights, state u,
    !! and supplied targets, the contract is
    !!
    !!   r_n = w_n (N u - g_n),   r_t = w_t (T u - g_t).
    !!
    !! The maps can come from scalar FEM/BEM/DtN/PML, vector traces, IGA, DG,
    !! or wall clients.  No constitutive law, equilibrium convention, or
    !! geometry construction is selected here.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortfem_kinds, only: dp
    use fortsparse, only: FORTSPARSE_INVALID_MATRIX, FORTSPARSE_OK, &
        fortsparse_status_t, status_set
    implicit none
    private

    public :: assemble_boundary_trace_residual
    public :: assemble_boundary_trace_residual_jvp
    public :: assemble_boundary_trace_residual_vjp

contains

    subroutine assemble_boundary_trace_residual( &
            normal_trace, tangential_trace, normal_weights, tangential_weights, &
            state, normal_target, tangential_target, normal_residual, &
            tangential_residual, status)
        real(dp), intent(in) :: normal_trace(:, :), tangential_trace(:, :)
        real(dp), intent(in) :: normal_weights(:), tangential_weights(:)
        real(dp), intent(in) :: state(:), normal_target(:), tangential_target(:)
        real(dp), intent(out) :: normal_residual(:), tangential_residual(:)
        type(fortsparse_status_t), intent(out) :: status

        normal_residual = 0.0_dp
        tangential_residual = 0.0_dp
        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            "boundary trace residual received incompatible arrays")
        if (.not. valid_inputs(normal_trace, tangential_trace, normal_weights, &
            tangential_weights, state, normal_target, tangential_target, &
            normal_residual, tangential_residual)) return
        normal_residual = normal_weights*(matmul(normal_trace, state) - normal_target)
        tangential_residual = tangential_weights*(matmul(tangential_trace, state) - &
            tangential_target)
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine assemble_boundary_trace_residual

    subroutine assemble_boundary_trace_residual_jvp( &
            normal_trace, tangential_trace, normal_weights, tangential_weights, state, &
            normal_target, tangential_target, normal_trace_dot, tangential_trace_dot, &
            normal_weights_dot, tangential_weights_dot, state_dot, normal_target_dot, &
            tangential_target_dot, normal_residual_dot, tangential_residual_dot, status)
        real(dp), intent(in) :: normal_trace(:, :), tangential_trace(:, :)
        real(dp), intent(in) :: normal_weights(:), tangential_weights(:)
        real(dp), intent(in) :: state(:), normal_target(:), tangential_target(:)
        real(dp), intent(in) :: normal_trace_dot(:, :), tangential_trace_dot(:, :)
        real(dp), intent(in) :: normal_weights_dot(:), tangential_weights_dot(:)
        real(dp), intent(in) :: state_dot(:), normal_target_dot(:), tangential_target_dot(:)
        real(dp), intent(out) :: normal_residual_dot(:), tangential_residual_dot(:)
        type(fortsparse_status_t), intent(out) :: status

        normal_residual_dot = 0.0_dp
        tangential_residual_dot = 0.0_dp
        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            "boundary trace residual JVP received incompatible arrays")
        if (.not. valid_inputs(normal_trace, tangential_trace, normal_weights, &
            tangential_weights, state, normal_target, tangential_target, &
            normal_residual_dot, tangential_residual_dot)) return
        if (.not. valid_direction(normal_trace, tangential_trace, state, normal_target, &
            tangential_target, normal_trace_dot, tangential_trace_dot, normal_weights_dot, &
            tangential_weights_dot, state_dot, normal_target_dot, tangential_target_dot)) return
        normal_residual_dot = normal_weights_dot*(matmul(normal_trace, state) - &
            normal_target) + normal_weights*(matmul(normal_trace_dot, state) + &
            matmul(normal_trace, state_dot) - normal_target_dot)
        tangential_residual_dot = tangential_weights_dot*(matmul(tangential_trace, state) - &
            tangential_target) + tangential_weights*(matmul(tangential_trace_dot, state) + &
            matmul(tangential_trace, state_dot) - tangential_target_dot)
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine assemble_boundary_trace_residual_jvp

    subroutine assemble_boundary_trace_residual_vjp( &
            normal_trace, tangential_trace, normal_weights, tangential_weights, state, &
            normal_target, tangential_target, normal_residual_bar, tangential_residual_bar, &
            normal_trace_bar, tangential_trace_bar, normal_weights_bar, &
            tangential_weights_bar, state_bar, normal_target_bar, tangential_target_bar, status)
        real(dp), intent(in) :: normal_trace(:, :), tangential_trace(:, :)
        real(dp), intent(in) :: normal_weights(:), tangential_weights(:)
        real(dp), intent(in) :: state(:), normal_target(:), tangential_target(:)
        real(dp), intent(in) :: normal_residual_bar(:), tangential_residual_bar(:)
        real(dp), intent(out) :: normal_trace_bar(:, :), tangential_trace_bar(:, :)
        real(dp), intent(out) :: normal_weights_bar(:), tangential_weights_bar(:)
        real(dp), intent(out) :: state_bar(:), normal_target_bar(:), tangential_target_bar(:)
        type(fortsparse_status_t), intent(out) :: status
        integer :: row, column
        real(dp) :: normal_raw, tangential_raw
        real(dp) :: normal_weighted_bar, tangential_weighted_bar

        normal_trace_bar = 0.0_dp
        tangential_trace_bar = 0.0_dp
        normal_weights_bar = 0.0_dp
        tangential_weights_bar = 0.0_dp
        state_bar = 0.0_dp
        normal_target_bar = 0.0_dp
        tangential_target_bar = 0.0_dp
        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            "boundary trace residual VJP received incompatible arrays")
        if (.not. valid_inputs(normal_trace, tangential_trace, normal_weights, &
            tangential_weights, state, normal_target, tangential_target, &
            normal_residual_bar, tangential_residual_bar)) return
        if (.not. valid_adjoint_outputs(normal_trace, tangential_trace, state, &
            normal_target, tangential_target, normal_trace_bar, tangential_trace_bar, &
            normal_weights_bar, tangential_weights_bar, state_bar, normal_target_bar, &
            tangential_target_bar)) return
        if (.not. all(ieee_is_finite(normal_residual_bar)) .or. &
            .not. all(ieee_is_finite(tangential_residual_bar))) return

        do row = 1, size(normal_weights)
            normal_raw = sum(normal_trace(row, :)*state) - normal_target(row)
            normal_weighted_bar = normal_weights(row)*normal_residual_bar(row)
            normal_weights_bar(row) = normal_residual_bar(row)*normal_raw
            normal_target_bar(row) = -normal_weighted_bar
            state_bar = state_bar + normal_trace(row, :)*normal_weighted_bar
            do column = 1, size(state)
                normal_trace_bar(row, column) = normal_weighted_bar*state(column)
            end do
        end do
        do row = 1, size(tangential_weights)
            tangential_raw = sum(tangential_trace(row, :)*state) - &
                tangential_target(row)
            tangential_weighted_bar = tangential_weights(row)* &
                tangential_residual_bar(row)
            tangential_weights_bar(row) = tangential_residual_bar(row)*tangential_raw
            tangential_target_bar(row) = -tangential_weighted_bar
            state_bar = state_bar + tangential_trace(row, :)*tangential_weighted_bar
            do column = 1, size(state)
                tangential_trace_bar(row, column) = tangential_weighted_bar*state(column)
            end do
        end do
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine assemble_boundary_trace_residual_vjp

    logical function valid_inputs( &
            normal_trace, tangential_trace, normal_weights, tangential_weights, state, &
            normal_target, tangential_target, normal_residual, tangential_residual) &
            result(valid)
        real(dp), intent(in) :: normal_trace(:, :), tangential_trace(:, :)
        real(dp), intent(in) :: normal_weights(:), tangential_weights(:)
        real(dp), intent(in) :: state(:), normal_target(:), tangential_target(:)
        real(dp), intent(in) :: normal_residual(:), tangential_residual(:)

        valid = .false.
        if (size(normal_trace, 1) < 1 .or. size(tangential_trace, 1) < 1 .or. &
            size(normal_trace, 2) < 1 .or. &
            size(tangential_trace, 2) /= size(normal_trace, 2)) return
        if (size(state) /= size(normal_trace, 2) .or. &
            size(normal_weights) /= size(normal_trace, 1) .or. &
            size(tangential_weights) /= size(tangential_trace, 1) .or. &
            size(normal_target) /= size(normal_weights) .or. &
            size(tangential_target) /= size(tangential_weights) .or. &
            size(normal_residual) /= size(normal_weights) .or. &
            size(tangential_residual) /= size(tangential_weights)) return
        if (any(normal_weights <= 0.0_dp) .or. any(tangential_weights <= 0.0_dp)) return
        if (.not. all(ieee_is_finite(normal_weights)) .or. &
            .not. all(ieee_is_finite(tangential_weights))) return
        if (.not. all(ieee_is_finite(normal_trace)) .or. &
            .not. all(ieee_is_finite(tangential_trace)) .or. &
            .not. all(ieee_is_finite(state)) .or. &
            .not. all(ieee_is_finite(normal_target)) .or. &
            .not. all(ieee_is_finite(tangential_target))) return
        valid = .true.
    end function valid_inputs

    logical function valid_direction( &
            normal_trace, tangential_trace, state, normal_target, tangential_target, &
            normal_trace_dot, tangential_trace_dot, normal_weights_dot, &
            tangential_weights_dot, state_dot, normal_target_dot, tangential_target_dot) &
            result(valid)
        real(dp), intent(in) :: normal_trace(:, :), tangential_trace(:, :)
        real(dp), intent(in) :: state(:), normal_target(:), tangential_target(:)
        real(dp), intent(in) :: normal_trace_dot(:, :), tangential_trace_dot(:, :)
        real(dp), intent(in) :: normal_weights_dot(:), tangential_weights_dot(:)
        real(dp), intent(in) :: state_dot(:), normal_target_dot(:), tangential_target_dot(:)

        valid = .false.
        if (any(shape(normal_trace_dot) /= shape(normal_trace)) .or. &
            any(shape(tangential_trace_dot) /= shape(tangential_trace)) .or. &
            size(normal_weights_dot) /= size(normal_target) .or. &
            size(tangential_weights_dot) /= size(tangential_target) .or. &
            size(state_dot) /= size(state) .or. &
            size(normal_target_dot) /= size(normal_target) .or. &
            size(tangential_target_dot) /= size(tangential_target)) return
        if (.not. all(ieee_is_finite(normal_trace_dot)) .or. &
            .not. all(ieee_is_finite(tangential_trace_dot)) .or. &
            .not. all(ieee_is_finite(normal_weights_dot)) .or. &
            .not. all(ieee_is_finite(tangential_weights_dot)) .or. &
            .not. all(ieee_is_finite(state_dot)) .or. &
            .not. all(ieee_is_finite(normal_target_dot)) .or. &
            .not. all(ieee_is_finite(tangential_target_dot))) return
        valid = .true.
    end function valid_direction

    logical function valid_adjoint_outputs( &
            normal_trace, tangential_trace, state, normal_target, tangential_target, &
            normal_trace_bar, tangential_trace_bar, normal_weights_bar, &
            tangential_weights_bar, state_bar, normal_target_bar, tangential_target_bar) &
            result(valid)
        real(dp), intent(in) :: normal_trace(:, :), tangential_trace(:, :)
        real(dp), intent(in) :: state(:), normal_target(:), tangential_target(:)
        real(dp), intent(in) :: normal_trace_bar(:, :), tangential_trace_bar(:, :)
        real(dp), intent(in) :: normal_weights_bar(:), tangential_weights_bar(:)
        real(dp), intent(in) :: state_bar(:), normal_target_bar(:), tangential_target_bar(:)

        valid = all(shape(normal_trace_bar) == shape(normal_trace)) .and. &
            all(shape(tangential_trace_bar) == shape(tangential_trace)) .and. &
            size(normal_weights_bar) == size(normal_target) .and. &
            size(tangential_weights_bar) == size(tangential_target) .and. &
            size(state_bar) == size(state) .and. &
            size(normal_target_bar) == size(normal_target) .and. &
            size(tangential_target_bar) == size(tangential_target)
    end function valid_adjoint_outputs

end module fortfem_boundary_trace_residual
