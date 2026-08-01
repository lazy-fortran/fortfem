module fortfem_coupled_field_residual
    !! Generic rectangular field residual plus caller-owned constraints.
    !!
    !! For a field operator A and constraint operator C, the contract is
    !!
    !!   r_f = A u - f,       r_c = C u - g.
    !!
    !! The operators are deliberately dense local blocks.  A caller can form
    !! them from finite-element, boundary, Fourier, tensor, or interface
    !! blocks and retain its own global sparse ownership.  This module only
    !! supplies the composable residual and its exact real JVP/VJP actions.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortfem_kinds, only: dp
    use fortsparse, only: FORTSPARSE_INVALID_MATRIX, FORTSPARSE_OK, &
        fortsparse_status_t, status_set
    implicit none
    private

    public :: assemble_coupled_field_residual
    public :: assemble_coupled_field_residual_jvp
    public :: assemble_coupled_field_residual_vjp

contains

    subroutine assemble_coupled_field_residual( &
            field_operator, state, source, constraint_operator, constraint_target, &
            field_residual, constraint_residual, status)
        real(dp), intent(in) :: field_operator(:, :), state(:), source(:)
        real(dp), intent(in) :: constraint_operator(:, :), constraint_target(:)
        real(dp), intent(out) :: field_residual(:), constraint_residual(:)
        type(fortsparse_status_t), intent(out) :: status

        field_residual = 0.0_dp
        constraint_residual = 0.0_dp
        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            "coupled field residual received incompatible arrays")
        if (.not. valid_shapes(field_operator, state, source, constraint_operator, &
            constraint_target, field_residual, constraint_residual)) return
        if (.not. finite_inputs(field_operator, state, source, constraint_operator, &
            constraint_target)) return
        field_residual = matmul(field_operator, state) - source
        constraint_residual = matmul(constraint_operator, state) - constraint_target
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine assemble_coupled_field_residual

    subroutine assemble_coupled_field_residual_jvp( &
            field_operator, state, source, constraint_operator, constraint_target, &
            field_operator_dot, state_dot, source_dot, constraint_operator_dot, &
            constraint_target_dot, field_residual_dot, constraint_residual_dot, &
            status)
        real(dp), intent(in) :: field_operator(:, :), state(:), source(:)
        real(dp), intent(in) :: constraint_operator(:, :), constraint_target(:)
        real(dp), intent(in) :: field_operator_dot(:, :), state_dot(:), source_dot(:)
        real(dp), intent(in) :: constraint_operator_dot(:, :)
        real(dp), intent(in) :: constraint_target_dot(:)
        real(dp), intent(out) :: field_residual_dot(:), constraint_residual_dot(:)
        type(fortsparse_status_t), intent(out) :: status

        field_residual_dot = 0.0_dp
        constraint_residual_dot = 0.0_dp
        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            "coupled field residual JVP received incompatible arrays")
        if (.not. valid_shapes(field_operator, state, source, constraint_operator, &
            constraint_target, field_residual_dot, constraint_residual_dot)) return
        if (.not. valid_shapes(field_operator_dot, state_dot, source_dot, &
            constraint_operator_dot, constraint_target_dot, field_residual_dot, &
            constraint_residual_dot)) return
        if (.not. finite_inputs(field_operator, state, source, constraint_operator, &
            constraint_target) .or. &
            .not. finite_inputs(field_operator_dot, state_dot, source_dot, &
            constraint_operator_dot, constraint_target_dot)) return
        field_residual_dot = matmul(field_operator_dot, state) + &
            matmul(field_operator, state_dot) - source_dot
        constraint_residual_dot = matmul(constraint_operator_dot, state) + &
            matmul(constraint_operator, state_dot) - constraint_target_dot
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine assemble_coupled_field_residual_jvp

    subroutine assemble_coupled_field_residual_vjp( &
            field_operator, state, source, constraint_operator, constraint_target, &
            field_residual_bar, constraint_residual_bar, field_operator_bar, &
            state_bar, source_bar, constraint_operator_bar, constraint_target_bar, &
            status)
        real(dp), intent(in) :: field_operator(:, :), state(:), source(:)
        real(dp), intent(in) :: constraint_operator(:, :), constraint_target(:)
        real(dp), intent(in) :: field_residual_bar(:), constraint_residual_bar(:)
        real(dp), intent(out) :: field_operator_bar(:, :), state_bar(:)
        real(dp), intent(out) :: source_bar(:), constraint_operator_bar(:, :)
        real(dp), intent(out) :: constraint_target_bar(:)
        type(fortsparse_status_t), intent(out) :: status
        integer :: row, column

        field_operator_bar = 0.0_dp
        state_bar = 0.0_dp
        source_bar = 0.0_dp
        constraint_operator_bar = 0.0_dp
        constraint_target_bar = 0.0_dp
        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            "coupled field residual VJP received incompatible arrays")
        if (.not. valid_shapes(field_operator, state, source, constraint_operator, &
            constraint_target, field_residual_bar, constraint_residual_bar)) return
        if (any(shape(field_operator_bar) /= shape(field_operator)) .or. &
            size(state_bar) /= size(state) .or. size(source_bar) /= size(source) .or. &
            any(shape(constraint_operator_bar) /= shape(constraint_operator)) .or. &
            size(constraint_target_bar) /= size(constraint_target)) return
        if (.not. finite_inputs(field_operator, state, source, constraint_operator, &
            constraint_target) .or. any(.not. ieee_is_finite(field_residual_bar)) .or. &
            any(.not. ieee_is_finite(constraint_residual_bar))) return

        do row = 1, size(field_operator, 1)
            do column = 1, size(field_operator, 2)
                field_operator_bar(row, column) = field_residual_bar(row)*state(column)
            end do
        end do
        do row = 1, size(constraint_operator, 1)
            do column = 1, size(constraint_operator, 2)
                constraint_operator_bar(row, column) = &
                    constraint_residual_bar(row)*state(column)
            end do
        end do
        state_bar = matmul(transpose(field_operator), field_residual_bar) + &
            matmul(transpose(constraint_operator), constraint_residual_bar)
        source_bar = -field_residual_bar
        constraint_target_bar = -constraint_residual_bar
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine assemble_coupled_field_residual_vjp

    logical function valid_shapes( &
            field_operator, state, source, constraint_operator, constraint_target, &
            field_residual, constraint_residual) result(valid)
        real(dp), intent(in) :: field_operator(:, :), state(:), source(:)
        real(dp), intent(in) :: constraint_operator(:, :), constraint_target(:)
        real(dp), intent(in) :: field_residual(:), constraint_residual(:)

        valid = .false.
        if (size(field_operator, 1) < 1 .or. size(field_operator, 2) < 1) return
        if (size(field_operator, 2) /= size(state) .or. &
            size(field_operator, 1) /= size(source) .or. &
            size(field_operator, 2) /= size(constraint_operator, 2) .or. &
            size(constraint_operator, 1) /= size(constraint_target) .or. &
            size(field_residual) /= size(source) .or. &
            size(constraint_residual) /= size(constraint_target)) return
        valid = .true.
    end function valid_shapes

    logical function finite_inputs( &
            field_operator, state, source, constraint_operator, constraint_target) &
            result(valid)
        real(dp), intent(in) :: field_operator(:, :), state(:), source(:)
        real(dp), intent(in) :: constraint_operator(:, :), constraint_target(:)

        valid = all(ieee_is_finite(field_operator)) .and. &
            all(ieee_is_finite(state)) .and. all(ieee_is_finite(source)) .and. &
            all(ieee_is_finite(constraint_operator)) .and. &
            all(ieee_is_finite(constraint_target))
    end function finite_inputs

end module fortfem_coupled_field_residual
