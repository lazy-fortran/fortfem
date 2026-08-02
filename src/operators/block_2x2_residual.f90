module fortfem_block_2x2_residual
    !! Neutral two-field rectangular block residual and exact derivatives.
    !!
    !! For independent states x and y, the caller-owned block system is
    !!
    !!   r_1 = A_11 x + A_12 y - f_1,
    !!   r_2 = A_21 x + A_22 y - f_2.
    !!
    !! The blocks may be local FEM, BEM, DtN, PML, interface, tensor, or
    !! Fourier actions.  This module deliberately does not assemble global
    !! sparse storage, choose a Schur complement, or impose a physical closure.
    !! It is the small composition contract used by those higher-level paths.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortfem_kinds, only: dp
    use fortsparse, only: FORTSPARSE_INVALID_MATRIX, FORTSPARSE_OK, &
        fortsparse_status_t, status_set
    implicit none
    private

    public :: assemble_block_2x2_residual
    public :: assemble_block_2x2_residual_jvp
    public :: assemble_block_2x2_residual_vjp

contains

    subroutine assemble_block_2x2_residual( &
            block_11, block_12, block_21, block_22, state_1, state_2, &
            rhs_1, rhs_2, residual_1, residual_2, status)
        real(dp), intent(in) :: block_11(:, :), block_12(:, :)
        real(dp), intent(in) :: block_21(:, :), block_22(:, :)
        real(dp), intent(in) :: state_1(:), state_2(:), rhs_1(:), rhs_2(:)
        real(dp), intent(out) :: residual_1(:), residual_2(:)
        type(fortsparse_status_t), intent(out) :: status

        residual_1 = 0.0_dp
        residual_2 = 0.0_dp
        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            "2x2 block residual received incompatible arrays")
        if (.not. valid_shapes(block_11, block_12, block_21, block_22, state_1, &
            state_2, rhs_1, rhs_2, residual_1, residual_2)) return
        if (.not. finite_inputs(block_11, block_12, block_21, block_22, state_1, &
            state_2, rhs_1, rhs_2)) return
        residual_1 = matmul(block_11, state_1) + matmul(block_12, state_2) - rhs_1
        residual_2 = matmul(block_21, state_1) + matmul(block_22, state_2) - rhs_2
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine assemble_block_2x2_residual

    subroutine assemble_block_2x2_residual_jvp( &
            block_11, block_12, block_21, block_22, state_1, state_2, rhs_1, rhs_2, &
            block_11_dot, block_12_dot, block_21_dot, block_22_dot, state_1_dot, &
            state_2_dot, rhs_1_dot, rhs_2_dot, residual_1_dot, residual_2_dot, status)
        real(dp), intent(in) :: block_11(:, :), block_12(:, :)
        real(dp), intent(in) :: block_21(:, :), block_22(:, :)
        real(dp), intent(in) :: state_1(:), state_2(:), rhs_1(:), rhs_2(:)
        real(dp), intent(in) :: block_11_dot(:, :), block_12_dot(:, :)
        real(dp), intent(in) :: block_21_dot(:, :), block_22_dot(:, :)
        real(dp), intent(in) :: state_1_dot(:), state_2_dot(:)
        real(dp), intent(in) :: rhs_1_dot(:), rhs_2_dot(:)
        real(dp), intent(out) :: residual_1_dot(:), residual_2_dot(:)
        type(fortsparse_status_t), intent(out) :: status

        residual_1_dot = 0.0_dp
        residual_2_dot = 0.0_dp
        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            "2x2 block residual JVP received incompatible arrays")
        if (.not. valid_shapes(block_11, block_12, block_21, block_22, state_1, &
            state_2, rhs_1, rhs_2, residual_1_dot, residual_2_dot)) return
        if (.not. valid_shapes(block_11_dot, block_12_dot, block_21_dot, &
            block_22_dot, state_1_dot, state_2_dot, rhs_1_dot, rhs_2_dot, &
            residual_1_dot, residual_2_dot)) return
        if (.not. same_shapes(block_11_dot, block_11) .or. &
            .not. same_shapes(block_12_dot, block_12) .or. &
            .not. same_shapes(block_21_dot, block_21) .or. &
            .not. same_shapes(block_22_dot, block_22) .or. &
            size(state_1_dot) /= size(state_1) .or. size(state_2_dot) /= size(state_2) .or. &
            size(rhs_1_dot) /= size(rhs_1) .or. size(rhs_2_dot) /= size(rhs_2)) return
        if (.not. finite_inputs(block_11, block_12, block_21, block_22, state_1, &
            state_2, rhs_1, rhs_2) .or. &
            .not. finite_inputs(block_11_dot, block_12_dot, block_21_dot, block_22_dot, &
            state_1_dot, state_2_dot, rhs_1_dot, rhs_2_dot)) return
        residual_1_dot = matmul(block_11_dot, state_1) + matmul(block_11, state_1_dot) + &
            matmul(block_12_dot, state_2) + matmul(block_12, state_2_dot) - rhs_1_dot
        residual_2_dot = matmul(block_21_dot, state_1) + matmul(block_21, state_1_dot) + &
            matmul(block_22_dot, state_2) + matmul(block_22, state_2_dot) - rhs_2_dot
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine assemble_block_2x2_residual_jvp

    subroutine assemble_block_2x2_residual_vjp( &
            block_11, block_12, block_21, block_22, state_1, state_2, rhs_1, rhs_2, &
            residual_1_bar, residual_2_bar, block_11_bar, block_12_bar, &
            block_21_bar, block_22_bar, state_1_bar, state_2_bar, rhs_1_bar, &
            rhs_2_bar, status)
        real(dp), intent(in) :: block_11(:, :), block_12(:, :)
        real(dp), intent(in) :: block_21(:, :), block_22(:, :)
        real(dp), intent(in) :: state_1(:), state_2(:), rhs_1(:), rhs_2(:)
        real(dp), intent(in) :: residual_1_bar(:), residual_2_bar(:)
        real(dp), intent(out) :: block_11_bar(:, :), block_12_bar(:, :)
        real(dp), intent(out) :: block_21_bar(:, :), block_22_bar(:, :)
        real(dp), intent(out) :: state_1_bar(:), state_2_bar(:)
        real(dp), intent(out) :: rhs_1_bar(:), rhs_2_bar(:)
        type(fortsparse_status_t), intent(out) :: status
        integer :: row, column

        block_11_bar = 0.0_dp
        block_12_bar = 0.0_dp
        block_21_bar = 0.0_dp
        block_22_bar = 0.0_dp
        state_1_bar = 0.0_dp
        state_2_bar = 0.0_dp
        rhs_1_bar = 0.0_dp
        rhs_2_bar = 0.0_dp
        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            "2x2 block residual VJP received incompatible arrays")
        if (.not. valid_shapes(block_11, block_12, block_21, block_22, state_1, &
            state_2, rhs_1, rhs_2, residual_1_bar, residual_2_bar)) return
        if (.not. same_shapes(block_11_bar, block_11) .or. &
            .not. same_shapes(block_12_bar, block_12) .or. &
            .not. same_shapes(block_21_bar, block_21) .or. &
            .not. same_shapes(block_22_bar, block_22) .or. &
            size(state_1_bar) /= size(state_1) .or. size(state_2_bar) /= size(state_2) .or. &
            size(rhs_1_bar) /= size(rhs_1) .or. size(rhs_2_bar) /= size(rhs_2)) return
        if (.not. finite_inputs(block_11, block_12, block_21, block_22, state_1, &
            state_2, rhs_1, rhs_2) .or. any(.not. ieee_is_finite(residual_1_bar)) .or. &
            any(.not. ieee_is_finite(residual_2_bar))) return

        do row = 1, size(block_11, 1)
            do column = 1, size(block_11, 2)
                block_11_bar(row, column) = residual_1_bar(row)*state_1(column)
            end do
        end do
        do row = 1, size(block_12, 1)
            do column = 1, size(block_12, 2)
                block_12_bar(row, column) = residual_1_bar(row)*state_2(column)
            end do
        end do
        do row = 1, size(block_21, 1)
            do column = 1, size(block_21, 2)
                block_21_bar(row, column) = residual_2_bar(row)*state_1(column)
            end do
        end do
        do row = 1, size(block_22, 1)
            do column = 1, size(block_22, 2)
                block_22_bar(row, column) = residual_2_bar(row)*state_2(column)
            end do
        end do
        state_1_bar = matmul(transpose(block_11), residual_1_bar) + &
            matmul(transpose(block_21), residual_2_bar)
        state_2_bar = matmul(transpose(block_12), residual_1_bar) + &
            matmul(transpose(block_22), residual_2_bar)
        rhs_1_bar = -residual_1_bar
        rhs_2_bar = -residual_2_bar
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine assemble_block_2x2_residual_vjp

    logical function valid_shapes( &
            block_11, block_12, block_21, block_22, state_1, state_2, rhs_1, rhs_2, &
            residual_1, residual_2) result(valid)
        real(dp), intent(in) :: block_11(:, :), block_12(:, :)
        real(dp), intent(in) :: block_21(:, :), block_22(:, :)
        real(dp), intent(in) :: state_1(:), state_2(:), rhs_1(:), rhs_2(:)
        real(dp), intent(in) :: residual_1(:), residual_2(:)

        valid = .false.
        if (size(block_11, 1) < 1 .or. size(block_11, 2) < 1 .or. &
            size(block_12, 1) < 1 .or. size(block_12, 2) < 1 .or. &
            size(block_21, 1) < 1 .or. size(block_21, 2) < 1 .or. &
            size(block_22, 1) < 1 .or. size(block_22, 2) < 1) return
        if (size(block_11, 1) /= size(block_12, 1) .or. &
            size(block_21, 1) /= size(block_22, 1) .or. &
            size(block_11, 2) /= size(block_21, 2) .or. &
            size(block_12, 2) /= size(block_22, 2) .or. &
            size(block_11, 2) /= size(state_1) .or. &
            size(block_12, 2) /= size(state_2) .or. &
            size(block_11, 1) /= size(rhs_1) .or. &
            size(block_21, 1) /= size(rhs_2) .or. &
            size(residual_1) /= size(rhs_1) .or. size(residual_2) /= size(rhs_2)) return
        valid = .true.
    end function valid_shapes

    logical function same_shapes(left, right) result(same)
        real(dp), intent(in) :: left(:, :), right(:, :)

        same = all(shape(left) == shape(right))
    end function same_shapes

    logical function finite_inputs( &
            block_11, block_12, block_21, block_22, state_1, state_2, rhs_1, rhs_2) &
            result(valid)
        real(dp), intent(in) :: block_11(:, :), block_12(:, :)
        real(dp), intent(in) :: block_21(:, :), block_22(:, :)
        real(dp), intent(in) :: state_1(:), state_2(:), rhs_1(:), rhs_2(:)

        valid = all(ieee_is_finite(block_11)) .and. all(ieee_is_finite(block_12)) .and. &
            all(ieee_is_finite(block_21)) .and. all(ieee_is_finite(block_22)) .and. &
            all(ieee_is_finite(state_1)) .and. all(ieee_is_finite(state_2)) .and. &
            all(ieee_is_finite(rhs_1)) .and. all(ieee_is_finite(rhs_2))
    end function finite_inputs

end module fortfem_block_2x2_residual
