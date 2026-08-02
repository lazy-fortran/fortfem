module fortfem_complex_coupled_field_residual
    !! Complex rectangular field residual plus caller-owned constraints.
    !!
    !! For a field operator A and constraint operator C, the contract is
    !!
    !!   r_f = A u - f,       r_c = C u - g.
    !!
    !! The blocks may be rectangular and may represent FEM, BEM, DtN, PML,
    !! wall, Fourier, or interface actions.  This module only supplies the
    !! composable complex residual and its exact real-part adjoint products;
    !! ownership of global sparse storage and physical closure stays with the
    !! caller.  No dense matrix other than the caller-provided local blocks is
    !! formed.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortfem_kinds, only: dp
    use fortsparse, only: FORTSPARSE_INVALID_MATRIX, FORTSPARSE_OK, &
        fortsparse_status_t, status_set
    implicit none
    private

    public :: assemble_complex_coupled_field_residual
    public :: assemble_complex_coupled_field_residual_jvp
    public :: assemble_complex_coupled_field_residual_vjp

    interface finite_complex
        module procedure finite_complex_vector
        module procedure finite_complex_matrix
    end interface finite_complex

contains

    subroutine assemble_complex_coupled_field_residual( &
            field_operator, state, source, constraint_operator, constraint_target, &
            field_residual, constraint_residual, status)
        complex(dp), intent(in) :: field_operator(:, :), state(:), source(:)
        complex(dp), intent(in) :: constraint_operator(:, :), constraint_target(:)
        complex(dp), intent(out) :: field_residual(:), constraint_residual(:)
        type(fortsparse_status_t), intent(out) :: status

        field_residual = cmplx(0.0_dp, 0.0_dp, dp)
        constraint_residual = cmplx(0.0_dp, 0.0_dp, dp)
        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            "complex coupled field residual received incompatible arrays")
        if (.not. valid_shapes(field_operator, state, source, constraint_operator, &
            constraint_target, field_residual, constraint_residual)) return
        if (.not. finite_inputs(field_operator, state, source, constraint_operator, &
            constraint_target)) return
        field_residual = matmul(field_operator, state) - source
        constraint_residual = matmul(constraint_operator, state) - constraint_target
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine assemble_complex_coupled_field_residual

    subroutine assemble_complex_coupled_field_residual_jvp( &
            field_operator, state, source, constraint_operator, constraint_target, &
            field_operator_dot, state_dot, source_dot, constraint_operator_dot, &
            constraint_target_dot, field_residual_dot, constraint_residual_dot, &
            status)
        complex(dp), intent(in) :: field_operator(:, :), state(:), source(:)
        complex(dp), intent(in) :: constraint_operator(:, :), constraint_target(:)
        complex(dp), intent(in) :: field_operator_dot(:, :), state_dot(:)
        complex(dp), intent(in) :: source_dot(:), constraint_operator_dot(:, :)
        complex(dp), intent(in) :: constraint_target_dot(:)
        complex(dp), intent(out) :: field_residual_dot(:), constraint_residual_dot(:)
        type(fortsparse_status_t), intent(out) :: status

        field_residual_dot = cmplx(0.0_dp, 0.0_dp, dp)
        constraint_residual_dot = cmplx(0.0_dp, 0.0_dp, dp)
        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            "complex coupled field residual JVP received incompatible arrays")
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
    end subroutine assemble_complex_coupled_field_residual_jvp

    subroutine assemble_complex_coupled_field_residual_vjp( &
            field_operator, state, source, constraint_operator, constraint_target, &
            field_residual_bar, constraint_residual_bar, field_operator_bar, &
            state_bar, source_bar, constraint_operator_bar, constraint_target_bar, &
            status)
        complex(dp), intent(in) :: field_operator(:, :), state(:), source(:)
        complex(dp), intent(in) :: constraint_operator(:, :), constraint_target(:)
        complex(dp), intent(in) :: field_residual_bar(:), constraint_residual_bar(:)
        complex(dp), intent(out) :: field_operator_bar(:, :), state_bar(:)
        complex(dp), intent(out) :: source_bar(:), constraint_operator_bar(:, :)
        complex(dp), intent(out) :: constraint_target_bar(:)
        type(fortsparse_status_t), intent(out) :: status
        integer :: column, row

        field_operator_bar = cmplx(0.0_dp, 0.0_dp, dp)
        state_bar = cmplx(0.0_dp, 0.0_dp, dp)
        source_bar = cmplx(0.0_dp, 0.0_dp, dp)
        constraint_operator_bar = cmplx(0.0_dp, 0.0_dp, dp)
        constraint_target_bar = cmplx(0.0_dp, 0.0_dp, dp)
        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            "complex coupled field residual VJP received incompatible arrays")
        if (.not. valid_shapes(field_operator, state, source, constraint_operator, &
            constraint_target, field_residual_bar, constraint_residual_bar)) return
        if (any(shape(field_operator_bar) /= shape(field_operator)) .or. &
            size(state_bar) /= size(state) .or. size(source_bar) /= size(source) .or. &
            any(shape(constraint_operator_bar) /= shape(constraint_operator)) .or. &
            size(constraint_target_bar) /= size(constraint_target)) return
        if (.not. finite_inputs(field_operator, state, source, constraint_operator, &
            constraint_target) .or. .not. finite_complex(field_residual_bar) .or. &
            .not. finite_complex(constraint_residual_bar)) return

        do row = 1, size(field_operator, 1)
            do column = 1, size(field_operator, 2)
                field_operator_bar(row, column) = field_residual_bar(row)* &
                    conjg(state(column))
            end do
        end do
        do row = 1, size(constraint_operator, 1)
            do column = 1, size(constraint_operator, 2)
                constraint_operator_bar(row, column) = &
                    constraint_residual_bar(row)*conjg(state(column))
            end do
        end do
        state_bar = matmul(conjg(transpose(field_operator)), field_residual_bar) + &
            matmul(conjg(transpose(constraint_operator)), constraint_residual_bar)
        source_bar = -field_residual_bar
        constraint_target_bar = -constraint_residual_bar
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine assemble_complex_coupled_field_residual_vjp

    logical function valid_shapes( &
            field_operator, state, source, constraint_operator, constraint_target, &
            field_residual, constraint_residual) result(valid)
        complex(dp), intent(in) :: field_operator(:, :), state(:), source(:)
        complex(dp), intent(in) :: constraint_operator(:, :), constraint_target(:)
        complex(dp), intent(in) :: field_residual(:), constraint_residual(:)

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
        complex(dp), intent(in) :: field_operator(:, :), state(:), source(:)
        complex(dp), intent(in) :: constraint_operator(:, :), constraint_target(:)

        valid = finite_complex(field_operator) .and. finite_complex(state) .and. &
            finite_complex(source) .and. finite_complex(constraint_operator) .and. &
            finite_complex(constraint_target)
    end function finite_inputs

    logical function finite_complex_vector(values) result(valid)
        complex(dp), intent(in) :: values(:)

        valid = all(ieee_is_finite(real(values, dp))) .and. &
            all(ieee_is_finite(aimag(values)))
    end function finite_complex_vector

    logical function finite_complex_matrix(values) result(valid)
        complex(dp), intent(in) :: values(:, :)

        valid = all(ieee_is_finite(real(values, dp))) .and. &
            all(ieee_is_finite(aimag(values)))
    end function finite_complex_matrix

end module fortfem_complex_coupled_field_residual
