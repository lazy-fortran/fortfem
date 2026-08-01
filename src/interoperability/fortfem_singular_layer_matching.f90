module fortfem_singular_layer_matching
    !! Neutral complex trace matching for an inner and an outer layer.
    !!
    !! The residual is the weighted oriented jump
    !!
    !!   r = w * (T_outer u_outer - T_inner u_inner - jump).
    !!
    !! It is deliberately agnostic about the PDE, singular asymptotics, or
    !! physical normalization.  GPEC, MARS, GLISS, FEM/BEM, DtN, and PML
    !! clients can supply their own trace matrices and states.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortfem_kinds, only: dp
    implicit none
    private

    public :: assemble_singular_layer_matching
    public :: assemble_singular_layer_matching_jvp
    public :: assemble_singular_layer_matching_vjp

    interface finite_complex
        module procedure finite_complex_vector
        module procedure finite_complex_matrix
    end interface finite_complex

contains

    subroutine assemble_singular_layer_matching( &
            outer_trace, inner_trace, weights, outer_state, inner_state, jump, &
            residual, status)
        complex(dp), intent(in) :: outer_trace(:, :), inner_trace(:, :)
        real(dp), intent(in) :: weights(:)
        complex(dp), intent(in) :: outer_state(:), inner_state(:), jump(:)
        complex(dp), intent(out) :: residual(:)
        integer, intent(out) :: status

        residual = cmplx(0.0_dp, 0.0_dp, dp)
        status = 1
        if (.not. valid_inputs( &
            outer_trace, inner_trace, weights, outer_state, inner_state, jump, &
            residual)) return
        residual = weights*(matmul(outer_trace, outer_state) - &
            matmul(inner_trace, inner_state) - jump)
        status = 0
    end subroutine assemble_singular_layer_matching

    subroutine assemble_singular_layer_matching_jvp( &
            outer_trace, inner_trace, weights, outer_state, inner_state, jump, &
            outer_trace_dot, inner_trace_dot, weights_dot, outer_state_dot, &
            inner_state_dot, jump_dot, residual_dot, status)
        complex(dp), intent(in) :: outer_trace(:, :), inner_trace(:, :)
        real(dp), intent(in) :: weights(:)
        complex(dp), intent(in) :: outer_state(:), inner_state(:), jump(:)
        complex(dp), intent(in) :: outer_trace_dot(:, :), inner_trace_dot(:, :)
        real(dp), intent(in) :: weights_dot(:)
        complex(dp), intent(in) :: outer_state_dot(:), inner_state_dot(:), jump_dot(:)
        complex(dp), intent(out) :: residual_dot(:)
        integer, intent(out) :: status

        complex(dp), allocatable :: raw(:), raw_dot(:)

        residual_dot = cmplx(0.0_dp, 0.0_dp, dp)
        status = 1
        if (.not. valid_inputs( &
            outer_trace, inner_trace, weights, outer_state, inner_state, jump, &
            residual_dot)) return
        if (.not. same_shape(outer_trace, outer_trace_dot)) return
        if (.not. same_shape(inner_trace, inner_trace_dot)) return
        if (size(weights_dot) /= size(weights)) return
        if (size(outer_state_dot) /= size(outer_state)) return
        if (size(inner_state_dot) /= size(inner_state)) return
        if (size(jump_dot) /= size(jump)) return
        if (.not. finite_complex(outer_trace_dot) .or. &
            .not. finite_complex(inner_trace_dot) .or. &
            .not. all(ieee_is_finite(weights_dot)) .or. &
            .not. finite_complex(outer_state_dot) .or. &
            .not. finite_complex(inner_state_dot) .or. &
            .not. finite_complex(jump_dot)) return
        allocate(raw(size(weights)), raw_dot(size(weights)))
        raw = matmul(outer_trace, outer_state) - &
            matmul(inner_trace, inner_state) - jump
        raw_dot = matmul(outer_trace_dot, outer_state) + &
            matmul(outer_trace, outer_state_dot) - &
            matmul(inner_trace_dot, inner_state) - &
            matmul(inner_trace, inner_state_dot) - jump_dot
        residual_dot = weights_dot*raw + weights*raw_dot
        status = 0
    end subroutine assemble_singular_layer_matching_jvp

    subroutine assemble_singular_layer_matching_vjp( &
            outer_trace, inner_trace, weights, outer_state, inner_state, jump, &
            residual_bar, outer_trace_bar, inner_trace_bar, weights_bar, &
            outer_state_bar, inner_state_bar, jump_bar, status)
        complex(dp), intent(in) :: outer_trace(:, :), inner_trace(:, :)
        real(dp), intent(in) :: weights(:)
        complex(dp), intent(in) :: outer_state(:), inner_state(:), jump(:)
        complex(dp), intent(in) :: residual_bar(:)
        complex(dp), intent(out) :: outer_trace_bar(:, :), inner_trace_bar(:, :)
        real(dp), intent(out) :: weights_bar(:)
        complex(dp), intent(out) :: outer_state_bar(:), inner_state_bar(:), jump_bar(:)
        integer, intent(out) :: status

        complex(dp), allocatable :: raw(:), weighted_bar(:)

        outer_trace_bar = cmplx(0.0_dp, 0.0_dp, dp)
        inner_trace_bar = cmplx(0.0_dp, 0.0_dp, dp)
        weights_bar = 0.0_dp
        outer_state_bar = cmplx(0.0_dp, 0.0_dp, dp)
        inner_state_bar = cmplx(0.0_dp, 0.0_dp, dp)
        jump_bar = cmplx(0.0_dp, 0.0_dp, dp)
        status = 1
        if (.not. valid_inputs( &
            outer_trace, inner_trace, weights, outer_state, inner_state, jump, &
            residual_bar)) return
        if (size(weights_bar) /= size(weights)) return
        if (.not. same_shape(outer_trace, outer_trace_bar)) return
        if (.not. same_shape(inner_trace, inner_trace_bar)) return
        if (size(outer_state_bar) /= size(outer_state)) return
        if (size(inner_state_bar) /= size(inner_state)) return
        if (size(jump_bar) /= size(jump)) return
        allocate(raw(size(weights)), weighted_bar(size(weights)))
        raw = matmul(outer_trace, outer_state) - &
            matmul(inner_trace, inner_state) - jump
        weighted_bar = weights*residual_bar
        call rank_one_product(outer_trace_bar, weighted_bar, outer_state)
        call rank_one_product(inner_trace_bar, weighted_bar, inner_state)
        inner_trace_bar = -inner_trace_bar
        outer_state_bar = matmul(conjg(transpose(outer_trace)), weighted_bar)
        inner_state_bar = -matmul(conjg(transpose(inner_trace)), weighted_bar)
        jump_bar = -weighted_bar
        weights_bar = real(conjg(residual_bar)*raw, dp)
        status = 0
    end subroutine assemble_singular_layer_matching_vjp

    logical function valid_inputs( &
            outer_trace, inner_trace, weights, outer_state, inner_state, jump, &
            residual) result(valid)
        complex(dp), intent(in) :: outer_trace(:, :), inner_trace(:, :)
        real(dp), intent(in) :: weights(:)
        complex(dp), intent(in) :: outer_state(:), inner_state(:), jump(:)
        complex(dp), intent(in) :: residual(:)
        integer :: constraint_count

        valid = .false.
        constraint_count = size(weights)
        if (constraint_count < 1) return
        if (size(outer_trace, 1) /= constraint_count) return
        if (size(inner_trace, 1) /= constraint_count) return
        if (size(outer_trace, 2) < 1 .or. size(inner_trace, 2) < 1) return
        if (size(outer_state) /= size(outer_trace, 2)) return
        if (size(inner_state) /= size(inner_trace, 2)) return
        if (size(jump) /= constraint_count .or. &
            size(residual) /= constraint_count) return
        if (any(weights <= 0.0_dp)) return
        if (.not. all(ieee_is_finite(weights))) return
        if (.not. finite_complex(outer_trace) .or. &
            .not. finite_complex(inner_trace) .or. &
            .not. finite_complex(outer_state) .or. &
            .not. finite_complex(inner_state) .or. &
            .not. finite_complex(jump)) return
        valid = .true.
    end function valid_inputs

    logical function same_shape(first, second) result(same)
        complex(dp), intent(in) :: first(:, :), second(:, :)

        same = all(shape(first) == shape(second))
    end function same_shape

    subroutine rank_one_product(matrix, left, right)
        complex(dp), intent(out) :: matrix(:, :)
        complex(dp), intent(in) :: left(:), right(:)
        integer :: column, row

        do column = 1, size(right)
            do row = 1, size(left)
                matrix(row, column) = left(row)*conjg(right(column))
            end do
        end do
    end subroutine rank_one_product

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

end module fortfem_singular_layer_matching
