module fortfem_maxwell_fem_bem_state_ad
    !! Differentiable implicit solve for a mixed Maxwell FEM--BEM block.
    !!
    !! Assembly owns the physical volume and trace forms.  This small neutral
    !! layer owns only the reusable state map A x = b, so callers can compose
    !! assembled FEM--BEM, FEM--DtN, or PML matrices with analytical matrix and
    !! right-hand-side products without differentiating through a factorization.
    use fortfem_kinds, only: dp
    use fortnum_linalg, only: &
        dense_solve, linear_solve_complex_jvp, linear_solve_complex_vjp
    implicit none
    private

    public :: solve_maxwell_fem_bem_linear_state
    public :: solve_maxwell_fem_bem_linear_state_jvp
    public :: solve_maxwell_fem_bem_linear_state_vjp

contains

    subroutine solve_maxwell_fem_bem_linear_state(matrix, right_hand_side, &
            state, status)
        complex(dp), intent(in) :: matrix(:, :), right_hand_side(:)
        complex(dp), allocatable, intent(out) :: state(:)
        integer, intent(out) :: status

        integer :: info

        allocate(state(size(right_hand_side)))
        state = cmplx(0.0_dp, 0.0_dp, dp)
        status = 1
        if (.not. valid_system(matrix, right_hand_side)) return
        call dense_solve(matrix, right_hand_side, state, info)
        if (info /= 0) then
            status = 2
            return
        end if
        status = 0
    end subroutine solve_maxwell_fem_bem_linear_state

    subroutine solve_maxwell_fem_bem_linear_state_jvp( &
            matrix, right_hand_side, matrix_dot, right_hand_side_dot, &
            state_dot, status)
        complex(dp), intent(in) :: matrix(:, :), right_hand_side(:)
        complex(dp), intent(in) :: matrix_dot(:, :), right_hand_side_dot(:)
        complex(dp), allocatable, intent(out) :: state_dot(:)
        complex(dp), allocatable :: state(:)
        integer, intent(out) :: status
        integer :: info

        allocate(state_dot(size(right_hand_side)))
        state_dot = cmplx(0.0_dp, 0.0_dp, dp)
        status = 1
        if (.not. valid_system(matrix, right_hand_side)) return
        if (any(shape(matrix_dot) /= shape(matrix)) .or. &
            size(right_hand_side_dot) /= size(right_hand_side)) return
        allocate(state(size(right_hand_side)))
        call dense_solve(matrix, right_hand_side, state, info)
        if (info /= 0) then
            status = 2
            return
        end if
        call linear_solve_complex_jvp( &
            size(right_hand_side), matrix, state, matrix_dot, &
            right_hand_side_dot, state_dot, info)
        if (info /= 0) then
            status = 2
            return
        end if
        status = 0
    end subroutine solve_maxwell_fem_bem_linear_state_jvp

    subroutine solve_maxwell_fem_bem_linear_state_vjp( &
            matrix, right_hand_side, state_bar, matrix_bar, right_hand_side_bar, &
            status)
        complex(dp), intent(in) :: matrix(:, :), right_hand_side(:)
        complex(dp), intent(in) :: state_bar(:)
        complex(dp), allocatable, intent(out) :: matrix_bar(:, :)
        complex(dp), allocatable, intent(out) :: right_hand_side_bar(:)
        complex(dp), allocatable :: state(:)
        integer, intent(out) :: status
        integer :: info

        allocate(matrix_bar(size(matrix, 1), size(matrix, 2)))
        allocate(right_hand_side_bar(size(right_hand_side)))
        matrix_bar = cmplx(0.0_dp, 0.0_dp, dp)
        right_hand_side_bar = cmplx(0.0_dp, 0.0_dp, dp)
        status = 1
        if (.not. valid_system(matrix, right_hand_side)) return
        if (size(state_bar) /= size(right_hand_side)) return
        allocate(state(size(right_hand_side)))
        call dense_solve(matrix, right_hand_side, state, info)
        if (info /= 0) then
            status = 2
            return
        end if
        call linear_solve_complex_vjp( &
            size(right_hand_side), matrix, state, state_bar, matrix_bar, &
            right_hand_side_bar, info)
        if (info /= 0) then
            status = 2
            return
        end if
        status = 0
    end subroutine solve_maxwell_fem_bem_linear_state_vjp

    logical function valid_system(matrix, right_hand_side)
        complex(dp), intent(in) :: matrix(:, :), right_hand_side(:)

        valid_system = size(matrix, 1) > 0 .and. &
            size(matrix, 1) == size(matrix, 2) .and. &
            size(matrix, 1) == size(right_hand_side)
    end function valid_system

end module fortfem_maxwell_fem_bem_state_ad
