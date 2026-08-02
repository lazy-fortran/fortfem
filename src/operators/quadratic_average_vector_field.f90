module fortfem_quadratic_average_vector_field
    !! Average-vector-field step for a caller-owned quadratic Hamiltonian.
    !!
    !! For H(x) = 1/2 x^T K x + g^T x and a constant skew interconnection J,
    !! the AVF/discrete-gradient update is
    !!
    !!   x_{n+1} - x_n = h J [ K (x_n+x_{n+1})/2 + g ].
    !!
    !! The caller owns K, J, and g.  When K is symmetric and J is skew, the
    !! update preserves H exactly (up to the linear solve) and is reversible
    !! under a signed step reversal.  Dissipation is deliberately not folded
    !! into this ideal contract; clients should compose it with a separate
    !! dissipative block.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortfem_kinds, only: dp
    use fortnum_linalg, only: dense_solve
    use fortsparse, only: fortsparse_status_t, status_set, &
        FORTSPARSE_INVALID_MATRIX, FORTSPARSE_OK, FORTSPARSE_SINGULAR
    implicit none
    private

    public :: advance_quadratic_avf
    public :: advance_quadratic_avf_jvp
    public :: advance_quadratic_avf_vjp

contains

    subroutine advance_quadratic_avf( &
            hamiltonian, interconnection, linear_term, time_step, state, &
            state_next, status)
        !! Advance one quadratic average-vector-field step.
        real(dp), intent(in) :: hamiltonian(:, :), interconnection(:, :)
        real(dp), intent(in) :: linear_term(:), time_step, state(:)
        real(dp), intent(out) :: state_next(:)
        type(fortsparse_status_t), intent(out) :: status

        integer :: state_count, info
        real(dp), allocatable :: matrix_a(:, :), matrix_b(:, :), rhs(:)

        state_next = 0.0_dp
        if (.not. validate_inputs(hamiltonian, interconnection, linear_term, &
                time_step, state, state_next, status)) return
        state_count = size(state)
        allocate(matrix_a(state_count, state_count), &
            matrix_b(state_count, state_count), rhs(state_count))
        call assemble_system(hamiltonian, interconnection, linear_term, &
            time_step, state, matrix_a, matrix_b, rhs)
        call dense_solve(matrix_a, rhs, state_next, info)
        if (info /= 0) then
            call status_set(status, FORTSPARSE_SINGULAR, &
                "quadratic AVF system is singular")
            return
        end if
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine advance_quadratic_avf

    subroutine advance_quadratic_avf_jvp( &
            hamiltonian, interconnection, linear_term, time_step, state, &
            hamiltonian_dot, interconnection_dot, linear_term_dot, &
            time_step_dot, state_dot, state_next_dot, status)
        !! Apply the exact tangent of one quadratic AVF step.
        real(dp), intent(in) :: hamiltonian(:, :), interconnection(:, :)
        real(dp), intent(in) :: linear_term(:), time_step, state(:)
        real(dp), intent(in) :: hamiltonian_dot(:, :), interconnection_dot(:, :)
        real(dp), intent(in) :: linear_term_dot(:), time_step_dot, state_dot(:)
        real(dp), intent(out) :: state_next_dot(:)
        type(fortsparse_status_t), intent(out) :: status

        integer :: state_count, info
        real(dp), allocatable :: matrix_a(:, :), matrix_b(:, :), rhs(:)
        real(dp), allocatable :: matrix_a_dot(:, :), rhs_dot(:), state_next(:)

        state_next_dot = 0.0_dp
        if (.not. validate_inputs(hamiltonian, interconnection, linear_term, &
                time_step, state, state_next_dot, status)) return
        state_count = size(state)
        if (.not. validate_direction(hamiltonian, interconnection, linear_term, &
                hamiltonian_dot, interconnection_dot, linear_term_dot, &
                state_dot, state_count, time_step_dot)) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "quadratic AVF JVP has incompatible increments")
            return
        end if
        allocate(matrix_a(state_count, state_count), &
            matrix_b(state_count, state_count), rhs(state_count), &
            matrix_a_dot(state_count, state_count), rhs_dot(state_count), &
            state_next(state_count))
        call assemble_system(hamiltonian, interconnection, linear_term, &
            time_step, state, matrix_a, matrix_b, rhs)
        call dense_solve(matrix_a, rhs, state_next, info)
        if (info /= 0) then
            call status_set(status, FORTSPARSE_SINGULAR, &
                "quadratic AVF JVP base system is singular")
            return
        end if
        matrix_a_dot = -0.5_dp*(time_step_dot*matmul(interconnection, &
            hamiltonian) + time_step*matmul(interconnection_dot, hamiltonian) &
            + time_step*matmul(interconnection, hamiltonian_dot))
        rhs_dot = matmul(0.5_dp*(time_step_dot*matmul(interconnection, &
            hamiltonian) + time_step*matmul(interconnection_dot, hamiltonian) &
            + time_step*matmul(interconnection, hamiltonian_dot)), state) &
            + matmul(matrix_b, state_dot) + time_step_dot*matmul(interconnection, &
            linear_term) + time_step*matmul(interconnection_dot, linear_term) &
            + time_step*matmul(interconnection, linear_term_dot)
        rhs_dot = rhs_dot - matmul(matrix_a_dot, state_next)
        call dense_solve(matrix_a, rhs_dot, state_next_dot, info)
        if (info /= 0) then
            call status_set(status, FORTSPARSE_SINGULAR, &
                "quadratic AVF JVP system is singular")
            return
        end if
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine advance_quadratic_avf_jvp

    subroutine advance_quadratic_avf_vjp( &
            hamiltonian, interconnection, linear_term, time_step, state, &
            state_next, state_next_bar, hamiltonian_bar, interconnection_bar, &
            linear_term_bar, time_step_bar, state_bar, status)
        !! Apply the real adjoint of one quadratic AVF step.
        real(dp), intent(in) :: hamiltonian(:, :), interconnection(:, :)
        real(dp), intent(in) :: linear_term(:), time_step, state(:)
        real(dp), intent(in) :: state_next(:), state_next_bar(:)
        real(dp), intent(out) :: hamiltonian_bar(:, :), interconnection_bar(:, :)
        real(dp), intent(out) :: linear_term_bar(:), time_step_bar, state_bar(:)
        type(fortsparse_status_t), intent(out) :: status

        integer :: state_count, info
        real(dp), allocatable :: matrix_a(:, :)
        real(dp), allocatable :: adjoint(:), half_state(:), hamiltonian_half(:)

        hamiltonian_bar = 0.0_dp
        interconnection_bar = 0.0_dp
        linear_term_bar = 0.0_dp
        time_step_bar = 0.0_dp
        state_bar = 0.0_dp
        if (.not. validate_inputs(hamiltonian, interconnection, linear_term, &
                time_step, state, state_next, status)) return
        state_count = size(state)
        if (size(state_next_bar) /= state_count .or. &
                size(hamiltonian_bar, 1) /= state_count .or. &
                size(hamiltonian_bar, 2) /= state_count .or. &
                size(interconnection_bar, 1) /= state_count .or. &
                size(interconnection_bar, 2) /= state_count .or. &
                size(linear_term_bar) /= state_count .or. &
                size(state_bar) /= state_count .or. &
                any(.not. ieee_is_finite(state_next_bar))) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "quadratic AVF VJP has incompatible cotangents")
            return
        end if
        allocate(matrix_a(state_count, state_count), adjoint(state_count), &
            half_state(state_count), hamiltonian_half(state_count))
        matrix_a = 0.0_dp
        matrix_a = identity(matrix_a) - 0.5_dp*time_step*matmul( &
            interconnection, hamiltonian)
        call dense_solve(transpose(matrix_a), state_next_bar, adjoint, info)
        if (info /= 0) then
            call status_set(status, FORTSPARSE_SINGULAR, &
                "quadratic AVF VJP transpose system is singular")
            return
        end if
        half_state = 0.5_dp*(state + state_next)
        hamiltonian_half = matmul(hamiltonian, half_state)
        hamiltonian_bar = outer_product(matmul(transpose(interconnection), &
            adjoint), time_step*half_state)
        interconnection_bar = outer_product(adjoint, time_step*(&
            hamiltonian_half + linear_term))
        linear_term_bar = time_step*matmul(transpose(interconnection), adjoint)
        time_step_bar = dot_product(adjoint, matmul(interconnection, &
            hamiltonian_half + linear_term))
        state_bar = matmul(identity(matrix_a) + 0.5_dp*time_step* &
            matmul(transpose(hamiltonian), transpose(interconnection)), adjoint)
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine advance_quadratic_avf_vjp

    subroutine assemble_system( &
            hamiltonian, interconnection, linear_term, time_step, state, &
            matrix_a, matrix_b, rhs)
        real(dp), intent(in) :: hamiltonian(:, :), interconnection(:, :)
        real(dp), intent(in) :: linear_term(:), time_step, state(:)
        real(dp), intent(out) :: matrix_a(:, :), matrix_b(:, :), rhs(:)

        integer :: state_count
        real(dp), allocatable :: identity_matrix(:, :), generator(:, :)

        state_count = size(state)
        allocate(identity_matrix(state_count, state_count), &
            generator(state_count, state_count))
        identity_matrix = 0.0_dp
        identity_matrix = identity(identity_matrix)
        generator = matmul(interconnection, hamiltonian)
        matrix_a = identity_matrix - 0.5_dp*time_step*generator
        matrix_b = identity_matrix + 0.5_dp*time_step*generator
        rhs = matmul(matrix_b, state) + time_step*matmul(interconnection, linear_term)
    end subroutine assemble_system

    logical function validate_inputs( &
            hamiltonian, interconnection, linear_term, time_step, state, target, &
            status) result(valid)
        real(dp), intent(in) :: hamiltonian(:, :), interconnection(:, :)
        real(dp), intent(in) :: linear_term(:), time_step, state(:), target(:)
        type(fortsparse_status_t), intent(out) :: status
        integer :: state_count

        valid = .false.
        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            "quadratic AVF has incompatible arrays")
        state_count = size(state)
        if (state_count < 1 .or. size(hamiltonian, 1) /= state_count .or. &
                size(hamiltonian, 2) /= state_count .or. &
                size(interconnection, 1) /= state_count .or. &
                size(interconnection, 2) /= state_count .or. &
                size(linear_term) /= state_count .or. size(target) /= state_count .or. &
                .not. ieee_is_finite(time_step)) return
        if (any(.not. ieee_is_finite(hamiltonian)) .or. &
                any(.not. ieee_is_finite(interconnection)) .or. &
                any(.not. ieee_is_finite(linear_term)) .or. &
                any(.not. ieee_is_finite(state)) .or. &
                any(.not. ieee_is_finite(target))) return
        valid = .true.
        call status_set(status, FORTSPARSE_OK, "")
    end function validate_inputs

    pure logical function validate_direction( &
            hamiltonian, interconnection, linear_term, hamiltonian_dot, &
            interconnection_dot, linear_term_dot, state_dot, state_count, &
            time_step_dot) result(valid)
        real(dp), intent(in) :: hamiltonian(:, :), interconnection(:, :)
        real(dp), intent(in) :: linear_term(:), hamiltonian_dot(:, :)
        real(dp), intent(in) :: interconnection_dot(:, :), linear_term_dot(:)
        real(dp), intent(in) :: state_dot(:), time_step_dot
        integer, intent(in) :: state_count

        valid = all(shape(hamiltonian_dot) == [state_count, state_count]) .and. &
            all(shape(interconnection_dot) == [state_count, state_count]) .and. &
            size(linear_term_dot) == state_count .and. size(state_dot) == state_count &
            .and. ieee_is_finite(time_step_dot) .and. &
            all(ieee_is_finite(hamiltonian_dot)) .and. &
            all(ieee_is_finite(interconnection_dot)) .and. &
            all(ieee_is_finite(linear_term_dot)) .and. all(ieee_is_finite(state_dot))
    end function validate_direction

    pure function identity(matrix) result(result_matrix)
        real(dp), intent(in) :: matrix(:, :)
        real(dp) :: result_matrix(size(matrix, 1), size(matrix, 2))
        integer :: index

        result_matrix = 0.0_dp
        do index = 1, min(size(matrix, 1), size(matrix, 2))
            result_matrix(index, index) = 1.0_dp
        end do
    end function identity

    pure function outer_product(left, right) result(product)
        real(dp), intent(in) :: left(:), right(:)
        real(dp) :: product(size(left), size(right))
        integer :: i, j

        do j = 1, size(right)
            do i = 1, size(left)
                product(i, j) = left(i)*right(j)
            end do
        end do
    end function outer_product

end module fortfem_quadratic_average_vector_field
