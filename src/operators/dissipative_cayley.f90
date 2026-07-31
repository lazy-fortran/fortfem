module fortfem_dissipative_cayley
    !! Energy-contracting Cayley step for a positive-time dissipative block.
    !!
    !! For a mass matrix M and caller-owned damping/resistivity D, the step is
    !!
    !!   x_next = (M + h D/2)^(-1) (M - h D/2) x.
    !!
    !! If M is symmetric positive definite and D is symmetric positive
    !! semidefinite, h >= 0 gives a non-increasing M-energy.  This block is
    !! intentionally separate from the Hamiltonian/symplectic wave steps:
    !! splitting clients can compose reversible ideal and dissipative pieces
    !! without calling the dissipative map symplectic.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortfem_kinds, only: dp
    use fortnum_linalg, only: dense_solve
    use fortsparse, only: fortsparse_status_t, status_set, &
        FORTSPARSE_INVALID_MATRIX, FORTSPARSE_OK, FORTSPARSE_SINGULAR
    implicit none
    private

    public :: advance_dissipative_cayley
    public :: advance_dissipative_cayley_jvp
    public :: advance_dissipative_cayley_vjp

contains

    subroutine advance_dissipative_cayley( &
            mass, damping, time_step, state, state_next, status)
        real(dp), intent(in) :: mass(:, :), damping(:, :), time_step, state(:)
        real(dp), intent(out) :: state_next(:)
        type(fortsparse_status_t), intent(out) :: status

        integer :: state_count, info
        real(dp), allocatable :: matrix_a(:, :), matrix_b(:, :), rhs(:)

        state_next = 0.0_dp
        if (.not. validate_inputs( &
            mass, damping, time_step, state, state_next, status)) return
        state_count = size(state)
        allocate(matrix_a(state_count, state_count), &
            matrix_b(state_count, state_count), rhs(state_count))
        matrix_a = mass + 0.5_dp*time_step*damping
        matrix_b = mass - 0.5_dp*time_step*damping
        rhs = matmul(matrix_b, state)
        call dense_solve(matrix_a, rhs, state_next, info)
        if (info /= 0) then
            call status_set(status, FORTSPARSE_SINGULAR, &
                "dissipative Cayley system is singular")
            return
        end if
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine advance_dissipative_cayley

    subroutine advance_dissipative_cayley_jvp( &
            mass, damping, time_step, state, mass_dot, damping_dot, &
            time_step_dot, state_dot, state_next_dot, status)
        real(dp), intent(in) :: mass(:, :), damping(:, :), time_step, state(:)
        real(dp), intent(in) :: mass_dot(:, :), damping_dot(:, :), time_step_dot
        real(dp), intent(in) :: state_dot(:)
        real(dp), intent(out) :: state_next_dot(:)
        type(fortsparse_status_t), intent(out) :: status

        integer :: state_count, info
        real(dp), allocatable :: matrix_a(:, :), matrix_b(:, :), rhs(:)
        real(dp), allocatable :: matrix_a_dot(:, :), matrix_b_dot(:, :)
        real(dp), allocatable :: rhs_dot(:), state_next(:)

        state_next_dot = 0.0_dp
        if (.not. validate_inputs( &
            mass, damping, time_step, state, state_next_dot, status)) return
        state_count = size(state)
        if (.not. validate_direction( &
            mass_dot, damping_dot, time_step_dot, state_dot, state_count)) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "dissipative Cayley JVP has incompatible increments")
            return
        end if
        allocate(matrix_a(state_count, state_count), &
            matrix_b(state_count, state_count), rhs(state_count), &
            matrix_a_dot(state_count, state_count), &
            matrix_b_dot(state_count, state_count), rhs_dot(state_count), &
            state_next(state_count))
        matrix_a = mass + 0.5_dp*time_step*damping
        matrix_b = mass - 0.5_dp*time_step*damping
        matrix_a_dot = mass_dot + 0.5_dp*(time_step_dot*damping + &
            time_step*damping_dot)
        matrix_b_dot = mass_dot - 0.5_dp*(time_step_dot*damping + &
            time_step*damping_dot)
        rhs = matmul(matrix_b, state)
        call dense_solve(matrix_a, rhs, state_next, info)
        if (info /= 0) then
            call status_set(status, FORTSPARSE_SINGULAR, &
                "dissipative Cayley JVP base system is singular")
            return
        end if
        rhs_dot = matmul(matrix_b_dot, state) + matmul(matrix_b, state_dot) - &
            matmul(matrix_a_dot, state_next)
        call dense_solve(matrix_a, rhs_dot, state_next_dot, info)
        if (info /= 0) then
            call status_set(status, FORTSPARSE_SINGULAR, &
                "dissipative Cayley JVP system is singular")
            return
        end if
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine advance_dissipative_cayley_jvp

    subroutine advance_dissipative_cayley_vjp( &
            mass, damping, time_step, state, state_next, state_next_bar, &
            mass_bar, damping_bar, time_step_bar, state_bar, status)
        real(dp), intent(in) :: mass(:, :), damping(:, :), time_step, state(:)
        real(dp), intent(in) :: state_next(:), state_next_bar(:)
        real(dp), intent(out) :: mass_bar(:, :), damping_bar(:, :)
        real(dp), intent(out) :: time_step_bar, state_bar(:)
        type(fortsparse_status_t), intent(out) :: status

        integer :: state_count, info
        real(dp), allocatable :: matrix_a(:, :), matrix_b(:, :), rhs(:)
        real(dp), allocatable :: rhs_bar(:), matrix_a_bar(:, :), matrix_b_bar(:, :)

        mass_bar = 0.0_dp
        damping_bar = 0.0_dp
        state_bar = 0.0_dp
        time_step_bar = 0.0_dp
        if (.not. validate_inputs( &
            mass, damping, time_step, state, state_next, status)) return
        if (size(state_next_bar) /= size(state) .or. &
            size(mass_bar, 1) /= size(state) .or. &
            size(mass_bar, 2) /= size(state) .or. &
            size(damping_bar, 1) /= size(state) .or. &
            size(damping_bar, 2) /= size(state) .or. &
            size(state_bar) /= size(state) .or. &
            any(.not. ieee_is_finite(state_next_bar))) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "dissipative Cayley VJP has incompatible cotangents")
            return
        end if
        state_count = size(state)
        allocate(matrix_a(state_count, state_count), &
            matrix_b(state_count, state_count), rhs(state_count), &
            rhs_bar(state_count), matrix_a_bar(state_count, state_count), &
            matrix_b_bar(state_count, state_count))
        matrix_a = mass + 0.5_dp*time_step*damping
        matrix_b = mass - 0.5_dp*time_step*damping
        call dense_solve(transpose(matrix_a), state_next_bar, rhs_bar, info)
        if (info /= 0) then
            call status_set(status, FORTSPARSE_SINGULAR, &
                "dissipative Cayley VJP transpose system is singular")
            return
        end if
        matrix_b_bar = outer_product(rhs_bar, state)
        matrix_a_bar = -outer_product(rhs_bar, state_next)
        state_bar = matmul(transpose(matrix_b), rhs_bar)
        mass_bar = matrix_a_bar + matrix_b_bar
        damping_bar = 0.5_dp*time_step*(matrix_a_bar - matrix_b_bar)
        time_step_bar = 0.5_dp*sum((matrix_a_bar - matrix_b_bar)*damping)
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine advance_dissipative_cayley_vjp

    logical function validate_inputs( &
            mass, damping, time_step, state, state_next, status) result(valid)
        real(dp), intent(in) :: mass(:, :), damping(:, :), time_step, state(:)
        real(dp), intent(in) :: state_next(:)
        type(fortsparse_status_t), intent(out) :: status
        integer :: state_count

        valid = .false.
        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            "dissipative Cayley has incompatible arrays")
        state_count = size(state)
        if (state_count < 1 .or. size(mass, 1) /= state_count .or. &
            size(mass, 2) /= state_count .or. size(damping, 1) /= state_count .or. &
            size(damping, 2) /= state_count .or. size(state_next) /= state_count .or. &
            .not. ieee_is_finite(time_step) .or. time_step < 0.0_dp) return
        if (any(.not. ieee_is_finite(mass)) .or. &
            any(.not. ieee_is_finite(damping)) .or. &
            any(.not. ieee_is_finite(state)) .or. &
            any(.not. ieee_is_finite(state_next))) return
        valid = .true.
        call status_set(status, FORTSPARSE_OK, "")
    end function validate_inputs

    pure logical function validate_direction( &
            mass_dot, damping_dot, time_step_dot, state_dot, state_count)
        real(dp), intent(in) :: mass_dot(:, :), damping_dot(:, :)
        real(dp), intent(in) :: time_step_dot, state_dot(:)
        integer, intent(in) :: state_count

        validate_direction = size(mass_dot, 1) == state_count .and. &
            size(mass_dot, 2) == state_count .and. &
            size(damping_dot, 1) == state_count .and. &
            size(damping_dot, 2) == state_count .and. &
            size(state_dot) == state_count .and. &
            ieee_is_finite(time_step_dot) .and. &
            all(ieee_is_finite(mass_dot)) .and. &
            all(ieee_is_finite(damping_dot)) .and. &
            all(ieee_is_finite(state_dot))
    end function validate_direction

    pure function outer_product(left, right) result(product)
        real(dp), intent(in) :: left(:), right(:)
        real(dp) :: product(size(left), size(right))

        product = spread(left, dim=2, ncopies=size(right))* &
            spread(right, dim=1, ncopies=size(left))
    end function outer_product

end module fortfem_dissipative_cayley
