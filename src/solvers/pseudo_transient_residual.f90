module fortfem_pseudo_transient_residual
    !! Neutral pseudo-transient continuation residual and derivatives.
    !!
    !! Given a caller-owned nonlinear residual R(u), this adds the stabilizing
    !! term M (u-u_old) / dt.  The wrapper does not select a mass matrix,
    !! stopping rule, or time integrator; symplectic updates remain separate.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortfem_kinds, only: dp
    use fortsparse, only: fortsparse_status_t, status_set, &
        FORTSPARSE_INVALID_MATRIX, FORTSPARSE_OK
    implicit none
    private

    public :: assemble_pseudo_transient_residual
    public :: assemble_pseudo_transient_residual_jvp
    public :: assemble_pseudo_transient_residual_vjp

contains

    subroutine assemble_pseudo_transient_residual( &
            residual, mass, state, previous_state, time_step, augmented_residual, status)
        real(dp), intent(in) :: residual(:), mass(:, :), state(:), previous_state(:)
        real(dp), intent(in) :: time_step
        real(dp), intent(out) :: augmented_residual(:)
        type(fortsparse_status_t), intent(out) :: status

        augmented_residual = 0.0_dp
        if (.not. valid_base_inputs( &
            residual, mass, state, previous_state, time_step, augmented_residual)) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "pseudo-transient residual has incompatible inputs")
            return
        end if
        augmented_residual = residual + matmul(mass, state - previous_state)/time_step
        if (.not. all(ieee_is_finite(augmented_residual))) then
            augmented_residual = 0.0_dp
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "pseudo-transient residual is non-finite")
            return
        end if
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine assemble_pseudo_transient_residual

    subroutine assemble_pseudo_transient_residual_jvp( &
            residual, mass, state, previous_state, time_step, residual_dot, mass_dot, &
            state_dot, previous_state_dot, time_step_dot, augmented_residual_dot, status)
        real(dp), intent(in) :: residual(:), mass(:, :), state(:), previous_state(:)
        real(dp), intent(in) :: time_step
        real(dp), intent(in) :: residual_dot(:), mass_dot(:, :), state_dot(:)
        real(dp), intent(in) :: previous_state_dot(:), time_step_dot
        real(dp), intent(out) :: augmented_residual_dot(:)
        type(fortsparse_status_t), intent(out) :: status

        augmented_residual_dot = 0.0_dp
        if (.not. valid_base_inputs( &
            residual, mass, state, previous_state, time_step, augmented_residual_dot)) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "pseudo-transient JVP has incompatible base inputs")
            return
        end if
        if (.not. valid_direction_inputs( &
            residual, mass, state, previous_state, residual_dot, mass_dot, state_dot, &
            previous_state_dot, time_step_dot)) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "pseudo-transient JVP has incompatible directions")
            return
        end if
        augmented_residual_dot = residual_dot + matmul(mass_dot, state - previous_state)/time_step + &
            matmul(mass, state_dot - previous_state_dot)/time_step - &
            matmul(mass, state - previous_state)*time_step_dot/time_step**2
        if (.not. all(ieee_is_finite(augmented_residual_dot))) then
            augmented_residual_dot = 0.0_dp
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "pseudo-transient JVP is non-finite")
            return
        end if
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine assemble_pseudo_transient_residual_jvp

    subroutine assemble_pseudo_transient_residual_vjp( &
            residual, mass, state, previous_state, time_step, augmented_residual_bar, &
            residual_bar, mass_bar, state_bar, previous_state_bar, time_step_bar, status)
        real(dp), intent(in) :: residual(:), mass(:, :), state(:), previous_state(:)
        real(dp), intent(in) :: time_step, augmented_residual_bar(:)
        real(dp), intent(out) :: residual_bar(:), mass_bar(:, :), state_bar(:)
        real(dp), intent(out) :: previous_state_bar(:), time_step_bar
        type(fortsparse_status_t), intent(out) :: status
        integer :: row, column

        residual_bar = 0.0_dp
        mass_bar = 0.0_dp
        state_bar = 0.0_dp
        previous_state_bar = 0.0_dp
        time_step_bar = 0.0_dp
        if (.not. valid_base_inputs( &
            residual, mass, state, previous_state, time_step, augmented_residual_bar)) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "pseudo-transient VJP has incompatible base inputs")
            return
        end if
        if (.not. valid_vjp_shapes( &
            residual, mass, state, previous_state, residual_bar, mass_bar, state_bar, &
            previous_state_bar)) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "pseudo-transient VJP has incompatible cotangent shapes")
            return
        end if
        residual_bar = augmented_residual_bar
        state_bar = matmul(transpose(mass), augmented_residual_bar)/time_step
        previous_state_bar = -state_bar
        do row = 1, size(mass, 1)
            do column = 1, size(mass, 2)
                mass_bar(row, column) = augmented_residual_bar(row)* &
                    (state(column) - previous_state(column))/time_step
            end do
        end do
        time_step_bar = -dot_product(state_bar, state - previous_state)/time_step
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine assemble_pseudo_transient_residual_vjp

    logical function valid_base_inputs( &
            residual, mass, state, previous_state, time_step, target) result(valid)
        real(dp), intent(in) :: residual(:), mass(:, :), state(:), previous_state(:)
        real(dp), intent(in) :: time_step, target(:)

        valid = size(residual) > 0 .and. size(mass, 1) == size(residual) .and. &
            size(mass, 2) == size(residual) .and. size(state) == size(residual) .and. &
            size(previous_state) == size(residual) .and. size(target) == size(residual)
        if (.not. valid) return
        if (.not. ieee_is_finite(time_step) .or. time_step <= 0.0_dp) then
            valid = .false.
            return
        end if
        valid = all(ieee_is_finite(residual)) .and. all(ieee_is_finite(mass)) .and. &
            all(ieee_is_finite(state)) .and. all(ieee_is_finite(previous_state)) .and. &
            all(ieee_is_finite(target))
    end function valid_base_inputs

    logical function valid_direction_inputs( &
            residual, mass, state, previous_state, residual_dot, mass_dot, state_dot, &
            previous_state_dot, time_step_dot) result(valid)
        real(dp), intent(in) :: residual(:), mass(:, :), state(:), previous_state(:)
        real(dp), intent(in) :: residual_dot(:), mass_dot(:, :), state_dot(:)
        real(dp), intent(in) :: previous_state_dot(:), time_step_dot

        valid = size(residual_dot) == size(residual) .and. &
            all(shape(mass_dot) == shape(mass)) .and. size(state_dot) == size(state) .and. &
            size(previous_state_dot) == size(previous_state)
        if (.not. valid) return
        valid = ieee_is_finite(time_step_dot) .and. all(ieee_is_finite(residual_dot)) .and. &
            all(ieee_is_finite(mass_dot)) .and. all(ieee_is_finite(state_dot)) .and. &
            all(ieee_is_finite(previous_state_dot))
    end function valid_direction_inputs

    logical function valid_vjp_shapes( &
            residual, mass, state, previous_state, residual_bar, mass_bar, state_bar, &
            previous_state_bar) result(valid)
        real(dp), intent(in) :: residual(:), mass(:, :), state(:), previous_state(:)
        real(dp), intent(in) :: residual_bar(:), mass_bar(:, :), state_bar(:)
        real(dp), intent(in) :: previous_state_bar(:)

        valid = size(residual_bar) == size(residual) .and. &
            all(shape(mass_bar) == shape(mass)) .and. size(state_bar) == size(state) .and. &
            size(previous_state_bar) == size(previous_state)
    end function valid_vjp_shapes

end module fortfem_pseudo_transient_residual
