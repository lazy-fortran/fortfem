module fortfem_pseudo_arclength_residual
    !! Neutral pseudo-arclength continuation residual and derivatives.
    !!
    !! The caller supplies a nonlinear residual R(state, parameter) and a
    !! frozen tangent predictor.  This layer adds the scalar arclength row
    !! without selecting a nonlinear solver, a physical parameter, or a
    !! continuation policy.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortfem_kinds, only: dp
    use fortsparse, only: fortsparse_status_t, status_set, &
        FORTSPARSE_INVALID_MATRIX, FORTSPARSE_OK
    implicit none
    private

    public :: assemble_pseudo_arclength_residual
    public :: assemble_pseudo_arclength_residual_jvp
    public :: assemble_pseudo_arclength_residual_vjp

contains

    subroutine assemble_pseudo_arclength_residual( &
            residual, state, parameter, previous_state, previous_parameter, &
            tangent_state, tangent_parameter, step, augmented_residual, status)
        real(dp), intent(in) :: residual(:), state(:), parameter(:)
        real(dp), intent(in) :: previous_state(:), previous_parameter(:)
        real(dp), intent(in) :: tangent_state(:), tangent_parameter(:), step
        real(dp), intent(out) :: augmented_residual(:)
        type(fortsparse_status_t), intent(out) :: status
        integer :: residual_size

        augmented_residual = 0.0_dp
        if (.not. valid_inputs( &
            residual, state, parameter, previous_state, previous_parameter, tangent_state, &
            tangent_parameter, step, augmented_residual)) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "pseudo-arclength residual has incompatible inputs")
            return
        end if
        residual_size = size(residual)
        augmented_residual(:residual_size) = residual
        augmented_residual(residual_size + 1) = dot_product( &
            tangent_state, state - previous_state) + dot_product( &
            tangent_parameter, parameter - previous_parameter) - step
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine assemble_pseudo_arclength_residual

    subroutine assemble_pseudo_arclength_residual_jvp( &
            residual, state, parameter, previous_state, previous_parameter, tangent_state, &
            tangent_parameter, step, residual_dot, state_dot, parameter_dot, previous_state_dot, &
            previous_parameter_dot, tangent_state_dot, tangent_parameter_dot, step_dot, &
            augmented_residual_dot, status)
        real(dp), intent(in) :: residual(:), state(:), parameter(:)
        real(dp), intent(in) :: previous_state(:), previous_parameter(:)
        real(dp), intent(in) :: tangent_state(:), tangent_parameter(:), step
        real(dp), intent(in) :: residual_dot(:), state_dot(:), parameter_dot(:)
        real(dp), intent(in) :: previous_state_dot(:), previous_parameter_dot(:)
        real(dp), intent(in) :: tangent_state_dot(:), tangent_parameter_dot(:), step_dot
        real(dp), intent(out) :: augmented_residual_dot(:)
        type(fortsparse_status_t), intent(out) :: status
        integer :: residual_size

        augmented_residual_dot = 0.0_dp
        if (.not. valid_inputs( &
            residual, state, parameter, previous_state, previous_parameter, tangent_state, &
            tangent_parameter, step, augmented_residual_dot)) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "pseudo-arclength JVP has incompatible base inputs")
            return
        end if
        if (.not. valid_directions( &
            residual, state, parameter, previous_state, previous_parameter, tangent_state, &
            tangent_parameter, residual_dot, state_dot, parameter_dot, previous_state_dot, &
            previous_parameter_dot, tangent_state_dot, tangent_parameter_dot, step_dot)) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "pseudo-arclength JVP has incompatible directions")
            return
        end if
        residual_size = size(residual)
        augmented_residual_dot(:residual_size) = residual_dot
        augmented_residual_dot(residual_size + 1) = dot_product( &
            tangent_state_dot, state - previous_state) + dot_product( &
            tangent_state, state_dot - previous_state_dot) + dot_product( &
            tangent_parameter_dot, parameter - previous_parameter) + dot_product( &
            tangent_parameter, parameter_dot - previous_parameter_dot) - step_dot
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine assemble_pseudo_arclength_residual_jvp

    subroutine assemble_pseudo_arclength_residual_vjp( &
            residual, state, parameter, previous_state, previous_parameter, tangent_state, &
            tangent_parameter, step, augmented_residual_bar, residual_bar, state_bar, &
            parameter_bar, previous_state_bar, previous_parameter_bar, tangent_state_bar, &
            tangent_parameter_bar, step_bar, status)
        real(dp), intent(in) :: residual(:), state(:), parameter(:)
        real(dp), intent(in) :: previous_state(:), previous_parameter(:)
        real(dp), intent(in) :: tangent_state(:), tangent_parameter(:), step
        real(dp), intent(in) :: augmented_residual_bar(:)
        real(dp), intent(out) :: residual_bar(:), state_bar(:), parameter_bar(:)
        real(dp), intent(out) :: previous_state_bar(:), previous_parameter_bar(:)
        real(dp), intent(out) :: tangent_state_bar(:), tangent_parameter_bar(:), step_bar
        type(fortsparse_status_t), intent(out) :: status
        real(dp) :: arc_bar
        integer :: residual_size

        residual_bar = 0.0_dp
        state_bar = 0.0_dp
        parameter_bar = 0.0_dp
        previous_state_bar = 0.0_dp
        previous_parameter_bar = 0.0_dp
        tangent_state_bar = 0.0_dp
        tangent_parameter_bar = 0.0_dp
        step_bar = 0.0_dp
        if (.not. valid_inputs( &
            residual, state, parameter, previous_state, previous_parameter, tangent_state, &
            tangent_parameter, step, augmented_residual_bar)) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "pseudo-arclength VJP has incompatible base inputs")
            return
        end if
        if (.not. valid_adjoint_shapes( &
            residual, state, parameter, previous_state, previous_parameter, tangent_state, &
            tangent_parameter, residual_bar, state_bar, parameter_bar, previous_state_bar, &
            previous_parameter_bar, tangent_state_bar, tangent_parameter_bar)) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "pseudo-arclength VJP has incompatible cotangent outputs")
            return
        end if
        residual_size = size(residual)
        residual_bar = augmented_residual_bar(:residual_size)
        arc_bar = augmented_residual_bar(residual_size + 1)
        state_bar = arc_bar*tangent_state
        previous_state_bar = -state_bar
        parameter_bar = arc_bar*tangent_parameter
        previous_parameter_bar = -parameter_bar
        tangent_state_bar = arc_bar*(state - previous_state)
        tangent_parameter_bar = arc_bar*(parameter - previous_parameter)
        step_bar = -arc_bar
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine assemble_pseudo_arclength_residual_vjp

    logical function valid_inputs( &
            residual, state, parameter, previous_state, previous_parameter, tangent_state, &
            tangent_parameter, step, target) result(valid)
        real(dp), intent(in) :: residual(:), state(:), parameter(:)
        real(dp), intent(in) :: previous_state(:), previous_parameter(:)
        real(dp), intent(in) :: tangent_state(:), tangent_parameter(:), step
        real(dp), intent(in) :: target(:)

        valid = size(residual) > 0 .and. size(state) > 0 .and. size(parameter) > 0 .and. &
            size(previous_state) == size(state) .and. size(previous_parameter) == size(parameter) .and. &
            size(tangent_state) == size(state) .and. size(tangent_parameter) == size(parameter) .and. &
            size(target) == size(residual) + 1 .and. ieee_is_finite(step) .and. &
            all(ieee_is_finite(residual)) .and. all(ieee_is_finite(state)) .and. &
            all(ieee_is_finite(parameter)) .and. all(ieee_is_finite(previous_state)) .and. &
            all(ieee_is_finite(previous_parameter)) .and. all(ieee_is_finite(tangent_state)) .and. &
            all(ieee_is_finite(tangent_parameter)) .and. all(ieee_is_finite(target))
    end function valid_inputs

    logical function valid_directions( &
            residual, state, parameter, previous_state, previous_parameter, tangent_state, &
            tangent_parameter, residual_dot, state_dot, parameter_dot, previous_state_dot, &
            previous_parameter_dot, tangent_state_dot, tangent_parameter_dot, step_dot) result(valid)
        real(dp), intent(in) :: residual(:), state(:), parameter(:)
        real(dp), intent(in) :: previous_state(:), previous_parameter(:)
        real(dp), intent(in) :: tangent_state(:), tangent_parameter(:)
        real(dp), intent(in) :: residual_dot(:), state_dot(:), parameter_dot(:)
        real(dp), intent(in) :: previous_state_dot(:), previous_parameter_dot(:)
        real(dp), intent(in) :: tangent_state_dot(:), tangent_parameter_dot(:), step_dot

        valid = size(residual_dot) == size(residual) .and. size(state_dot) == size(state) .and. &
            size(parameter_dot) == size(parameter) .and. size(previous_state_dot) == size(previous_state) .and. &
            size(previous_parameter_dot) == size(previous_parameter) .and. &
            size(tangent_state_dot) == size(tangent_state) .and. &
            size(tangent_parameter_dot) == size(tangent_parameter) .and. ieee_is_finite(step_dot) .and. &
            all(ieee_is_finite(residual_dot)) .and. all(ieee_is_finite(state_dot)) .and. &
            all(ieee_is_finite(parameter_dot)) .and. all(ieee_is_finite(previous_state_dot)) .and. &
            all(ieee_is_finite(previous_parameter_dot)) .and. all(ieee_is_finite(tangent_state_dot)) .and. &
            all(ieee_is_finite(tangent_parameter_dot))
    end function valid_directions

    logical function valid_adjoint_shapes( &
            residual, state, parameter, previous_state, previous_parameter, tangent_state, &
            tangent_parameter, residual_bar, state_bar, parameter_bar, previous_state_bar, &
            previous_parameter_bar, tangent_state_bar, tangent_parameter_bar) result(valid)
        real(dp), intent(in) :: residual(:), state(:), parameter(:)
        real(dp), intent(in) :: previous_state(:), previous_parameter(:)
        real(dp), intent(in) :: tangent_state(:), tangent_parameter(:)
        real(dp), intent(in) :: residual_bar(:), state_bar(:), parameter_bar(:)
        real(dp), intent(in) :: previous_state_bar(:), previous_parameter_bar(:)
        real(dp), intent(in) :: tangent_state_bar(:), tangent_parameter_bar(:)

        valid = size(residual_bar) == size(residual) .and. size(state_bar) == size(state) .and. &
            size(parameter_bar) == size(parameter) .and. size(previous_state_bar) == size(previous_state) .and. &
            size(previous_parameter_bar) == size(previous_parameter) .and. &
            size(tangent_state_bar) == size(tangent_state) .and. &
            size(tangent_parameter_bar) == size(tangent_parameter)
    end function valid_adjoint_shapes

end module fortfem_pseudo_arclength_residual
