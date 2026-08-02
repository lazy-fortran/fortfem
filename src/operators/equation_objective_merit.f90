module fortfem_equation_objective_merit
    !! Neutral normalized weighted merit for one metadata block.
    !!
    !! The metadata record owns targets, scales, weights, and activation.  This
    !! module supplies only the algebraic merit
    !!
    !!   1/2 sum_i active * w_i ((v_i - target_i)/scale_i)^2,
    !!
    !! together with exact value/JVP/VJP actions.  Bounds, constraints,
    !! profiles, and optimizer decisions remain caller-owned.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortfem_kinds, only: dp
    use fortfem_equation_objective_metadata, only: &
        equation_objective_metadata_t, validate_equation_objective_metadata
    use fortsparse, only: FORTSPARSE_INVALID_MATRIX, FORTSPARSE_OK, &
        fortsparse_status_t, status_set
    implicit none
    private

    public :: evaluate_equation_objective_merit
    public :: evaluate_equation_objective_merit_jvp
    public :: evaluate_equation_objective_merit_vjp

contains

    subroutine evaluate_equation_objective_merit( &
            metadata, values, merit, status)
        type(equation_objective_metadata_t), intent(in) :: metadata
        real(dp), intent(in) :: values(:)
        real(dp), intent(out) :: merit
        type(fortsparse_status_t), intent(out) :: status
        integer :: row
        real(dp) :: normalized

        merit = 0.0_dp
        if (.not. valid_values(metadata, values, status)) return
        if (.not. metadata%active) then
            call status_set(status, FORTSPARSE_OK, "")
            return
        end if
        do row = 1, metadata%value_count
            normalized = (values(row) - metadata%target(row))/ &
                metadata%scale(row)
            merit = merit + 0.5_dp*metadata%weight(row)*normalized**2
        end do
        if (.not. ieee_is_finite(merit)) then
            merit = 0.0_dp
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "equation/objective merit is non-finite")
            return
        end if
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine evaluate_equation_objective_merit

    subroutine evaluate_equation_objective_merit_jvp( &
            metadata, values, values_dot, merit_dot, status)
        type(equation_objective_metadata_t), intent(in) :: metadata
        real(dp), intent(in) :: values(:), values_dot(:)
        real(dp), intent(out) :: merit_dot
        type(fortsparse_status_t), intent(out) :: status
        integer :: row
        real(dp) :: normalized

        merit_dot = 0.0_dp
        if (.not. valid_values(metadata, values, status)) return
        if (size(values_dot) /= metadata%value_count .or. &
            .not. all(ieee_is_finite(values_dot))) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "equation/objective merit JVP has invalid values")
            return
        end if
        if (.not. metadata%active) then
            call status_set(status, FORTSPARSE_OK, "")
            return
        end if
        do row = 1, metadata%value_count
            normalized = (values(row) - metadata%target(row))/ &
                metadata%scale(row)
            merit_dot = merit_dot + metadata%weight(row)*normalized* &
                values_dot(row)/metadata%scale(row)
        end do
        if (.not. ieee_is_finite(merit_dot)) then
            merit_dot = 0.0_dp
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "equation/objective merit JVP is non-finite")
            return
        end if
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine evaluate_equation_objective_merit_jvp

    subroutine evaluate_equation_objective_merit_vjp( &
            metadata, values, merit_bar, values_bar, status)
        type(equation_objective_metadata_t), intent(in) :: metadata
        real(dp), intent(in) :: values(:), merit_bar
        real(dp), intent(out) :: values_bar(:)
        type(fortsparse_status_t), intent(out) :: status
        integer :: row
        real(dp) :: normalized

        values_bar = 0.0_dp
        if (.not. valid_values(metadata, values, status)) return
        if (size(values_bar) /= metadata%value_count .or. &
            .not. ieee_is_finite(merit_bar)) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "equation/objective merit VJP has invalid values")
            return
        end if
        if (.not. metadata%active) then
            call status_set(status, FORTSPARSE_OK, "")
            return
        end if
        do row = 1, metadata%value_count
            normalized = (values(row) - metadata%target(row))/ &
                metadata%scale(row)
            values_bar(row) = merit_bar*metadata%weight(row)*normalized/ &
                metadata%scale(row)
        end do
        if (.not. all(ieee_is_finite(values_bar))) then
            values_bar = 0.0_dp
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "equation/objective merit VJP is non-finite")
            return
        end if
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine evaluate_equation_objective_merit_vjp

    logical function valid_values(metadata, values, status) result(valid)
        type(equation_objective_metadata_t), intent(in) :: metadata
        real(dp), intent(in) :: values(:)
        type(fortsparse_status_t), intent(out) :: status

        valid = .false.
        if (.not. validate_equation_objective_metadata(metadata, status)) return
        if (size(values) /= metadata%value_count .or. &
            .not. all(ieee_is_finite(values))) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "equation/objective merit received invalid values")
            return
        end if
        valid = .true.
    end function valid_values

end module fortfem_equation_objective_merit
