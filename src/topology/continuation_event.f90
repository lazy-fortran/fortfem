module fortfem_continuation_event
    !! Deterministic signed-margin diagnostics for continuation events.
    !!
    !! A caller supplies one signed margin per fixed-topology constraint.  A
    !! sign crossing is a discrete event; near-zero margins are a warning that
    !! a cut, separatrix, resonance, or other branch may be about to change.
    !! The classifier reports the event only and never differentiates a branch.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortfem_kinds, only: dp
    use fortsparse, only: fortsparse_status_t, status_set, &
        FORTSPARSE_INVALID_MATRIX, FORTSPARSE_OK
    implicit none
    private

    integer, parameter, public :: CONTINUATION_EVENT_NONE = 0
    integer, parameter, public :: CONTINUATION_EVENT_SIGN_CROSSING = 1
    integer, parameter, public :: CONTINUATION_EVENT_NEAR_ZERO = 2

    public :: classify_continuation_event

contains

    subroutine classify_continuation_event( &
            previous_margin, current_margin, tolerance, event_code, event_index, &
            minimum_margin, status)
        real(dp), intent(in) :: previous_margin(:), current_margin(:), tolerance
        integer, intent(out) :: event_code, event_index
        real(dp), intent(out) :: minimum_margin
        type(fortsparse_status_t), intent(out) :: status
        integer :: index
        real(dp) :: candidate

        event_code = CONTINUATION_EVENT_NONE
        event_index = 0
        minimum_margin = 0.0_dp
        if (.not. valid_inputs(previous_margin, current_margin, tolerance)) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "continuation event margins are incompatible")
            return
        end if

        minimum_margin = huge(1.0_dp)
        do index = 1, size(previous_margin)
            candidate = min(abs(previous_margin(index)), abs(current_margin(index)))
            minimum_margin = min(minimum_margin, candidate)
            if (event_code == CONTINUATION_EVENT_NONE .and. &
                (previous_margin(index) == 0.0_dp .or. current_margin(index) == 0.0_dp)) then
                event_code = CONTINUATION_EVENT_SIGN_CROSSING
                event_index = index
            end if
            if (event_code == CONTINUATION_EVENT_NONE .and. &
                ((previous_margin(index) < 0.0_dp .and. current_margin(index) > 0.0_dp) .or. &
                (previous_margin(index) > 0.0_dp .and. current_margin(index) < 0.0_dp))) then
                event_code = CONTINUATION_EVENT_SIGN_CROSSING
                event_index = index
            end if
        end do
        if (event_code == CONTINUATION_EVENT_NONE) then
            do index = 1, size(previous_margin)
                candidate = min(abs(previous_margin(index)), abs(current_margin(index)))
                if (candidate <= tolerance) then
                    event_code = CONTINUATION_EVENT_NEAR_ZERO
                    event_index = index
                    exit
                end if
            end do
        end if
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine classify_continuation_event

    logical function valid_inputs(previous_margin, current_margin, tolerance) result(valid)
        real(dp), intent(in) :: previous_margin(:), current_margin(:), tolerance

        valid = size(previous_margin) > 0 .and. size(current_margin) == size(previous_margin)
        if (.not. valid) return
        if (.not. ieee_is_finite(tolerance) .or. tolerance < 0.0_dp) then
            valid = .false.
            return
        end if
        valid = all(ieee_is_finite(previous_margin)) .and. &
            all(ieee_is_finite(current_margin))
    end function valid_inputs

end module fortfem_continuation_event
