program test_continuation_event
    use check, only: check_condition, check_summary
    use fortfem_feec, only: CONTINUATION_EVENT_NONE
    use fortfem_time, only: &
        CONTINUATION_EVENT_NEAR_ZERO, CONTINUATION_EVENT_SIGN_CROSSING, &
        classify_continuation_event
    use fortfem_kinds, only: dp
    use fortsparse, only: FORTSPARSE_INVALID_MATRIX, FORTSPARSE_OK, &
        fortsparse_status_t
    implicit none

    logical :: all_passed

    all_passed = .true.
    call check_contract()
    call check_summary("Continuation event classifier")
    if (.not. all_passed) error stop 1

contains

    subroutine check_contract()
        real(dp) :: previous_margin(3), current_margin(3), minimum_margin
        integer :: event_code, event_index
        type(fortsparse_status_t) :: status

        previous_margin = [0.4_dp, -0.3_dp, 0.8_dp]
        current_margin = [0.2_dp, 0.1_dp, 0.7_dp]
        call classify_continuation_event( &
            previous_margin, current_margin, 0.05_dp, event_code, event_index, &
            minimum_margin, status)
        call record(status%code == FORTSPARSE_OK .and. &
            event_code == CONTINUATION_EVENT_SIGN_CROSSING .and. event_index == 2 .and. &
            abs(minimum_margin - 0.1_dp) < 1.0e-14_dp, &
            "signed crossing oracle")

        previous_margin = [0.4_dp, 0.03_dp, 0.8_dp]
        current_margin = [0.2_dp, 0.01_dp, 0.7_dp]
        call classify_continuation_event( &
            previous_margin, current_margin, 0.02_dp, event_code, event_index, &
            minimum_margin, status)
        call record(status%code == FORTSPARSE_OK .and. &
            event_code == CONTINUATION_EVENT_NEAR_ZERO .and. event_index == 2 .and. &
            abs(minimum_margin - 0.01_dp) < 1.0e-14_dp, &
            "near-zero margin oracle")

        previous_margin = [0.4_dp, -0.3_dp, 0.8_dp]
        current_margin = [0.2_dp, -0.1_dp, 0.7_dp]
        call classify_continuation_event( &
            previous_margin, current_margin, 0.05_dp, event_code, event_index, &
            minimum_margin, status)
        call record(status%code == FORTSPARSE_OK .and. &
            event_code == CONTINUATION_EVENT_NONE .and. event_index == 0 .and. &
            abs(minimum_margin - 0.1_dp) < 1.0e-14_dp, &
            "no-event oracle")

        call classify_continuation_event( &
            previous_margin, current_margin, -1.0_dp, event_code, event_index, &
            minimum_margin, status)
        call record(status%code == FORTSPARSE_INVALID_MATRIX, &
            "negative event tolerance is rejected")
    end subroutine check_contract

    subroutine record(condition, message)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: message

        call check_condition(condition, message)
        all_passed = all_passed .and. condition
    end subroutine record

end program test_continuation_event
