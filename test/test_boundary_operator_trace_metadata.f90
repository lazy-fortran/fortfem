program test_boundary_operator_trace_metadata
    use check, only: check_condition, check_summary
    use fortfem_boundary, only: &
        BOUNDARY_OPERATOR_BACKEND_BEM, &
        BOUNDARY_OPERATOR_TRACE_CHANNEL_NORMAL, &
        BOUNDARY_OPERATOR_TRACE_CHANNEL_TANGENTIAL, &
        boundary_operator_contract_t, &
        initialize_boundary_operator_contract, &
        initialize_boundary_operator_trace_metadata, &
        validate_boundary_operator_contract
    implicit none

    type(boundary_operator_contract_t) :: contract, invalid
    integer :: status
    logical :: all_passed

    all_passed = .true.
    call initialize_boundary_operator_contract(contract, &
        BOUNDARY_OPERATOR_BACKEND_BEM, "maxwell", "Hcurl-tangential", 6, 6, &
        .true., .true., .true., .true., .true., .true., "V/m", "SI", &
        "trace-channel oracle", "torus-fixed-1", status)
    call record_condition(status == 0, "base boundary contract initializes")
    call record_condition(contract%trace_channel == 1 .and. &
        trim(contract%work_pairing) == "scalar-l2", &
        "legacy initialization supplies neutral scalar metadata")

    call initialize_boundary_operator_trace_metadata(contract, &
        BOUNDARY_OPERATOR_TRACE_CHANNEL_TANGENTIAL, "surface-tangential-work", status)
    call record_condition(status == 0 .and. &
        validate_boundary_operator_contract(contract, status), &
        "trace metadata initializes and validates")
    call record_condition( &
        contract%trace_channel == BOUNDARY_OPERATOR_TRACE_CHANNEL_TANGENTIAL .and. &
        trim(contract%work_pairing) == "surface-tangential-work", &
        "trace metadata records channel and work pairing")

    invalid = contract
    call initialize_boundary_operator_trace_metadata( &
        invalid, BOUNDARY_OPERATOR_TRACE_CHANNEL_NORMAL, "surface-normal-work", status)
    call record_condition(status == 0 .and. &
        validate_boundary_operator_contract(invalid, status) .and. &
        invalid%trace_channel == BOUNDARY_OPERATOR_TRACE_CHANNEL_NORMAL .and. &
        invalid%trace_channel /= contract%trace_channel, &
        "normal and tangential ports remain distinguishable")

    invalid%trace_channel = BOUNDARY_OPERATOR_TRACE_CHANNEL_NORMAL
    invalid%work_pairing = ""
    call record_condition( &
        .not. validate_boundary_operator_contract(invalid, status) .and. &
        status /= 0, "validation rejects a missing work pairing label")

    call initialize_boundary_operator_trace_metadata( &
        contract, 99, "surface-work", status)
    call record_condition(status /= 0 .and. &
        contract%trace_channel == BOUNDARY_OPERATOR_TRACE_CHANNEL_TANGENTIAL, &
        "trace metadata rejects an unknown channel without mutation")

    call check_summary("boundary operator trace metadata")
    if (.not. all_passed) error stop 1

contains

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_boundary_operator_trace_metadata
