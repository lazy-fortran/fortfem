program test_boundary_backend_registry
    !! Direct-facade oracle for every neutral exterior-operator backend tag.
    !!
    !! The backend names are part of the serialized boundary contract.  This
    !! client deliberately imports only ``fortfem_boundary`` and checks the
    !! registry independently of any numerical BEM, DtN, NESTOR, BIEST, or
    !! virtual-casing implementation.
    use check, only: check_condition, check_summary
    use fortfem_boundary, only: &
        BOUNDARY_OPERATOR_BACKEND_BEM, &
        BOUNDARY_OPERATOR_BACKEND_BIEST, &
        BOUNDARY_OPERATOR_BACKEND_DTN, &
        BOUNDARY_OPERATOR_BACKEND_FEM, &
        BOUNDARY_OPERATOR_BACKEND_NESTOR, &
        BOUNDARY_OPERATOR_BACKEND_PML, &
        BOUNDARY_OPERATOR_BACKEND_USER, &
        BOUNDARY_OPERATOR_BACKEND_VIRTUAL_CASING, &
        BOUNDARY_OPERATOR_TRACE_CHANNEL_MIXED, &
        boundary_operator_contract_t, &
        initialize_boundary_operator_contract, &
        initialize_boundary_operator_trace_metadata, &
        validate_boundary_operator_contract
    implicit none

    integer, parameter :: backend_count = 8
    integer, parameter :: backend_kind(backend_count) = [ &
        BOUNDARY_OPERATOR_BACKEND_FEM, BOUNDARY_OPERATOR_BACKEND_BEM, &
        BOUNDARY_OPERATOR_BACKEND_DTN, BOUNDARY_OPERATOR_BACKEND_PML, &
        BOUNDARY_OPERATOR_BACKEND_NESTOR, BOUNDARY_OPERATOR_BACKEND_BIEST, &
        BOUNDARY_OPERATOR_BACKEND_VIRTUAL_CASING, BOUNDARY_OPERATOR_BACKEND_USER]
    character(len=32), parameter :: expected_name(backend_count) = [ &
        character(len=32) :: "FEM", "BEM", "DtN", "PML", "NESTOR", "BIEST", &
        "virtual-casing", "user"]
    type(boundary_operator_contract_t) :: contract, copy, invalid
    integer :: backend, status

    do backend = 1, backend_count
        call initialize_boundary_operator_contract( &
            contract, backend_kind(backend), "maxwell", "Hcurl-trace", 5, 7, &
            .true., .true., .true., .true., .true., .true., "SI", &
            "work-normalized", "backend-registry oracle", "torus-fixed-1", status)
        call check_condition(status == 0 .and. &
            validate_boundary_operator_contract(contract, status) .and. &
            contract%backend_kind == backend_kind(backend) .and. &
            trim(contract%backend_name) == trim(expected_name(backend)) .and. &
            contract%row_count == 5 .and. contract%column_count == 7, &
            "boundary facade preserves every backend registry entry")

        call initialize_boundary_operator_trace_metadata( &
            contract, BOUNDARY_OPERATOR_TRACE_CHANNEL_MIXED, &
            "complex-work", status)
        call check_condition(status == 0 .and. &
            validate_boundary_operator_contract(contract, status) .and. &
            contract%trace_channel == BOUNDARY_OPERATOR_TRACE_CHANNEL_MIXED .and. &
            trim(contract%work_pairing) == "complex-work", &
            "backend registry accepts the mixed work-conjugate trace")

        copy = contract
        call check_condition(validate_boundary_operator_contract(copy, status) .and. &
            trim(copy%backend_name) == trim(expected_name(backend)) .and. &
            copy%trace_channel == contract%trace_channel, &
            "backend metadata has independent value-copy semantics")
    end do

    call initialize_boundary_operator_contract(invalid, 0, "maxwell", "trace", &
        1, 1, .true., .true., .true., .true., .true., .true., "SI", "unit", &
        "invalid backend oracle", "invalid", status)
    call check_condition(status /= 0, &
        "boundary facade rejects an unknown backend registry entry")

    call check_summary("boundary backend registry")
end program test_boundary_backend_registry
