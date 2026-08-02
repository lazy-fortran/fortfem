program test_boundary_operator_contract
    use, intrinsic :: iso_fortran_env, only: real64
    use check, only: check_condition, check_summary
    use fortfem_api, only: &
        BOUNDARY_OPERATOR_BACKEND_BEM, &
        BOUNDARY_OPERATOR_BACKEND_DTN, &
        boundary_operator_contract_t, &
        initialize_boundary_operator_contract, &
        validate_boundary_operator_contract
    implicit none

    integer, parameter :: dp = real64
    type(boundary_operator_contract_t) :: contract, copy, invalid
    integer :: status
    logical :: all_passed

    all_passed = .true.
    call initialize_boundary_operator_contract(contract, &
        BOUNDARY_OPERATOR_BACKEND_BEM, "helmholtz", "H1-trace", 12, 8, &
        .true., .true., .true., .true., .true., .true., "m", "unit-normalized", &
        "manufactured circle oracle; revision test-1", "circle-fixed-1", status)
    call record_condition(status == 0, "boundary operator contract initializes")
    call record_condition(validate_boundary_operator_contract(contract, status), &
        "boundary operator contract validates")
    call record_condition(contract%backend_kind == BOUNDARY_OPERATOR_BACKEND_BEM .and. &
        contract%row_count == 12 .and. contract%column_count == 8 .and. &
        contract%complex_valued .and. contract%matrix_free .and. contract%assembled .and. &
        contract%has_jvp .and. contract%has_vjp .and. contract%has_residual .and. &
        trim(contract%equation) == "helmholtz" .and. &
        trim(contract%space) == "H1-trace" .and. &
        trim(contract%provenance) == "manufactured circle oracle; revision test-1", &
        "boundary operator metadata matches the supplied provenance")

    copy = contract
    call record_condition(validate_boundary_operator_contract(copy, status) .and. &
        copy%backend_kind == contract%backend_kind .and. &
        copy%row_count == contract%row_count .and. &
        copy%provenance == contract%provenance, &
        "boundary operator contract has value-copy semantics")

    call initialize_boundary_operator_contract(invalid, &
        BOUNDARY_OPERATOR_BACKEND_DTN, "laplace", "trace", 4, 4, .false., .false., &
        .true., .true., .true., .true., "1", "", "missing topology", "circle-fixed-1", status)
    call record_condition(status /= 0, &
        "boundary operator contract rejects an unavailable action")

    call initialize_boundary_operator_contract(invalid, &
        BOUNDARY_OPERATOR_BACKEND_DTN, "laplace", "trace", 4, 4, .false., .true., &
        .true., .true., .true., .true., "1", "unit", "", "circle-fixed-1", status)
    call record_condition(status /= 0, &
        "boundary operator contract rejects missing provenance")

    call check_summary("boundary operator contract")
    if (.not. all_passed) error stop 1

contains

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_boundary_operator_contract
