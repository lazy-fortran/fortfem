program test_cmake_boundary_consumer
    !! CMake consumer smoke for the boundary facade: the same `use
    !! fortfem_boundary` that fails when the CMake build does not compile and
    !! export the public Fortran API modules.
    use check, only: check_condition, check_summary
    use fortfem_boundary, only: &
        BOUNDARY_OPERATOR_BACKEND_DTN, &
        BOUNDARY_OPERATOR_TRACE_CHANNEL_NORMAL, &
        boundary_operator_contract_t, &
        initialize_boundary_operator_contract, &
        initialize_boundary_operator_trace_metadata, &
        validate_boundary_operator_contract
    implicit none

    type(boundary_operator_contract_t) :: contract
    integer :: status
    logical :: all_passed

    all_passed = .true.
    call initialize_boundary_operator_contract(contract, &
        BOUNDARY_OPERATOR_BACKEND_DTN, "elastic-acoustic", "H1-trace", 6, 6, &
        .true., .false., .true., .true., .true., .true., "m", "unit-normalized", &
        "cmake consumer smoke; planar dtn interface", "planar-fixed-1", status)
    call record_condition(status == 0, "boundary facade contract initializes")
    call initialize_boundary_operator_trace_metadata(contract, &
        BOUNDARY_OPERATOR_TRACE_CHANNEL_NORMAL, "Re<int u, v>", status)
    call record_condition(status == 0, "boundary facade accepts normal traces")
    call record_condition(validate_boundary_operator_contract(contract, status), &
        "boundary facade contract validates")
    call check_summary("CMake boundary facade consumer")
    if (.not. all_passed) error stop 1

contains

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: description

        if (.not. condition) all_passed = .false.
        call check_condition(condition, description)
    end subroutine record_condition

end program test_cmake_boundary_consumer
