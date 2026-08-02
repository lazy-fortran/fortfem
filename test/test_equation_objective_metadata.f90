program test_equation_objective_metadata
    use check, only: check_condition, check_summary
    use fortfem_equation_objective_metadata, only: &
        equation_objective_metadata_t, &
        initialize_equation_objective_metadata, &
        validate_equation_objective_metadata, &
        OBJECTIVE_METADATA_KIND_CONSTRAINT, &
        OBJECTIVE_METADATA_KIND_OBJECTIVE
    use fortfem_kinds, only: dp
    use fortsparse, only: fortsparse_status_t
    implicit none

    real(dp), parameter :: target(2) = [1.0_dp, -2.0_dp]
    real(dp), parameter :: lower_bound(2) = [-3.0_dp, -4.0_dp]
    real(dp), parameter :: upper_bound(2) = [3.0_dp, 4.0_dp]
    real(dp), parameter :: weight(2) = [2.0_dp, 0.5_dp]
    real(dp), parameter :: scale(2) = [10.0_dp, 20.0_dp]
    type(equation_objective_metadata_t) :: metadata, copied, invalid
    type(fortsparse_status_t) :: status
    logical :: all_passed

    all_passed = .true.
    call initialize_equation_objective_metadata(metadata, "flux", &
        OBJECTIVE_METADATA_KIND_OBJECTIVE, target, lower_bound, upper_bound, &
        weight, scale, .true., .false., 7, 0, status, units="Wb", &
        provenance="manufactured oracle", has_parameter_tangent=.true., &
        parameter_tangent_id=11, has_hvp=.true., hvp_id=13)
    call record_condition(status%code == 0 .and. &
        validate_equation_objective_metadata(metadata, status), &
        "typed objective metadata initializes and validates")
    call record_condition(metadata%value_count == 2 .and. &
        metadata%name == "flux" .and. metadata%units == "Wb" .and. &
        metadata%provenance == "manufactured oracle" .and. metadata%active .and. &
        .not. metadata%fixed .and. metadata%kkt_id == 7 .and. &
        metadata%parameter_tangent_id == 11 .and. metadata%hvp_id == 13, &
        "metadata retains numeric and derivative capability fields")

    copied = metadata
    copied%target(1) = -1.0_dp
    copied%name = "copied"
    call record_condition(metadata%target(1) == target(1) .and. &
        metadata%name == "flux" .and. &
        validate_equation_objective_metadata(copied, status), &
        "metadata assignment makes a deterministic deep copy")

    call initialize_equation_objective_metadata(invalid, "charge", &
        OBJECTIVE_METADATA_KIND_CONSTRAINT, target, lower_bound, upper_bound, &
        [0.0_dp, 1.0_dp], scale, .true., .false., 0, 4, status)
    call record_condition(status%code /= 0 .and. &
        .not. validate_equation_objective_metadata(invalid, status), &
        "metadata rejects a non-positive row weight")

    invalid = metadata
    invalid%target(2) = 9.0_dp
    call record_condition( &
        .not. validate_equation_objective_metadata(invalid, status) .and. &
        status%code /= 0, "metadata rejects a target outside its bounds")

    invalid = metadata
    invalid%has_hvp = .false.
    call record_condition( &
        .not. validate_equation_objective_metadata(invalid, status) .and. &
        status%code /= 0, "metadata rejects an orphan HVP capability identifier")

    call check_summary("equation/objective typed metadata")
    if (.not. all_passed) error stop 1

contains

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_equation_objective_metadata
