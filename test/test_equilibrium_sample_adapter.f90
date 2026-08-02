program test_equilibrium_sample_adapter
    use check, only: check_condition, check_summary
    use fortfem_api, only: &
        build_equilibrium_interchange_sample_set, &
        equilibrium_interchange_t, equilibrium_normalization_t, &
        initialize_equilibrium_interchange, interchange_sample_set_t, &
        validate_interchange_samples
    use fortfem_kinds, only: dp
    implicit none

    integer, parameter :: sample_count = 4
    real(dp) :: mapped(2, sample_count), physical(2, sample_count)
    real(dp) :: coefficient_values(2, sample_count)
    real(dp) :: profile_coordinates(2), profile_values(1, 2)
    real(dp) :: boundary_coordinates(2, 4), weights(sample_count)
    character(len=32) :: coefficient_names(2), profile_names(1), boundary_names(1)
    type(equilibrium_normalization_t) :: normalization
    type(equilibrium_interchange_t) :: interchange
    type(interchange_sample_set_t) :: samples
    integer :: status
    logical :: all_passed

    all_passed = .true.
    mapped = reshape([ &
        0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, &
        1.0_dp, 1.0_dp, 0.0_dp, 1.0_dp], shape(mapped))
    physical = reshape([ &
        2.0_dp, -1.0_dp, 3.0_dp, -1.0_dp, &
        3.0_dp, 0.0_dp, 2.0_dp, 0.0_dp], shape(physical))
    coefficient_values(1, :) = [1.0_dp, 2.0_dp, 4.0_dp, 3.0_dp]
    coefficient_values(2, :) = [0.5_dp, 0.25_dp, -0.25_dp, -0.5_dp]
    profile_coordinates = [0.0_dp, 1.0_dp]
    profile_values(1, :) = [2.0_dp, 1.0_dp]
    boundary_coordinates = physical
    weights = [1.0_dp, 2.0_dp, 1.5_dp, 0.5_dp]
    coefficient_names = [character(len=32) :: "pressure", "Bz"]
    profile_names = [character(len=32) :: "density"]
    boundary_names = [character(len=32) :: "plasma"]

    call initialize_equilibrium_interchange( &
        interchange, mapped, physical, coefficient_names, coefficient_values, &
        profile_coordinates, profile_names, profile_values, [1, 5], &
        boundary_names, boundary_coordinates, normalization, status)
    call record_condition(status == 0, &
        "sample adapter accepts a valid equilibrium interchange record")

    call build_equilibrium_interchange_sample_set( &
        interchange, [1, 2], weights, "manufactured-adapter", &
        "analytic-physical-grid", samples, status)
    call record_condition(status == 0, &
        "sample adapter projects selected coefficients")
    call record_condition(validate_interchange_samples(samples, status) .and. &
        status == 0, "projected sample set validates")
    call record_condition(samples%spatial_dimension == 2 .and. &
        samples%component_count == 2 .and. samples%sample_count == sample_count, &
        "projected sample metadata matches the physical grid")
    call record_condition(all(samples%coordinates == physical) .and. &
        all(samples%weights == weights), &
        "projected samples retain physical coordinates and weights")
    call record_condition(all(samples%values == coefficient_values), &
        "selected component values match the independent coefficient oracle")
    call record_condition(samples%producer == "manufactured-adapter" .and. &
        samples%provenance == "analytic-physical-grid", &
        "sample provenance is explicit and retained")

    call build_equilibrium_interchange_sample_set( &
        interchange, [2], weights, "manufactured-adapter", &
        "single-component-oracle", samples, status)
    call record_condition(status == 0 .and. samples%component_count == 1 .and. &
        all(samples%values(1, :) == coefficient_values(2, :)), &
        "single-component projection preserves component ordering")

    call build_equilibrium_interchange_sample_set( &
        interchange, [1, 1], weights, "manufactured-adapter", &
        "duplicate-selector", samples, status)
    call record_condition(status /= 0, "duplicate component selectors are rejected")
    call build_equilibrium_interchange_sample_set( &
        interchange, [3], weights, "manufactured-adapter", &
        "out-of-range-selector", samples, status)
    call record_condition(status /= 0, "out-of-range component selectors are rejected")

    call check_summary("equilibrium sample adapter")
    if (.not. all_passed) error stop 1

contains

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_equilibrium_sample_adapter
