program test_equilibrium_interchange
    use, intrinsic :: iso_fortran_env, only: real64
    use, intrinsic :: ieee_arithmetic, only: ieee_quiet_nan, ieee_value
    use check, only: check_condition, check_summary
    use fortfem_api, only: &
        equilibrium_interchange_t, equilibrium_normalization_t, &
        initialize_equilibrium_interchange, validate_equilibrium_interchange
    implicit none

    integer, parameter :: dp = real64
    integer, parameter :: sample_count = 4, profile_sample_count = 3
    integer, parameter :: coefficient_count = 2, profile_count = 2
    integer, parameter :: boundary_point_count = 6, boundary_count = 2
    real(dp) :: mapped_coordinates(3, sample_count)
    real(dp) :: physical_coordinates(3, sample_count)
    real(dp) :: coefficient_values(coefficient_count, sample_count)
    real(dp) :: profile_coordinates(profile_sample_count)
    real(dp) :: profile_values(profile_count, profile_sample_count)
    real(dp) :: boundary_coordinates(3, boundary_point_count)
    integer :: boundary_offsets(boundary_count + 1)
    character(len=32) :: coefficient_names(coefficient_count)
    character(len=32) :: profile_names(profile_count)
    character(len=32) :: boundary_names(boundary_count)
    type(equilibrium_normalization_t) :: normalization
    type(equilibrium_interchange_t) :: data, copy, invalid
    integer :: status
    logical :: all_passed

    all_passed = .true.
    call build_manufactured_data( &
        mapped_coordinates, physical_coordinates, coefficient_values, &
        profile_coordinates, profile_values, boundary_coordinates, &
        boundary_offsets, coefficient_names, profile_names, boundary_names)
    normalization%length_unit = "m"
    normalization%magnetic_field_unit = "T"
    normalization%pressure_unit = "Pa"
    normalization%current_unit = "A/m2"
    normalization%length_scale = 2.0_dp
    normalization%magnetic_field_scale = 3.0_dp
    normalization%pressure_scale = 4.0_dp
    normalization%current_scale = 5.0_dp

    call initialize_equilibrium_interchange( &
        data, mapped_coordinates, physical_coordinates, coefficient_names, &
        coefficient_values, profile_coordinates, profile_names, profile_values, &
        boundary_offsets, boundary_names, boundary_coordinates, normalization, &
        status)
    call record_condition(status == 0, "manufactured interchange data initializes")
    call record_condition(validate_equilibrium_interchange(data, status), &
        "initialized interchange data validates")
    call record_condition(status == 0, "validation status is successful")
    call record_condition(data%spatial_dimension == 3, "spatial dimension is retained")
    call record_condition(data%sample_count == sample_count, "sample count is retained")
    call record_condition(data%coefficient_count == coefficient_count, &
        "coefficient count is retained")
    call record_condition(data%profile_count == profile_count, &
        "profile count is retained")
    call record_condition(data%boundary_count == boundary_count, &
        "boundary count is retained")
    call record_condition(all(abs(data%physical_coordinates - &
        physical_coordinates) < 1.0e-14_dp), &
        "physical mapped samples are retained")
    call record_condition(all(abs(data%coefficient_values - &
        coefficient_values) < 1.0e-14_dp), &
        "named coefficient samples are retained")
    call record_condition(data%boundary_offsets(1) == 1 .and. &
        data%boundary_offsets(boundary_count + 1) == boundary_point_count + 1, &
        "boundary offsets retain the segmented geometry")
    call record_condition(data%normalization%length_unit == "m" .and. &
        abs(data%normalization%pressure_scale - 4.0_dp) < 1.0e-14_dp, &
        "units and normalization are retained")
    data%producer = "manufactured-torus"
    data%provenance = "analytic-coordinate-oracle"

    copy = data
    call record_condition(validate_equilibrium_interchange(copy, status), &
        "interchange assignment preserves a valid deep copy")
    call record_condition(all(copy%mapped_coordinates == data%mapped_coordinates), &
        "interchange assignment copies mapped coordinates")
    call record_condition(copy%producer == data%producer .and. &
        copy%provenance == data%provenance, &
        "interchange assignment copies provenance")

    invalid = data
    invalid%boundary_offsets(2) = invalid%boundary_offsets(3) + 1
    call record_condition(.not. validate_equilibrium_interchange(invalid, status), &
        "nonmonotone boundary offsets are rejected")
    invalid = data
    invalid%coefficient_values(1, 1) = ieee_value(1.0_dp, ieee_quiet_nan)
    call record_condition(.not. validate_equilibrium_interchange(invalid, status), &
        "nonfinite coefficient samples are rejected")
    invalid = data
    invalid%normalization%pressure_scale = 0.0_dp
    call record_condition(.not. validate_equilibrium_interchange(invalid, status), &
        "nonpositive normalization scales are rejected")

    call check_summary("equilibrium interchange")
    if (.not. all_passed) error stop 1

contains

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

    subroutine build_manufactured_data( &
            mapped, physical, coefficients, profile_x, profiles, boundary, &
            offsets, coefficient_labels, profile_labels, boundary_labels)
        real(dp), intent(out) :: mapped(:, :), physical(:, :)
        real(dp), intent(out) :: coefficients(:, :), profile_x(:), profiles(:, :)
        real(dp), intent(out) :: boundary(:, :)
        integer, intent(out) :: offsets(:)
        character(len=*), intent(out) :: coefficient_labels(:), profile_labels(:)
        character(len=*), intent(out) :: boundary_labels(:)
        integer :: point

        mapped = reshape([ &
            0.1_dp, 0.2_dp, 0.3_dp, 0.4_dp, &
            0.0_dp, 0.5_dp, 1.0_dp, 1.5_dp, &
            0.0_dp, 0.2_dp, 0.4_dp, 0.6_dp], shape(mapped))
        do point = 1, sample_count
            physical(1, point) = (2.0_dp + mapped(1, point)*cos(mapped(2, point)))* &
                cos(mapped(3, point))
            physical(2, point) = (2.0_dp + mapped(1, point)*cos(mapped(2, point)))* &
                sin(mapped(3, point))
            physical(3, point) = mapped(1, point)*sin(mapped(2, point))
        end do
        coefficients(1, :) = mapped(1, :)**2 + cos(mapped(2, :))
        coefficients(2, :) = sin(mapped(2, :)) + mapped(3, :)**2
        profile_x = [0.0_dp, 0.5_dp, 1.0_dp]
        profiles(1, :) = profile_x**2
        profiles(2, :) = 1.0_dp - profile_x
        boundary(:, 1:3) = physical(:, 1:3)
        boundary(:, 4:6) = physical(:, 2:4)
        offsets = [1, 4, 7]
        coefficient_labels = [character(len=32) :: "pressure", "magnetic_flux"]
        profile_labels = [character(len=32) :: "density", "temperature"]
        boundary_labels = [character(len=32) :: "plasma", "wall"]
    end subroutine build_manufactured_data

end program test_equilibrium_interchange
