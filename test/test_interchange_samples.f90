program test_interchange_samples
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use check, only: check_condition, check_summary
    use fortfem_interop, only: &
        compare_interchange_samples, initialize_interchange_samples, &
        interchange_sample_set_t, validate_interchange_samples
    implicit none

    real(dp) :: coordinates(2, 3), values(2, 3), weights(3)
    real(dp) :: candidate_values(2, 3), candidate_coordinates(2, 3)
    real(dp) :: absolute_error, relative_error, maximum_error
    type(interchange_sample_set_t) :: reference, candidate, copy
    integer :: status
    logical :: all_passed

    all_passed = .true.
    coordinates = reshape([ &
        0.0_dp, 0.0_dp, &
        1.0_dp, 0.0_dp, &
        0.0_dp, 1.0_dp], shape(coordinates))
    values = reshape([ &
        1.0_dp, 2.0_dp, &
        2.0_dp, 4.0_dp, &
        3.0_dp, 6.0_dp], shape(values))
    weights = [1.0_dp, 2.0_dp, 3.0_dp]

    call initialize_interchange_samples( &
        reference, coordinates, values, weights, "reference", "analytic", &
        status)
    call record_condition(status == 0, "valid physical samples initialize")
    call record_condition(validate_interchange_samples(reference, status), &
        "initialized physical samples validate")
    call record_condition(status == 0, "valid physical samples report status")

    copy = reference
    copy%values(1, 1) = -7.0_dp
    call record_condition(reference%values(1, 1) == 1.0_dp, &
        "sample assignment is a deep copy")

    candidate_values = values
    candidate_values(1, 2) = candidate_values(1, 2) + 0.5_dp
    call initialize_interchange_samples( &
        candidate, coordinates, candidate_values, weights, "candidate", &
        "independent", status)
    call record_condition(status == 0, "candidate samples initialize")
    call compare_interchange_samples(reference, candidate, 1.0e-14_dp, &
        absolute_error, relative_error, maximum_error, status)
    call record_condition(status == 0, "matching sample sets compare")
    call record_condition( &
        abs(absolute_error - sqrt(2.0_dp*0.5_dp**2)) < 1.0e-14_dp, &
        "weighted sample L2 error is independent")
    call record_condition( &
        abs(relative_error - absolute_error/sqrt(1.0_dp + 2.0_dp*4.0_dp + &
        3.0_dp*9.0_dp + 1.0_dp*4.0_dp + 2.0_dp*16.0_dp + &
        3.0_dp*36.0_dp)) < 1.0e-14_dp, &
        "relative sample error uses the reference norm")
    call record_condition(maximum_error == 0.5_dp, &
        "sample Linf error is independent")

    candidate_coordinates = coordinates
    candidate_coordinates(1, 2) = candidate_coordinates(1, 2) + 1.0e-3_dp
    call initialize_interchange_samples( &
        candidate, candidate_coordinates, candidate_values, weights, &
        "candidate", "mismatched", status)
    call record_condition(status == 0, "mismatched sample set initializes")
    call compare_interchange_samples(reference, candidate, 1.0e-6_dp, &
        absolute_error, relative_error, maximum_error, status)
    call record_condition(status /= 0, "different sample coordinates reject")
    call record_condition(absolute_error == huge(1.0_dp), &
        "rejected coordinate comparison clears the error")

    weights(2) = -1.0_dp
    call initialize_interchange_samples( &
        candidate, coordinates, candidate_values, weights, "candidate", &
        "invalid-weight", status)
    call record_condition(status /= 0, "non-positive sample weights reject")

    call check_summary("interchange physical sample contract")
    if (.not. all_passed) error stop 1

contains

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_interchange_samples
