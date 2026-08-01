program test_property_helpers
    use, intrinsic :: iso_fortran_env, only: int32, int64
    use check, only: check_condition, check_property, check_summary, &
        property_random_integer, property_random_unit, property_rng_initialize, &
        property_rng_t, property_shrink_integer
    implicit none

    type(property_rng_t) :: first_rng, second_rng
    real(kind(0.0d0)) :: first_unit, second_unit
    integer :: first_integer, second_integer, index
    integer(int32) :: first_failed, shrunk, failure_seed, failure_shrunk
    logical :: all_passed, expected_failure_passed

    call property_rng_initialize(first_rng, 123456_int32)
    call property_rng_initialize(second_rng, 123456_int32)
    do index = 1, 12
        first_unit = property_random_unit(first_rng)
        second_unit = property_random_unit(second_rng)
        call check_condition(first_unit == second_unit .and. &
            first_unit >= 0.0d0 .and. first_unit <= 1.0d0, &
            "property RNG is deterministic and bounded")
        first_integer = property_random_integer(first_rng, -7, 9)
        second_integer = property_random_integer(second_rng, -7, 9)
        call check_condition(first_integer == second_integer .and. &
            first_integer >= -7 .and. first_integer <= 9, &
            "property RNG integer samples are reproducible and in range")
    end do
    call check_condition(property_shrink_integer(17_int32) == 8_int32 .and. &
        property_shrink_integer(-17_int32) == -8_int32 .and. &
        property_shrink_integer(1_int32) == 0_int32, &
        "property shrinker moves signed witnesses toward zero")

    call check_property( &
        "seeded property runner completes reproducible cases", 98765_int32, 24, &
        bounded_case, all_passed, first_failed, shrunk)
    call check_condition(all_passed .and. first_failed == 0_int32 .and. &
        shrunk == 0_int32, "seeded property runner reports no failure seed")
    call check_property( &
        "intentional shrinking probe", 24680_int32, 1, failing_case, &
        expected_failure_passed, failure_seed, failure_shrunk, &
        report_result=.false.)
    call check_condition(.not. expected_failure_passed .and. &
        failure_seed /= 0_int32 .and. &
        abs(int(failure_shrunk, int64)) <= 1_int64, &
        "property runner shrinks a failing seed deterministically")
    call check_summary("deterministic property helpers")
    if (.not. all_passed) error stop 1

contains

    logical function bounded_case(case_seed)
        integer(int32), intent(in) :: case_seed
        type(property_rng_t) :: rng
        real(kind(0.0d0)) :: value
        integer :: sample, trial

        bounded_case = .true.
        call property_rng_initialize(rng, case_seed)
        do trial = 1, 8
            value = property_random_unit(rng)
            sample = property_random_integer(rng, -100, 100)
            bounded_case = bounded_case .and. value >= 0.0d0 .and. &
                value <= 1.0d0 .and. sample >= -100 .and. sample <= 100
        end do
    end function bounded_case

    logical function failing_case(case_seed)
        integer(int32), intent(in) :: case_seed

        failing_case = .false.
    end function failing_case

end program test_property_helpers
