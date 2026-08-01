module check
    use, intrinsic :: iso_fortran_env, only: int32, int64
    implicit none
    private

    public :: check_condition
    public :: check_summary
    public :: check_property
    public :: property_rng_t
    public :: property_rng_initialize
    public :: property_random_unit
    public :: property_random_integer
    public :: property_shrink_integer

    type :: property_rng_t
        private
        integer(int32) :: state = 1_int32
    end type property_rng_t

    abstract interface
        logical function property_case_interface(case_seed)
            import :: int32
            integer(int32), intent(in) :: case_seed
        end function property_case_interface
    end interface

    integer :: n_tests = 0
    integer :: n_passed = 0

contains

    subroutine property_rng_initialize(rng, seed)
        type(property_rng_t), intent(out) :: rng
        integer(int32), intent(in) :: seed

        rng%state = seed
        if (rng%state == 0_int32) rng%state = int(z'6D2B79F5', int32)
    end subroutine property_rng_initialize

    real(kind(0.0d0)) function property_random_unit(rng)
        type(property_rng_t), intent(inout) :: rng
        integer(int64) :: nonnegative_state

        rng%state = ieor(rng%state, shiftl(rng%state, 13))
        rng%state = ieor(rng%state, shiftr(rng%state, 17))
        rng%state = ieor(rng%state, shiftl(rng%state, 5))
        nonnegative_state = iand(int(rng%state, int64), &
            int(z'7FFFFFFF', int64))
        property_random_unit = real(nonnegative_state, kind(0.0d0))/ &
            real(int(z'7FFFFFFF', int64), kind(0.0d0))
    end function property_random_unit

    integer function property_random_integer(rng, lower, upper)
        type(property_rng_t), intent(inout) :: rng
        integer, intent(in) :: lower, upper
        integer(int64) :: span, sample

        if (upper <= lower) then
            property_random_integer = lower
            return
        end if
        span = int(upper, int64) - int(lower, int64) + 1_int64
        sample = int(property_random_unit(rng)*real(span, kind(0.0d0)), int64)
        property_random_integer = lower + int(modulo(sample, span))
    end function property_random_integer

    integer(int32) function property_shrink_integer(value)
        integer(int32), intent(in) :: value

        property_shrink_integer = value/2_int32
    end function property_shrink_integer

    subroutine check_property(test_name, seed, case_count, callback, all_passed, &
            first_failed_seed, shrunk_seed, report_result)
        character(len=*), intent(in) :: test_name
        integer(int32), intent(in) :: seed
        integer, intent(in) :: case_count
        procedure(property_case_interface) :: callback
        logical, intent(out) :: all_passed
        integer(int32), intent(out), optional :: first_failed_seed, shrunk_seed
        logical, intent(in), optional :: report_result
        type(property_rng_t) :: rng
        integer(int32) :: case_seed, candidate, trial
        integer :: case_index
        logical :: passed

        all_passed = case_count > 0
        if (present(first_failed_seed)) first_failed_seed = 0_int32
        if (present(shrunk_seed)) shrunk_seed = 0_int32
        call property_rng_initialize(rng, seed)
        do case_index = 1, max(case_count, 0)
            case_seed = int(property_random_integer( &
                rng, -1000000000, 1000000000), int32)
            passed = callback(case_seed)
            if (.not. passed) then
                all_passed = .false.
                candidate = case_seed
                do while (abs(int(candidate, int64)) > 1_int64)
                    trial = property_shrink_integer(candidate)
                    if (callback(trial)) exit
                    candidate = trial
                end do
                if (present(first_failed_seed)) first_failed_seed = case_seed
                if (present(shrunk_seed)) shrunk_seed = candidate
                write (*, '(a,i0,a,i0)') &
                    "[PROPERTY FAIL] seed=", case_seed, " shrunk=", candidate
                exit
            end if
        end do
        if (present(report_result)) then
            if (report_result) call check_condition(all_passed, test_name)
        else
            call check_condition(all_passed, test_name)
        end if
    end subroutine check_property

    subroutine check_condition(condition, description)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: description

        n_tests = n_tests + 1
        if (condition) then
            n_passed = n_passed + 1
            print *, "[PASS] ", description
        else
            print *, "[FAIL] ", description
        end if
    end subroutine check_condition

    subroutine check_summary(test_name)
        character(len=*), intent(in) :: test_name

        print *, ""
        print *, "=== ", test_name, " Summary ==="
        print *, "Tests passed: ", n_passed, " / ", n_tests

        if (n_passed < n_tests) then
            stop 1
        end if
    end subroutine check_summary

end module check
