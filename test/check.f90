module check
    implicit none
    private
    
    public :: check_condition
    public :: check_summary
    
    integer :: n_tests = 0
    integer :: n_passed = 0
    
contains
    
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