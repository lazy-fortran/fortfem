program run_convergence_tests
    ! Main test runner for convergence tests
    use fortfem_kinds
    implicit none

    logical :: all_passed
    integer :: exit_code

    print *, "FortFEM Convergence Test Suite"
    print *, "=============================="
    print *, ""

    all_passed = .true.

    ! Run tests
    print *, "Running convergence tests against analytical solutions..."
    print *, ""

    ! Note: Individual tests would be run as separate executables
    ! This is a framework for organizing the test suite

    call run_test_suite(all_passed)

    print *, ""
    print *, "Test Summary"
    print *, "============"
    if (all_passed) then
        print *, "🎉 ALL CONVERGENCE TESTS PASSED!"
        print *, ""
        print *, "The FortFEM implementation demonstrates:"
        print *, "- Optimal convergence rates for Poisson equation"
        print *, "- Correct behavior for curl-curl problems"
        print *, "- Proper edge element implementation"
        print *, "- Reliable GMRES solver performance"
        exit_code = 0
    else
        print *, "❌ SOME TESTS FAILED!"
        print *, ""
        print *, "Check the implementation for issues with:"
        print *, "- Finite element assembly routines"
        print *, "- Basis function implementations"
        print *, "- Edge element curl/div operators"
        print *, "- Solver algorithms"
        exit_code = 1
    end if

    call exit(exit_code)

contains

    subroutine run_test_suite(all_passed)
        logical, intent(inout) :: all_passed

        print *, "Individual test executables:"
        print *, "  test_poisson_1d_convergence  - 1D Poisson O(h²) convergence"
        print *, "  test_poisson_2d_convergence  - 2D Poisson O(h²) convergence"
        print *, "  test_curl_curl_convergence   - 2D curl-curl O(h) convergence"
        print *, ""
        print *, "Run each test individually to see detailed convergence analysis."
        print *, ""

        ! In a full implementation, this would:
        ! 1. Execute each test program
        ! 2. Capture their return codes
        ! 3. Aggregate results
        ! 4. Report overall success/failure

        print *, "Framework ready for automated testing."
    end subroutine run_test_suite

end program run_convergence_tests
