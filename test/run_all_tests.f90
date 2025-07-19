program run_all_tests
    ! Master test runner for FortFEM
    ! Runs all unit tests and reports results
    
    implicit none
    
    integer :: n_tests, n_passed, n_failed
    integer :: status
    character(len=100) :: test_name
    logical :: test_passed
    
    print *, "======================================"
    print *, "       FortFEM Test Suite"
    print *, "======================================"
    print *, ""
    
    n_tests = 0
    n_passed = 0
    n_failed = 0
    
    ! New API tests
    call run_test("Simple API", "test_simple_api", status)
    call update_counts(status, n_tests, n_passed, n_failed)
    
    call run_test("Forms API Structure", "test_forms_api", status)
    call update_counts(status, n_tests, n_passed, n_failed)
    
    ! Core module tests
    call run_test("Kinds", "test_kinds", status)
    call update_counts(status, n_tests, n_passed, n_failed)
    
    call run_test("Mesh 1D", "test_mesh_1d", status)
    call update_counts(status, n_tests, n_passed, n_failed)
    
    call run_test("Mesh 2D", "test_mesh_2d", status)
    call update_counts(status, n_tests, n_passed, n_failed)
    
    call run_test("Assembly 1D", "test_assembly_1d", status)
    call update_counts(status, n_tests, n_passed, n_failed)
    
    call run_test("Assembly 2D", "test_assembly_2d", status)
    call update_counts(status, n_tests, n_passed, n_failed)
    
    call run_test("Basis 1D", "test_basis_1d", status)
    call update_counts(status, n_tests, n_passed, n_failed)
    
    call run_test("Basis 2D", "test_basis_2d", status)
    call update_counts(status, n_tests, n_passed, n_failed)
    
    call run_test("P1 Basis 2D", "test_basis_p1_2d", status)
    call update_counts(status, n_tests, n_passed, n_failed)
    
    call run_test("P2 Basis 2D", "test_basis_p2_2d", status)
    call update_counts(status, n_tests, n_passed, n_failed)
    
    call run_test("Sparse Matrix", "test_sparse_matrix", status)
    call update_counts(status, n_tests, n_passed, n_failed)
    
    call run_test("Solver Interface", "test_solver_interface", status)
    call update_counts(status, n_tests, n_passed, n_failed)
    
    ! Application tests
    call run_test("Poisson 1D", "test_poisson_1d", status)
    call update_counts(status, n_tests, n_passed, n_failed)
    
    call run_test("Poisson 1D Benchmark", "test_poisson_1d_benchmark", status)
    call update_counts(status, n_tests, n_passed, n_failed)
    
    call run_test("Poisson 2D Benchmark", "test_poisson_2d_benchmark", status)
    call update_counts(status, n_tests, n_passed, n_failed)
    
    call run_test("Weak Forms", "test_weak_forms", status)
    call update_counts(status, n_tests, n_passed, n_failed)
    
    ! Print summary
    print *, ""
    print *, "======================================"
    print *, "         Test Summary"
    print *, "======================================"
    print '(A,I3)', " Total tests run: ", n_tests
    print '(A,I3)', " Tests passed:    ", n_passed
    print '(A,I3)', " Tests failed:    ", n_failed
    print *, ""
    
    if (n_failed == 0) then
        print *, "✅ ALL TESTS PASSED! 🎉"
        stop 0
    else
        print *, "❌ SOME TESTS FAILED!"
        stop 1
    end if
    
contains

    subroutine run_test(description, test_program, status)
        character(len=*), intent(in) :: description, test_program
        integer, intent(out) :: status
        character(len=256) :: command
        
        print '(A,A,A)', "Running ", trim(description), "..."
        
        ! Build command to run test
        write(command, '(A,A,A)') "fpm test ", trim(test_program), " 2>&1 > test_output.tmp"
        
        ! Execute test
        call execute_command_line(trim(command), exitstat=status)
        
        if (status == 0) then
            print '(A,A,A)', "  ✅ ", trim(description), " passed"
        else
            print '(A,A,A)', "  ❌ ", trim(description), " failed"
            ! Show output for failed tests
            call execute_command_line("cat test_output.tmp")
        end if
        
        ! Clean up
        call execute_command_line("rm -f test_output.tmp")
    end subroutine run_test
    
    subroutine update_counts(status, n_tests, n_passed, n_failed)
        integer, intent(in) :: status
        integer, intent(inout) :: n_tests, n_passed, n_failed
        
        n_tests = n_tests + 1
        if (status == 0) then
            n_passed = n_passed + 1
        else
            n_failed = n_failed + 1
        end if
    end subroutine update_counts

end program run_all_tests