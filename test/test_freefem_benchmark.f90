program test_freefem_benchmark
    ! TDD test for FreeFEM benchmark comparison
    use fortfem_api
    use fortfem_kinds
    implicit none

    integer :: test_count = 0, passed_tests = 0
    
    write(*,*) "=== FreeFEM Benchmark TDD Tests ==="
    write(*,*) ""
    
    call test_benchmark_setup()
    call test_mesh_comparison()
    call test_timing_measurement()
    
    ! Summary
    write(*,*) ""
    write(*,'(A,I0,A,I0)') "Tests passed: ", passed_tests, "/", test_count
    if (passed_tests == test_count) then
        write(*,*) "✓ All benchmark tests passed!"
    else
        write(*,*) "✗ Some tests failed!"
        stop 1
    end if
    
contains

    subroutine test_benchmark_setup()
        character(len=*), parameter :: test_name = "Benchmark Setup"
        type(mesh_t) :: mesh
        real(dp) :: start_time, end_time
        
        call start_test(test_name)
        
        ! Test that we can measure timing for mesh generation
        call cpu_time(start_time)
        mesh = unit_square_mesh(10)  ! 10x10 mesh for benchmark
        call cpu_time(end_time)
        
        write(*,'(A,ES10.3,A)') "  FortFEM mesh time: ", end_time - start_time, " seconds"
        write(*,'(A,I0)') "  Generated vertices: ", mesh%data%n_vertices
        write(*,'(A,I0)') "  Generated triangles: ", mesh%data%n_triangles
        
        ! Validation for benchmark
        call assert_equal(mesh%data%n_vertices, 100, "10x10 = 100 vertices")
        call assert_true(mesh%data%n_triangles > 0, "Has triangles")
        call assert_true(end_time > start_time, "Time measurement works")
        
        call end_test()
    end subroutine
    
    subroutine test_mesh_comparison()
        character(len=*), parameter :: test_name = "Mesh Quality Comparison"
        type(mesh_t) :: mesh_coarse, mesh_fine
        real(dp) :: coarse_time, fine_time, start_time, end_time
        
        call start_test(test_name)
        
        ! Coarse mesh timing
        call cpu_time(start_time)
        mesh_coarse = unit_square_mesh(5)
        call cpu_time(end_time)
        coarse_time = end_time - start_time
        
        ! Fine mesh timing  
        call cpu_time(start_time)
        mesh_fine = unit_square_mesh(20)
        call cpu_time(end_time)
        fine_time = end_time - start_time
        
        write(*,'(A,I0,A,ES10.3,A)') "  Coarse (5x5): ", mesh_coarse%data%n_triangles, " triangles, ", coarse_time, "s"
        write(*,'(A,I0,A,ES10.3,A)') "  Fine (20x20): ", mesh_fine%data%n_triangles, " triangles, ", fine_time, "s"
        
        ! Scaling validation
        call assert_true(mesh_fine%data%n_triangles > mesh_coarse%data%n_triangles, "Fine mesh has more triangles")
        call assert_true(fine_time >= coarse_time, "Fine mesh takes longer (expected)")
        
        call end_test()
    end subroutine
    
    subroutine test_timing_measurement()
        character(len=*), parameter :: test_name = "Timing Measurement Accuracy"
        real(dp) :: times(5), mean_time, std_dev
        integer :: i
        
        call start_test(test_name)
        
        ! Run multiple timing measurements for statistical accuracy
        do i = 1, 5
            times(i) = time_mesh_generation()
        end do
        
        mean_time = sum(times) / 5.0_dp
        std_dev = sqrt(sum((times - mean_time)**2) / 4.0_dp)
        
        write(*,'(A,ES10.3,A,ES10.3,A)') "  Mean time: ", mean_time, " ± ", std_dev, " seconds"
        write(*,'(A,ES8.1,A)') "  Relative std dev: ", 100.0_dp * std_dev / mean_time, "%"
        
        ! Validation for benchmark reliability
        call assert_true(mean_time > 0.0_dp, "Mean time positive")
        call assert_true(std_dev / mean_time < 0.5_dp, "Timing reasonably consistent")
        
        call end_test()
    end subroutine
    
    ! Helper functions
    real(dp) function time_mesh_generation() result(elapsed)
        real(dp) :: start_time, end_time
        type(mesh_t) :: mesh
        
        call cpu_time(start_time)
        mesh = unit_square_mesh(15)  ! Standard benchmark size
        call cpu_time(end_time)
        
        elapsed = end_time - start_time
    end function
    
    ! Test framework helpers
    subroutine start_test(test_name)
        character(len=*), intent(in) :: test_name
        test_count = test_count + 1
        write(*,'(A,I0,A,A)') "Test ", test_count, ": ", test_name
    end subroutine
    
    subroutine end_test()
        passed_tests = passed_tests + 1
        write(*,*) "  ✓ PASSED"
    end subroutine
    
    subroutine assert_equal(actual, expected, description)
        integer, intent(in) :: actual, expected
        character(len=*), intent(in) :: description
        if (actual /= expected) then
            write(*,'(A,A,A,I0,A,I0)') "  ✗ FAILED: ", description, " - got ", actual, ", expected ", expected
            stop 1
        end if
    end subroutine
    
    subroutine assert_true(condition, description)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: description
        if (.not. condition) then
            write(*,'(A,A)') "  ✗ FAILED: ", description
            stop 1
        end if
    end subroutine

end program test_freefem_benchmark