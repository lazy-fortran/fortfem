program test_mesh_plotting
    ! TDD test for mesh plotting functionality
    use fortfem_api
    use fortfem_kinds
    implicit none

    integer :: test_count = 0, passed_tests = 0
    
    write(*,*) "=== Mesh Plotting TDD Tests ==="
    write(*,*) ""
    
    call test_plot_function_exists()
    call test_simple_mesh_plot()
    call test_plot_with_options()
    
    ! Summary
    write(*,*) ""
    write(*,'(A,I0,A,I0)') "Tests passed: ", passed_tests, "/", test_count
    if (passed_tests == test_count) then
        write(*,*) "✓ All mesh plotting tests passed!"
    else
        write(*,*) "✗ Some tests failed!"
        stop 1
    end if
    
contains

    subroutine test_plot_function_exists()
        character(len=*), parameter :: test_name = "Plot Function Interface"
        type(mesh_t) :: mesh
        
        call start_test(test_name)
        
        mesh = unit_square_mesh(5)
        
        ! Test that plot function can be called
        ! Expected: Simple wireframe plot showing triangulation
        write(*,*) "  Testing basic plot interface..."
        call plot(mesh)  ! Should generate mesh plot
        write(*,*) "  ✓ plot(mesh) interface works"
        
        call end_test()
    end subroutine
    
    subroutine test_simple_mesh_plot()
        character(len=*), parameter :: test_name = "Simple Mesh Plot"
        type(mesh_t) :: mesh
        
        call start_test(test_name)
        
        mesh = unit_square_mesh(4)
        
        ! Test basic mesh plotting with filename
        write(*,*) "  Testing plot with filename..."
        call plot(mesh, filename="test_mesh.png")
        write(*,*) "  ✓ plot(mesh, filename) works"
        
        call end_test()
    end subroutine
    
    subroutine test_plot_with_options()
        character(len=*), parameter :: test_name = "Plot with Options"
        type(mesh_t) :: mesh
        
        call start_test(test_name)
        
        mesh = rectangle_mesh(3, 4, [0.0_dp, 2.0_dp, 0.0_dp, 1.0_dp])
        
        ! Test mesh plotting with various options
        write(*,*) "  Testing plot with title and options..."
        call plot(mesh, filename="test_mesh_titled.png", title="Rectangle Mesh", &
                  show_vertices=.true., show_edges=.true.)
        write(*,*) "  ✓ plot(mesh, options) works"
        
        call end_test()
    end subroutine
    
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

end program test_mesh_plotting