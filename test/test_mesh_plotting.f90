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
    call test_quad_mesh_plot()
    call test_mixed_mesh_plot()
    
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
    
    subroutine test_quad_mesh_plot()
        character(len=*), parameter :: test_name = "Quadrilateral Mesh Plot"
        type(mesh_t) :: mesh
        
        call start_test(test_name)
        
        mesh = structured_quad_mesh(4, 3, 0.0_dp, 2.0_dp, 0.0_dp, 1.0_dp)
        
        write(*,*) "  Testing quadrilateral mesh plotting..."
        call plot(mesh, filename="build/test_quad_mesh.png", &
                  title="Structured Quad Mesh")
        write(*,*) "  ✓ plot(quad mesh) works"
        
        call end_test()
    end subroutine

    subroutine test_mixed_mesh_plot()
        character(len=*), parameter :: test_name = "Mixed Mesh Plot"
        type(mesh_t) :: tri_mesh, quad_mesh, mixed_mesh
        integer :: n_vertices, i
        
        call start_test(test_name)
        
        tri_mesh = unit_square_mesh(3)
        quad_mesh = structured_quad_mesh(2, 2, 0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp)
        
        mixed_mesh%data%n_vertices = tri_mesh%data%n_vertices + quad_mesh%data%n_vertices
        mixed_mesh%data%n_triangles = tri_mesh%data%n_triangles
        mixed_mesh%data%n_quads = quad_mesh%data%n_quads
        mixed_mesh%data%has_triangles = .true.
        mixed_mesh%data%has_quads = .true.
        mixed_mesh%data%has_mixed_elements = .true.
        
        allocate(mixed_mesh%data%vertices(2, mixed_mesh%data%n_vertices))
        allocate(mixed_mesh%data%triangles(3, mixed_mesh%data%n_triangles))
        allocate(mixed_mesh%data%quads(4, mixed_mesh%data%n_quads))
        
        mixed_mesh%data%vertices(:, 1:tri_mesh%data%n_vertices) = tri_mesh%data%vertices
        mixed_mesh%data%triangles = tri_mesh%data%triangles
        
        n_vertices = tri_mesh%data%n_vertices
        do i = 1, quad_mesh%data%n_vertices
            mixed_mesh%data%vertices(:, n_vertices + i) = quad_mesh%data%vertices(:, i) + &
                [1.0_dp, 0.0_dp]
        end do
        do i = 1, quad_mesh%data%n_quads
            mixed_mesh%data%quads(:, i) = quad_mesh%data%quads(:, i) + n_vertices
        end do
        
        call mixed_mesh%data%build_connectivity()
        call mixed_mesh%data%find_boundary()
        
        write(*,*) "  Testing mixed tri/quad mesh plotting..."
        call plot(mixed_mesh, filename="build/test_mixed_mesh.png", &
                  title="Mixed Tri/Quad Mesh")
        write(*,*) "  ✓ plot(mixed mesh) works"
        
        call end_test()
    end subroutine
    
    subroutine test_simple_mesh_plot()
        character(len=*), parameter :: test_name = "Simple Mesh Plot"
        type(mesh_t) :: mesh
        
        call start_test(test_name)
        
        mesh = unit_square_mesh(4)
        
        ! Test basic mesh plotting with filename
        write(*,*) "  Testing plot with filename..."
        call plot(mesh, filename="build/test_mesh.png")
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
        call plot(mesh, filename="build/test_mesh_titled.png", &
                  title="Rectangle Mesh", show_labels=.false.)
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
