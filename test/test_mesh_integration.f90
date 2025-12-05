program test_mesh_integration
    ! Test integration between frontend mesh API and Delaunay backend
    use fortfem_api
    use fortfem_kinds
    implicit none

    integer :: test_count = 0, passed_tests = 0
    
    write(*,*) "=== FortFEM Mesh Integration Tests ==="
    write(*,*) ""
    
    call test_unit_disk_triangulation()
    call test_circle_boundary_mesh()
    call test_square_boundary_mesh()
    call test_mesh_from_domain_square()
    
    ! Summary
    write(*,*) ""
    write(*,'(A,I0,A,I0)') "Tests passed: ", passed_tests, "/", test_count
    if (passed_tests == test_count) then
        write(*,*) "✓ All integration tests passed!"
    else
        write(*,*) "✗ Some tests failed!"
        stop 1
    end if
    
contains

    subroutine test_unit_disk_triangulation()
        character(len=*), parameter :: test_name = "Unit Disk Triangulation"
        type(mesh_t) :: mesh
        
        call start_test(test_name)
        
        ! Test unit disk mesh generation using Delaunay backend
        mesh = unit_disk_mesh(resolution=0.2_dp)
        
        write(*,'(A,I0)') "  Generated vertices: ", mesh%data%n_vertices
        write(*,'(A,I0)') "  Generated triangles: ", mesh%data%n_triangles
        
        ! Basic validation
        call assert_true(mesh%data%n_vertices >= 3, "At least 3 vertices")
        call assert_true(mesh%data%n_triangles >= 1, "At least 1 triangle")
        call assert_true(allocated(mesh%data%vertices), "Vertices allocated")
        call assert_true(allocated(mesh%data%triangles), "Triangles allocated")
        
        ! Visual inspection output
        call plot(mesh, filename="build/test_unit_disk_mesh.png")
        
        call end_test()
    end subroutine
    
    subroutine test_circle_boundary_mesh()
        character(len=*), parameter :: test_name = "Circle Boundary Mesh"
        type(boundary_t) :: boundary
        type(mesh_t) :: mesh
        
        call start_test(test_name)
        
        ! Create circle boundary
        boundary = circle_boundary([0.0_dp, 0.0_dp], 1.0_dp, 12)
        
        ! Generate mesh from boundary
        mesh = mesh_from_boundary(boundary, resolution=0.3_dp)
        
        write(*,'(A,I0)') "  Boundary points: ", boundary%n_points
        write(*,'(A,I0)') "  Generated vertices: ", mesh%data%n_vertices
        write(*,'(A,I0)') "  Generated triangles: ", mesh%data%n_triangles
        
        ! Validation
        call assert_true(mesh%data%n_vertices >= boundary%n_points, "At least boundary points")
        call assert_true(mesh%data%n_triangles >= 1, "At least 1 triangle")
        call assert_true(boundary%is_closed, "Boundary is closed")
        
        ! Visual inspection output
        call plot(mesh, filename="build/test_circle_boundary_mesh.png")
        
        call end_test()
    end subroutine
    
    subroutine test_square_boundary_mesh()
        character(len=*), parameter :: test_name = "Square Boundary Mesh"
        type(boundary_t) :: boundary
        type(mesh_t) :: mesh
        integer :: n
        integer :: n_label_1, n_label_2, n_label_3, n_label_4
        
        call start_test(test_name)
        
        ! Create rectangle boundary
        n = 5
        boundary = rectangle_boundary([0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp], n)
        
        ! Generate mesh from boundary  
        mesh = mesh_from_boundary(boundary, resolution=0.2_dp)
        
        write(*,'(A,I0)') "  Boundary points: ", boundary%n_points
        write(*,'(A,I0)') "  Generated vertices: ", mesh%data%n_vertices
        write(*,'(A,I0)') "  Generated triangles: ", mesh%data%n_triangles
        
        ! Validation
        call assert_true(boundary%n_points == 4 * n, &
            "Square boundary: correct number of points")
        call assert_true(size(boundary%labels) == 4 * n - 1, &
            "Square boundary: correct number of edge labels")
        call assert_true(mesh%data%n_vertices >= 4, &
            "At least 4 vertices for square")
        call assert_true(mesh%data%n_triangles >= 2, &
            "At least 2 triangles for square")
        call assert_true(boundary%is_closed, "Boundary is closed")

        n_label_1 = count(boundary%labels == 1)
        n_label_2 = count(boundary%labels == 2)
        n_label_3 = count(boundary%labels == 3)
        n_label_4 = count(boundary%labels == 4)

        call assert_true(n_label_1 == n - 1, &
            "Square boundary: correct count for label 1")
        call assert_true(n_label_2 == n - 1, &
            "Square boundary: correct count for label 2")
        call assert_true(n_label_3 == n - 1, &
            "Square boundary: correct count for label 3")
        call assert_true(n_label_4 == n + 2, &
            "Square boundary: correct count for label 4")
        
        ! Visual inspection output
        call plot(mesh, filename="build/test_square_boundary_mesh.png")
        
        call end_test()
    end subroutine

    subroutine test_mesh_from_domain_square()
        character(len=*), parameter :: test_name = "Mesh From Domain: Unit Square"
        real(dp), parameter :: vertices(2,4) = reshape([&
            0.0_dp, 0.0_dp, &
            1.0_dp, 0.0_dp, &
            1.0_dp, 1.0_dp, &
            0.0_dp, 1.0_dp], [2, 4])
        integer, parameter :: segments(2,4) = reshape([&
            1, 2, &
            2, 3, &
            3, 4, &
            4, 1], [2, 4])
        type(mesh_t) :: mesh

        call start_test(test_name)

        mesh = mesh_from_domain(vertices, segments)

        write(*,'(A,I0)') "  Domain vertices: ", size(vertices, 2)
        write(*,'(A,I0)') "  Mesh vertices:   ", mesh%data%n_vertices
        write(*,'(A,I0)') "  Mesh triangles:  ", mesh%data%n_triangles

        call assert_true(mesh%data%n_vertices >= 4, &
            "Mesh from domain: at least 4 vertices")
        call assert_true(mesh%data%n_triangles >= 2, &
            "Mesh from domain: at least 2 triangles")
        call assert_true(allocated(mesh%data%vertices), &
            "Mesh from domain: vertices allocated")
        call assert_true(allocated(mesh%data%triangles), &
            "Mesh from domain: triangles allocated")

        call end_test()
    end subroutine test_mesh_from_domain_square
    
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
    
    subroutine assert_true(condition, description)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: description
        if (.not. condition) then
            write(*,'(A,A)') "  ✗ FAILED: ", description
            stop 1
        end if
    end subroutine

end program test_mesh_integration
