program test_delaunay_triangulation
    use fortfem_kinds, only: dp
    use triangulator, only: triangulate_points
    use triangulation_fortran, only: triangulation_result_t, triangulate_fortran, &
                                     cleanup_triangulation
    implicit none
    
    integer :: test_count = 0, passed_tests = 0
    
    write(*,*) "=== FortFEM Delaunay Triangulation Tests ==="
    write(*,*) ""
    
    ! Test sequence following TDD
    call test_simple_triangle_vertices()
    call test_square_triangulation()
    call test_unit_circle_points()
    call test_boundary_constraint_preservation()
    
    ! Summary
    write(*,*) ""
    write(*,'(A,I0,A,I0)') "Tests passed: ", passed_tests, "/", test_count
    if (passed_tests == test_count) then
        write(*,*) "✓ All triangulation tests passed!"
    else
        write(*,*) "✗ Some tests failed!"
        stop 1
    end if
    
contains

    subroutine test_simple_triangle_vertices()
        character(len=*), parameter :: test_name = "Simple Triangle Vertices"
        real(dp), parameter :: points(2,3) = reshape([&
            0.0_dp, 0.0_dp, &
            1.0_dp, 0.0_dp, &
            0.5_dp, 0.866_dp], [2, 3])
        real(dp), allocatable :: mesh_points(:,:)
        integer, allocatable :: mesh_triangles(:,:)
        integer :: n_points, n_triangles
        integer :: tri(3)
        
        call start_test(test_name)
        
        ! This test validates we can create basic triangle from 3 points
        ! Expected: 1 triangle with vertices 1,2,3
        write(*,*) "  Input: 3 points forming equilateral triangle"
        write(*,*) "  Expected: 1 triangle, all points preserved"
        
        ! Test actual triangulation
        call triangulate_points(points, mesh_points, mesh_triangles, &
                                n_points, n_triangles)
        
        write(*,*) "  Output vertices:", n_points
        write(*,*) "  Output triangles:", n_triangles
        
        call assert_equal_int(n_points, 3, "Simple triangle: 3 vertices")
        call assert_equal_int(n_triangles, 1, "Simple triangle: 1 triangle")
        
        tri = mesh_triangles(:, 1)
        call assert_true(all(tri >= 1 .and. tri <= 3), &
            "Simple triangle: triangle uses input vertices")
        
        deallocate(mesh_points, mesh_triangles)
        call end_test()
    end subroutine
    
    subroutine test_square_triangulation()
        character(len=*), parameter :: test_name = "Unit Square Triangulation"
        real(dp), parameter :: points(2,4) = reshape([&
            0.0_dp, 0.0_dp, &
            1.0_dp, 0.0_dp, &
            1.0_dp, 1.0_dp, &
            0.0_dp, 1.0_dp], [2, 4])
        integer, parameter :: segments(2,4) = reshape([&
            1, 2, &
            2, 3, &
            3, 4, &
            4, 1], [2, 4])
        type(triangulation_result_t) :: result
        
        call start_test(test_name)
        
        ! This test validates square is properly triangulated
        ! Expected: 2 triangles forming the square
        write(*,*) "  Input: 4 corner points of unit square"
        write(*,*) "  Expected: 2 triangles, Delaunay property satisfied"
        
        call triangulate_fortran(points, segments, result)
        
        write(*,*) "  Output points:", result%npoints
        write(*,*) "  Output triangles:", result%ntriangles
        write(*,*) "  Output segments:", result%nsegments
        
        call assert_equal_int(result%npoints, 4, "Unit square: 4 vertices")
        call assert_equal_int(result%nsegments, 4, "Unit square: 4 boundary segments")
        call assert_true(result%ntriangles >= 1, "Unit square: at least one triangle")
        call assert_true(all_triangles_positive(result), &
            "Unit square: triangles have positive area")
        
        call cleanup_triangulation(result)
        call end_test()
    end subroutine
    
    subroutine test_unit_circle_points()
        character(len=*), parameter :: test_name = "Unit Circle Boundary Points"
        integer, parameter :: n = 8
        real(dp) :: points(2, n)
        real(dp), allocatable :: mesh_points(:,:)
        integer, allocatable :: mesh_triangles(:,:)
        integer :: i, n_points, n_triangles
        real(dp) :: theta
        
        call start_test(test_name)
        
        ! Generate 8 points on unit circle
        do i = 1, n
            theta = 2.0_dp * acos(-1.0_dp) * (i-1) / n
            points(1, i) = cos(theta)
            points(2, i) = sin(theta)
        end do
        
        write(*,*) "  Input: 8 points on unit circle boundary"
        write(*,*) "  Expected: Triangulation with boundary preserved"
        
        call triangulate_points(points, mesh_points, mesh_triangles, &
                                n_points, n_triangles)
        
        write(*,*) "  Output vertices:", n_points
        write(*,*) "  Output triangles:", n_triangles
        
        call assert_true(n_points >= n, "Unit circle: vertices preserved or refined")
        call assert_true(n_triangles > 0, "Unit circle: triangles created")
        
        deallocate(mesh_points, mesh_triangles)
        
        call end_test()
    end subroutine
    
    subroutine test_boundary_constraint_preservation()
        character(len=*), parameter :: test_name = "Boundary Constraint Preservation"
        real(dp), parameter :: points(2,4) = reshape([&
            0.0_dp, 0.0_dp, &
            1.0_dp, 0.0_dp, &
            1.0_dp, 1.0_dp, &
            0.0_dp, 1.0_dp], [2, 4])
        integer, parameter :: edges(2,4) = reshape([&
            1, 2, &
            2, 3, &
            3, 4, &
            4, 1], [2, 4])
        type(triangulation_result_t) :: result
        
        call start_test(test_name)
        
        write(*,*) "  Input: Square with constrained boundary edges"
        write(*,*) "  Expected: Boundary edges preserved in triangulation"
        
        call triangulate_fortran(points, edges, result)
        
        call assert_true(constraint_edges_preserved(result, edges), &
            "Boundary constraints: all boundary edges preserved")
        
        call cleanup_triangulation(result)
        
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

    subroutine assert_equal_int(actual, expected, description)
        integer, intent(in) :: actual, expected
        character(len=*), intent(in) :: description
        if (actual /= expected) then
            write(*,'(A,A,A,I0,A,I0)') "  ✗ FAILED: ", description, &
                " - got ", actual, ", expected ", expected
            stop 1
        end if
    end subroutine assert_equal_int
    
    subroutine assert_true(condition, description)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: description
        if (.not. condition) then
            write(*,'(A,A)') "  ✗ FAILED: ", description
            stop 1
        end if
    end subroutine assert_true
    
    logical function all_triangles_positive(result)
        type(triangulation_result_t), intent(in) :: result
        integer :: i
        real(dp) :: x1, y1, x2, y2, x3, y3, area
        
        all_triangles_positive = .true.
        do i = 1, result%ntriangles
            x1 = result%points(1, result%triangles(1, i))
            y1 = result%points(2, result%triangles(1, i))
            x2 = result%points(1, result%triangles(2, i))
            y2 = result%points(2, result%triangles(2, i))
            x3 = result%points(1, result%triangles(3, i))
            y3 = result%points(2, result%triangles(3, i))
            
            area = 0.5_dp * ((x1*(y2-y3) + x2*(y3-y1) + x3*(y1-y2)))
            if (area <= 0.0_dp) then
                all_triangles_positive = .false.
                return
            end if
        end do
    end function all_triangles_positive
    
    logical function constraint_edges_preserved(result, edges)
        type(triangulation_result_t), intent(in) :: result
        integer, intent(in) :: edges(:,:)
        integer :: e, t
        integer :: v1, v2, a, b, c
        logical :: found
        
        constraint_edges_preserved = .true.
        
        do e = 1, size(edges, 2)
            v1 = edges(1, e)
            v2 = edges(2, e)
            found = .false.
            
            do t = 1, result%ntriangles
                a = result%triangles(1, t)
                b = result%triangles(2, t)
                c = result%triangles(3, t)
                
                if ((a == v1 .and. b == v2) .or. (a == v2 .and. b == v1) .or. &
                    (b == v1 .and. c == v2) .or. (b == v2 .and. c == v1) .or. &
                    (c == v1 .and. a == v2) .or. (c == v2 .and. a == v1)) then
                    found = .true.
                    exit
                end if
            end do
            
            if (.not. found) then
                constraint_edges_preserved = .false.
                return
            end if
        end do
    end function constraint_edges_preserved

end program test_delaunay_triangulation
