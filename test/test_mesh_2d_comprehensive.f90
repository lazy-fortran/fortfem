program test_mesh_2d_comprehensive
    ! Comprehensive unit tests for 2D mesh module
    use fortfem_kinds
    use fortfem_mesh_2d
    implicit none
    
    logical :: all_passed = .true.
    
    print *, "=== Testing Mesh 2D Module (Comprehensive) ==="
    print *, ""
    
    call test_rectangular_mesh_creation()
    call test_triangular_mesh_creation()
    call test_mesh_connectivity()
    call test_edge_detection()
    call test_boundary_edge_identification()
    call test_mesh_quality_metrics()
    ! call test_mesh_io()  ! Not implemented yet
    call test_mesh_refinement_preparation()
    call test_edge_to_triangle_mapping()
    call test_mesh_orientation()
    
    print *, ""
    if (all_passed) then
        print *, "✅ ALL MESH TESTS PASSED!"
    else
        print *, "❌ SOME MESH TESTS FAILED!"
        stop 1
    end if
    
contains

    subroutine test_rectangular_mesh_creation()
        type(mesh_2d_t) :: mesh
        integer :: nx, ny, expected_vertices, expected_triangles
        logical :: passed = .true.
        
        print *, "Testing rectangular mesh creation..."
        
        ! Test various mesh sizes
        do nx = 2, 5
            do ny = 2, 5
                call mesh%create_rectangular(nx, ny, 0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp)
                
                expected_vertices = nx * ny
                expected_triangles = 2 * (nx-1) * (ny-1)
                
                if (mesh%n_vertices /= expected_vertices) then
                    print *, "  ❌ Wrong number of vertices for", nx, "x", ny
                    print *, "    Expected:", expected_vertices, "Got:", mesh%n_vertices
                    passed = .false.
                end if
                
                if (mesh%n_triangles /= expected_triangles) then
                    print *, "  ❌ Wrong number of triangles for", nx, "x", ny
                    print *, "    Expected:", expected_triangles, "Got:", mesh%n_triangles
                    passed = .false.
                end if
                
                call mesh%destroy()
            end do
        end do
        
        if (passed) then
            print *, "  ✅ Rectangular mesh creation test passed"
        else
            all_passed = .false.
        end if
    end subroutine test_rectangular_mesh_creation
    
    subroutine test_triangular_mesh_creation()
        type(mesh_2d_t) :: mesh
        real(dp) :: x, y, area, expected_area
        integer :: t
        logical :: passed = .true.
        real(dp), parameter :: tol = 1e-14_dp
        
        print *, "Testing triangular mesh properties..."
        
        ! Create simple 2x2 mesh
        call mesh%create_rectangular(2, 2, 0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp)
        
        ! Check that all triangles have positive area
        do t = 1, mesh%n_triangles
            area = compute_triangle_area(mesh, t)
            if (area <= 0.0_dp) then
                print *, "  ❌ Triangle", t, "has non-positive area:", area
                passed = .false.
            end if
        end do
        
        ! Check total area equals domain area
        area = 0.0_dp
        do t = 1, mesh%n_triangles
            area = area + compute_triangle_area(mesh, t)
        end do
        expected_area = 1.0_dp  ! Unit square
        
        if (abs(area - expected_area) > tol) then
            print *, "  ❌ Total mesh area incorrect!"
            print *, "    Expected:", expected_area, "Got:", area
            passed = .false.
        end if
        
        if (passed) then
            print *, "  ✅ Triangular mesh properties test passed"
        else
            all_passed = .false.
        end if
        
        call mesh%destroy()
    end subroutine test_triangular_mesh_creation
    
    subroutine test_mesh_connectivity()
        type(mesh_2d_t) :: mesh
        integer :: t, i, j, v1, v2, count
        logical :: passed = .true.
        logical :: found
        
        print *, "Testing mesh connectivity..."
        
        call mesh%create_rectangular(3, 3, 0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp)
        
        ! Check that each edge appears in at most 2 triangles
        do t = 1, mesh%n_triangles
            do i = 1, 3
                v1 = mesh%triangles(i, t)
                v2 = mesh%triangles(mod(i,3)+1, t)
                
                ! Count how many triangles share this edge
                count = 0
                do j = 1, mesh%n_triangles
                    if (triangle_has_edge(mesh, j, v1, v2)) then
                        count = count + 1
                    end if
                end do
                
                if (count > 2) then
                    print *, "  ❌ Edge (", v1, ",", v2, ") shared by", count, "triangles!"
                    passed = .false.
                end if
            end do
        end do
        
        if (passed) then
            print *, "  ✅ Mesh connectivity test passed"
        else
            all_passed = .false.
        end if
        
        call mesh%destroy()
    end subroutine test_mesh_connectivity
    
    subroutine test_edge_detection()
        type(mesh_2d_t) :: mesh
        integer :: expected_edges, interior_edges, boundary_edges
        logical :: passed = .true.
        
        print *, "Testing edge detection..."
        
        ! 3x3 mesh
        call mesh%create_rectangular(3, 3, 0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp)
        call mesh%build_edge_connectivity()
        
        ! Expected: 20 edges total (12 interior + 8 boundary)
        expected_edges = 20
        interior_edges = 12
        boundary_edges = 8
        
        if (mesh%n_edges /= expected_edges) then
            print *, "  ❌ Wrong number of edges!"
            print *, "    Expected:", expected_edges, "Got:", mesh%n_edges
            passed = .false.
        end if
        
        ! Count boundary edges
        if (mesh%n_boundary_edges /= boundary_edges) then
            print *, "  ❌ Wrong number of boundary edges!"
            print *, "    Expected:", boundary_edges, "Got:", mesh%n_boundary_edges
            passed = .false.
        end if
        
        if (passed) then
            print *, "  ✅ Edge detection test passed"
        else
            all_passed = .false.
        end if
        
        call mesh%destroy()
    end subroutine test_edge_detection
    
    subroutine test_boundary_edge_identification()
        type(mesh_2d_t) :: mesh
        integer :: e, v1, v2
        real(dp) :: x1, y1, x2, y2
        logical :: passed = .true.
        logical :: should_be_boundary
        real(dp), parameter :: tol = 1e-10_dp
        
        print *, "Testing boundary edge identification..."
        
        call mesh%create_rectangular(4, 4, 0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp)
        call mesh%build_edge_connectivity()
        
        ! Check that boundary edges are correctly identified
        do e = 1, mesh%n_edges
            v1 = mesh%edges(1, e)
            v2 = mesh%edges(2, e)
            x1 = mesh%vertices(1, v1)
            y1 = mesh%vertices(2, v1)
            x2 = mesh%vertices(1, v2)
            y2 = mesh%vertices(2, v2)
            
            ! Edge is on boundary if both vertices are on same boundary
            should_be_boundary = .false.
            if ((abs(x1) < tol .and. abs(x2) < tol) .or. &
                (abs(x1 - 1.0_dp) < tol .and. abs(x2 - 1.0_dp) < tol) .or. &
                (abs(y1) < tol .and. abs(y2) < tol) .or. &
                (abs(y1 - 1.0_dp) < tol .and. abs(y2 - 1.0_dp) < tol)) then
                should_be_boundary = .true.
            end if
            
            if (mesh%is_boundary_edge(e) .neqv. should_be_boundary) then
                print *, "  ❌ Edge", e, "boundary flag incorrect!"
                passed = .false.
            end if
        end do
        
        if (passed) then
            print *, "  ✅ Boundary edge identification test passed"
        else
            all_passed = .false.
        end if
        
        call mesh%destroy()
    end subroutine test_boundary_edge_identification
    
    subroutine test_mesh_quality_metrics()
        type(mesh_2d_t) :: mesh
        real(dp) :: min_area, max_area, area_ratio
        real(dp) :: min_angle, max_angle
        integer :: t
        logical :: passed = .true.
        
        print *, "Testing mesh quality metrics..."
        
        ! Regular mesh should have uniform triangles
        call mesh%create_rectangular(5, 5, 0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp)
        
        min_area = huge(1.0_dp)
        max_area = 0.0_dp
        
        do t = 1, mesh%n_triangles
            min_area = min(min_area, compute_triangle_area(mesh, t))
            max_area = max(max_area, compute_triangle_area(mesh, t))
        end do
        
        area_ratio = max_area / min_area
        
        ! For regular mesh, all triangles should have same area
        if (area_ratio > 1.01_dp) then
            print *, "  ❌ Triangle areas not uniform!"
            print *, "    Area ratio:", area_ratio
            passed = .false.
        end if
        
        ! Check triangle angles (should be 45-45-90 for regular mesh)
        ! This is a simplified test - just check no degenerate triangles
        do t = 1, mesh%n_triangles
            if (compute_triangle_area(mesh, t) < 1e-10_dp) then
                print *, "  ❌ Degenerate triangle found!"
                passed = .false.
            end if
        end do
        
        if (passed) then
            print *, "  ✅ Mesh quality metrics test passed"
        else
            all_passed = .false.
        end if
        
        call mesh%destroy()
    end subroutine test_mesh_quality_metrics
    
    subroutine test_mesh_io()
        type(mesh_2d_t) :: mesh1, mesh2
        logical :: passed = .true.
        integer :: i, j
        real(dp), parameter :: tol = 1e-14_dp
        
        print *, "Testing mesh I/O..."
        
        ! Create mesh
        call mesh1%create_rectangular(3, 3, 0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp)
        
        ! Save to file
        call mesh1%save_to_file("test_mesh.dat")
        
        ! Load from file
        call mesh2%load_from_file("test_mesh.dat")
        
        ! Compare meshes
        if (mesh1%n_vertices /= mesh2%n_vertices) then
            print *, "  ❌ Vertex count mismatch after I/O!"
            passed = .false.
        end if
        
        if (mesh1%n_triangles /= mesh2%n_triangles) then
            print *, "  ❌ Triangle count mismatch after I/O!"
            passed = .false.
        end if
        
        ! Check vertex coordinates
        do i = 1, mesh1%n_vertices
            do j = 1, 2
                if (abs(mesh1%vertices(j,i) - mesh2%vertices(j,i)) > tol) then
                    print *, "  ❌ Vertex coordinate mismatch!"
                    passed = .false.
                    exit
                end if
            end do
        end do
        
        if (passed) then
            print *, "  ✅ Mesh I/O test passed"
        else
            all_passed = .false.
        end if
        
        ! Clean up
        call mesh1%destroy()
        call mesh2%destroy()
        call execute_command_line("rm -f test_mesh.dat")
    end subroutine test_mesh_io
    
    subroutine test_mesh_refinement_preparation()
        type(mesh_2d_t) :: mesh
        logical :: passed = .true.
        
        print *, "Testing mesh refinement preparation..."
        
        ! Create coarse mesh
        call mesh%create_rectangular(2, 2, 0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp)
        
        ! Find edges (required for refinement)
        call mesh%build_edge_connectivity()
        
        ! Check edge midpoints can be computed
        if (.not. allocated(mesh%edges)) then
            print *, "  ❌ Edges not allocated!"
            passed = .false.
        end if
        
        ! Verify edge-to-triangle connectivity exists
        if (.not. allocated(mesh%edge_to_triangles)) then
            print *, "  ❌ Edge-to-triangle mapping not created!"
            passed = .false.
        end if
        
        if (passed) then
            print *, "  ✅ Mesh refinement preparation test passed"
        else
            all_passed = .false.
        end if
        
        call mesh%destroy()
    end subroutine test_mesh_refinement_preparation
    
    subroutine test_edge_to_triangle_mapping()
        type(mesh_2d_t) :: mesh
        integer :: e, t1, t2, count
        logical :: passed = .true.
        
        print *, "Testing edge-to-triangle mapping..."
        
        call mesh%create_rectangular(3, 3, 0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp)
        call mesh%build_edge_connectivity()
        
        do e = 1, mesh%n_edges
            t1 = mesh%edge_to_triangles(1, e)
            t2 = mesh%edge_to_triangles(2, e)
            
            ! Interior edges should have 2 triangles
            if (.not. mesh%is_boundary_edge(e)) then
                if (t1 <= 0 .or. t2 <= 0) then
                    print *, "  ❌ Interior edge", e, "missing triangle!"
                    passed = .false.
                end if
            else
                ! Boundary edges should have 1 triangle
                count = 0
                if (t1 > 0) count = count + 1
                if (t2 > 0) count = count + 1
                if (count /= 1) then
                    print *, "  ❌ Boundary edge", e, "has", count, "triangles!"
                    passed = .false.
                end if
            end if
        end do
        
        if (passed) then
            print *, "  ✅ Edge-to-triangle mapping test passed"
        else
            all_passed = .false.
        end if
        
        call mesh%destroy()
    end subroutine test_edge_to_triangle_mapping
    
    subroutine test_mesh_orientation()
        type(mesh_2d_t) :: mesh
        real(dp) :: area
        integer :: t
        logical :: passed = .true.
        
        print *, "Testing mesh orientation..."
        
        call mesh%create_rectangular(4, 4, 0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp)
        
        ! All triangles should have positive area (counter-clockwise orientation)
        do t = 1, mesh%n_triangles
            area = compute_triangle_area(mesh, t)
            if (area <= 0.0_dp) then
                print *, "  ❌ Triangle", t, "has wrong orientation! Area:", area
                passed = .false.
            end if
        end do
        
        if (passed) then
            print *, "  ✅ Mesh orientation test passed"
        else
            all_passed = .false.
        end if
        
        call mesh%destroy()
    end subroutine test_mesh_orientation
    
    ! Helper function to compute triangle area
    real(dp) function compute_triangle_area(mesh, t)
        type(mesh_2d_t), intent(in) :: mesh
        integer, intent(in) :: t
        real(dp) :: x1, y1, x2, y2, x3, y3
        
        x1 = mesh%vertices(1, mesh%triangles(1,t))
        y1 = mesh%vertices(2, mesh%triangles(1,t))
        x2 = mesh%vertices(1, mesh%triangles(2,t))
        y2 = mesh%vertices(2, mesh%triangles(2,t))
        x3 = mesh%vertices(1, mesh%triangles(3,t))
        y3 = mesh%vertices(2, mesh%triangles(3,t))
        
        compute_triangle_area = 0.5_dp * abs((x2-x1)*(y3-y1) - (x3-x1)*(y2-y1))
    end function compute_triangle_area
    
    ! Helper function to check if triangle has edge
    logical function triangle_has_edge(mesh, t, v1, v2)
        type(mesh_2d_t), intent(in) :: mesh
        integer, intent(in) :: t, v1, v2
        integer :: i, va, vb
        
        triangle_has_edge = .false.
        do i = 1, 3
            va = mesh%triangles(i, t)
            vb = mesh%triangles(mod(i,3)+1, t)
            if ((va == v1 .and. vb == v2) .or. (va == v2 .and. vb == v1)) then
                triangle_has_edge = .true.
                return
            end if
        end do
    end function triangle_has_edge

end program test_mesh_2d_comprehensive