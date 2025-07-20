program test_edge_debug
    ! Debug edge existence in triangulated mesh
    use fortfem_kinds
    use delaunay_types  
    use constrained_delaunay
    implicit none

    type(mesh_t) :: mesh
    real(dp) :: points(2, 4)
    integer, allocatable :: constraint_segments(:,:)
    integer :: i, j, v1, v2, v3
    
    write(*,*) "=== Edge Existence Debug ==="
    
    ! Create square vertices
    points(:, 1) = [0.0_dp, 0.0_dp]
    points(:, 2) = [1.0_dp, 0.0_dp]
    points(:, 3) = [1.0_dp, 1.0_dp]
    points(:, 4) = [0.0_dp, 1.0_dp]
    
    allocate(constraint_segments(2, 0))
    call constrained_delaunay_triangulate(points, constraint_segments, mesh)
    
    write(*,*) "Total vertices:", mesh%npoints
    write(*,*) "Valid triangles:"
    
    do i = 1, mesh%ntriangles
        if (mesh%triangles(i)%valid) then
            v1 = mesh%triangles(i)%vertices(1)
            v2 = mesh%triangles(i)%vertices(2) 
            v3 = mesh%triangles(i)%vertices(3)
            
            write(*,'(A,I0,A,3I0,A,3I0)') "  Triangle ", i, ": vertices=", v1, v2, v3, &
                " IDs=", mesh%points(v1)%id, mesh%points(v2)%id, mesh%points(v3)%id
                
            ! Check each edge of this triangle
            write(*,'(A,I0,A,I0,A,L)') "    Edge ", v1, "-", v2, ": ", constraint_edge_exists(mesh, v1, v2)
            write(*,'(A,I0,A,I0,A,L)') "    Edge ", v2, "-", v3, ": ", constraint_edge_exists(mesh, v2, v3)
            write(*,'(A,I0,A,I0,A,L)') "    Edge ", v3, "-", v1, ": ", constraint_edge_exists(mesh, v3, v1)
        end if
    end do

end program test_edge_debug