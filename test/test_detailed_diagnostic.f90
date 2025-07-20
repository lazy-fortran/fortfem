program test_detailed_diagnostic
    ! Detailed diagnostic test to understand Bowyer-Watson behavior
    use fortfem_kinds
    use delaunay_types
    use bowyer_watson
    implicit none

    type(mesh_t) :: mesh
    real(dp) :: points(2, 3)
    integer :: i, super_tri_idx, point_idx
    
    write(*,*) "=== Detailed Bowyer-Watson Diagnostic ==="
    
    ! Create simple triangle: (0,0), (1,0), (0.5,1)
    points(:, 1) = [0.0_dp, 0.0_dp]
    points(:, 2) = [1.0_dp, 0.0_dp]
    points(:, 3) = [0.5_dp, 1.0_dp]
    
    ! Initialize mesh with appropriate size
    call create_mesh(mesh, size(points, 2) + 3, 6 * size(points, 2), 3 * size(points, 2))
    
    ! Create super-triangle containing all points
    call create_super_triangle(points, mesh, super_tri_idx)
    
    ! Add input points to mesh
    do i = 1, size(points, 2)
        point_idx = add_point(mesh, points(1, i), points(2, i), i)
    end do
    
    ! Insert each point using Bowyer-Watson algorithm
    do i = 4, mesh%npoints  ! Start after super-triangle vertices
        call insert_point(mesh, i)
    end do
    
    ! DON'T remove super-triangle for debugging
    ! call remove_super_triangle(mesh)
    
    write(*,*) "After Delaunay triangulation:"
    write(*,*) "  Vertices:", mesh%npoints
    write(*,*) "  Triangles:", mesh%ntriangles
    
    write(*,*) "All triangles (before super-triangle removal):"
    do i = 1, mesh%ntriangles
        if (mesh%triangles(i)%valid) then
            write(*,'(A,I0,A,3I0,A,3I0)') "  Triangle ", i, ": vertices=", &
                mesh%triangles(i)%vertices, " vertex_IDs=", &
                mesh%points(mesh%triangles(i)%vertices(1))%id, &
                mesh%points(mesh%triangles(i)%vertices(2))%id, &
                mesh%points(mesh%triangles(i)%vertices(3))%id
        else
            write(*,'(A,I0,A)') "  Triangle ", i, ": [invalid]"
        end if
    end do
    
    write(*,*) "Super-triangle vertices:", mesh%super_vertices

end program test_detailed_diagnostic