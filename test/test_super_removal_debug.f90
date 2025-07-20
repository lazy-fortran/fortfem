program test_super_removal_debug
    ! Debug super-triangle removal
    use fortfem_kinds
    use delaunay_types
    use bowyer_watson
    implicit none

    type(mesh_t) :: mesh
    real(dp) :: points(2, 3)
    integer :: i, super_tri_idx, point_idx
    logical :: contains_super_vertex
    integer :: j, v, valid_count
    
    write(*,*) "=== Super-Triangle Removal Debug ==="
    
    ! Create simple triangle and run full Bowyer-Watson
    points(:, 1) = [0.0_dp, 0.0_dp]
    points(:, 2) = [1.0_dp, 0.0_dp]
    points(:, 3) = [0.5_dp, 1.0_dp]
    
    call create_mesh(mesh, 6, 18, 9)
    call create_super_triangle(points, mesh, super_tri_idx)
    
    do i = 1, 3
        point_idx = add_point(mesh, points(1, i), points(2, i), i)
    end do
    
    do i = 4, mesh%npoints
        call insert_point(mesh, i)
    end do
    
    write(*,*) "Super-triangle vertices:", mesh%super_vertices
    write(*,*) "Total triangles before removal:", mesh%ntriangles
    
    ! Analyze each triangle
    valid_count = 0
    do i = 1, mesh%ntriangles
        if (.not. mesh%triangles(i)%valid) cycle
        
        contains_super_vertex = .false.
        write(*,'(A,I0,A,3I0)', advance='no') "Triangle ", i, ": vertices=", mesh%triangles(i)%vertices
        
        do j = 1, 3
            v = mesh%triangles(i)%vertices(j)
            if (any(mesh%super_vertices == v)) then
                contains_super_vertex = .true.
            end if
        end do
        
        write(*,'(A,3I0)', advance='no') " IDs=", &
            mesh%points(mesh%triangles(i)%vertices(1))%id, &
            mesh%points(mesh%triangles(i)%vertices(2))%id, &
            mesh%points(mesh%triangles(i)%vertices(3))%id
            
        if (contains_super_vertex) then
            write(*,*) " [USES SUPER-TRIANGLE - REMOVE]"
        else
            write(*,*) " [KEEP]"
            valid_count = valid_count + 1
        end if
    end do
    
    write(*,*) "Triangles to keep:", valid_count
    
    ! Now do actual removal
    call remove_super_triangle(mesh)
    
    ! Count final valid triangles
    valid_count = 0
    do i = 1, mesh%ntriangles
        if (mesh%triangles(i)%valid) valid_count = valid_count + 1
    end do
    write(*,*) "Final valid triangles:", valid_count

end program test_super_removal_debug