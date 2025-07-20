program test_constraint_diagnostic
    ! Detailed diagnostic for constraint edge insertion
    use fortfem_kinds
    use delaunay_types
    use constrained_delaunay
    implicit none

    type(mesh_t) :: mesh
    real(dp) :: points(2, 4)
    integer, allocatable :: constraint_segments(:,:)
    integer :: i, j, valid_before, valid_after

    write(*,*) "=== Constraint Edge Insertion Diagnostic ==="
    
    ! Create square vertices
    points(:, 1) = [0.0_dp, 0.0_dp]
    points(:, 2) = [1.0_dp, 0.0_dp]
    points(:, 3) = [1.0_dp, 1.0_dp]
    points(:, 4) = [0.0_dp, 1.0_dp]
    
    ! First triangulate without constraints
    allocate(constraint_segments(2, 0))
    call constrained_delaunay_triangulate(points, constraint_segments, mesh)
    
    valid_before = 0
    do i = 1, mesh%ntriangles
        if (mesh%triangles(i)%valid) valid_before = valid_before + 1
    end do
    
    write(*,*) "Triangles before constraint:", valid_before
    
    ! Show triangles before constraint
    do i = 1, mesh%ntriangles
        if (mesh%triangles(i)%valid) then
            write(*,'(A,I0,A,3I0)') "  Triangle ", i, ": vertices=", mesh%triangles(i)%vertices
        end if
    end do
    
    ! Now add diagonal constraint (vertex 4 to vertex 6 - real vertices)
    deallocate(constraint_segments)
    allocate(constraint_segments(2, 1))
    constraint_segments(:, 1) = [4, 6]  ! Constraint from first to third real vertex
    call constrained_delaunay_triangulate(points, constraint_segments, mesh)
    
    valid_after = 0
    do i = 1, mesh%ntriangles
        if (mesh%triangles(i)%valid) valid_after = valid_after + 1
    end do
    
    write(*,*) "Triangles after constraint:", valid_after
    
    ! Show triangles after constraint
    do i = 1, mesh%ntriangles
        if (mesh%triangles(i)%valid) then
            write(*,'(A,I0,A,3I0)') "  Triangle ", i, ": vertices=", mesh%triangles(i)%vertices
        end if
    end do
    
    ! Check if constraint edge exists
    write(*,*) "Constraint edge 4-6 exists:", constraint_edge_exists(mesh, 4, 6)

end program test_constraint_diagnostic