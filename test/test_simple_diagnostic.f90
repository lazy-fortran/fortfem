program test_simple_diagnostic
    ! Diagnostic test to understand triangulation behavior
    use fortfem_kinds
    use delaunay_types
    use constrained_delaunay
    implicit none

    type(mesh_t) :: mesh
    real(dp) :: points(2, 3)
    integer, allocatable :: constraint_segments(:,:)
    integer :: i
    
    write(*,*) "=== Simple Triangulation Diagnostic ==="
    
    ! Create simple triangle: (0,0), (1,0), (0.5,1)
    points(:, 1) = [0.0_dp, 0.0_dp]
    points(:, 2) = [1.0_dp, 0.0_dp]
    points(:, 3) = [0.5_dp, 1.0_dp]
    
    allocate(constraint_segments(2, 0))
    call constrained_delaunay_triangulate(points, constraint_segments, mesh)
    
    write(*,*) "Input points: 3"
    write(*,*) "Output vertices:", mesh%npoints
    write(*,*) "Output triangles:", mesh%ntriangles
    
    write(*,*) "Triangles:"
    do i = 1, mesh%ntriangles
        if (mesh%triangles(i)%valid) then
            write(*,'(A,I0,A,3I0)') "  Triangle ", i, ": ", mesh%triangles(i)%vertices
        else
            write(*,'(A,I0,A)') "  Triangle ", i, ": [invalid]"
        end if
    end do
    
    write(*,*) "Vertices:"
    do i = 1, mesh%npoints
        write(*,'(A,I0,A,2F8.3)') "  Vertex ", i, ": ", mesh%points(i)%x, mesh%points(i)%y
    end do

end program test_simple_diagnostic