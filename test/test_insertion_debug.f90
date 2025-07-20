program test_insertion_debug
    ! Debug point insertion step by step
    use fortfem_kinds
    use delaunay_types
    use bowyer_watson
    implicit none

    type(mesh_t) :: mesh
    real(dp) :: points(2, 3)
    integer :: i, super_tri_idx, point_idx
    
    write(*,*) "=== Point Insertion Debug ==="
    
    ! Create simple triangle: (0,0), (1,0), (0.5,1)
    points(:, 1) = [0.0_dp, 0.0_dp]
    points(:, 2) = [1.0_dp, 0.0_dp]
    points(:, 3) = [0.5_dp, 1.0_dp]
    
    ! Initialize mesh
    call create_mesh(mesh, 6, 18, 9)
    
    ! Create super-triangle
    write(*,*) "Creating super-triangle..."
    call create_super_triangle(points, mesh, super_tri_idx)
    write(*,*) "Super-triangle created. Triangle index:", super_tri_idx
    write(*,*) "Super-triangle vertices:", mesh%super_vertices
    
    ! Add input points
    write(*,*) "Adding input points..."
    do i = 1, 3
        point_idx = add_point(mesh, points(1, i), points(2, i), i)
        write(*,'(A,I0,A,2F8.3,A,I0)') "  Point ", i, " at (", points(:, i), ") -> vertex ", point_idx
    end do
    
    write(*,*) "After adding points:"
    write(*,*) "  Total vertices:", mesh%npoints
    write(*,*) "  Total triangles:", mesh%ntriangles
    
    ! Insert points one by one with detailed output
    do i = 4, mesh%npoints
        write(*,*) ""
        write(*,'(A,I0,A,2F8.3,A,I0)') "Inserting vertex ", i, " at (", &
            mesh%points(i)%x, mesh%points(i)%y, ") with ID ", mesh%points(i)%id
        
        ! Count valid triangles before insertion
        write(*,*) "  Valid triangles before insertion:", count_valid_triangles(mesh)
        
        call insert_point(mesh, i)
        
        ! Count valid triangles after insertion
        write(*,*) "  Valid triangles after insertion:", count_valid_triangles(mesh)
        
        ! Show all triangles
        call show_all_triangles(mesh)
    end do

contains

    function count_valid_triangles(mesh) result(count)
        type(mesh_t), intent(in) :: mesh
        integer :: count, j
        count = 0
        do j = 1, mesh%ntriangles
            if (mesh%triangles(j)%valid) count = count + 1
        end do
    end function

    subroutine show_all_triangles(mesh)
        type(mesh_t), intent(in) :: mesh
        integer :: j
        
        write(*,*) "  Current triangles:"
        do j = 1, mesh%ntriangles
            if (mesh%triangles(j)%valid) then
                write(*,'(A,I0,A,3I0,A,3I0)') "    Triangle ", j, ": vertices=", &
                    mesh%triangles(j)%vertices, " IDs=", &
                    mesh%points(mesh%triangles(j)%vertices(1))%id, &
                    mesh%points(mesh%triangles(j)%vertices(2))%id, &
                    mesh%points(mesh%triangles(j)%vertices(3))%id
            end if
        end do
    end subroutine

end program test_insertion_debug