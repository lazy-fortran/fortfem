program test_debug_dof_mapping
    use fortfem_kinds, only: dp
    use fortfem_mesh_2d
    use fortfem_hcurl_space
    implicit none

    type(mesh_2d_t) :: mesh
    type(hcurl_space_t) :: space
    
    ! Create simple 2x2 mesh
    call mesh%create_rectangular(3, 3, 0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp)
    call mesh%build_edge_connectivity()
    call mesh%build_edge_dof_numbering()
    
    call space%init(mesh)
    
    print *, "=== Mesh Debug Info ==="
    print *, "Triangles:", mesh%n_triangles
    print *, "Edges:", mesh%n_edges
    print *, "Vertices:", mesh%n_vertices
    print *, "DOFs:", space%get_n_dofs()
    
    call debug_triangles(mesh)
    call debug_edges(mesh)
    call debug_triangle_edge_mapping(mesh)
    
    call space%destroy()
    call mesh%destroy()

contains

    subroutine debug_triangles(mesh)
        type(mesh_2d_t), intent(in) :: mesh
        integer :: t
        
        print *, ""
        print *, "=== Triangle Vertices ==="
        do t = 1, min(4, mesh%n_triangles)
            print '("Triangle", I2, ": vertices", 3I3)', t, mesh%triangles(:, t)
            print '("  Coords: (", F4.1, ",", F4.1, ") (", F4.1, ",", F4.1, ") (", F4.1, ",", F4.1, ")")', &
                mesh%vertices(1, mesh%triangles(1, t)), mesh%vertices(2, mesh%triangles(1, t)), &
                mesh%vertices(1, mesh%triangles(2, t)), mesh%vertices(2, mesh%triangles(2, t)), &
                mesh%vertices(1, mesh%triangles(3, t)), mesh%vertices(2, mesh%triangles(3, t))
        end do
    end subroutine
    
    subroutine debug_edges(mesh)
        type(mesh_2d_t), intent(in) :: mesh
        integer :: e, vertices(2)
        real(dp) :: length, tangent(2)
        
        print *, ""
        print *, "=== Edge Information ==="
        do e = 1, min(8, mesh%n_edges)
            call mesh%get_edge_vertices(e, vertices)
            call mesh%get_edge_length_tangent(e, length, tangent)
            print '("Edge", I2, ": vertices", 2I3, " DOF", I3, " len=", F4.2, " tan=(", F5.2, ",", F5.2, ")")', &
                e, vertices, mesh%edge_to_dof(e), length, tangent
        end do
    end subroutine
    
    subroutine debug_triangle_edge_mapping(mesh)
        type(mesh_2d_t), intent(in) :: mesh
        integer :: t, triangle_dofs(3), e
        
        print *, ""
        print *, "=== Triangle-Edge DOF Mapping ==="
        do t = 1, min(4, mesh%n_triangles)
            call mesh%get_triangle_edge_dofs(t, triangle_dofs)
            print '("Triangle", I2, ": edge DOFs", 3I3)', t, triangle_dofs
            
            do e = 1, 3
                print '("  Local edge", I1, " -> DOF", I2)', e, triangle_dofs(e)
            end do
        end do
    end subroutine

end program test_debug_dof_mapping