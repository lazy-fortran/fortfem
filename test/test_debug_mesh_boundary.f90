program test_debug_mesh_boundary
    ! Debug mesh boundary edge identification
    use fortfem_kinds, only: dp
    use fortfem_mesh_2d
    use fortfem_hcurl_space
    implicit none
    
    type(mesh_2d_t) :: mesh
    type(hcurl_space_t) :: space
    integer :: n, nx, ny, edge_idx, n_dofs, n_boundary
    
    ! Test with the failing mesh size: n=4 (level=1)
    n = 4
    nx = n + 1  ! 5
    ny = n + 1  ! 5
    
    call mesh%create_rectangular(nx, ny, 0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp)
    call mesh%build_edge_connectivity()
    call mesh%build_edge_dof_numbering()
    
    call space%init(mesh)
    n_dofs = space%get_n_dofs()
    
    print *, "Mesh info for n=", n
    print *, "  vertices =", mesh%n_vertices, "triangles =", mesh%n_triangles
    print *, "  edges =", mesh%n_edges, "DOFs =", n_dofs
    
    ! Count and list boundary edges
    n_boundary = 0
    print *, ""
    print *, "Boundary edge analysis:"
    do edge_idx = 1, n_dofs
        if (mesh%is_boundary_edge(edge_idx)) then
            n_boundary = n_boundary + 1
            write(*, '(A,I0,A,I0)') "  Edge ", edge_idx, " is boundary (DOF ", edge_idx, ")"
        end if
    end do
    
    print *, ""
    print *, "Total boundary edges:", n_boundary
    print *, "Expected for 5x5 grid: 16 boundary edges"
    
    call space%destroy()
    call mesh%destroy()
end program test_debug_mesh_boundary