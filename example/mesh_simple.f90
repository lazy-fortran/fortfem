program mesh_simple
    ! Demonstrate mesh creation and visualization
    use fortfem
    implicit none
    
    type(mesh_2d_t) :: mesh
    real(dp), allocatable :: dummy(:)
    
    print *, "Simple Mesh Demo"
    print *, "================"
    
    ! Create unit square mesh
    call create_unit_square_mesh(mesh, n=10)
    
    ! Print info
    print '(a,i0)', "Vertices:  ", mesh%n_vertices
    print '(a,i0)', "Triangles: ", mesh%n_triangles
    print '(a,i0)', "Edges:     ", mesh%n_edges
    
    ! Save mesh (with dummy data for VTK format)
    allocate(dummy(mesh%n_vertices))
    dummy = 0.0_dp
    call write_vtk("mesh.vtk", mesh, dummy)
    print *, "Mesh saved to mesh.vtk"
    
    ! Clean up
    call mesh%destroy()
    deallocate(dummy)
    
end program mesh_simple