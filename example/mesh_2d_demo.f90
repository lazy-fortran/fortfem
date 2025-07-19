program mesh_2d_demo
    ! Create and analyze 2D meshes
    use fortfem
    use fortplot
    implicit none
    
    type(mesh_2d_t) :: mesh
    real(dp), allocatable :: areas(:)
    
    ! Create mesh
    call create_unit_square_mesh(mesh, n=10)
    
    ! Print statistics
    print '(a)', "Mesh statistics:"
    print '(a,i0)', "  Vertices:  ", mesh%n_vertices
    print '(a,i0)', "  Triangles: ", mesh%n_triangles
    print '(a,i0)', "  Edges:     ", mesh%n_edges
    
    ! Analyze quality
    areas = mesh%compute_areas()
    print '(a,f8.6)', "  Min area: ", minval(areas)
    print '(a,f8.6)', "  Max area: ", maxval(areas)
    
    ! Save mesh
    call write_vtk("mesh.vtk", mesh)
    print *, "Mesh saved to mesh.vtk"
    
    ! Simple plot
    call figure()
    call plot_mesh_edges(mesh)
    call xlabel('x')
    call ylabel('y')
    call title('Unit Square Mesh')
    call savefig('mesh.png')
    
end program mesh_2d_demo