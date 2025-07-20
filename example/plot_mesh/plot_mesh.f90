program plot_mesh_example
    ! Example demonstrating mesh plotting functionality in FortFEM
    
    use fortfem_api
    implicit none
    
    type(mesh_t) :: mesh
    integer :: n_refinements
    
    ! Create unit square mesh with different refinement levels
    n_refinements = 5
    mesh = unit_square_mesh(n_refinements)
    
    ! Plot basic mesh
    call plot(mesh, filename="mesh_basic.png", title="Unit Square Mesh")
    
    ! Create finer mesh
    n_refinements = 10
    mesh = unit_square_mesh(n_refinements)
    
    ! Plot finer mesh
    call plot(mesh, filename="mesh_fine.png", title="Refined Unit Square Mesh")
    
    ! Create coarse mesh for clarity
    n_refinements = 3
    mesh = unit_square_mesh(n_refinements)
    
    ! Plot with custom title
    call plot(mesh, filename="mesh_coarse.png", title="Coarse Mesh (3x3)")
    
    ! Clean up
    call mesh%destroy()
    
    print *, "Mesh plotting examples completed!"
    print *, "Generated files:"
    print *, "  - mesh_basic.png"
    print *, "  - mesh_fine.png"
    print *, "  - mesh_coarse.png"
    
end program plot_mesh_example