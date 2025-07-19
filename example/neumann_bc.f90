program neumann_bc
    ! Solve -∆u = f with mixed boundary conditions
    ! u = 0 on left/right, ∂u/∂n = g on top/bottom
    use fortfem
    implicit none
    
    type(mesh_2d_t) :: mesh
    type(sparse_matrix_t) :: A
    real(dp), allocatable :: u(:), f(:)
    
    ! Create mesh
    call create_unit_square_mesh(mesh, n=20)
    
    ! Assemble Poisson system
    call assemble_poisson_2d(mesh, A, f)
    
    ! Add Neumann BC contribution: ∫gv ds on top/bottom
    call add_neumann_bc(mesh, f, top=1.0_dp, bottom=-1.0_dp)
    
    ! Apply Dirichlet BC on left/right
    call apply_dirichlet_bc(mesh, A, f, left=0.0_dp, right=0.0_dp)
    
    ! Solve
    call solve_sparse(A, f, u)
    
    ! Output
    print '(a,f8.5)', "Solution range: [", minval(u), ",", maxval(u), "]"
    
end program neumann_bc