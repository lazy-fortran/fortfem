program p2_poisson
    ! Solve Poisson equation with P2 (quadratic) elements
    use fortfem
    implicit none
    
    type(mesh_2d_t) :: mesh
    type(sparse_matrix_t) :: A
    real(dp), allocatable :: u(:), f(:)
    
    ! Create mesh
    call create_unit_square_mesh(mesh, n=10)
    
    ! Assemble with P2 elements (6 DOFs per triangle)
    call assemble_poisson_p2(mesh, A, f)
    
    ! Apply boundary conditions
    call apply_zero_bc_p2(mesh, A, f)
    
    ! Solve
    call solve_sparse(A, f, u)
    
    ! Output
    print '(a,i0)', "P2 DOFs: ", size(u)
    print '(a,f8.5)', "Max solution: ", maxval(u)
    
end program p2_poisson