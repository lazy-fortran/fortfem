program poisson_simple
    ! Solve -∆u = f on unit square with u = 0 on boundary
    use fortfem
    implicit none
    
    type(mesh_2d_t) :: mesh
    type(sparse_matrix_t) :: A
    real(dp), allocatable :: u(:), f(:)
    integer :: info
    
    print *, "Simple Poisson Solver"
    print *, "===================="
    
    ! Create mesh
    call create_unit_square_mesh(mesh, n=20)
    print '(a,i0)', "Mesh vertices: ", mesh%n_vertices
    
    ! Assemble system
    call assemble_poisson_2d(mesh, A, f)
    print *, "System assembled"
    
    ! Apply boundary conditions
    call apply_zero_bc(mesh, A, f)
    print *, "Boundary conditions applied"
    
    ! Solve (simple LAPACK solver for dense)
    allocate(u(size(f)))
    u = f
    call solve_lapack_dense(A, u, info)
    
    if (info == 0) then
        print *, "Solution computed successfully"
        print '(a,f8.5)', "Max solution value: ", maxval(abs(u))
        call write_vtk("poisson_solution.vtk", mesh, u)
        print *, "Solution saved to poisson_solution.vtk"
    else
        print *, "Solver failed with info = ", info
    end if
    
    ! Clean up
    call mesh%destroy()
    call A%destroy()
    
end program poisson_simple