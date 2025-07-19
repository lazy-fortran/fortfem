program test_simple_api
    ! Test the new simplified API
    use fortfem
    use check
    implicit none
    
    type(mesh_2d_t) :: mesh
    type(sparse_matrix_t) :: A
    real(dp), allocatable :: u(:), f(:)
    integer :: info
    
    print *, "Testing simplified FortFEM API"
    print *, "=============================="
    
    ! Test mesh creation
    call create_unit_square_mesh(mesh, n=4)
    call check_condition(mesh%n_vertices > 0, "Mesh has vertices")
    call check_condition(mesh%n_triangles > 0, "Mesh has triangles")
    
    ! Test assembly
    call assemble_poisson_2d(mesh, A, f)
    call check_condition(A%n > 0, "Matrix assembled")
    call check_condition(size(f) > 0, "RHS assembled")
    
    ! Test boundary conditions
    call apply_zero_bc(mesh, A, f)
    call check_condition(.true., "BC applied")
    
    ! Test simple solver
    allocate(u(size(f)))
    u = f
    call solve_lapack_dense(A, u, info)
    call check_condition(info == 0, "Solver converged")
    call check_condition(maxval(abs(u)) > 0.0_dp, "Solution non-trivial")
    
    ! Test VTK output
    call write_vtk("test_solution.vtk", mesh, u)
    call check_condition(.true., "VTK written")
    
    ! Clean up
    call mesh%destroy()
    call A%destroy()
    deallocate(u, f)
    
    call check_summary("Simple API Test")
    
end program test_simple_api