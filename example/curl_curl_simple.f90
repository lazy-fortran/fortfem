program curl_curl_simple
    ! Solve curl(curl(E)) + E = J with E×n = 0 on boundary
    use fortfem
    implicit none
    
    type(mesh_2d_t) :: mesh
    type(sparse_matrix_t) :: A
    real(dp), allocatable :: E(:), J(:)
    
    ! Create mesh
    call create_unit_square_mesh(mesh, n=20)
    
    ! Assemble curl-curl system
    call assemble_curl_curl(mesh, A, J)
    
    ! Apply essential boundary conditions E×n = 0
    call apply_tangential_bc(mesh, A, J)
    
    ! Solve
    call solve_sparse(A, J, E)
    
    ! Output
    print '(a,f8.5)', "Max |E|: ", maxval(abs(E))
    call write_vtk("electric_field.vtk", mesh, E, vector=.true.)
    
end program curl_curl_simple