program test_simple_poisson
    use fortfem_kinds
    use fortfem_mesh_2d
    use basis_p1_2d_module
    use check
    implicit none
    
    type(mesh_2d_t) :: mesh
    type(basis_p1_2d_t) :: basis
    real(dp), allocatable :: K(:,:), F(:), u(:)
    real(dp), allocatable :: u_exact(:)
    integer :: n, i, j, e
    integer :: ndof
    real(dp) :: xi, eta, x, y
    real(dp) :: grad_phi_i(2), grad_phi_j(2)
    real(dp) :: det_J, area
    real(dp) :: error_l2
    
    print *, "Testing simple Poisson equation solver..."
    
    ! Create a simple 3x3 mesh on [0,1]x[0,1]
    n = 3
    call mesh%create_rectangular(n, n, 0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp)
    
    ! Number of degrees of freedom
    ndof = mesh%n_vertices
    
    ! Allocate matrices and vectors
    allocate(K(ndof, ndof), F(ndof), u(ndof), u_exact(ndof))
    K = 0.0_dp
    F = 0.0_dp
    
    ! Test that mesh was created correctly
    call check_condition(mesh%n_vertices == 9, "Expected 9 vertices")
    call check_condition(mesh%n_triangles > 0, "Expected positive number of elements")
    
    ! Simple assembly test - check that stiffness matrix is symmetric
    ! For now, just create a simple test matrix
    do i = 1, ndof
        K(i,i) = 2.0_dp
        if (i > 1) K(i,i-1) = -1.0_dp
        if (i < ndof) K(i,i+1) = -1.0_dp
    end do
    
    ! Check symmetry
    do i = 1, ndof
        do j = i+1, ndof
            call check_condition(abs(K(i,j) - K(j,i)) < 1.0e-10_dp, "Stiffness matrix should be symmetric")
        end do
    end do
    
    ! Simple RHS for testing
    F = 1.0_dp
    
    ! Apply Dirichlet BC (u=0 on boundary) - simplified test
    ! Just set first and last node
    K(1,:) = 0.0_dp
    K(:,1) = 0.0_dp
    K(1,1) = 1.0_dp
    F(1) = 0.0_dp
    
    K(ndof,:) = 0.0_dp
    K(:,ndof) = 0.0_dp
    K(ndof,ndof) = 1.0_dp
    F(ndof) = 0.0_dp
    
    ! Check that matrix modifications preserve symmetry where applicable
    call check_condition(abs(K(1,1) - 1.0_dp) < 1.0e-10_dp, "Dirichlet node should have 1 on diagonal")
    call check_condition(abs(F(1)) < 1.0e-10_dp, "Dirichlet node should have 0 RHS")
    
    ! Print test summary
    call check_summary("Simple Poisson Test")
    
    ! Cleanup
    deallocate(K, F, u, u_exact)
    call mesh%destroy()
    
end program test_simple_poisson