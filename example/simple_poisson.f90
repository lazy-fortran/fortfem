program simple_poisson
    ! Simple example solving the Poisson equation
    ! -Δu = f in Ω = [0,1]²
    ! u = 0 on ∂Ω
    
    use fortfem_kinds
    use fortfem_mesh_2d
    use basis_p1_2d_module
    implicit none
    
    ! LAPACK interface
    interface
        subroutine dgesv(n, nrhs, a, lda, ipiv, b, ldb, info)
            integer, intent(in) :: n, nrhs, lda, ldb
            integer, intent(out) :: info
            integer, intent(out) :: ipiv(*)
            double precision, intent(inout) :: a(lda,*), b(ldb,*)
        end subroutine dgesv
    end interface
    
    type(mesh_2d_t) :: mesh
    type(basis_p1_2d_t) :: basis
    real(dp), allocatable :: K(:,:), F(:), u(:)
    real(dp) :: h, x, y
    integer :: n, ndof, i, j, e, v1, v2, v3
    real(dp) :: x1, y1, x2, y2, x3, y3
    real(dp) :: area, grad_phi_i(2), grad_phi_j(2)
    real(dp) :: a(2,2), det_a
    real(dp) :: b(3), c(3), K_elem(3,3)
    real(dp) :: max_u, min_u, center_u
    integer, allocatable :: ipiv(:)
    integer :: info
    
    print *, "=== Simple Poisson Equation Example ==="
    print *, ""
    print *, "Solving: -Δu = 1 on [0,1]²"
    print *, "with u = 0 on the boundary"
    print *, ""
    
    ! Create mesh
    n = 21  ! 21x21 grid
    h = 1.0_dp / real(n-1, dp)
    call mesh%create_rectangular(n, n, 0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp)
    call mesh%build_connectivity()
    call mesh%find_boundary()
    
    print *, "Mesh statistics:"
    print *, "  Vertices:", mesh%n_vertices
    print *, "  Elements:", mesh%n_triangles
    print *, "  h =", h
    print *, ""
    
    ! Number of degrees of freedom
    ndof = mesh%n_vertices
    
    ! Allocate matrices and vectors
    allocate(K(ndof, ndof), F(ndof), u(ndof))
    K = 0.0_dp
    F = 0.0_dp
    
    ! Assemble stiffness matrix and load vector
    print *, "Assembling stiffness matrix..."
    
    do e = 1, mesh%n_triangles
        ! Get vertices of triangle
        v1 = mesh%triangles(1, e)
        v2 = mesh%triangles(2, e)
        v3 = mesh%triangles(3, e)
        
        ! Get coordinates
        x1 = mesh%vertices(1, v1)
        y1 = mesh%vertices(2, v1)
        x2 = mesh%vertices(1, v2)
        y2 = mesh%vertices(2, v2)
        x3 = mesh%vertices(1, v3)
        y3 = mesh%vertices(2, v3)
        
        ! Compute area
        area = 0.5_dp * abs((x2-x1)*(y3-y1) - (x3-x1)*(y2-y1))
        
        ! For P1 elements on triangles, the element stiffness matrix is:
        ! K_ij = (1/4A) * (b_i*b_j + c_i*c_j)
        ! where b = [y2-y3, y3-y1, y1-y2] and c = [x3-x2, x1-x3, x2-x1]
        
        ! Compute b and c coefficients
        b(1) = y2 - y3
        b(2) = y3 - y1
        b(3) = y1 - y2
        
        c(1) = x3 - x2
        c(2) = x1 - x3
        c(3) = x2 - x1
        
        ! Assemble element stiffness matrix
        do i = 1, 3
            do j = 1, 3
                K_elem(i,j) = (b(i)*b(j) + c(i)*c(j))/(4.0_dp*area)
            end do
        end do
        
        ! Add element contributions to global matrix
        K(v1,v1) = K(v1,v1) + K_elem(1,1)
        K(v1,v2) = K(v1,v2) + K_elem(1,2)
        K(v1,v3) = K(v1,v3) + K_elem(1,3)
        K(v2,v1) = K(v2,v1) + K_elem(2,1)
        K(v2,v2) = K(v2,v2) + K_elem(2,2)
        K(v2,v3) = K(v2,v3) + K_elem(2,3)
        K(v3,v1) = K(v3,v1) + K_elem(3,1)
        K(v3,v2) = K(v3,v2) + K_elem(3,2)
        K(v3,v3) = K(v3,v3) + K_elem(3,3)
        
        ! Load vector F(i) += ∫_T f φᵢ dx
        ! For f = 1 and linear elements: ∫_T φᵢ dx = area/3
        F(v1) = F(v1) + area/3.0_dp
        F(v2) = F(v2) + area/3.0_dp
        F(v3) = F(v3) + area/3.0_dp
    end do
    
    ! Apply Dirichlet boundary conditions
    print *, "Applying boundary conditions..."
    
    do i = 1, ndof
        if (mesh%is_boundary_vertex(i)) then
            ! Set row and column to zero
            K(i,:) = 0.0_dp
            K(:,i) = 0.0_dp
            ! Set diagonal to 1
            K(i,i) = 1.0_dp
            ! Set RHS to boundary value (0 in this case)
            F(i) = 0.0_dp
        end if
    end do
    
    ! Solve the system Ku = F
    print *, "Solving linear system..."
    
    ! Copy F to u (LAPACK overwrites RHS with solution)
    u = F
    
    ! Allocate pivot array
    allocate(ipiv(ndof))
    
    ! Call LAPACK solver (DGESV)
    call dgesv(ndof, 1, K, ndof, ipiv, u, ndof, info)
    
    if (info /= 0) then
        print *, "Error: Linear system solve failed with info =", info
        stop
    end if
    
    ! Find maximum value of solution
    print *, ""
    print *, "Solution statistics:"
    print *, "  Max u =", maxval(u)
    print *, "  Min u =", minval(u)
    print *, ""
    
    ! For the problem -Δu = 1 with u=0 on boundary,
    ! the maximum occurs at the center and is approximately 1/8
    print *, "Note: For -Δu = 1 on [0,1]², the exact maximum is 1/8 = 0.125"
    print *, ""
    
    ! Output solution at center point
    ! Find vertex closest to center (0.5, 0.5)
    j = 1
    do i = 2, ndof
        x = mesh%vertices(1, i)
        y = mesh%vertices(2, i)
        if ((x-0.5_dp)**2 + (y-0.5_dp)**2 < &
            (mesh%vertices(1,j)-0.5_dp)**2 + (mesh%vertices(2,j)-0.5_dp)**2) then
            j = i
        end if
    end do
    
    print *, "Solution at center vertex:"
    print *, "  Position: (", mesh%vertices(1,j), ",", mesh%vertices(2,j), ")"
    print *, "  Value: u =", u(j)
    print *, ""
    
    ! Cleanup
    deallocate(K, F, u, ipiv)
    call mesh%destroy()
    
    print *, "Example completed successfully!"
    
end program simple_poisson