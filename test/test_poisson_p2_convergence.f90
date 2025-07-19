program test_poisson_p2_convergence
    ! Test convergence of P2 elements for Poisson equation
    use fortfem_kinds
    use fortfem_mesh_2d
    use basis_p2_2d_module
    use fortfem_sparse_matrix
    use fortfem_solver_interface
    use fortfem_gauss_quadrature_2d
    implicit none
    
    integer :: n, i, j, k, t, info
    integer :: n_meshes
    real(dp), allocatable :: h_values(:), l2_errors(:), h1_errors(:)
    type(mesh_2d_t) :: mesh
    type(basis_p2_2d_t) :: basis
    type(gauss_quadrature_triangle_t) :: quad
    type(triplet_matrix_t) :: A_triplet
    type(csr_matrix_t) :: A_csr
    type(lapack_dense_solver_t) :: solver
    real(dp), allocatable :: b(:), u(:), u_exact(:)
    real(dp) :: pi_val
    
    print *, "=== Testing P2 Convergence for Poisson Equation ==="
    print *, ""
    
    pi_val = 4.0_dp * atan(1.0_dp)
    
    ! Test on sequence of refined meshes
    n_meshes = 4
    allocate(h_values(n_meshes), l2_errors(n_meshes), h1_errors(n_meshes))
    
    do i = 1, n_meshes
        n = 2**(i+1)  ! 4, 8, 16, 32
        h_values(i) = 1.0_dp / real(n-1, dp)
        
        ! Create mesh
        call mesh%create_rectangular(n, n, 0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp)
        call mesh%build_edge_connectivity()
        
        ! Solve Poisson problem with P2 elements
        call solve_poisson_p2(mesh, basis, l2_errors(i), h1_errors(i))
        
        print '(A,I3,A,F8.5,A,E12.5,A,E12.5)', &
            "n = ", n, ", h = ", h_values(i), &
            ", L2 error = ", l2_errors(i), &
            ", H1 error = ", h1_errors(i)
        
        call mesh%destroy()
    end do
    
    ! Check convergence rates
    print *, ""
    print *, "Convergence rates:"
    do i = 2, n_meshes
        print '(A,I2,A,F6.3,A,F6.3)', &
            "Level ", i, ": L2 rate = ", &
            log(l2_errors(i-1)/l2_errors(i))/log(2.0_dp), &
            ", H1 rate = ", &
            log(h1_errors(i-1)/h1_errors(i))/log(2.0_dp)
    end do
    
    ! Expected rates: L2 ~ O(h^3), H1 ~ O(h^2) for P2 elements
    print *, ""
    print *, "Expected rates: L2 ~ 3.0, H1 ~ 2.0"
    
    deallocate(h_values, l2_errors, h1_errors)
    
contains

    subroutine solve_poisson_p2(mesh, basis, l2_error, h1_error)
        type(mesh_2d_t), intent(in) :: mesh
        type(basis_p2_2d_t), intent(in) :: basis
        real(dp), intent(out) :: l2_error, h1_error
        
        integer :: n_dofs, n_interior_dofs
        integer :: v, e, t, i, j, k, g
        integer, allocatable :: dof_map(:), interior_dof_map(:)
        real(dp), allocatable :: x(:), rhs(:)
        real(dp) :: vertices(2,3), xi, eta, w, det_j, jac(2,2)
        real(dp) :: phi_i, phi_j, grad_i(2), grad_j(2)
        real(dp) :: x_phys, y_phys, f_val
        real(dp) :: Ae(6,6), be(6)
        integer :: local_dofs(6), global_i, global_j
        logical, allocatable :: is_boundary_dof(:)
        
        ! Count DOFs: vertices + edge midpoints
        n_dofs = mesh%n_vertices + mesh%n_edges
        allocate(dof_map(n_dofs), is_boundary_dof(n_dofs))
        
        ! Mark boundary DOFs
        is_boundary_dof = .false.
        
        ! Boundary vertices
        do v = 1, mesh%n_vertices
            if (is_boundary_vertex(mesh, v)) then
                is_boundary_dof(v) = .true.
            end if
        end do
        
        ! Boundary edge midpoints
        do e = 1, mesh%n_edges
            if (mesh%is_boundary_edge(e)) then
                is_boundary_dof(mesh%n_vertices + e) = .true.
            end if
        end do
        
        ! Create interior DOF mapping
        n_interior_dofs = count(.not. is_boundary_dof)
        allocate(interior_dof_map(n_interior_dofs))
        
        k = 0
        do i = 1, n_dofs
            if (.not. is_boundary_dof(i)) then
                k = k + 1
                interior_dof_map(k) = i
                dof_map(i) = k
            else
                dof_map(i) = -1
            end if
        end do
        
        ! Initialize system
        allocate(x(n_interior_dofs), rhs(n_interior_dofs))
        x = 0.0_dp
        rhs = 0.0_dp
        
        ! Initialize quadrature for P2 (need order 4 for exact integration)
        call quad%init(4)
        
        ! Initialize sparse matrix
        call A_triplet%init(n_interior_dofs, n_interior_dofs * 36)
        
        ! Assemble system
        do t = 1, mesh%n_triangles
            ! Get element vertices
            do i = 1, 3
                vertices(:,i) = mesh%vertices(:, mesh%triangles(i,t))
            end do
            
            ! Local DOF numbers for P2 element
            ! Vertices
            local_dofs(1:3) = mesh%triangles(:,t)
            
            ! Edge midpoints
            do i = 1, 3
                e = get_triangle_edge(mesh, t, i)
                local_dofs(3+i) = mesh%n_vertices + e
            end do
            
            ! Compute element matrices
            Ae = 0.0_dp
            be = 0.0_dp
            
            call basis%compute_jacobian(vertices, jac, det_j)
            
            ! Quadrature loop
            do g = 1, quad%n_points
                xi = quad%xi(g)
                eta = quad%eta(g)
                w = quad%weights(g)
                
                call basis%transform_to_physical(xi, eta, vertices, x_phys, y_phys)
                
                ! Element stiffness matrix
                do i = 1, 6
                    grad_i = basis%grad(i, xi, eta)
                    ! Transform gradient to physical coordinates
                    call transform_gradient(grad_i, jac, det_j)
                    
                    do j = 1, 6
                        grad_j = basis%grad(j, xi, eta)
                        call transform_gradient(grad_j, jac, det_j)
                        
                        Ae(i,j) = Ae(i,j) + w * det_j * dot_product(grad_i, grad_j)
                    end do
                    
                    ! Element load vector
                    phi_i = basis%eval(i, xi, eta)
                    f_val = source_function(x_phys, y_phys)
                    be(i) = be(i) + w * det_j * f_val * phi_i
                end do
            end do
            
            ! Assemble into global system
            do i = 1, 6
                global_i = dof_map(local_dofs(i))
                if (global_i > 0) then
                    rhs(global_i) = rhs(global_i) + be(i)
                    
                    do j = 1, 6
                        global_j = dof_map(local_dofs(j))
                        if (global_j > 0) then
                            call A_triplet%add(global_i, global_j, Ae(i,j))
                        end if
                    end do
                end if
            end do
        end do
        
        ! Convert to CSR and solve
        call A_triplet%to_csr(A_csr)
        call solver%init()
        call solver%solve(A_csr, rhs, x, info)
        
        if (info /= 0) then
            print *, "Solver failed with info =", info
            stop 1
        end if
        
        ! Compute errors
        call compute_errors(mesh, basis, x, interior_dof_map, dof_map, &
                           l2_error, h1_error)
        
        ! Clean up
        deallocate(dof_map, is_boundary_dof, interior_dof_map, x, rhs)
        call A_triplet%destroy()
        call A_csr%destroy()
        call solver%destroy()
        call quad%destroy()
        
    end subroutine solve_poisson_p2
    
    subroutine compute_errors(mesh, basis, u_h, interior_dof_map, dof_map, &
                             l2_error, h1_error)
        type(mesh_2d_t), intent(in) :: mesh
        type(basis_p2_2d_t), intent(in) :: basis
        real(dp), intent(in) :: u_h(:)
        integer, intent(in) :: interior_dof_map(:), dof_map(:)
        real(dp), intent(out) :: l2_error, h1_error
        
        integer :: t, i, g, e
        real(dp) :: vertices(2,3), xi, eta, w, det_j, jac(2,2)
        real(dp) :: x_phys, y_phys
        real(dp) :: u_exact_val, u_h_val, grad_exact(2), grad_h(2)
        real(dp) :: l2_norm2, h1_norm2
        integer :: local_dofs(6), global_dof
        type(gauss_quadrature_triangle_t) :: quad
        
        ! Initialize high-order quadrature for accurate error computation
        call quad%init(5)
        
        l2_norm2 = 0.0_dp
        h1_norm2 = 0.0_dp
        
        do t = 1, mesh%n_triangles
            ! Get element vertices
            do i = 1, 3
                vertices(:,i) = mesh%vertices(:, mesh%triangles(i,t))
            end do
            
            ! Local DOF numbers
            local_dofs(1:3) = mesh%triangles(:,t)
            do i = 1, 3
                e = get_triangle_edge(mesh, t, i)
                local_dofs(3+i) = mesh%n_vertices + e
            end do
            
            call basis%compute_jacobian(vertices, jac, det_j)
            
            ! Quadrature loop
            do g = 1, quad%n_points
                xi = quad%xi(g)
                eta = quad%eta(g)
                w = quad%weights(g)
                
                call basis%transform_to_physical(xi, eta, vertices, x_phys, y_phys)
                
                ! Exact solution and gradient
                u_exact_val = exact_solution(x_phys, y_phys)
                grad_exact = exact_gradient(x_phys, y_phys)
                
                ! Numerical solution
                u_h_val = 0.0_dp
                grad_h = 0.0_dp
                
                do i = 1, 6
                    global_dof = dof_map(local_dofs(i))
                    if (global_dof > 0) then
                        u_h_val = u_h_val + u_h(global_dof) * basis%eval(i, xi, eta)
                        
                        grad_h = grad_h + u_h(global_dof) * basis%grad(i, xi, eta)
                    end if
                end do
                
                ! Transform gradient to physical coordinates
                call transform_gradient(grad_h, jac, det_j)
                
                ! Accumulate errors
                l2_norm2 = l2_norm2 + w * det_j * (u_exact_val - u_h_val)**2
                h1_norm2 = h1_norm2 + w * det_j * sum((grad_exact - grad_h)**2)
            end do
        end do
        
        l2_error = sqrt(l2_norm2)
        h1_error = sqrt(h1_norm2)
        
        call quad%destroy()
        
    end subroutine compute_errors
    
    logical function is_boundary_vertex(mesh, v)
        type(mesh_2d_t), intent(in) :: mesh
        integer, intent(in) :: v
        real(dp) :: x, y
        real(dp), parameter :: tol = 1e-10_dp
        
        x = mesh%vertices(1, v)
        y = mesh%vertices(2, v)
        
        is_boundary_vertex = (abs(x) < tol .or. abs(x - 1.0_dp) < tol .or. &
                             abs(y) < tol .or. abs(y - 1.0_dp) < tol)
    end function is_boundary_vertex
    
    integer function get_triangle_edge(mesh, t, local_edge)
        type(mesh_2d_t), intent(in) :: mesh
        integer, intent(in) :: t, local_edge
        integer :: v1, v2, e
        
        ! Get vertices of the edge
        v1 = mesh%triangles(local_edge, t)
        v2 = mesh%triangles(mod(local_edge, 3) + 1, t)
        
        ! Find global edge number
        do e = 1, mesh%n_edges
            if ((mesh%edges(1,e) == v1 .and. mesh%edges(2,e) == v2) .or. &
                (mesh%edges(1,e) == v2 .and. mesh%edges(2,e) == v1)) then
                get_triangle_edge = e
                return
            end if
        end do
        
        print *, "Error: Edge not found!"
        stop 1
    end function get_triangle_edge
    
    subroutine transform_gradient(grad_ref, jac, det_j)
        real(dp), intent(inout) :: grad_ref(2)
        real(dp), intent(in) :: jac(2,2), det_j
        real(dp) :: jac_inv(2,2), grad_phys(2)
        
        ! Compute inverse Jacobian
        jac_inv(1,1) = jac(2,2) / det_j
        jac_inv(1,2) = -jac(1,2) / det_j
        jac_inv(2,1) = -jac(2,1) / det_j
        jac_inv(2,2) = jac(1,1) / det_j
        
        ! Transform gradient: grad_phys = J^(-T) * grad_ref
        grad_phys(1) = jac_inv(1,1) * grad_ref(1) + jac_inv(2,1) * grad_ref(2)
        grad_phys(2) = jac_inv(1,2) * grad_ref(1) + jac_inv(2,2) * grad_ref(2)
        
        grad_ref = grad_phys
    end subroutine transform_gradient
    
    ! Test problem: -Δu = f with u = sin(πx)sin(πy) on ∂Ω
    pure function exact_solution(x, y) result(u)
        real(dp), intent(in) :: x, y
        real(dp) :: u
        real(dp), parameter :: pi_const = 3.14159265358979323846_dp
        u = sin(pi_const * x) * sin(pi_const * y)
    end function exact_solution
    
    pure function exact_gradient(x, y) result(grad)
        real(dp), intent(in) :: x, y
        real(dp) :: grad(2)
        real(dp), parameter :: pi_const = 3.14159265358979323846_dp
        grad(1) = pi_const * cos(pi_const * x) * sin(pi_const * y)
        grad(2) = pi_const * sin(pi_const * x) * cos(pi_const * y)
    end function exact_gradient
    
    pure function source_function(x, y) result(f)
        real(dp), intent(in) :: x, y
        real(dp) :: f
        real(dp), parameter :: pi_const = 3.14159265358979323846_dp
        f = 2.0_dp * pi_const**2 * sin(pi_const * x) * sin(pi_const * y)
    end function source_function

end program test_poisson_p2_convergence