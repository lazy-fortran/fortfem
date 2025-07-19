program test_poisson_2d_benchmark
    ! Poisson 2D benchmark with analytical solution
    ! Problem: -Δu = f in (0,1)×(0,1), u = 0 on ∂Ω
    ! Analytical solution: u(x,y) = sin(πx)sin(πy)
    ! Source term: f(x,y) = 2π²sin(πx)sin(πy)
    use fortfem_kinds, only: dp
    use poisson_2d_module
    use fortfem_mesh_2d
    implicit none
    
    integer, parameter :: n_levels = 5
    integer :: level, nx, ny, n_vertices
    real(dp) :: h, L2_error, H1_error
    real(dp), parameter :: pi = 4.0_dp * atan(1.0_dp)
    type(poisson_2d_t) :: solver
    type(mesh_2d_t) :: mesh
    real(dp), allocatable :: u_numerical(:), u_exact(:)
    integer, allocatable :: boundary_nodes(:)
    real(dp), allocatable :: boundary_values(:)
    
    print *, "FortFEM Poisson 2D Benchmark"
    print *, "============================"
    print *, "Problem: -Δu = 2π²sin(πx)sin(πy), u = 0 on ∂Ω"
    print *, "Analytical: u(x,y) = sin(πx)sin(πy)"
    print *, ""
    print *, "h        L2_error        H1_error        DOFs"
    print *, "--------------------------------------------"
    
    do level = 1, n_levels
        ! Create mesh: start with 4×4, double each level
        nx = 4 * 2**(level-1)
        ny = nx
        n_vertices = (nx + 1) * (ny + 1)
        h = 1.0_dp / real(nx, dp)
        
        ! Create rectangular mesh
        call mesh%create_rectangular(nx + 1, ny + 1, 0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp)
        
        ! Initialize solver
        call solver%init("lapack")
        call solver%set_mesh(mesh)
        
        ! Set homogeneous Dirichlet boundary conditions u = 0 on all boundaries
        call get_boundary_nodes(mesh, boundary_nodes)
        allocate(boundary_values(size(boundary_nodes)))
        boundary_values = 0.0_dp
        call solver%set_dirichlet_bc(boundary_nodes, boundary_values)
        
        ! Solve system
        call solver%solve(poisson_source_2d)
        u_numerical = solver%get_solution()
        
        ! Compute exact solution
        call compute_exact_solution_2d(mesh, u_exact)
        
        ! Compute errors
        call compute_errors_2d(mesh, u_numerical, u_exact, L2_error, H1_error)
        
        ! Print results with solution statistics
        write(*, '(F8.6, 4X, E12.6, 4X, E12.6, 4X, I0)') h, L2_error, H1_error, n_vertices
        if (level == 1) then
            print *, "Debug: max |u_num| =", maxval(abs(u_numerical))
            print *, "Debug: max |u_exact| =", maxval(abs(u_exact))
            print *, "Debug: max |u_exact - u_num| =", maxval(abs(u_exact - u_numerical))
            print *, "Debug: interior solution (5x5 node) =", u_numerical(13) ! should be ~sin(π*0.5)*sin(π*0.5) ≈ 1.0
            print *, "Debug: expected at center =", sin(pi * 0.5_dp) * sin(pi * 0.5_dp)
            print *, "Debug: scaling factor =", u_numerical(13) / (sin(pi * 0.5_dp) * sin(pi * 0.5_dp))
        end if
        
        ! Clean up
        deallocate(u_exact, boundary_nodes, boundary_values)
        call solver%destroy()
        call mesh%destroy()
    end do
    
    print *, ""
    print *, "Expected convergence: L2 O(h²), H1 O(h)"

contains

    ! Source function: f(x,y) = 2π²sin(πx)sin(πy)
    pure function poisson_source_2d(x, y) result(f)
        real(dp), intent(in) :: x, y
        real(dp) :: f
        f = 2.0_dp * pi**2 * sin(pi * x) * sin(pi * y)
        ! Debug: print first few calls
        !if (abs(x - 0.5_dp) < 0.01_dp .and. abs(y - 0.5_dp) < 0.01_dp) then
        !    print *, "Source at center: f(0.5,0.5) =", f, "should be ≈", 2.0_dp * pi**2
        !end if
    end function poisson_source_2d
    
    ! Get boundary nodes for rectangular domain
    subroutine get_boundary_nodes(mesh, boundary_nodes)
        type(mesh_2d_t), intent(in) :: mesh
        integer, allocatable, intent(out) :: boundary_nodes(:)
        
        integer :: i, n_boundary, count
        real(dp) :: x, y
        real(dp), parameter :: tol = 1e-12_dp
        
        ! Count boundary nodes
        n_boundary = 0
        do i = 1, mesh%n_vertices
            x = mesh%vertices(1, i)
            y = mesh%vertices(2, i)
            if (abs(x) < tol .or. abs(x - 1.0_dp) < tol .or. &
                abs(y) < tol .or. abs(y - 1.0_dp) < tol) then
                n_boundary = n_boundary + 1
            end if
        end do
        
        ! Store boundary node indices
        allocate(boundary_nodes(n_boundary))
        count = 0
        do i = 1, mesh%n_vertices
            x = mesh%vertices(1, i)
            y = mesh%vertices(2, i)
            if (abs(x) < tol .or. abs(x - 1.0_dp) < tol .or. &
                abs(y) < tol .or. abs(y - 1.0_dp) < tol) then
                count = count + 1
                boundary_nodes(count) = i
            end if
        end do
    end subroutine get_boundary_nodes
    
    ! Compute exact solution at mesh vertices
    subroutine compute_exact_solution_2d(mesh, u_exact)
        type(mesh_2d_t), intent(in) :: mesh
        real(dp), allocatable, intent(out) :: u_exact(:)
        integer :: i
        real(dp) :: x, y
        
        allocate(u_exact(mesh%n_vertices))
        do i = 1, mesh%n_vertices
            x = mesh%vertices(1, i)
            y = mesh%vertices(2, i)
            u_exact(i) = sin(pi * x) * sin(pi * y)
        end do
    end subroutine compute_exact_solution_2d
    
    ! Compute L2 and H1 errors using Gaussian quadrature
    subroutine compute_errors_2d(mesh, u_numerical, u_exact, L2_error, H1_error)
        type(mesh_2d_t), intent(in) :: mesh
        real(dp), intent(in) :: u_numerical(:), u_exact(:)
        real(dp), intent(out) :: L2_error, H1_error
        
        ! 3-point edge midpoint rule for triangles (same as assembly)
        real(dp), parameter :: xi_quad(3) = [0.5_dp, 0.5_dp, 0.0_dp]
        real(dp), parameter :: eta_quad(3) = [0.0_dp, 0.5_dp, 0.5_dp]
        real(dp), parameter :: w_quad(3) = [1.0_dp/6.0_dp, 1.0_dp/6.0_dp, 1.0_dp/6.0_dp]
        
        integer :: t, q, n1, n2, n3
        real(dp) :: xi, eta, w, triangle_area
        real(dp) :: x_phys, y_phys
        real(dp) :: phi1, phi2, phi3
        real(dp) :: dphidx1, dphidx2, dphidx3, dphidy1, dphidy2, dphidy3
        real(dp) :: u_num, u_ex, dudx_num, dudx_ex, dudy_num, dudy_ex
        real(dp) :: error_u, error_du
        real(dp) :: x1, y1, x2, y2, x3, y3, det_J
        
        L2_error = 0.0_dp
        H1_error = 0.0_dp
        
        ! Integrate over all triangles
        do t = 1, mesh%n_triangles
            n1 = mesh%triangles(1, t)
            n2 = mesh%triangles(2, t)
            n3 = mesh%triangles(3, t)
            
            ! Triangle vertices
            x1 = mesh%vertices(1, n1); y1 = mesh%vertices(2, n1)
            x2 = mesh%vertices(1, n2); y2 = mesh%vertices(2, n2)
            x3 = mesh%vertices(1, n3); y3 = mesh%vertices(2, n3)
            
            ! Triangle area
            triangle_area = 0.5_dp * abs((x1-x3)*(y2-y3) - (x2-x3)*(y1-y3))
            det_J = 2.0_dp * triangle_area
            
            ! Loop over quadrature points
            do q = 1, 3
                xi = xi_quad(q)
                eta = eta_quad(q)
                w = w_quad(q)
                
                ! Linear basis functions on reference triangle
                phi1 = 1.0_dp - xi - eta
                phi2 = xi
                phi3 = eta
                
                ! Physical coordinates using consistent basis functions  
                ! x = x1*phi1 + x2*phi2 + x3*phi3 = x1*(1-xi-eta) + x2*xi + x3*eta
                x_phys = x1 * phi1 + x2 * phi2 + x3 * phi3
                y_phys = y1 * phi1 + y2 * phi2 + y3 * phi3
                
                ! Basis function gradients using same approach as assembly module
                ! Jacobian J = [[x2-x1, x3-x1], [y2-y1, y3-y1]]
                ! J^{-1} = 1/det * [[y3-y1, x1-x3], [y1-y2, x2-x1]]
                ! grad_phys = J^{-T} * grad_ref 
                ! For basis functions: grad_phi1 = [-1,-1], grad_phi2 = [1,0], grad_phi3 = [0,1]
                
                ! Compute J^{-T} explicitly
                ! J^{-T} = 1/det * [[y3-y1, y1-y2], [x1-x3, x2-x1]]
                dphidx1 = (-(y3-y1) - (y1-y2)) / det_J   ! -1*(y3-y1)/det + (-1)*(y1-y2)/det
                dphidx2 = (y3-y1) / det_J                 ! 1*(y3-y1)/det + 0*(y1-y2)/det
                dphidx3 = (y1-y2) / det_J                 ! 0*(y3-y1)/det + 1*(y1-y2)/det
                
                dphidy1 = (-(x1-x3) - (x2-x1)) / det_J   ! -1*(x1-x3)/det + (-1)*(x2-x1)/det  
                dphidy2 = (x1-x3) / det_J                 ! 1*(x1-x3)/det + 0*(x2-x1)/det
                dphidy3 = (x2-x1) / det_J                 ! 0*(x1-x3)/det + 1*(x2-x1)/det
                
                ! Numerical solution and gradient
                u_num = u_numerical(n1) * phi1 + u_numerical(n2) * phi2 + u_numerical(n3) * phi3
                dudx_num = u_numerical(n1) * dphidx1 + u_numerical(n2) * dphidx2 + u_numerical(n3) * dphidx3
                dudy_num = u_numerical(n1) * dphidy1 + u_numerical(n2) * dphidy2 + u_numerical(n3) * dphidy3
                
                ! Exact solution and gradient
                u_ex = sin(pi * x_phys) * sin(pi * y_phys)
                dudx_ex = pi * cos(pi * x_phys) * sin(pi * y_phys)
                dudy_ex = pi * sin(pi * x_phys) * cos(pi * y_phys)
                
                ! Errors
                error_u = u_num - u_ex
                error_du = (dudx_num - dudx_ex)**2 + (dudy_num - dudy_ex)**2
                
                ! Integrate (use abs(det_J) for consistency)
                L2_error = L2_error + error_u**2 * w * abs(det_J)
                H1_error = H1_error + (error_u**2 + error_du) * w * abs(det_J)
            end do
        end do
        
        L2_error = sqrt(L2_error)
        H1_error = sqrt(H1_error)
    end subroutine compute_errors_2d

end program test_poisson_2d_benchmark