program test_poisson_1d_benchmark
    ! Poisson 1D benchmark with analytical solution
    ! Problem: -u'' = f in (0,1), u(0) = u(1) = 0
    ! Analytical solution: u(x) = sin(π*x)
    ! Source term: f(x) = π²*sin(π*x)
    use fortfem_kinds, only: dp
    use fortfem_poisson_1d
    use fortfem_mesh_1d
    implicit none
    
    ! LAPACK interface
    interface
        subroutine dgesv(n, nrhs, a, lda, ipiv, b, ldb, info)
            integer, intent(in) :: n, nrhs, lda, ldb
            double precision, intent(inout) :: a(lda,*), b(ldb,*)
            integer, intent(out) :: ipiv(*), info
        end subroutine dgesv
    end interface
    
    integer, parameter :: n_levels = 6
    integer :: level, n_nodes, n
    real(dp) :: h, L2_error, H1_error
    real(dp), parameter :: pi = 4.0_dp * atan(1.0_dp)
    type(poisson_1d_solver_t) :: solver
    real(dp), allocatable :: u_numerical(:), u_exact(:)
    
    print *, "FortFEM Poisson 1D Benchmark"
    print *, "============================"
    print *, "Problem: -u'' = π²sin(πx), u(0) = u(1) = 0"
    print *, "Analytical: u(x) = sin(πx)"
    print *, ""
    print *, "h        L2_error        H1_error        DOFs"
    print *, "--------------------------------------------"
    
    do level = 1, n_levels
        ! Create mesh: start with 8 elements, double each level
        n = 8 * 2**(level-1)
        n_nodes = n + 1
        h = 1.0_dp / real(n, dp)
        
        ! Initialize solver
        call solver%init(n_nodes, 0.0_dp, 1.0_dp)
        
        ! Solve system
        call solver%solve(poisson_source, u_numerical)
        
        ! Compute exact solution
        allocate(u_exact(n_nodes))
        call compute_exact_solution(solver%mesh, u_exact)
        
        ! Compute errors
        call compute_errors(solver%mesh, u_numerical, u_exact, L2_error, H1_error)
        
        ! Print results
        write(*, '(F8.6, 4X, E12.6, 4X, E12.6, 4X, I0)') h, L2_error, H1_error, n_nodes
        
        ! Clean up
        deallocate(u_exact)
        call solver%deallocate()
    end do
    
    print *, ""
    print *, "Expected convergence: L2 O(h²), H1 O(h)"

contains

    ! Source function: f(x) = π²*sin(π*x)
    function poisson_source(x) result(f)
        real(dp), intent(in) :: x
        real(dp) :: f
        f = pi**2 * sin(pi * x)
    end function poisson_source
    
    ! Compute exact solution at mesh nodes
    subroutine compute_exact_solution(mesh, u_exact)
        type(mesh_1d_t), intent(in) :: mesh
        real(dp), intent(out) :: u_exact(:)
        integer :: i
        
        do i = 1, mesh%n_nodes
            u_exact(i) = sin(pi * mesh%nodes(i))
        end do
    end subroutine compute_exact_solution
    
    ! Compute L2 and H1 errors using Gaussian quadrature
    subroutine compute_errors(mesh, u_numerical, u_exact, L2_error, H1_error)
        type(mesh_1d_t), intent(in) :: mesh
        real(dp), intent(in) :: u_numerical(:), u_exact(:)
        real(dp), intent(out) :: L2_error, H1_error
        
        ! 2-point Gaussian quadrature
        real(dp), parameter :: xi_quad(2) = [-0.5773502692_dp, 0.5773502692_dp]
        real(dp), parameter :: w_quad(2) = [1.0_dp, 1.0_dp]
        
        integer :: e, q, n1, n2
        real(dp) :: h_elem, xi, w, x_phys
        real(dp) :: phi1, phi2, dphi1, dphi2
        real(dp) :: u_num, u_ex, du_num, du_ex
        real(dp) :: error_u, error_du
        
        L2_error = 0.0_dp
        H1_error = 0.0_dp
        
        ! Integrate over all elements
        do e = 1, mesh%n_elements
            n1 = mesh%connectivity(1, e)
            n2 = mesh%connectivity(2, e)
            h_elem = mesh%element_size(e)
            
            ! Loop over quadrature points
            do q = 1, 2
                xi = xi_quad(q)
                w = w_quad(q)
                
                ! Physical coordinate
                x_phys = mesh%nodes(n1) + 0.5_dp * (1.0_dp + xi) * h_elem
                
                ! Linear basis functions on reference element [-1,1]
                phi1 = 0.5_dp * (1.0_dp - xi)
                phi2 = 0.5_dp * (1.0_dp + xi)
                
                ! Derivatives (chain rule: d/dx = d/dxi * dxi/dx)
                dphi1 = -0.5_dp * (2.0_dp / h_elem)  ! -1 * (2/h)
                dphi2 = 0.5_dp * (2.0_dp / h_elem)   ! +1 * (2/h)
                
                ! Numerical solution and derivative
                u_num = u_numerical(n1) * phi1 + u_numerical(n2) * phi2
                du_num = u_numerical(n1) * dphi1 + u_numerical(n2) * dphi2
                
                ! Exact solution and derivative
                u_ex = sin(pi * x_phys)
                du_ex = pi * cos(pi * x_phys)
                
                ! Errors
                error_u = u_num - u_ex
                error_du = du_num - du_ex
                
                ! Integrate (Jacobian = h/2 for mapping from [-1,1] to element)
                L2_error = L2_error + error_u**2 * w * (h_elem / 2.0_dp)
                H1_error = H1_error + (error_u**2 + error_du**2) * w * (h_elem / 2.0_dp)
            end do
        end do
        
        L2_error = sqrt(L2_error)
        H1_error = sqrt(H1_error)
    end subroutine compute_errors

end program test_poisson_1d_benchmark