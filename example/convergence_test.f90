program convergence_test
    ! Test convergence rates for Poisson equation
    use fortfem
    implicit none
    
    type(mesh_2d_t) :: mesh
    type(sparse_matrix_t) :: A
    real(dp), allocatable :: u(:), f(:), exact(:)
    real(dp) :: h, error, error_old
    integer :: n, i
    
    print '(a)', "Convergence test for -∆u = f"
    print '(a)', "n     h        L2 error   Rate"
    print '(a)', "--------------------------------"
    
    error_old = 1.0_dp
    
    do i = 1, 4
        n = 10 * 2**(i-1)
        h = 1.0_dp / real(n, dp)
        
        ! Create mesh and solve
        call create_unit_square_mesh(mesh, n=n)
        call assemble_poisson_2d(mesh, A, f)
        call apply_zero_bc(mesh, A, f)
        call solve_sparse(A, f, u)
        
        ! Compute error (simplified)
        error = compute_l2_error(mesh, u)
        
        ! Print results
        if (i > 1) then
            print '(i3,3f10.5)', n, h, error, log(error_old/error)/log(2.0_dp)
        else
            print '(i3,2f10.5,a)', n, h, error, "    -"
        end if
        
        error_old = error
    end do
    
end program convergence_test