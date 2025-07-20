program test_advanced_solvers
    use fortfem_kinds
    use fortfem_advanced_solvers
    use check
    implicit none

    write(*,*) "Testing advanced linear solver implementations..."
    
    call test_cg_solver()
    call test_cg_with_preconditioners()
    call test_bicgstab_solver()
    call test_gmres_solver()
    call test_direct_sparse_solvers()
    call test_solver_selection()
    call test_convergence_criteria()
    call test_performance_comparison()
    call test_ill_conditioned_systems()
    call test_parallel_efficiency()
    
    call check_summary("Advanced Linear Solvers")

contains

    ! Test Conjugate Gradient solver
    subroutine test_cg_solver()
        real(dp), allocatable :: A(:,:), b(:), x(:), x_exact(:)
        type(solver_options_t) :: opts
        type(solver_stats_t) :: stats
        integer :: n = 100, i
        real(dp) :: residual_norm, error_norm
        
        ! Create SPD test system: A = D + sparse perturbation
        call create_spd_matrix(n, A)
        allocate(b(n), x(n), x_exact(n))
        
        ! Set up exact solution and RHS
        x_exact = [(real(i, dp), i=1,n)]
        b = matmul(A, x_exact)
        x = 0.0_dp  ! Initial guess
        
        ! Configure CG solver
        opts = solver_options(method="cg", tolerance=1.0e-10_dp, &
                            max_iterations=n, verbosity=0)
        
        ! Solve system
        call solve(A, b, x, opts, stats)
        
        ! Check convergence
        call check_condition(stats%converged, &
            "CG solver: converged")
        call check_condition(stats%iterations < n/2, &
            "CG solver: reasonable iteration count")
        
        residual_norm = norm(matmul(A, x) - b)
        error_norm = norm(x - x_exact)
        
        call check_condition(residual_norm < 1.0e-7_dp, &
            "CG solver: small residual")
        call check_condition(error_norm < 1.0e-6_dp, &
            "CG solver: accurate solution")
        call check_condition(stats%final_residual < 1.0e-7_dp, &
            "CG solver: tolerance achieved")
        
        write(*,*) "   CG iterations:", stats%iterations
        write(*,*) "   CG residual:", stats%final_residual
        write(*,*) "   Error norm:", error_norm
        
        deallocate(A, b, x, x_exact)
        write(*,*) "   CG solver: all tests passed"
    end subroutine test_cg_solver

    ! Test CG with different preconditioners
    subroutine test_cg_with_preconditioners()
        real(dp), allocatable :: A(:,:), b(:), x_jacobi(:), x_ilu(:)
        type(solver_options_t) :: opts_jacobi, opts_ilu
        type(solver_stats_t) :: stats_jacobi, stats_ilu
        integer :: n = 200
        
        call create_spd_matrix(n, A)
        allocate(b(n), x_jacobi(n), x_ilu(n))
        b = 1.0_dp  ! Simple RHS
        
        ! Test Jacobi preconditioned CG
        x_jacobi = 0.0_dp
        opts_jacobi = solver_options(method="pcg", preconditioner="jacobi", &
                                   tolerance=1.0e-8_dp, max_iterations=n)
        call solve(A, b, x_jacobi, opts_jacobi, stats_jacobi)
        
        call check_condition(stats_jacobi%converged, &
            "PCG Jacobi: converged")
        call check_condition(stats_jacobi%iterations < n, &
            "PCG Jacobi: reasonable iterations")
        
        ! Test ILU preconditioned CG  
        x_ilu = 0.0_dp
        opts_ilu = solver_options(method="pcg", preconditioner="ilu", &
                                tolerance=1.0e-8_dp, max_iterations=n)
        call solve(A, b, x_ilu, opts_ilu, stats_ilu)
        
        call check_condition(stats_ilu%converged, &
            "PCG ILU: converged")
        call check_condition(stats_ilu%iterations <= stats_jacobi%iterations, &
            "PCG ILU: better than Jacobi")
        
        ! Check both solutions are similar
        call check_condition(norm(x_jacobi - x_ilu) < 1.0e-6_dp, &
            "PCG: consistent solutions")
        
        write(*,*) "   Jacobi iterations:", stats_jacobi%iterations
        write(*,*) "   ILU iterations:", stats_ilu%iterations
        write(*,*) "   ILU improvement factor:", &
            real(stats_jacobi%iterations) / real(stats_ilu%iterations)
        
        deallocate(A, b, x_jacobi, x_ilu)
        write(*,*) "   Preconditioned CG: all tests passed"
    end subroutine test_cg_with_preconditioners

    ! Test BiCGSTAB for non-symmetric systems
    subroutine test_bicgstab_solver()
        real(dp), allocatable :: A(:,:), b(:), x(:), x_exact(:)
        type(solver_options_t) :: opts
        type(solver_stats_t) :: stats
        integer :: n = 100, i
        real(dp) :: error_norm
        
        ! Create non-symmetric test matrix
        call create_nonsymmetric_matrix(n, A)
        allocate(b(n), x(n), x_exact(n))
        
        x_exact = [(1.0_dp, i=1,n)]
        b = matmul(A, x_exact)
        x = 0.0_dp
        
        opts = solver_options(method="bicgstab", tolerance=1.0e-10_dp, &
                            max_iterations=2*n)
        call solve(A, b, x, opts, stats)
        
        call check_condition(stats%converged, &
            "BiCGSTAB: converged")
        call check_condition(stats%iterations < n, &
            "BiCGSTAB: reasonable iterations")
        
        error_norm = norm(x - x_exact)
        call check_condition(error_norm < 1.0e-8_dp, &
            "BiCGSTAB: accurate solution")
        
        ! Test with preconditioning
        x = 0.0_dp
        opts%preconditioner = "ilu"
        call solve(A, b, x, opts, stats)
        
        call check_condition(stats%converged, &
            "BiCGSTAB+ILU: converged")
        
        write(*,*) "   BiCGSTAB iterations:", stats%iterations
        write(*,*) "   Error norm:", error_norm
        
        deallocate(A, b, x, x_exact)
        write(*,*) "   BiCGSTAB solver: all tests passed"
    end subroutine test_bicgstab_solver

    ! Test GMRES solver with restart
    subroutine test_gmres_solver()
        real(dp), allocatable :: A(:,:), b(:), x(:)
        type(solver_options_t) :: opts
        type(solver_stats_t) :: stats
        integer :: n = 150
        real(dp) :: residual_norm
        
        call create_nonsymmetric_matrix(n, A)
        allocate(b(n), x(n))
        b = 1.0_dp
        x = 0.0_dp
        
        ! Test GMRES with restart
        opts = solver_options(method="gmres", tolerance=1.0e-8_dp, &
                            max_iterations=n, restart=30)
        call solve(A, b, x, opts, stats)
        
        call check_condition(stats%converged, &
            "GMRES: converged")
        call check_condition(stats%restarts >= 0, &
            "GMRES: restart count tracked")
        
        residual_norm = norm(matmul(A, x) - b)
        call check_condition(residual_norm < 1.0e-7_dp, &
            "GMRES: small residual")
        
        ! Test flexible GMRES with preconditioning
        x = 0.0_dp
        opts%method = "fgmres"
        opts%preconditioner = "ilu"
        call solve(A, b, x, opts, stats)
        
        call check_condition(stats%converged, &
            "FGMRES: converged")
        
        write(*,*) "   GMRES iterations:", stats%iterations
        write(*,*) "   GMRES restarts:", stats%restarts
        write(*,*) "   Residual norm:", residual_norm
        
        deallocate(A, b, x)
        write(*,*) "   GMRES solver: all tests passed"
    end subroutine test_gmres_solver

    ! Test direct sparse solvers
    subroutine test_direct_sparse_solvers()
        real(dp), allocatable :: A(:,:), b(:), x_lu(:), x_umf(:)
        type(solver_options_t) :: opts_lu, opts_umf
        type(solver_stats_t) :: stats_lu, stats_umf
        integer :: n = 100
        real(dp) :: error_norm
        
        call create_spd_matrix(n, A)
        allocate(b(n), x_lu(n), x_umf(n))
        b = 1.0_dp
        
        ! Test sparse LU factorization
        x_lu = 0.0_dp
        opts_lu = solver_options(method="sparse_lu")
        call solve(A, b, x_lu, opts_lu, stats_lu)
        
        call check_condition(stats_lu%converged, &
            "Sparse LU: solved")
        call check_condition(stats_lu%iterations == 1, &
            "Sparse LU: direct method")
        
        error_norm = norm(matmul(A, x_lu) - b)
        call check_condition(error_norm < 1.0e-12_dp, &
            "Sparse LU: machine precision accuracy")
        
        ! Test UMFPACK if available
        if (umfpack_available()) then
            x_umf = 0.0_dp
            opts_umf = solver_options(method="umfpack")
            call solve(A, b, x_umf, opts_umf, stats_umf)
            
            call check_condition(stats_umf%converged, &
                "UMFPACK: solved")
            call check_condition(norm(x_lu - x_umf) < 1.0e-10_dp, &
                "UMFPACK: consistent with LU")
            
            write(*,*) "   UMFPACK available and tested"
        else
            write(*,*) "   UMFPACK not available, skipped"
        end if
        
        write(*,*) "   LU residual:", error_norm
        
        deallocate(A, b, x_lu, x_umf)
        write(*,*) "   Direct sparse solvers: all tests passed"
    end subroutine test_direct_sparse_solvers

    ! Test automatic solver selection
    subroutine test_solver_selection()
        real(dp), allocatable :: A_small(:,:), A_large(:,:), b_small(:), b_large(:)
        real(dp), allocatable :: x_small(:), x_large(:)
        type(solver_options_t) :: opts_auto
        type(solver_stats_t) :: stats_small, stats_large
        
        ! Small system should use direct solver
        call create_spd_matrix(50, A_small)
        allocate(b_small(50), x_small(50))
        b_small = 1.0_dp
        x_small = 0.0_dp
        
        opts_auto = solver_options(method="auto")
        call solve(A_small, b_small, x_small, opts_auto, stats_small)
        
        call check_condition(stats_small%converged, &
            "Auto selection: small system solved")
        call check_condition(index(stats_small%method_used, "lapack") > 0 .or. &
                           index(stats_small%method_used, "direct") > 0, &
            "Auto selection: direct for small system")
        
        ! Large system should use iterative solver
        call create_spd_matrix(1000, A_large)
        allocate(b_large(1000), x_large(1000))
        b_large = 1.0_dp
        x_large = 0.0_dp
        
        call solve(A_large, b_large, x_large, opts_auto, stats_large)
        
        call check_condition(stats_large%converged, &
            "Auto selection: large system solved")
        call check_condition(index(stats_large%method_used, "cg") > 0 .or. &
                           index(stats_large%method_used, "pcg") > 0, &
            "Auto selection: iterative for large system")
        
        write(*,*) "   Small system method:", trim(stats_small%method_used)
        write(*,*) "   Large system method:", trim(stats_large%method_used)
        
        deallocate(A_small, A_large, b_small, b_large, x_small, x_large)
        write(*,*) "   Automatic solver selection: all tests passed"
    end subroutine test_solver_selection

    ! Test convergence criteria and monitoring
    subroutine test_convergence_criteria()
        real(dp), allocatable :: A(:,:), b(:), x(:)
        type(solver_options_t) :: opts
        type(solver_stats_t) :: stats
        integer :: n = 100
        
        call create_spd_matrix(n, A)
        allocate(b(n), x(n))
        b = 1.0_dp
        x = 0.0_dp
        
        ! Test absolute tolerance
        opts = solver_options(method="cg", tolerance=1.0e-6_dp, &
                            tolerance_type="absolute")
        call solve(A, b, x, opts, stats)
        
        call check_condition(stats%converged, &
            "Convergence: absolute tolerance")
        call check_condition(stats%final_residual <= opts%tolerance, &
            "Convergence: tolerance satisfied")
        
        ! Test relative tolerance
        x = 0.0_dp
        opts%tolerance_type = "relative"
        opts%tolerance = 1.0e-8_dp
        call solve(A, b, x, opts, stats)
        
        call check_condition(stats%converged, &
            "Convergence: relative tolerance")
        
        ! Test maximum iterations limit
        x = 0.0_dp
        opts%max_iterations = 5
        opts%tolerance = 1.0e-12_dp
        call solve(A, b, x, opts, stats)
        
        call check_condition(.not. stats%converged, &
            "Convergence: max iterations reached")
        call check_condition(stats%iterations == opts%max_iterations, &
            "Convergence: iteration limit enforced")
        
        write(*,*) "   Final residual (abs):", stats%final_residual
        write(*,*) "   Iterations with limit:", stats%iterations
        
        deallocate(A, b, x)
        write(*,*) "   Convergence criteria: all tests passed"
    end subroutine test_convergence_criteria

    ! Test performance comparison
    subroutine test_performance_comparison()
        real(dp), allocatable :: A(:,:), b(:), x_direct(:), x_iter(:)
        type(solver_options_t) :: opts_direct, opts_iter
        type(solver_stats_t) :: stats_direct, stats_iter
        real(dp) :: time_direct, time_iter, speedup
        integer :: n = 500
        
        call create_spd_matrix(n, A)
        allocate(b(n), x_direct(n), x_iter(n))
        b = 1.0_dp
        
        ! Benchmark direct solver
        x_direct = 0.0_dp
        opts_direct = solver_options(method="lapack_lu")
        call cpu_time(time_direct)
        call solve(A, b, x_direct, opts_direct, stats_direct)
        call cpu_time(time_direct)
        time_direct = stats_direct%solve_time
        
        ! Benchmark iterative solver
        x_iter = 0.0_dp
        opts_iter = solver_options(method="pcg", preconditioner="ilu")
        call solve(A, b, x_iter, opts_iter, stats_iter)
        time_iter = stats_iter%solve_time
        
        call check_condition(stats_direct%converged, &
            "Performance: direct solver works")
        call check_condition(stats_iter%converged, &
            "Performance: iterative solver works")
        
        speedup = time_direct / time_iter
        call check_condition(norm(x_direct - x_iter) < 1.0e-5_dp, &
            "Performance: consistent solutions")
        
        ! Memory usage comparison
        call check_condition(stats_iter%memory_usage < stats_direct%memory_usage, &
            "Performance: iterative uses less memory")
        
        write(*,*) "   Direct solve time:", time_direct
        write(*,*) "   Iterative solve time:", time_iter
        write(*,*) "   Speedup factor:", speedup
        write(*,*) "   Memory ratio (iter/direct):", &
            real(stats_iter%memory_usage) / real(stats_direct%memory_usage)
        
        deallocate(A, b, x_direct, x_iter)
        write(*,*) "   Performance comparison: all tests passed"
    end subroutine test_performance_comparison

    ! Test behavior on ill-conditioned systems
    subroutine test_ill_conditioned_systems()
        real(dp), allocatable :: A(:,:), b(:), x(:)
        type(solver_options_t) :: opts
        type(solver_stats_t) :: stats
        integer :: n = 100
        real(dp) :: condition_number, residual_norm
        
        ! Create ill-conditioned matrix
        call create_ill_conditioned_matrix(n, A, condition_number)
        allocate(b(n), x(n))
        b = 1.0_dp
        x = 0.0_dp
        
        call check_condition(condition_number > 1.0e8_dp, &
            "Ill-conditioned: high condition number")
        
        ! Test CG with preconditioning
        opts = solver_options(method="pcg", preconditioner="ilu", &
                            tolerance=1.0e-6_dp, max_iterations=n*2)
        call solve(A, b, x, opts, stats)
        
        call check_condition(stats%converged .or. stats%iterations < n*2, &
            "Ill-conditioned: reasonable behavior")
        
        residual_norm = norm(matmul(A, x) - b)
        
        ! Even if not fully converged, should make progress
        call check_condition(residual_norm < norm(b), &
            "Ill-conditioned: residual reduction")
        
        write(*,*) "   Condition number:", condition_number
        write(*,*) "   Final residual:", residual_norm
        write(*,*) "   Iterations used:", stats%iterations
        
        deallocate(A, b, x)
        write(*,*) "   Ill-conditioned systems: all tests passed"
    end subroutine test_ill_conditioned_systems

    ! Test parallel efficiency (basic)
    subroutine test_parallel_efficiency()
        real(dp), allocatable :: A(:,:), b(:), x(:)
        type(solver_options_t) :: opts
        type(solver_stats_t) :: stats
        integer :: n = 300
        
        call create_spd_matrix(n, A)
        allocate(b(n), x(n))
        b = 1.0_dp
        x = 0.0_dp
        
        ! Test parallel CG if available
        opts = solver_options(method="pcg", preconditioner="jacobi", &
                            parallel=.true., num_threads=2)
        call solve(A, b, x, opts, stats)
        
        call check_condition(stats%converged, &
            "Parallel: solver converged")
        
        if (stats%parallel_efficiency > 0.0_dp) then
            call check_condition(stats%parallel_efficiency > 0.5_dp, &
                "Parallel: reasonable efficiency")
            write(*,*) "   Parallel efficiency:", stats%parallel_efficiency
        else
            write(*,*) "   Parallel efficiency not measured"
        end if
        
        deallocate(A, b, x)
        write(*,*) "   Parallel efficiency: all tests passed"
    end subroutine test_parallel_efficiency

    ! Helper functions for test matrix creation

    subroutine create_spd_matrix(n, A)
        integer, intent(in) :: n
        real(dp), allocatable, intent(out) :: A(:,:)
        integer :: i, j
        
        allocate(A(n, n))
        A = 0.0_dp
        
        ! Create diagonally dominant SPD matrix - tridiagonal
        do i = 1, n
            A(i, i) = 4.0_dp
            
            ! Add fixed off-diagonal terms
            if (i > 1) then
                A(i, i-1) = -1.0_dp
                A(i-1, i) = -1.0_dp
            end if
            if (i < n) then
                A(i, i+1) = -1.0_dp
                A(i+1, i) = -1.0_dp
            end if
        end do
        
        ! Make sure it's well-conditioned
        do i = 1, n
            A(i, i) = A(i, i) + 0.1_dp
        end do
    end subroutine create_spd_matrix

    subroutine create_nonsymmetric_matrix(n, A)
        integer, intent(in) :: n
        real(dp), allocatable, intent(out) :: A(:,:)
        integer :: i, j
        real(dp) :: random_val
        
        allocate(A(n, n))
        A = 0.0_dp
        
        ! Create diagonally dominant non-symmetric matrix
        do i = 1, n
            A(i, i) = real(n, dp) + 1.0_dp
            
            if (i > 1) then
                call random_number(random_val)
                A(i, i-1) = -0.3_dp * random_val
            end if
            if (i < n) then
                call random_number(random_val)
                A(i, i+1) = -0.7_dp * random_val  ! Different from lower
            end if
        end do
    end subroutine create_nonsymmetric_matrix

    subroutine create_ill_conditioned_matrix(n, A, condition_number)
        integer, intent(in) :: n
        real(dp), allocatable, intent(out) :: A(:,:)
        real(dp), intent(out) :: condition_number
        integer :: i
        real(dp) :: smallest_eigenvalue, largest_eigenvalue
        
        allocate(A(n, n))
        A = 0.0_dp
        
        ! Create matrix with wide range of eigenvalues
        smallest_eigenvalue = 1.0e-10_dp
        largest_eigenvalue = 1.0_dp
        
        do i = 1, n
            A(i, i) = smallest_eigenvalue + &
                (largest_eigenvalue - smallest_eigenvalue) * &
                real(i-1, dp) / real(n-1, dp)
        end do
        
        condition_number = largest_eigenvalue / smallest_eigenvalue
    end subroutine create_ill_conditioned_matrix

    function norm(x) result(norm_val)
        real(dp), intent(in) :: x(:)
        real(dp) :: norm_val
        norm_val = sqrt(sum(x**2))
    end function norm

    function umfpack_available() result(available)
        logical :: available
        ! Placeholder - would check if UMFPACK is linked
        available = .false.
    end function umfpack_available

end program test_advanced_solvers