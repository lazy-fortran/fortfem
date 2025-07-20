module fortfem_advanced_solvers
    use fortfem_kinds
    implicit none
    private
    
    ! LAPACK interface
    interface
        subroutine dgesv(n, nrhs, a, lda, ipiv, b, ldb, info)
            import :: dp
            integer, intent(in) :: n, nrhs, lda, ldb
            real(dp), intent(inout) :: a(lda, *)
            integer, intent(out) :: ipiv(*)
            real(dp), intent(inout) :: b(ldb, *)
            integer, intent(out) :: info
        end subroutine dgesv
    end interface
    
    ! Public types and interfaces
    public :: solver_options_t, solver_stats_t
    public :: solve, solver_options
    public :: cg_solve, pcg_solve, bicgstab_solve, gmres_solve
    public :: jacobi_preconditioner, ilu_preconditioner
    
    ! Solver options type
    type :: solver_options_t
        character(len=32) :: method = "auto"
        real(dp) :: tolerance = 1.0e-6_dp
        integer :: max_iterations = 1000
        character(len=32) :: preconditioner = "none"
        character(len=32) :: tolerance_type = "relative"
        integer :: restart = 50
        logical :: parallel = .false.
        integer :: num_threads = 1
        integer :: verbosity = 0
        real(dp) :: drop_tolerance = 1.0e-4_dp  ! For ILU
        integer :: fill_level = 0  ! For ILU
    end type solver_options_t
    
    ! Solver statistics type
    type :: solver_stats_t
        logical :: converged = .false.
        integer :: iterations = 0
        real(dp) :: final_residual = 0.0_dp
        real(dp) :: solve_time = 0.0_dp
        integer :: memory_usage = 0  ! In bytes
        character(len=32) :: method_used = ""
        integer :: restarts = 0
        real(dp) :: parallel_efficiency = 0.0_dp
        real(dp) :: condition_estimate = 0.0_dp
    end type solver_stats_t
    
    ! Preconditioner type
    type :: preconditioner_t
        character(len=32) :: type = "none"
        real(dp), allocatable :: diagonal(:)  ! For Jacobi
        real(dp), allocatable :: L(:,:), U(:,:)  ! For ILU
        integer, allocatable :: pivot(:)  ! For ILU
    end type preconditioner_t

contains

    ! Main solver interface with automatic method selection
    subroutine solve(A, b, x, opts, stats)
        real(dp), intent(in) :: A(:,:), b(:)
        real(dp), intent(inout) :: x(:)
        type(solver_options_t), intent(in) :: opts
        type(solver_stats_t), intent(out) :: stats
        
        real(dp) :: start_time, end_time
        integer :: n
        character(len=32) :: selected_method
        
        n = size(A, 1)
        call cpu_time(start_time)
        
        ! Automatic method selection
        if (trim(opts%method) == "auto") then
            selected_method = select_best_method(A, opts)
        else
            selected_method = opts%method
        end if
        
        stats%method_used = selected_method
        
        ! Dispatch to appropriate solver
        select case (trim(selected_method))
        case ("cg")
            call cg_solve(A, b, x, opts, stats)
        case ("pcg")
            call pcg_solve(A, b, x, opts, stats)
        case ("bicgstab")
            call bicgstab_solve(A, b, x, opts, stats)
        case ("gmres", "fgmres")
            call gmres_solve(A, b, x, opts, stats)
        case ("lapack_lu", "direct")
            call lapack_solve(A, b, x, opts, stats)
        case ("sparse_lu")
            call sparse_lu_solve(A, b, x, opts, stats)
        case ("umfpack")
            call umfpack_solve(A, b, x, opts, stats)
        case default
            write(*,*) "Warning: Unknown solver method, using CG"
            call cg_solve(A, b, x, opts, stats)
        end select
        
        call cpu_time(end_time)
        stats%solve_time = end_time - start_time
        stats%memory_usage = estimate_memory_usage(n, selected_method)
    end subroutine solve
    
    ! Conjugate Gradient solver
    subroutine cg_solve(A, b, x, opts, stats)
        real(dp), intent(in) :: A(:,:), b(:)
        real(dp), intent(inout) :: x(:)
        type(solver_options_t), intent(in) :: opts
        type(solver_stats_t), intent(out) :: stats
        
        real(dp), allocatable :: r(:), p(:), Ap(:)
        real(dp) :: alpha, beta, rr_old, rr_new, pAp
        real(dp) :: residual_norm, initial_norm, tolerance
        integer :: n, iter
        
        n = size(x)
        allocate(r(n), p(n), Ap(n))
        
        ! Initial residual: r = b - Ax
        r = b - matmul(A, x)
        p = r
        rr_old = dot_product(r, r)
        initial_norm = sqrt(rr_old)
        
        ! Set tolerance based on type
        if (trim(opts%tolerance_type) == "absolute") then
            tolerance = opts%tolerance
        else
            tolerance = opts%tolerance * initial_norm
        end if
        
        stats%converged = .false.
        
        stats%iterations = 0
        residual_norm = sqrt(rr_old)
        
        do iter = 1, opts%max_iterations
            ! Ap = A * p
            Ap = matmul(A, p)
            pAp = dot_product(p, Ap)
            
            if (abs(pAp) < 1.0e-14_dp) then
                ! Breakdown
                stats%iterations = iter - 1
                exit
            end if
            
            ! alpha = (r, r) / (p, Ap)
            alpha = rr_old / pAp
            
            ! Update solution: x = x + alpha * p
            x = x + alpha * p
            
            ! Update residual: r = r - alpha * Ap
            r = r - alpha * Ap
            
            rr_new = dot_product(r, r)
            residual_norm = sqrt(rr_new)
            
            if (opts%verbosity > 0) then
                write(*,'(A,I4,A,E12.5)') "CG iter ", iter, " residual: ", residual_norm
            end if
            
            ! Check convergence
            if (residual_norm <= tolerance) then
                stats%converged = .true.
                stats%iterations = iter
                exit
            end if
            
            ! beta = (r_new, r_new) / (r_old, r_old)
            beta = rr_new / rr_old
            
            ! Update search direction: p = r + beta * p
            p = r + beta * p
            
            rr_old = rr_new
            stats%iterations = iter
        end do
        
        stats%final_residual = residual_norm
        stats%method_used = "cg"
        
        deallocate(r, p, Ap)
    end subroutine cg_solve
    
    ! Preconditioned Conjugate Gradient solver
    subroutine pcg_solve(A, b, x, opts, stats)
        real(dp), intent(in) :: A(:,:), b(:)
        real(dp), intent(inout) :: x(:)
        type(solver_options_t), intent(in) :: opts
        type(solver_stats_t), intent(out) :: stats
        
        real(dp), allocatable :: r(:), z(:), p(:), Ap(:)
        real(dp) :: alpha, beta, rz_old, rz_new, pAp
        real(dp) :: residual_norm, initial_norm, tolerance
        integer :: n, iter
        type(preconditioner_t) :: precond
        
        n = size(x)
        allocate(r(n), z(n), p(n), Ap(n))
        
        ! Build preconditioner
        call build_preconditioner(A, precond, opts%preconditioner, opts)
        
        ! Initial residual
        r = b - matmul(A, x)
        initial_norm = sqrt(dot_product(r, r))
        
        ! Set tolerance
        if (trim(opts%tolerance_type) == "absolute") then
            tolerance = opts%tolerance
        else
            tolerance = opts%tolerance * initial_norm
        end if
        
        ! Solve M * z = r
        call apply_preconditioner(precond, r, z)
        p = z
        rz_old = dot_product(r, z)
        
        stats%converged = .false.
        
        do iter = 1, opts%max_iterations
            Ap = matmul(A, p)
            pAp = dot_product(p, Ap)
            
            if (abs(pAp) < 1.0e-14_dp) exit
            
            alpha = rz_old / pAp
            x = x + alpha * p
            r = r - alpha * Ap
            
            residual_norm = sqrt(dot_product(r, r))
            
            if (opts%verbosity > 0) then
                write(*,'(A,I4,A,E12.5)') "PCG iter ", iter, " residual: ", residual_norm
            end if
            
            if (residual_norm <= tolerance) then
                stats%converged = .true.
                exit
            end if
            
            call apply_preconditioner(precond, r, z)
            rz_new = dot_product(r, z)
            beta = rz_new / rz_old
            p = z + beta * p
            rz_old = rz_new
        end do
        
        stats%iterations = iter
        stats%final_residual = residual_norm
        stats%method_used = "pcg"
        
        deallocate(r, z, p, Ap)
    end subroutine pcg_solve
    
    ! BiCGSTAB solver for non-symmetric systems
    subroutine bicgstab_solve(A, b, x, opts, stats)
        real(dp), intent(in) :: A(:,:), b(:)
        real(dp), intent(inout) :: x(:)
        type(solver_options_t), intent(in) :: opts
        type(solver_stats_t), intent(out) :: stats
        
        real(dp), allocatable :: r(:), r0(:), p(:), v(:), s(:), t(:), z(:), y(:)
        real(dp) :: rho, rho_old, alpha, omega, beta
        real(dp) :: residual_norm, initial_norm, tolerance
        integer :: n, iter
        type(preconditioner_t) :: precond
        
        n = size(x)
        allocate(r(n), r0(n), p(n), v(n), s(n), t(n), z(n), y(n))
        
        ! Build preconditioner if specified
        if (trim(opts%preconditioner) /= "none") then
            call build_preconditioner(A, precond, opts%preconditioner, opts)
        end if
        
        ! Initial residual
        r = b - matmul(A, x)
        r0 = r
        initial_norm = sqrt(dot_product(r, r))
        
        if (trim(opts%tolerance_type) == "absolute") then
            tolerance = opts%tolerance
        else
            tolerance = opts%tolerance * initial_norm
        end if
        
        p = r
        rho = 1.0_dp
        alpha = 1.0_dp
        omega = 1.0_dp
        
        stats%converged = .false.
        
        do iter = 1, opts%max_iterations
            rho_old = rho
            rho = dot_product(r0, r)
            
            if (abs(rho) < 1.0e-14_dp) exit  ! Breakdown
            
            beta = (rho / rho_old) * (alpha / omega)
            p = r + beta * (p - omega * v)
            
            ! Apply preconditioner to p if available
            if (trim(opts%preconditioner) /= "none") then
                call apply_preconditioner(precond, p, z)
                v = matmul(A, z)
            else
                v = matmul(A, p)
            end if
            
            alpha = rho / dot_product(r0, v)
            s = r - alpha * v
            
            ! Check for early convergence
            residual_norm = sqrt(dot_product(s, s))
            if (residual_norm <= tolerance) then
                if (trim(opts%preconditioner) /= "none") then
                    x = x + alpha * z
                else
                    x = x + alpha * p
                end if
                stats%converged = .true.
                exit
            end if
            
            ! Apply preconditioner to s if available
            if (trim(opts%preconditioner) /= "none") then
                call apply_preconditioner(precond, s, y)
                t = matmul(A, y)
            else
                t = matmul(A, s)
            end if
            
            omega = dot_product(t, s) / dot_product(t, t)
            
            if (trim(opts%preconditioner) /= "none") then
                x = x + alpha * z + omega * y
            else
                x = x + alpha * p + omega * s
            end if
            
            r = s - omega * t
            
            residual_norm = sqrt(dot_product(r, r))
            
            if (opts%verbosity > 0) then
                write(*,'(A,I4,A,E12.5)') "BiCGSTAB iter ", iter, " residual: ", residual_norm
            end if
            
            if (residual_norm <= tolerance) then
                stats%converged = .true.
                exit
            end if
            
            if (abs(omega) < 1.0e-14_dp) exit  ! Breakdown
        end do
        
        stats%iterations = iter
        stats%final_residual = residual_norm
        stats%method_used = "bicgstab"
        
        deallocate(r, r0, p, v, s, t, z, y)
    end subroutine bicgstab_solve
    
    ! GMRES solver with restart
    subroutine gmres_solve(A, b, x, opts, stats)
        real(dp), intent(in) :: A(:,:), b(:)
        real(dp), intent(inout) :: x(:)
        type(solver_options_t), intent(in) :: opts
        type(solver_stats_t), intent(out) :: stats
        
        real(dp), allocatable :: V(:,:), H(:,:), r(:), w(:), y(:), c(:), s(:), g(:)
        real(dp) :: beta, residual_norm, initial_norm, tolerance, temp
        integer :: n, m, iter, total_iter, restart_count, i, j, k
        
        n = size(x)
        m = min(opts%restart, n)
        
        allocate(V(n, m+1), H(m+1, m), r(n), w(n), y(m), c(m), s(m), g(m+1))
        
        ! Initial residual
        r = b - matmul(A, x)
        initial_norm = sqrt(dot_product(r, r))
        
        if (trim(opts%tolerance_type) == "absolute") then
            tolerance = opts%tolerance
        else
            tolerance = opts%tolerance * initial_norm
        end if
        
        total_iter = 0
        restart_count = 0
        stats%converged = .false.
        
        do while (total_iter < opts%max_iterations .and. .not. stats%converged)
            ! Start of restart cycle
            beta = sqrt(dot_product(r, r))
            
            if (beta <= tolerance) then
                stats%converged = .true.
                exit
            end if
            
            V(:, 1) = r / beta
            g = 0.0_dp
            g(1) = beta
            
            ! Arnoldi process
            do j = 1, m
                total_iter = total_iter + 1
                if (total_iter > opts%max_iterations) exit
                
                w = matmul(A, V(:, j))
                
                ! Modified Gram-Schmidt orthogonalization
                do i = 1, j
                    H(i, j) = dot_product(w, V(:, i))
                    w = w - H(i, j) * V(:, i)
                end do
                
                H(j+1, j) = sqrt(dot_product(w, w))
                
                if (H(j+1, j) < 1.0e-14_dp) then
                    m = j  ! Lucky breakdown
                    exit
                end if
                
                V(:, j+1) = w / H(j+1, j)
                
                ! Apply previous Givens rotations
                do k = 1, j-1
                    temp = c(k) * H(k, j) + s(k) * H(k+1, j)
                    H(k+1, j) = -s(k) * H(k, j) + c(k) * H(k+1, j)
                    H(k, j) = temp
                end do
                
                ! Compute new Givens rotation
                if (abs(H(j+1, j)) < 1.0e-14_dp) then
                    c(j) = 1.0_dp
                    s(j) = 0.0_dp
                else
                    temp = sqrt(H(j, j)**2 + H(j+1, j)**2)
                    c(j) = H(j, j) / temp
                    s(j) = H(j+1, j) / temp
                end if
                
                ! Apply new Givens rotation
                H(j, j) = c(j) * H(j, j) + s(j) * H(j+1, j)
                H(j+1, j) = 0.0_dp
                
                ! Update residual norm estimate
                g(j+1) = -s(j) * g(j)
                g(j) = c(j) * g(j)
                
                residual_norm = abs(g(j+1))
                
                if (opts%verbosity > 0) then
                    write(*,'(A,I4,A,E12.5)') "GMRES iter ", total_iter, " residual: ", residual_norm
                end if
                
                if (residual_norm <= tolerance) then
                    m = j
                    stats%converged = .true.
                    exit
                end if
            end do
            
            ! Solve upper triangular system H * y = g
            do i = m, 1, -1
                y(i) = g(i)
                do j = i+1, m
                    y(i) = y(i) - H(i, j) * y(j)
                end do
                y(i) = y(i) / H(i, i)
            end do
            
            ! Update solution
            do j = 1, m
                x = x + y(j) * V(:, j)
            end do
            
            if (.not. stats%converged) then
                ! Compute new residual for restart
                r = b - matmul(A, x)
                restart_count = restart_count + 1
            end if
        end do
        
        stats%iterations = total_iter
        stats%restarts = restart_count
        stats%final_residual = residual_norm
        stats%method_used = "gmres"
        
        deallocate(V, H, r, w, y, c, s, g)
    end subroutine gmres_solve
    
    ! Direct solver using LAPACK
    subroutine lapack_solve(A, b, x, opts, stats)
        real(dp), intent(in) :: A(:,:), b(:)
        real(dp), intent(inout) :: x(:)
        type(solver_options_t), intent(in) :: opts
        type(solver_stats_t), intent(out) :: stats
        
        real(dp), allocatable :: A_copy(:,:)
        integer, allocatable :: ipiv(:)
        integer :: n, info
        
        n = size(x)
        allocate(A_copy(n, n), ipiv(n))
        
        A_copy = A
        x = b
        
        ! LU factorization and solve
        call dgesv(n, 1, A_copy, n, ipiv, x, n, info)
        
        stats%converged = (info == 0)
        stats%iterations = 1
        
        if (stats%converged) then
            stats%final_residual = sqrt(sum((matmul(A, x) - b)**2))
        else
            stats%final_residual = huge(1.0_dp)
        end if
        
        stats%method_used = "lapack_lu"
        
        deallocate(A_copy, ipiv)
    end subroutine lapack_solve
    
    ! Placeholder for sparse LU solver
    subroutine sparse_lu_solve(A, b, x, opts, stats)
        real(dp), intent(in) :: A(:,:), b(:)
        real(dp), intent(inout) :: x(:)
        type(solver_options_t), intent(in) :: opts
        type(solver_stats_t), intent(out) :: stats
        
        ! Fallback to LAPACK for now
        call lapack_solve(A, b, x, opts, stats)
        stats%method_used = "sparse_lu"
    end subroutine sparse_lu_solve
    
    ! Placeholder for UMFPACK solver
    subroutine umfpack_solve(A, b, x, opts, stats)
        real(dp), intent(in) :: A(:,:), b(:)
        real(dp), intent(inout) :: x(:)
        type(solver_options_t), intent(in) :: opts
        type(solver_stats_t), intent(out) :: stats
        
        ! Fallback to LAPACK for now
        call lapack_solve(A, b, x, opts, stats)
        stats%method_used = "umfpack"
    end subroutine umfpack_solve
    
    ! Preconditioner construction
    subroutine build_preconditioner(A, precond, precond_type, opts)
        real(dp), intent(in) :: A(:,:)
        type(preconditioner_t), intent(out) :: precond
        character(len=*), intent(in) :: precond_type
        type(solver_options_t), intent(in) :: opts
        
        integer :: n, i
        
        n = size(A, 1)
        precond%type = precond_type
        
        select case (trim(precond_type))
        case ("jacobi")
            allocate(precond%diagonal(n))
            do i = 1, n
                if (abs(A(i, i)) > 1.0e-14_dp) then
                    precond%diagonal(i) = 1.0_dp / A(i, i)
                else
                    precond%diagonal(i) = 1.0_dp
                end if
            end do
            
        case ("ilu")
            call build_ilu_preconditioner(A, precond, opts)
            
        case default
            ! No preconditioning
            precond%type = "none"
        end select
    end subroutine build_preconditioner
    
    ! Apply preconditioner: solve M * z = r
    subroutine apply_preconditioner(precond, r, z)
        type(preconditioner_t), intent(in) :: precond
        real(dp), intent(in) :: r(:)
        real(dp), intent(out) :: z(:)
        
        integer :: n, i
        
        n = size(r)
        
        select case (trim(precond%type))
        case ("jacobi")
            do i = 1, n
                z(i) = precond%diagonal(i) * r(i)
            end do
            
        case ("ilu")
            call solve_ilu(precond, r, z)
            
        case default
            z = r  ! Identity preconditioner
        end select
    end subroutine apply_preconditioner
    
    ! Build ILU preconditioner
    subroutine build_ilu_preconditioner(A, precond, opts)
        real(dp), intent(in) :: A(:,:)
        type(preconditioner_t), intent(out) :: precond
        type(solver_options_t), intent(in) :: opts
        
        integer :: n, i, j, k
        real(dp) :: factor
        
        n = size(A, 1)
        allocate(precond%L(n, n), precond%U(n, n))
        
        ! Simple ILU(0) factorization
        precond%L = 0.0_dp
        precond%U = 0.0_dp
        
        ! Copy matrix structure
        do i = 1, n
            do j = 1, n
                if (abs(A(i, j)) > 1.0e-14_dp) then
                    if (i > j) then
                        precond%L(i, j) = A(i, j)
                    else
                        precond%U(i, j) = A(i, j)
                    end if
                end if
            end do
            precond%L(i, i) = 1.0_dp
        end do
        
        ! ILU factorization
        do k = 1, n-1
            if (abs(precond%U(k, k)) < 1.0e-14_dp) then
                precond%U(k, k) = 1.0e-14_dp  ! Avoid zero pivot
            end if
            
            do i = k+1, n
                if (abs(precond%L(i, k)) > 1.0e-14_dp) then
                    factor = precond%L(i, k) / precond%U(k, k)
                    precond%L(i, k) = factor
                    
                    do j = k+1, n
                        if (abs(precond%U(i, j)) > 1.0e-14_dp .or. &
                            abs(precond%U(k, j)) > 1.0e-14_dp) then
                            precond%U(i, j) = precond%U(i, j) - factor * precond%U(k, j)
                        end if
                    end do
                end if
            end do
        end do
    end subroutine build_ilu_preconditioner
    
    ! Solve ILU system: (L*U) * z = r
    subroutine solve_ilu(precond, r, z)
        type(preconditioner_t), intent(in) :: precond
        real(dp), intent(in) :: r(:)
        real(dp), intent(out) :: z(:)
        
        real(dp), allocatable :: y(:)
        integer :: n, i, j
        
        n = size(r)
        allocate(y(n))
        
        ! Forward solve: L * y = r
        do i = 1, n
            y(i) = r(i)
            do j = 1, i-1
                y(i) = y(i) - precond%L(i, j) * y(j)
            end do
        end do
        
        ! Backward solve: U * z = y
        do i = n, 1, -1
            z(i) = y(i)
            do j = i+1, n
                z(i) = z(i) - precond%U(i, j) * z(j)
            end do
            z(i) = z(i) / precond%U(i, i)
        end do
        
        deallocate(y)
    end subroutine solve_ilu
    
    ! Method selection heuristics
    function select_best_method(A, opts) result(method)
        real(dp), intent(in) :: A(:,:)
        type(solver_options_t), intent(in) :: opts
        character(len=32) :: method
        
        integer :: n
        logical :: is_symmetric
        real(dp) :: density
        
        n = size(A, 1)
        is_symmetric = check_symmetry(A)
        density = count(abs(A) > 1.0e-14_dp) / real(n*n, dp)
        
        ! Selection heuristics
        if (n <= 500) then
            method = "lapack_lu"
        else if (n <= 5000 .and. density > 0.1_dp) then
            method = "sparse_lu"
        else if (is_symmetric) then
            method = "pcg"
        else
            method = "bicgstab"
        end if
    end function select_best_method
    
    ! Helper functions
    function check_symmetry(A) result(is_symmetric)
        real(dp), intent(in) :: A(:,:)
        logical :: is_symmetric
        integer :: i, j, n
        real(dp) :: tol = 1.0e-12_dp
        
        n = size(A, 1)
        is_symmetric = .true.
        
        do i = 1, n
            do j = i+1, n
                if (abs(A(i, j) - A(j, i)) > tol * max(abs(A(i, j)), abs(A(j, i)))) then
                    is_symmetric = .false.
                    return
                end if
            end do
        end do
    end function check_symmetry
    
    function estimate_memory_usage(n, method) result(memory_bytes)
        integer, intent(in) :: n
        character(len=*), intent(in) :: method
        integer :: memory_bytes
        
        select case (trim(method))
        case ("lapack_lu", "direct")
            memory_bytes = n * n * 8  ! Dense storage
        case ("cg", "pcg", "bicgstab", "gmres")
            memory_bytes = n * 8 * 10  ! Several vectors
        case default
            memory_bytes = n * 8 * 5   ! Conservative estimate
        end select
    end function estimate_memory_usage
    
    ! Constructor function for solver options
    function solver_options(method, tolerance, max_iterations, preconditioner, &
                          tolerance_type, restart, parallel, num_threads, &
                          verbosity, drop_tolerance, fill_level) result(opts)
        character(len=*), intent(in), optional :: method, preconditioner, tolerance_type
        real(dp), intent(in), optional :: tolerance, drop_tolerance
        integer, intent(in), optional :: max_iterations, restart, num_threads, verbosity, fill_level
        logical, intent(in), optional :: parallel
        type(solver_options_t) :: opts
        
        ! Set defaults
        opts%method = "auto"
        opts%tolerance = 1.0e-6_dp
        opts%max_iterations = 1000
        opts%preconditioner = "none"
        opts%tolerance_type = "relative"
        opts%restart = 50
        opts%parallel = .false.
        opts%num_threads = 1
        opts%verbosity = 0
        opts%drop_tolerance = 1.0e-4_dp
        opts%fill_level = 0
        
        ! Override with provided values
        if (present(method)) opts%method = method
        if (present(tolerance)) opts%tolerance = tolerance
        if (present(max_iterations)) opts%max_iterations = max_iterations
        if (present(preconditioner)) opts%preconditioner = preconditioner
        if (present(tolerance_type)) opts%tolerance_type = tolerance_type
        if (present(restart)) opts%restart = restart
        if (present(parallel)) opts%parallel = parallel
        if (present(num_threads)) opts%num_threads = num_threads
        if (present(verbosity)) opts%verbosity = verbosity
        if (present(drop_tolerance)) opts%drop_tolerance = drop_tolerance
        if (present(fill_level)) opts%fill_level = fill_level
    end function solver_options

    ! Specific solver interfaces
    
    function jacobi_preconditioner() result(precond_name)
        character(len=32) :: precond_name
        precond_name = "jacobi"
    end function jacobi_preconditioner
    
    function ilu_preconditioner(drop_tol, fill_level) result(precond_name)
        real(dp), intent(in), optional :: drop_tol
        integer, intent(in), optional :: fill_level
        character(len=32) :: precond_name
        precond_name = "ilu"
    end function ilu_preconditioner

end module fortfem_advanced_solvers