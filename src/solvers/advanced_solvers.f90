module fortfem_advanced_solvers
    use fortfem_kinds, only: dp
    use fortfem_sparse_matrix, only: sparse_matrix_t, spmv
    use fortfem_krylov_solvers, only: gmres_impl, bicgstab_impl
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
    public :: solve, solve_sparse, solver_options
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
        real(dp), allocatable :: L(:, :), U(:, :)  ! For ILU
        integer, allocatable :: pivot(:)  ! For ILU
    end type preconditioner_t

    abstract interface
        subroutine matvec_proc(v_in, v_out)
            import :: dp
            real(dp), intent(in) :: v_in(:)
            real(dp), intent(out) :: v_out(:)
        end subroutine matvec_proc

        subroutine precond_proc(r_in, z_out)
            import :: dp
            real(dp), intent(in) :: r_in(:)
            real(dp), intent(out) :: z_out(:)
        end subroutine precond_proc
    end interface

contains

    ! Main solver interface with automatic method selection
    subroutine solve(A, b, x, opts, stats)
        real(dp), intent(in) :: A(:, :), b(:)
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
            write (*, *) "Warning: Unknown solver method, using CG"
            call cg_solve(A, b, x, opts, stats)
        end select

        call cpu_time(end_time)
        stats%solve_time = end_time - start_time
        stats%memory_usage = estimate_memory_usage(n, selected_method)
    end subroutine solve

    subroutine solve_sparse(A_sparse, b, x, opts, stats)
        type(sparse_matrix_t), intent(in) :: A_sparse
        real(dp), intent(in) :: b(:)
        real(dp), intent(inout) :: x(:)
        type(solver_options_t), intent(in) :: opts
        type(solver_stats_t), intent(out) :: stats

        type(solver_options_t) :: local_opts
        type(preconditioner_t) :: precond
        real(dp) :: start_time, end_time
        integer :: n

        n = A_sparse%nrows

        if (A_sparse%ncols /= n) then
            error stop "solve_sparse: only square matrices are supported"
        end if

        if (size(b) /= n .or. size(x) /= n) then
            error stop "solve_sparse: incompatible vector size"
        end if

        local_opts = opts

        if (trim(local_opts%method) == "auto") then
            local_opts%method = "pcg"
        end if

        if (trim(local_opts%preconditioner) /= "none") then
            call build_sparse_preconditioner(A_sparse, precond, &
                                             local_opts%preconditioner, local_opts)
        else
            precond%type = "none"
        end if

        call cpu_time(start_time)

        select case (trim(local_opts%method))
        case ("cg")
            call cg_solve_sparse_impl(A_sparse, b, x, local_opts, stats)
        case ("pcg")
            call pcg_solve_sparse_impl(A_sparse, b, x, local_opts, stats, &
                                       precond)
        case default
            local_opts%method = "pcg"
            call pcg_solve_sparse_impl(A_sparse, b, x, local_opts, stats, &
                                       precond)
        end select

        call cpu_time(end_time)

        stats%solve_time = end_time - start_time
        stats%memory_usage = estimate_sparse_memory_usage(A_sparse, &
                                                          trim(stats%method_used))
    end subroutine solve_sparse

    subroutine cg_solve_sparse_impl(A_sparse, b, x, opts, stats)
        type(sparse_matrix_t), intent(in) :: A_sparse
        real(dp), intent(in) :: b(:)
        real(dp), intent(inout) :: x(:)
        type(solver_options_t), intent(in) :: opts
        type(solver_stats_t), intent(out) :: stats

        call cg_solve_operator(sparse_matvec, b, x, opts, stats)

    contains

        subroutine sparse_matvec(v_in, v_out)
            real(dp), intent(in) :: v_in(:)
            real(dp), intent(out) :: v_out(:)

            call spmv(A_sparse, v_in, v_out)
        end subroutine sparse_matvec

    end subroutine cg_solve_sparse_impl

    subroutine pcg_solve_sparse_impl(A_sparse, b, x, opts, stats, precond)
        type(sparse_matrix_t), intent(in) :: A_sparse
        real(dp), intent(in) :: b(:)
        real(dp), intent(inout) :: x(:)
        type(solver_options_t), intent(in) :: opts
        type(solver_stats_t), intent(out) :: stats
        type(preconditioner_t), intent(in) :: precond

        call pcg_solve_operator(sparse_matvec, sparse_apply_precond, b, x, &
                                opts, stats)

    contains

        subroutine sparse_matvec(v_in, v_out)
            real(dp), intent(in) :: v_in(:)
            real(dp), intent(out) :: v_out(:)

            call spmv(A_sparse, v_in, v_out)
        end subroutine sparse_matvec

        subroutine sparse_apply_precond(r_in, z_out)
            real(dp), intent(in) :: r_in(:)
            real(dp), intent(out) :: z_out(:)

            call apply_preconditioner(precond, r_in, z_out)
        end subroutine sparse_apply_precond

    end subroutine pcg_solve_sparse_impl

    ! Conjugate Gradient solver (dense interface)
    subroutine cg_solve(A, b, x, opts, stats)
        real(dp), intent(in) :: A(:, :), b(:)
        real(dp), intent(inout) :: x(:)
        type(solver_options_t), intent(in) :: opts
        type(solver_stats_t), intent(out) :: stats

        call cg_solve_operator(dense_matvec, b, x, opts, stats)

    contains

        subroutine dense_matvec(v_in, v_out)
            real(dp), intent(in) :: v_in(:)
            real(dp), intent(out) :: v_out(:)

            v_out = matmul(A, v_in)
        end subroutine dense_matvec

    end subroutine cg_solve

    subroutine cg_solve_operator(matvec, b, x, opts, stats)
        procedure(matvec_proc) :: matvec
        real(dp), intent(in) :: b(:)
        real(dp), intent(inout) :: x(:)
        type(solver_options_t), intent(in) :: opts
        type(solver_stats_t), intent(out) :: stats

        real(dp), allocatable :: r(:), p(:), Ap(:), Ax(:)
        real(dp) :: alpha, beta, rr_old, rr_new, pAp
        real(dp) :: residual_norm, initial_norm, tolerance
        integer :: n, iter

        n = size(x)
        allocate (r(n), p(n), Ap(n), Ax(n))

        call matvec(x, Ax)
        r = b - Ax
        p = r
        rr_old = dot_product(r, r)
        initial_norm = sqrt(rr_old)

        if (trim(opts%tolerance_type) == "absolute") then
            tolerance = opts%tolerance
        else
            tolerance = opts%tolerance*initial_norm
        end if

        stats%converged = .false.
        stats%iterations = 0
        residual_norm = sqrt(rr_old)

        do iter = 1, opts%max_iterations
            call matvec(p, Ap)
            pAp = dot_product(p, Ap)

            if (abs(pAp) < 1.0e-14_dp) then
                stats%iterations = iter - 1
                exit
            end if

            alpha = rr_old/pAp
            x = x + alpha*p
            r = r - alpha*Ap

            rr_new = dot_product(r, r)
            residual_norm = sqrt(rr_new)

            if (opts%verbosity > 0) then
                write (*, '(A,I4,A,E12.5)') "CG iter ", iter, " residual: ", &
                    residual_norm
            end if

            if (residual_norm <= tolerance) then
                stats%converged = .true.
                stats%iterations = iter
                exit
            end if

            beta = rr_new/rr_old
            p = r + beta*p

            rr_old = rr_new
            stats%iterations = iter
        end do

        if (.not. stats%converged) then
            if (residual_norm <= tolerance) then
                stats%converged = .true.
            end if
        end if

        stats%final_residual = residual_norm
        stats%method_used = "cg"

        deallocate (r, p, Ap, Ax)
    end subroutine cg_solve_operator

    ! Preconditioned Conjugate Gradient solver (dense interface)
    subroutine pcg_solve(A, b, x, opts, stats)
        real(dp), intent(in) :: A(:, :), b(:)
        real(dp), intent(inout) :: x(:)
        type(solver_options_t), intent(in) :: opts
        type(solver_stats_t), intent(out) :: stats

        type(preconditioner_t) :: precond

        call build_preconditioner(A, precond, opts%preconditioner, opts)

        call pcg_solve_operator(dense_matvec, dense_apply_precond, b, x, &
                                opts, stats)

    contains

        subroutine dense_matvec(v_in, v_out)
            real(dp), intent(in) :: v_in(:)
            real(dp), intent(out) :: v_out(:)

            v_out = matmul(A, v_in)
        end subroutine dense_matvec

        subroutine dense_apply_precond(r_in, z_out)
            real(dp), intent(in) :: r_in(:)
            real(dp), intent(out) :: z_out(:)

            call apply_preconditioner(precond, r_in, z_out)
        end subroutine dense_apply_precond

    end subroutine pcg_solve

    subroutine pcg_solve_operator(matvec, apply_precond, b, x, opts, stats)
        procedure(matvec_proc) :: matvec
        procedure(precond_proc) :: apply_precond
        real(dp), intent(in) :: b(:)
        real(dp), intent(inout) :: x(:)
        type(solver_options_t), intent(in) :: opts
        type(solver_stats_t), intent(out) :: stats

        real(dp), allocatable :: r(:), z(:), p(:), Ap(:), Ax(:)
        real(dp) :: alpha, beta, rz_old, rz_new, pAp
        real(dp) :: residual_norm, initial_norm, tolerance
        integer :: n, iter

        n = size(x)
        allocate (r(n), z(n), p(n), Ap(n), Ax(n))

        call matvec(x, Ax)
        r = b - Ax
        initial_norm = sqrt(dot_product(r, r))

        if (trim(opts%tolerance_type) == "absolute") then
            tolerance = opts%tolerance
        else
            tolerance = opts%tolerance*initial_norm
        end if

        call apply_precond(r, z)
        p = z
        rz_old = dot_product(r, z)

        stats%converged = .false.
        residual_norm = initial_norm
        stats%iterations = 0

        do iter = 1, opts%max_iterations
            call matvec(p, Ap)
            pAp = dot_product(p, Ap)

            if (abs(pAp) < 1.0e-14_dp) then
                exit
            end if

            alpha = rz_old/pAp
            x = x + alpha*p
            r = r - alpha*Ap

            residual_norm = sqrt(dot_product(r, r))

            if (opts%verbosity > 0) then
                write (*, '(A,I4,A,E12.5)') "PCG iter ", iter, " residual: ", &
                    residual_norm
            end if

            if (residual_norm <= tolerance) then
                stats%converged = .true.
                stats%iterations = iter
                exit
            end if

            call apply_precond(r, z)
            rz_new = dot_product(r, z)
            beta = rz_new/rz_old
            p = z + beta*p
            rz_old = rz_new
            stats%iterations = iter
        end do

        if (.not. stats%converged) then
            if (residual_norm <= tolerance) then
                stats%converged = .true.
            end if
        end if

        stats%final_residual = residual_norm
        stats%method_used = "pcg"

        deallocate (r, z, p, Ap, Ax)
    end subroutine pcg_solve_operator

    ! BiCGSTAB solver for non-symmetric systems
    subroutine bicgstab_solve(A, b, x, opts, stats)
        real(dp), intent(in) :: A(:, :), b(:)
        real(dp), intent(inout) :: x(:)
        type(solver_options_t), intent(in) :: opts
        type(solver_stats_t), intent(out) :: stats

        type(preconditioner_t) :: precond
        logical :: use_precond

        use_precond = trim(opts%preconditioner) /= "none"

        if (use_precond) then
            call build_preconditioner(A, precond, opts%preconditioner, opts)
            call bicgstab_impl(A, b, x, precond%diagonal, precond%L, precond%U, &
                               use_precond, opts%tolerance, opts%max_iterations, &
                               opts%tolerance_type, opts%verbosity, stats%converged, &
                               stats%iterations, stats%final_residual)
        else
            call bicgstab_impl(A, b, x, use_precond=.false., tol=opts%tolerance, &
                               max_iter=opts%max_iterations, &
                                   tol_type=opts%tolerance_type, &
                               verbosity=opts%verbosity, converged=stats%converged, &
                               iterations=stats%iterations, &
                                   final_resid=stats%final_residual)
        end if

        stats%method_used = "bicgstab"
    end subroutine bicgstab_solve

    ! GMRES solver with restart
    subroutine gmres_solve(A, b, x, opts, stats)
        real(dp), intent(in) :: A(:, :), b(:)
        real(dp), intent(inout) :: x(:)
        type(solver_options_t), intent(in) :: opts
        type(solver_stats_t), intent(out) :: stats

        call gmres_impl(A, b, x, opts%tolerance, opts%max_iterations, &
                        opts%restart, opts%tolerance_type, opts%verbosity, &
                        stats%converged, stats%iterations, stats%restarts, &
                        stats%final_residual)

        stats%method_used = "gmres"
    end subroutine gmres_solve

    ! Direct solver using LAPACK
    subroutine lapack_solve(A, b, x, opts, stats)
        real(dp), intent(in) :: A(:, :), b(:)
        real(dp), intent(inout) :: x(:)
        type(solver_options_t), intent(in) :: opts
        type(solver_stats_t), intent(out) :: stats

        real(dp), allocatable :: A_copy(:, :)
        integer, allocatable :: ipiv(:)
        integer :: n, info

        n = size(x)
        allocate (A_copy(n, n), ipiv(n))

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

        deallocate (A_copy, ipiv)
    end subroutine lapack_solve

    ! Sparse LU solver entry using dense LU factors
    subroutine sparse_lu_solve(A, b, x, opts, stats)
        real(dp), intent(in) :: A(:, :), b(:)
        real(dp), intent(inout) :: x(:)
        type(solver_options_t), intent(in) :: opts
        type(solver_stats_t), intent(out) :: stats

        call lapack_solve(A, b, x, opts, stats)
        stats%method_used = "sparse_lu"
    end subroutine sparse_lu_solve

    ! UMFPACK-compatible solver entry delegating to LAPACK
    subroutine umfpack_solve(A, b, x, opts, stats)
        real(dp), intent(in) :: A(:, :), b(:)
        real(dp), intent(inout) :: x(:)
        type(solver_options_t), intent(in) :: opts
        type(solver_stats_t), intent(out) :: stats

        call lapack_solve(A, b, x, opts, stats)
        stats%method_used = "umfpack"
    end subroutine umfpack_solve

    ! Preconditioner construction
    subroutine build_preconditioner(A, precond, precond_type, opts)
        real(dp), intent(in) :: A(:, :)
        type(preconditioner_t), intent(out) :: precond
        character(len=*), intent(in) :: precond_type
        type(solver_options_t), intent(in) :: opts

        integer :: n, i

        n = size(A, 1)
        precond%type = precond_type

        select case (trim(precond_type))
        case ("jacobi")
            allocate (precond%diagonal(n))
            do i = 1, n
                if (abs(A(i, i)) > 1.0e-14_dp) then
                    precond%diagonal(i) = 1.0_dp/A(i, i)
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

    subroutine build_sparse_preconditioner(A_sparse, precond, precond_type, &
                                           opts)
        type(sparse_matrix_t), intent(in) :: A_sparse
        type(preconditioner_t), intent(out) :: precond
        character(len=*), intent(in) :: precond_type
        type(solver_options_t), intent(in) :: opts

        integer :: n, i, j, k
        integer :: row_start, row_end
        real(dp) :: a_ii
        real(dp), allocatable :: A_dense(:, :)

        n = A_sparse%nrows
        precond%type = precond_type

        select case (trim(precond_type))
        case ("jacobi")
            allocate (precond%diagonal(n))
            do i = 1, n
                a_ii = 0.0_dp
                row_start = A_sparse%row_ptr(i)
                row_end = A_sparse%row_ptr(i + 1) - 1
                do k = row_start, row_end
                    if (A_sparse%col_ind(k) == i) then
                        a_ii = A_sparse%values(k)
                        exit
                    end if
                end do
                if (abs(a_ii) > 1.0e-14_dp) then
                    precond%diagonal(i) = 1.0_dp/a_ii
                else
                    precond%diagonal(i) = 1.0_dp
                end if
            end do

        case ("ilu")
            allocate (A_dense(n, n))
            A_dense = 0.0_dp

            do i = 1, n
                row_start = A_sparse%row_ptr(i)
                row_end = A_sparse%row_ptr(i + 1) - 1
                do k = row_start, row_end
                    j = A_sparse%col_ind(k)
                    A_dense(i, j) = A_sparse%values(k)
                end do
            end do

            call build_ilu_preconditioner(A_dense, precond, opts)

            deallocate (A_dense)

        case default
            precond%type = "none"
        end select
    end subroutine build_sparse_preconditioner

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
                z(i) = precond%diagonal(i)*r(i)
            end do

        case ("ilu")
            call solve_ilu(precond, r, z)

        case default
            z = r  ! Identity preconditioner
        end select
    end subroutine apply_preconditioner

    ! Build ILU preconditioner
    subroutine build_ilu_preconditioner(A, precond, opts)
        real(dp), intent(in) :: A(:, :)
        type(preconditioner_t), intent(out) :: precond
        type(solver_options_t), intent(in) :: opts

        integer :: n, i, j, k
        real(dp) :: factor

        n = size(A, 1)
        allocate (precond%L(n, n), precond%U(n, n))

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
        do k = 1, n - 1
            if (abs(precond%U(k, k)) < 1.0e-14_dp) then
                if (opts%verbosity > 0) then
                    write (*, *) "ILU warning: near-zero pivot at ", k
                end if
                precond%U(k, k) = sign(1.0e-12_dp, precond%U(k, k))
            end if

            do i = k + 1, n
                if (abs(precond%L(i, k)) > 1.0e-14_dp) then
                    factor = precond%L(i, k)/precond%U(k, k)
                    precond%L(i, k) = factor

                    do j = k + 1, n
                        if (abs(precond%U(i, j)) > 1.0e-14_dp .or. &
                            abs(precond%U(k, j)) > 1.0e-14_dp) then
                            precond%U(i, j) = precond%U(i, j) - factor*precond%U(k, j)
                        end if
                    end do
                end if
            end do
        end do

        ! Check last diagonal element
        if (abs(precond%U(n, n)) < 1.0e-14_dp) then
            if (opts%verbosity > 0) then
                write (*, *) "ILU warning: near-zero pivot at ", n
            end if
            precond%U(n, n) = sign(1.0e-12_dp, precond%U(n, n))
        end if
    end subroutine build_ilu_preconditioner

    ! Solve ILU system: (L*U) * z = r
    subroutine solve_ilu(precond, r, z)
        type(preconditioner_t), intent(in) :: precond
        real(dp), intent(in) :: r(:)
        real(dp), intent(out) :: z(:)

        real(dp), allocatable :: y(:)
        integer :: n, i, j

        n = size(r)
        allocate (y(n))

        ! Forward solve: L * y = r
        do i = 1, n
            y(i) = r(i)
            do j = 1, i - 1
                y(i) = y(i) - precond%L(i, j)*y(j)
            end do
        end do

        ! Backward solve: U * z = y
        do i = n, 1, -1
            z(i) = y(i)
            do j = i + 1, n
                z(i) = z(i) - precond%U(i, j)*z(j)
            end do
            z(i) = z(i)/precond%U(i, i)
        end do

        deallocate (y)
    end subroutine solve_ilu

    ! Method selection heuristics
    function select_best_method(A, opts) result(method)
        real(dp), intent(in) :: A(:, :)
        type(solver_options_t), intent(in) :: opts
        character(len=32) :: method

        integer :: n
        logical :: is_symmetric
        real(dp) :: density

        n = size(A, 1)
        is_symmetric = check_symmetry(A)
        density = count(abs(A) > 1.0e-14_dp)/real(n*n, dp)

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
        real(dp), intent(in) :: A(:, :)
        logical :: is_symmetric
        integer :: i, j, n
        real(dp) :: tol = 1.0e-12_dp

        n = size(A, 1)
        is_symmetric = .true.

        do i = 1, n
            do j = i + 1, n
                if (abs(A(i, j) - A(j, i)) > tol*max(abs(A(i, j)), abs(A(j, i)))) then
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
            memory_bytes = n*n*8  ! Dense storage
        case ("cg", "pcg", "bicgstab", "gmres")
            memory_bytes = n*8*10  ! Several vectors
        case default
            memory_bytes = n*8*5   ! Conservative estimate
        end select
    end function estimate_memory_usage

    function estimate_sparse_memory_usage(A_sparse, method) result(memory_bytes)
        type(sparse_matrix_t), intent(in) :: A_sparse
        character(len=*), intent(in) :: method
        integer :: memory_bytes
        integer :: n, nnz

        n = A_sparse%nrows
        nnz = size(A_sparse%values)

        select case (trim(method))
        case ("lapack_lu", "direct")
            memory_bytes = n*n*8
        case ("cg", "pcg", "bicgstab", "gmres")
            memory_bytes = nnz*8 + n*8*10
        case default
            memory_bytes = nnz*8 + n*8*5
        end select
    end function estimate_sparse_memory_usage

    ! Constructor function for solver options
    function solver_options(method, tolerance, max_iterations, preconditioner, &
                            tolerance_type, restart, parallel, num_threads, &
                            verbosity, drop_tolerance, fill_level) result(opts)
        character(len=*), intent(in), optional :: method, preconditioner, tolerance_type
        real(dp), intent(in), optional :: tolerance, drop_tolerance
        integer, intent(in), optional :: max_iterations, restart, num_threads, &
                                         verbosity, fill_level
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
