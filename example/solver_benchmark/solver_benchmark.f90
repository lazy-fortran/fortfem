program solver_benchmark
    use fortfem_kinds, only: dp
    use fortfem_api, only: mesh_t, function_space_t, dirichlet_bc_t, &
                           unit_square_mesh, function_space, dirichlet_bc, &
                           assemble_laplacian_system
    use fortfem_advanced_solvers, only: solver_options_t, solver_stats_t, &
                                        solver_options, solve
    implicit none

    call run_solver_benchmark()

contains

    subroutine run_solver_benchmark()
        integer, parameter :: n_cases = 3
        integer, dimension(n_cases) :: mesh_sizes = [24, 40, 64]
        integer :: i
        integer, allocatable :: dofs(:), pcg_iters(:)
        real(dp), allocatable :: direct_times(:), pcg_times(:)
        real(dp), allocatable :: direct_residuals(:), pcg_residuals(:)
        real(dp) :: speedup
        real(dp), allocatable :: K(:, :), F(:)
        type(mesh_t) :: mesh
        type(function_space_t) :: Vh
        type(dirichlet_bc_t) :: bc

        write (*, *) "=== FortFEM Solver Benchmark (Poisson) ==="

        allocate (dofs(n_cases), pcg_iters(n_cases))
        allocate (direct_times(n_cases), pcg_times(n_cases))
        allocate (direct_residuals(n_cases), pcg_residuals(n_cases))

        do i = 1, n_cases
            call build_laplacian_system(mesh_sizes(i), mesh, Vh, bc, K, F, &
                                        dofs(i))
            call benchmark_solvers(K, F, direct_times(i), pcg_times(i), &
                                   pcg_iters(i), direct_residuals(i), &
                                   pcg_residuals(i))

            speedup = max(direct_times(i), 1.0e-12_dp)/ &
                      max(pcg_times(i), 1.0e-12_dp)

            write (*, '(A,I4,A,I8,A,ES12.4,A,ES12.4,A,I6)') &
                " case ", i, ": DOFs=", dofs(i), "  t_direct=", &
                direct_times(i), "  t_pcg=", pcg_times(i), &
                "  it_pcg=", pcg_iters(i)
            write (*, '(A,ES12.4)') "    speedup_pcg_over_direct=", speedup
        end do

        call write_benchmark_report(mesh_sizes, dofs, direct_times, &
                                    pcg_times, pcg_iters, direct_residuals, &
                                    pcg_residuals)
    end subroutine run_solver_benchmark

    subroutine build_laplacian_system(n_mesh, mesh, Vh, bc, K, F, ndof)
        integer, intent(in) :: n_mesh
        type(mesh_t), intent(out) :: mesh
        type(function_space_t), intent(out) :: Vh
        type(dirichlet_bc_t), intent(out) :: bc
        real(dp), allocatable, intent(out) :: K(:, :), F(:)
        integer, intent(out) :: ndof

        mesh = unit_square_mesh(n_mesh)
        Vh = function_space(mesh, "Lagrange", 1)
        bc = dirichlet_bc(Vh, 0.0_dp)

        call assemble_laplacian_system(Vh, bc, K, F)
        ndof = size(F)
    end subroutine build_laplacian_system

    subroutine benchmark_solvers(K, F, direct_time, pcg_time, pcg_iters, &
                                 direct_residual, pcg_residual)
        real(dp), intent(in) :: K(:, :), F(:)
        real(dp), intent(out) :: direct_time, pcg_time
        integer, intent(out) :: pcg_iters
        real(dp), intent(out) :: direct_residual, pcg_residual

        real(dp), allocatable :: x_direct(:), x_pcg(:)
        type(solver_options_t) :: opts_direct, opts_pcg
        type(solver_stats_t) :: stats_direct, stats_pcg
        integer :: ndof
        real(dp) :: target_tolerance

        ndof = size(F)
        allocate (x_direct(ndof), x_pcg(ndof))

        x_direct = 0.0_dp
        x_pcg = 0.0_dp

        target_tolerance = 1.0e-6_dp

        opts_direct = solver_options(method="lapack_lu", &
                                     tolerance=1.0e-10_dp, &
                                     max_iterations=ndof)
        call solve(K, F, x_direct, opts_direct, stats_direct)

        opts_pcg = solver_options(method="pcg", preconditioner="ilu", &
                                  tolerance=target_tolerance, &
                                  tolerance_type="absolute", &
                                  max_iterations=5*ndof)
        call solve(K, F, x_pcg, opts_pcg, stats_pcg)

        direct_residual = norm2(K, F, x_direct)
        pcg_residual = norm2(K, F, x_pcg)

        direct_time = stats_direct%solve_time
        pcg_time = stats_pcg%solve_time
        pcg_iters = stats_pcg%iterations
    end subroutine benchmark_solvers

    pure function norm2(A, b, x) result(residual_norm)
        real(dp), intent(in) :: A(:, :), b(:), x(:)
        real(dp) :: residual_norm
        real(dp), allocatable :: r(:)

        allocate (r(size(b)))
        r = matmul(A, x) - b
        residual_norm = sqrt(sum(r*r))
    end function norm2

    subroutine write_benchmark_report(mesh_sizes, dofs, direct_times, &
                                      pcg_times, pcg_iters, &
                                      direct_residuals, pcg_residuals)
        integer, intent(in) :: mesh_sizes(:), dofs(:), pcg_iters(:)
        real(dp), intent(in) :: direct_times(:), pcg_times(:)
        real(dp), intent(in) :: direct_residuals(:), pcg_residuals(:)

        integer :: unit_num, ios, i, n_cases
        character(len=*), parameter :: filename = &
                                       "artifacts/solver_benchmarks/"// &
                                       "poisson_solver_benchmark.txt"

        n_cases = size(mesh_sizes)

        call ensure_benchmark_directory()

        open (newunit=unit_num, file=filename, status="replace", &
              action="write", iostat=ios)

        if (ios /= 0) then
            write (*, *) "Warning: Could not write solver benchmark file"
            return
        end if

        write (unit_num, '(A)') "# FortFEM Poisson solver benchmark"
        write (unit_num, '(A)') "# mesh_size  dofs  t_direct  t_pcg  it_pcg "// &
            " residual_direct  residual_pcg"

        do i = 1, n_cases
            write (unit_num, &
                   '(I6,1X,I8,1X,ES15.8,1X,ES15.8,1X,I6,1X,ES15.8,1X,ES15.8)') &
                mesh_sizes(i), dofs(i), direct_times(i), pcg_times(i), &
                pcg_iters(i), direct_residuals(i), pcg_residuals(i)
        end do

        close (unit_num)

        write (*, *) "Benchmark report written to", trim(filename)
    end subroutine write_benchmark_report

    subroutine ensure_benchmark_directory()
        logical :: exists

        inquire (file="artifacts/solver_benchmarks", exist=exists)
        if (.not. exists) then
            call execute_command_line("mkdir -p artifacts/solver_benchmarks")
        end if
    end subroutine ensure_benchmark_directory

end program solver_benchmark
