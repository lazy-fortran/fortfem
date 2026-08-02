program solver_benchmark
    use fortfem_kinds, only: dp
    use fortfem_core, only: mesh_t, function_space_t, dirichlet_bc_t, &
        unit_square_mesh
    use fortfem_feec, only: function_space, dirichlet_bc, &
        assemble_laplacian_system, sparse_from_dense, sparse_matrix_t, &
        build_sparse_ilut_row, build_sparse_ichol_row, &
        sparse_ilut_factor_t, sparse_incomplete_cholesky_factor_t, &
        apply_sparse_ilut, apply_sparse_incomplete_cholesky
    use fortfem_advanced_solvers, only: solver_options_t, solver_stats_t, &
        solver_options, solve
    use fortsparse, only: fortsparse_status_t
    use fortplot, only: colorbar, figure, legend, pcolormesh, plot, savefig, &
        set_xscale, set_yscale, title, xlabel, ylabel
    implicit none

    character(*), parameter :: output_directory = &
        "output/example/solver_benchmark"

    call run_solver_benchmark()

contains

    subroutine run_solver_benchmark()
        integer, parameter :: n_cases = 3
        integer, dimension(n_cases) :: mesh_sizes = [24, 40, 64]
        integer :: i
        integer, allocatable :: dofs(:), pcg_iters(:)
        real(dp), allocatable :: direct_times(:), pcg_times(:)
        real(dp), allocatable :: direct_residuals(:), pcg_residuals(:)
        real(dp), allocatable :: row_ilut_times(:), row_ichol_times(:)
        real(dp) :: speedup
        real(dp), allocatable :: K(:, :), F(:)
        type(mesh_t) :: mesh
        type(function_space_t) :: Vh
        type(dirichlet_bc_t) :: bc

        write (*, *) "=== FortFEM Solver Benchmark (Poisson) ==="

        allocate (dofs(n_cases), pcg_iters(n_cases))
        allocate (direct_times(n_cases), pcg_times(n_cases))
        allocate (direct_residuals(n_cases), pcg_residuals(n_cases))
        allocate (row_ilut_times(n_cases), row_ichol_times(n_cases))

        do i = 1, n_cases
            call build_laplacian_system(mesh_sizes(i), mesh, Vh, bc, K, F, &
                dofs(i))
            call benchmark_solvers(K, F, direct_times(i), pcg_times(i), &
                pcg_iters(i), direct_residuals(i), &
                pcg_residuals(i))
            call benchmark_sparse_factors(K, F, row_ilut_times(i), &
                row_ichol_times(i))

            speedup = max(direct_times(i), 1.0e-12_dp)/ &
                max(pcg_times(i), 1.0e-12_dp)

            write (*, '(A,I4,A,I8,A,ES12.4,A,ES12.4,A,I6)') &
                " case ", i, ": DOFs=", dofs(i), "  t_direct=", &
                direct_times(i), "  t_pcg=", pcg_times(i), &
                "  it_pcg=", pcg_iters(i)
            write (*, '(A,ES12.4)') "    speedup_pcg_over_direct=", speedup
            write (*, '(A,2(ES12.4,1X))') &
                "    t_row_ilut/t_row_ichol=", row_ilut_times(i), &
                row_ichol_times(i)
        end do

        call write_benchmark_report(mesh_sizes, dofs, direct_times, &
            pcg_times, pcg_iters, direct_residuals, &
            pcg_residuals, row_ilut_times, row_ichol_times)

        call ensure_example_output_directory()
        call render_reference_solution()
        call plot_solver_times(mesh_sizes, direct_times, pcg_times)
        call plot_solver_residuals(mesh_sizes, direct_residuals, pcg_residuals)
        call plot_sparse_factor_times(mesh_sizes, row_ilut_times, &
            row_ichol_times)
    end subroutine run_solver_benchmark

    subroutine render_reference_solution()
        integer, parameter :: nx = 32, ny = 32
        real(dp) :: x_edges(nx + 1), y_edges(ny + 1), values(nx, ny)
        real(dp) :: x, y
        integer :: i, j

        do i = 1, nx + 1
            x_edges(i) = real(i - 1, dp)/real(nx, dp)
        end do
        do j = 1, ny + 1
            y_edges(j) = real(j - 1, dp)/real(ny, dp)
        end do
        do j = 1, ny
            y = 0.5_dp*(y_edges(j) + y_edges(j + 1))
            do i = 1, nx
                x = 0.5_dp*(x_edges(i) + x_edges(i + 1))
                values(i, j) = sin(acos(-1.0_dp)*x)* &
                    sin(acos(-1.0_dp)*y)
            end do
        end do

        call figure(figsize=[8.5_dp, 5.5_dp])
        call pcolormesh(x_edges, y_edges, values, cmap="viridis")
        call colorbar(label="manufactured Poisson solution u")
        call xlabel("x")
        call ylabel("y")
        call title("Poisson solution used by the solver benchmark")
        call savefig(output_directory//"/poisson_solution_2d.png")
    end subroutine render_reference_solution

    subroutine plot_solver_times(mesh_sizes, direct_times, pcg_times)
        integer, intent(in) :: mesh_sizes(:)
        real(dp), intent(in) :: direct_times(:), pcg_times(:)

        real(dp) :: sizes(size(mesh_sizes)), direct_plot(size(mesh_sizes))
        real(dp) :: pcg_plot(size(mesh_sizes))

        sizes = real(mesh_sizes, dp)
        direct_plot = max(direct_times, tiny(1.0_dp))
        pcg_plot = max(pcg_times, tiny(1.0_dp))
        call figure(figsize=[9.0_dp, 5.5_dp])
        call plot(sizes, direct_plot, &
            label="dense direct", marker="o")
        call plot(sizes, pcg_plot, &
            label="PCG + ILU", marker="s")
        call set_xscale("log")
        call set_yscale("log")
        call xlabel("mesh resolution")
        call ylabel("solve time [s]")
        call title("Poisson solver timing benchmark")
        call legend()
        call savefig(output_directory//"/solver_times_1d.png")
    end subroutine plot_solver_times

    subroutine plot_solver_residuals(mesh_sizes, direct_residuals, pcg_residuals)
        integer, intent(in) :: mesh_sizes(:)
        real(dp), intent(in) :: direct_residuals(:), pcg_residuals(:)

        real(dp) :: sizes(size(mesh_sizes)), direct_plot(size(mesh_sizes))
        real(dp) :: pcg_plot(size(mesh_sizes))

        sizes = real(mesh_sizes, dp)
        direct_plot = max(direct_residuals, tiny(1.0_dp))
        pcg_plot = max(pcg_residuals, tiny(1.0_dp))
        call figure(figsize=[9.0_dp, 5.5_dp])
        call plot(sizes, direct_plot, &
            label="dense direct", marker="o")
        call plot(sizes, pcg_plot, &
            label="PCG + ILU", marker="s")
        call set_xscale("log")
        call set_yscale("log")
        call xlabel("mesh resolution")
        call ylabel("residual norm")
        call title("Poisson solver residual benchmark")
        call legend()
        call savefig(output_directory//"/solver_residuals_1d.png")
    end subroutine plot_solver_residuals

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

    subroutine benchmark_sparse_factors(K, F, ilut_time, ichol_time)
        real(dp), intent(in) :: K(:, :), F(:)
        real(dp), intent(out) :: ilut_time, ichol_time

        type(sparse_matrix_t) :: sparse_K
        type(sparse_ilut_factor_t) :: ilut_factor
        type(sparse_incomplete_cholesky_factor_t) :: ichol_factor
        type(fortsparse_status_t) :: status
        real(dp), allocatable :: factor_apply(:), symmetric_K(:, :)
        real(dp) :: start_time, end_time
        real(dp) :: diagonal_shift
        integer :: max_fill, row

        ! The assembled Dirichlet matrix is not guaranteed to retain the
        ! eliminated columns symmetrically.  Use a shifted symmetric copy for
        ! the SPD ICHOL comparison; the direct/PCG numbers above use the
        ! original Poisson system.
        allocate(symmetric_K(size(K, 1), size(K, 2)))
        symmetric_K = 0.5_dp*(K + transpose(K))
        diagonal_shift = 1.0e-10_dp*max(1.0_dp, maxval(abs(K)))
        do row = 1, size(symmetric_K, 1)
            symmetric_K(row, row) = symmetric_K(row, row) + diagonal_shift
        end do
        call sparse_from_dense(symmetric_K, sparse_K)
        deallocate(symmetric_K)
        max_fill = min(32, size(F))

        call cpu_time(start_time)
        call build_sparse_ilut_row( &
            sparse_K, 0.0_dp, max_fill, ilut_factor, status)
        call cpu_time(end_time)
        if (status%code /= 0) error stop "row ILUT benchmark factorization failed"
        ilut_time = max(end_time - start_time, tiny(1.0_dp))
        call apply_sparse_ilut(ilut_factor, F, factor_apply, status)
        if (status%code /= 0) error stop "row ILUT benchmark apply failed"
        deallocate(factor_apply)

        call cpu_time(start_time)
        call build_sparse_ichol_row( &
            sparse_K, 0.0_dp, max_fill, ichol_factor, status)
        call cpu_time(end_time)
        if (status%code /= 0) then
            write (*, '(A,1X,I0,1X,A)') "row ICHOL status:", status%code, &
                trim(status%msg)
            error stop "row ICHOL benchmark factorization failed"
        end if
        ichol_time = max(end_time - start_time, tiny(1.0_dp))
        call apply_sparse_incomplete_cholesky( &
            ichol_factor, F, factor_apply, status)
        if (status%code /= 0) error stop "row ICHOL benchmark apply failed"
        deallocate(factor_apply)
    end subroutine benchmark_sparse_factors

    subroutine plot_sparse_factor_times(mesh_sizes, ilut_times, ichol_times)
        integer, intent(in) :: mesh_sizes(:)
        real(dp), intent(in) :: ilut_times(:), ichol_times(:)

        real(dp) :: sizes(size(mesh_sizes))
        real(dp) :: ilut_plot(size(mesh_sizes)), ichol_plot(size(mesh_sizes))

        sizes = real(mesh_sizes, dp)
        ilut_plot = max(ilut_times, tiny(1.0_dp))
        ichol_plot = max(ichol_times, tiny(1.0_dp))
        call figure(figsize=[9.0_dp, 5.5_dp])
        call plot(sizes, ilut_plot, label="row ILUT construction", marker="o")
        call plot(sizes, ichol_plot, label="row ICHOL construction", marker="s")
        call set_xscale("log")
        call set_yscale("log")
        call xlabel("mesh resolution")
        call ylabel("factorization time [s]")
        call title("Memory-scalable sparse factor construction")
        call legend()
        call savefig(output_directory//"/sparse_factor_times_1d.png")
    end subroutine plot_sparse_factor_times

    pure function norm2(A, b, x) result(residual_norm)
        real(dp), intent(in) :: A(:, :), b(:), x(:)
        real(dp) :: residual_norm
        real(dp), allocatable :: r(:)

        allocate (r(size(b)))
        r = matmul(A, x) - b
        residual_norm = sqrt(sum(r*r))
    end function norm2

    subroutine write_benchmark_report(mesh_sizes, dofs, direct_times, &
            pcg_times, pcg_iters, direct_residuals, pcg_residuals, &
            row_ilut_times, row_ichol_times)
        integer, intent(in) :: mesh_sizes(:), dofs(:), pcg_iters(:)
        real(dp), intent(in) :: direct_times(:), pcg_times(:)
        real(dp), intent(in) :: direct_residuals(:), pcg_residuals(:)
        real(dp), intent(in) :: row_ilut_times(:), row_ichol_times(:)

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
        write (unit_num, '(A)') &
            "# row factors use a shifted symmetric copy for ICHOL"
        write (unit_num, '(A)') "# mesh_size  dofs  t_direct  t_pcg  it_pcg "// &
            " residual_direct  residual_pcg  t_row_ilut  t_row_ichol"

        do i = 1, n_cases
            write (unit_num, &
                '(I6,1X,I8,1X,ES15.8,1X,ES15.8,1X,I6,1X,ES15.8,1X,ES15.8,'// &
                '1X,ES15.8,1X,ES15.8)') &
                mesh_sizes(i), dofs(i), direct_times(i), pcg_times(i), &
                pcg_iters(i), direct_residuals(i), pcg_residuals(i), &
                row_ilut_times(i), row_ichol_times(i)
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

    subroutine ensure_example_output_directory()
        logical :: exists

        inquire (file=output_directory, exist=exists)
        if (.not. exists) then
            call execute_command_line("mkdir -p "//output_directory)
        end if
    end subroutine ensure_example_output_directory

end program solver_benchmark
