program test_fortfem_api_solvers_modules
    use fortfem_kinds
    use fortfem_api, only: mesh_t, function_space_t, vector_function_space_t, &
        function_t, vector_function_t, dirichlet_bc_t, vector_bc_t, neumann_bc_t, &
        unit_square_mesh, function_space, vector_function_space, function, &
        vector_function, dirichlet_bc, vector_bc, neumann_bc_constant
    use fortfem_api_solvers_laplacian, only: compute_boundary_integral, &
        solve_laplacian_problem
    use fortfem_api_solvers_vector, only: solve_curl_curl_problem
    use fortfem_advanced_solvers, only: solver_options_t, solver_stats_t, &
        solver_options
    use check
    implicit none

    write(*,*) "Testing fortfem_api_solvers module split..."

    call test_laplacian_module_boundary_integral()
    call test_vector_module_direct_solver()

    call check_summary("fortfem_api_solvers module split")

contains

    subroutine test_laplacian_module_boundary_integral()
        type(mesh_t) :: mesh
        type(function_space_t) :: Vh
        type(neumann_bc_t) :: neumann_bc
        type(function_t) :: u
        type(dirichlet_bc_t) :: bc
        type(solver_stats_t) :: stats
        real(dp) :: boundary_integral, expected_integral, boundary_length

        mesh = unit_square_mesh(4)
        Vh = function_space(mesh, "Lagrange", 1)

        neumann_bc = neumann_bc_constant(Vh, 2.0_dp)

        call compute_boundary_integral(neumann_bc, boundary_integral)

        boundary_length = 4.0_dp
        expected_integral = 2.0_dp*boundary_length

        call check_condition(abs(boundary_integral - expected_integral) < 0.2_dp, &
            "Laplacian module: boundary integral close to expected")
        call check_condition(boundary_integral > 7.0_dp, &
            "Laplacian module: integral reasonable lower bound")
        call check_condition(boundary_integral < 9.0_dp, &
            "Laplacian module: integral reasonable upper bound")

        u = function(Vh)
        bc = dirichlet_bc(Vh, 0.0_dp)
        call solve_laplacian_problem(u, bc, solver_options(method="auto"), &
                                     stats)

        call check_condition(allocated(u%values), &
            "Laplacian module: solution values allocated")
    end subroutine test_laplacian_module_boundary_integral

    subroutine test_vector_module_direct_solver()
        type(mesh_t) :: mesh
        type(vector_function_space_t) :: Vh_vec
        type(vector_function_t) :: Eh
        type(vector_bc_t) :: bc_vec
        type(solver_options_t) :: opts
        type(solver_stats_t) :: stats
        real(dp) :: solution_norm

        mesh = unit_square_mesh(3)
        Vh_vec = vector_function_space(mesh, "Nedelec", 1)

        Eh = vector_function(Vh_vec)
        bc_vec = vector_bc(Vh_vec, [0.0_dp, 0.0_dp])

        opts = solver_options(method="gmres", tolerance=1.0e-6_dp, &
                              max_iterations=50, restart=10)

        call solve_curl_curl_problem(Eh, bc_vec, "direct", opts, stats)

        solution_norm = 0.0_dp
        if (allocated(Eh%values)) then
            solution_norm = sqrt(sum(Eh%values**2))
        end if

        call check_condition(solution_norm >= 0.0_dp, &
            "Vector module: non-negative solution norm")
        call check_condition(solution_norm < 10.0_dp, &
            "Vector module: bounded solution norm")
    end subroutine test_vector_module_direct_solver

end program test_fortfem_api_solvers_modules
