module fortfem_api_solvers
    use fortfem_api_solvers_laplacian, only: assemble_laplacian_system, &
        solve_scalar, solve_laplacian_problem, solve_laplacian_problem_p2, &
        solve_mixed_bc, solve_neumann, compute_boundary_integral, &
        solve_laplacian_with_neumann, solve_pure_neumann_problem, &
        solve_generic_problem
    use fortfem_api_solvers_vector, only: solve_vector, &
        solve_curl_curl_problem, solve_generic_vector_problem
    use fortfem_advanced_solvers, only: solver_options_t, solver_stats_t, &
        solver_options, cg_solve, pcg_solve, bicgstab_solve, gmres_solve, &
        jacobi_preconditioner, ilu_preconditioner
    implicit none

    private

    public :: solver_options_t, solver_stats_t
    public :: solver_options
    public :: cg_solve, pcg_solve, bicgstab_solve, gmres_solve
    public :: jacobi_preconditioner, ilu_preconditioner

    public :: assemble_laplacian_system
    public :: solve
    public :: solve_scalar
    public :: solve_vector
    public :: solve_laplacian_problem
    public :: solve_laplacian_problem_p2
    public :: solve_mixed_bc
    public :: solve_neumann
    public :: compute_boundary_integral
    public :: solve_laplacian_with_neumann
    public :: solve_pure_neumann_problem
    public :: solve_curl_curl_problem
    public :: solve_generic_problem
    public :: solve_generic_vector_problem

    interface solve
        module procedure solve_scalar
        module procedure solve_vector
    end interface solve

end module fortfem_api_solvers
