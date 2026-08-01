module fortfem_api_solvers
    use fortfem_api_solvers_laplacian, only: assemble_laplacian_system, &
        solve_scalar, solve_laplacian_problem, solve_laplacian_problem_p2, &
        solve_mixed_bc, solve_neumann, compute_boundary_integral, &
        solve_laplacian_with_neumann, solve_pure_neumann_problem, &
        solve_generic_problem
    use fortfem_api_solvers_vector, only: solve_vector, solve_vector_mixed, &
        solve_curl_curl_problem, solve_generic_vector_problem
    use fortfem_advanced_solvers, only: solver_options_t, solver_stats_t, &
        solve_sparse, solver_options, cg_solve, pcg_solve, bicgstab_solve, &
        gmres_solve, &
        cg_solve_jvp, cg_solve_vjp, pcg_solve_jvp, pcg_solve_vjp, &
        bicgstab_solve_jvp, bicgstab_solve_vjp, gmres_solve_jvp, &
        gmres_solve_vjp, &
        jacobi_preconditioner, ilu_preconditioner, ichol_controlled_preconditioner
    use fortfem_sparse_matrix, only: sparse_matrix_t, sparse_from_dense, &
        spmv, spmv_jvp, spmv_vjp
    use fortfem_sparse_direct, only: sparse_direct_factor_t, &
        sparse_direct_factor_csc, sparse_direct_factor_transpose_csc, &
        sparse_direct_solve_csc, sparse_direct_solve_factored, &
        sparse_direct_solve_factored_jvp, sparse_direct_solve_factored_vjp, &
        sparse_direct_free
    implicit none

    private

    public :: solver_options_t, solver_stats_t
    public :: solver_options
    public :: solve_sparse
    public :: cg_solve, pcg_solve, bicgstab_solve, gmres_solve
    public :: cg_solve_jvp, cg_solve_vjp
    public :: pcg_solve_jvp, pcg_solve_vjp
    public :: bicgstab_solve_jvp, bicgstab_solve_vjp
    public :: gmres_solve_jvp, gmres_solve_vjp
    public :: jacobi_preconditioner, ilu_preconditioner, &
        ichol_controlled_preconditioner
    public :: sparse_matrix_t, sparse_from_dense, spmv, spmv_jvp, spmv_vjp
    public :: sparse_direct_factor_t
    public :: sparse_direct_factor_csc, sparse_direct_factor_transpose_csc
    public :: sparse_direct_solve_csc
    public :: sparse_direct_solve_factored, sparse_direct_solve_factored_jvp
    public :: sparse_direct_solve_factored_vjp, sparse_direct_free

    public :: assemble_laplacian_system
    public :: solve
    public :: solve_scalar
    public :: solve_vector
    public :: solve_vector_mixed
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
        module procedure solve_vector_mixed
    end interface solve

end module fortfem_api_solvers
