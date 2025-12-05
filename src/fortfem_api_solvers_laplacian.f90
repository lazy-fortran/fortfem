module fortfem_api_solvers_laplacian
    use fortfem_kinds, only: dp
    use fortfem_api_types, only: function_space_t, function_t, &
        dirichlet_bc_t, neumann_bc_t
    use fortfem_api_forms, only: form_equation_t
    use fortfem_solvers_laplacian, only: assemble_laplacian_system, &
        solve_scalar, solve_laplacian_problem, solve_generic_problem
    use fortfem_solvers_p2, only: solve_laplacian_problem_p2
    use fortfem_solvers_neumann, only: solve_mixed_bc, solve_neumann, &
        compute_boundary_integral, solve_laplacian_with_neumann, &
        solve_pure_neumann_problem
    implicit none

    private

    public :: assemble_laplacian_system
    public :: solve_scalar
    public :: solve_laplacian_problem
    public :: solve_laplacian_problem_p2
    public :: solve_mixed_bc
    public :: solve_neumann
    public :: compute_boundary_integral
    public :: solve_laplacian_with_neumann
    public :: solve_pure_neumann_problem
    public :: solve_generic_problem

end module fortfem_api_solvers_laplacian

