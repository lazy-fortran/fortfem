module fortfem_api
    use fortfem_kinds
    use fortfem_boundary, only: boundary_t
    use fortfem_api_types
    use fortfem_api_mesh
    use fortfem_api_spaces
    use fortfem_api_forms
    use fortfem_api_solvers
    use fortfem_api_plot
    use fortfem_advanced_solvers, only: solver_options_t, solver_stats_t,   &
        solver_options, cg_solve, pcg_solve, bicgstab_solve, gmres_solve,   &
        jacobi_preconditioner, ilu_preconditioner
    implicit none

    private

    ! Public types
    public :: mesh_t
    public :: function_space_t
    public :: vector_function_space_t
    public :: function_t
    public :: vector_function_t
    public :: trial_function_t
    public :: test_function_t
    public :: vector_trial_function_t
    public :: vector_test_function_t
    public :: dirichlet_bc_t
    public :: vector_bc_t
    public :: neumann_bc_t
    public :: boundary_t
    public :: simple_expression_t
    public :: form_expr_t
    public :: form_equation_t

    ! Public mesh constructors and refinement
    public :: unit_square_mesh
    public :: rectangle_mesh
    public :: unit_disk_mesh
    public :: circle_boundary
    public :: rectangle_boundary
    public :: line_segment
    public :: arc_segment
    public :: l_shape_boundary
    public :: mesh_from_boundary
    public :: mesh_from_arrays
    public :: mesh_from_triangle_files
    public :: mesh_from_domain
    public :: structured_quad_mesh
    public :: refine_uniform
    public :: refine_adaptive

    ! Function space and field constructors
    public :: function_space
    public :: vector_function_space
    public :: function
    public :: vector_function
    public :: trial_function
    public :: test_function
    public :: vector_trial_function
    public :: vector_test_function
    public :: constant
    public :: dirichlet_bc
    public :: dirichlet_bc_on_boundary
    public :: vector_bc
    public :: neumann_bc_constant
    public :: neumann_bc_on_boundary

    ! Form operations
    public :: inner
    public :: grad
    public :: curl
    public :: dx
    public :: compile_form
    public :: operator(*)
    public :: operator(+)
    public :: operator(==)

    ! Solver interface
    public :: assemble_laplacian_system
    public :: solve
    public :: solve_mixed_bc
    public :: solve_neumann
    public :: compute_boundary_integral

    ! Advanced solver types and functions
    public :: solver_options_t, solver_stats_t
    public :: solver_options
    public :: cg_solve, pcg_solve, bicgstab_solve, gmres_solve
    public :: jacobi_preconditioner, ilu_preconditioner

    ! Plotting interface
    public :: plot

end module fortfem_api

