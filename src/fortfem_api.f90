module fortfem_api
    use fortfem_kinds
    use fortfem_boundary, only: boundary_t
    use fortfem_api_types
    use fortfem_api_mesh
    use fortfem_api_spaces
    use fortfem_api_forms
    use fortfem_api_solvers
    use fortfem_api_plot
    use fortfem_planar_helmholtz_dtn, only: apply_planar_helmholtz_dtn
    use fortfem_circular_dtn_2d, only: apply_circular_helmholtz_dtn, &
        circular_helmholtz_dtn_eigenvalue
    use fortfem_laplace_boundary_operators_2d, only: &
        assemble_laplace_adjoint_double_layer_constant, &
        assemble_laplace_double_layer_constant, &
        assemble_laplace_hypersingular_linear, &
        assemble_laplace_single_layer_constant
    use fortfem_helmholtz_boundary_operators_2d, only: &
        assemble_helmholtz_adjoint_double_layer_constant, &
        assemble_helmholtz_double_layer_constant, &
        assemble_helmholtz_hypersingular_linear, &
        assemble_helmholtz_single_layer_constant, &
        assemble_helmholtz_single_layer_linear
    use fortfem_helmholtz_exterior_2d, only: &
        evaluate_helmholtz_combined_potential_adaptive_constant, &
        evaluate_helmholtz_combined_potential_constant, &
        solve_helmholtz_cfie_constant
    use fortfem_triangle_duffy_quadrature, only: triangle_duffy_quadrature
    use fortfem_triangle_lagrange_arbitrary_order, only: &
        assignment(=), evaluate_triangle_lagrange_basis, &
        initialize_triangle_lagrange_basis, triangle_lagrange_basis_t, &
        triangle_lagrange_nodes
    use fortfem_triangle_nedelec_arbitrary_order, only: &
        assignment(=), evaluate_triangle_nedelec_first_kind, &
        initialize_triangle_nedelec_first_kind, triangle_nedelec_dof_count, &
        triangle_nedelec_first_kind_t
    use fortfem_triangle_rt_arbitrary_order, only: &
        assignment(=), evaluate_triangle_raviart_thomas, &
        initialize_triangle_raviart_thomas, triangle_rt_basis_t, &
        triangle_rt_dof_count
    use fortfem_edge_moment_orientation, only: apply_edge_moment_orientation
    use fortfem_triangle_piola_maps, only: &
        map_triangle_nedelec_covariant, map_triangle_rt_contravariant
    use fortfem_edge_interpolation_2d, only: &
        interpolate_axisymmetric_rt_edge_dofs, &
        interpolate_nedelec_edge_dofs, interpolate_rt_edge_dofs
    use fortfem_rt_field_2d, only: evaluate_rt_field_2d, &
        reconstruct_axisymmetric_fourier_toroidal, rt_l2_norm
    use fortfem_sparse_direct, only: sparse_direct_factor_t, &
        sparse_direct_factor_csc, sparse_direct_free, &
        sparse_direct_solve_factored
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
    public :: apply_planar_helmholtz_dtn
    public :: apply_circular_helmholtz_dtn
    public :: circular_helmholtz_dtn_eigenvalue
    public :: evaluate_helmholtz_combined_potential_adaptive_constant
    public :: evaluate_helmholtz_combined_potential_constant
    public :: solve_helmholtz_cfie_constant
    public :: triangle_duffy_quadrature
    public :: assignment(=)
    public :: evaluate_triangle_lagrange_basis
    public :: initialize_triangle_lagrange_basis
    public :: triangle_lagrange_basis_t
    public :: triangle_lagrange_nodes
    public :: evaluate_triangle_nedelec_first_kind
    public :: initialize_triangle_nedelec_first_kind
    public :: triangle_nedelec_dof_count
    public :: triangle_nedelec_first_kind_t
    public :: evaluate_triangle_raviart_thomas
    public :: initialize_triangle_raviart_thomas
    public :: triangle_rt_basis_t
    public :: triangle_rt_dof_count
    public :: apply_edge_moment_orientation
    public :: map_triangle_nedelec_covariant
    public :: map_triangle_rt_contravariant
    public :: assemble_helmholtz_adjoint_double_layer_constant
    public :: assemble_helmholtz_double_layer_constant
    public :: assemble_helmholtz_hypersingular_linear
    public :: assemble_helmholtz_single_layer_constant
    public :: assemble_helmholtz_single_layer_linear
    public :: assemble_laplace_adjoint_double_layer_constant
    public :: assemble_laplace_double_layer_constant
    public :: assemble_laplace_hypersingular_linear
    public :: assemble_laplace_single_layer_constant
    public :: interpolate_nedelec_edge_dofs, interpolate_rt_edge_dofs
    public :: interpolate_axisymmetric_rt_edge_dofs
    public :: evaluate_rt_field_2d
    public :: reconstruct_axisymmetric_fourier_toroidal, rt_l2_norm
    public :: sparse_direct_factor_t, sparse_direct_factor_csc
    public :: sparse_direct_solve_factored, sparse_direct_free

    ! Advanced solver types and functions
    public :: solver_options_t, solver_stats_t
    public :: solver_options
    public :: cg_solve, pcg_solve, bicgstab_solve, gmres_solve
    public :: jacobi_preconditioner, ilu_preconditioner

    ! Plotting interface
    public :: plot

end module fortfem_api
