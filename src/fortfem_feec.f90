module fortfem_feec
    !! Canonical facade for metric-independent FEEC diagnostics and gauges.
    !!
    !! The facade owns no implementation.  It exposes the de Rham
    !! composition/commuting-diagram contracts and the tree-cotree reduction
    !! used by direct H(curl) solves, while leaving element and assembly
    !! implementations in their existing lower-layer modules.  In
    !! particular, this module deliberately does not import ``fortfem_api``.
    use fortfem_feec_exact_sequence, only: &
        assemble_feec_exact_sequence, &
        assemble_feec_exact_sequence_jvp, &
        assemble_feec_exact_sequence_vjp
    use fortfem_feec_commuting_projection, only: &
        assemble_feec_commuting_projection, &
        assemble_feec_commuting_projection_jvp, &
        assemble_feec_commuting_projection_vjp
    use fortfem_tetra_nedelec_arbitrary_order, only: &
        evaluate_tetra_nedelec_first_kind, &
        evaluate_tetra_nedelec_first_kind_jvp, &
        evaluate_tetra_nedelec_first_kind_vjp, &
        initialize_tetra_nedelec_first_kind, &
        tetra_nedelec_dof_count, tetra_nedelec_first_kind_t
    use fortfem_tetra_nedelec_interpolation, only: &
        interpolate_reference_tetra_nedelec
    use fortfem_tetra_vector_evaluation, only: &
        evaluate_tetra_nedelec_interpolant_at_point
    use fortfem_tetra_edge_dof_map, only: build_tetra_edge_dof_map
    use fortfem_tetra_duffy_quadrature, only: tetra_duffy_quadrature
    use fortfem_tetra_nedelec_global_dof_map, only: build_tetra_nedelec_dof_map
    use fortfem_tetra_piola_maps, only: map_tetra_nedelec_covariant
    use fortfem_assembly_bspline_polar_2d, only: &
        assemble_bspline_polar_h1_operator_csc, &
        assemble_bspline_polar_hcurl_operator_csc, &
        assemble_bspline_polar_l2_mass_csc, restrict_bspline_polar_operator_csc
    use fortfem_assembly_bspline_2d, only: &
        assemble_bspline_h1_operator_csc, &
        assemble_bspline_h1_operator_csc_jvp, &
        assemble_bspline_h1_operator_csc_vjp
    use fortfem_bspline_feec, only: evaluate_bspline_basis
    use fortfem_bspline_polar, only: &
        build_bspline_polar_feec_2d_operators, &
        build_bspline_polar_feec_2d_extractions, &
        build_bspline_polar_h1_extraction, evaluate_periodic_bspline_basis
    use fortfem_field_aligned_constitutive_tensor, only: &
        evaluate_field_aligned_constitutive_tensor, &
        evaluate_field_aligned_constitutive_tensor_jvp, &
        evaluate_field_aligned_constitutive_tensor_vjp
    use fortfem_tensor_diffusion_matrix, only: &
        assemble_tensor_diffusion_matrix, &
        assemble_tensor_diffusion_matrix_jvp, &
        assemble_tensor_diffusion_matrix_vjp, &
        assemble_tensor_diffusion_matrix_3d, &
        assemble_tensor_diffusion_matrix_3d_jvp, &
        assemble_tensor_diffusion_matrix_3d_vjp
    use fortfem_field_aligned_flux, only: &
        evaluate_field_aligned_flux, evaluate_field_aligned_flux_jvp, &
        evaluate_field_aligned_flux_vjp
    use fortfem_fci_parallel_operator, only: &
        apply_fci_parallel_diffusion, apply_fci_parallel_diffusion_jvp, &
        apply_fci_parallel_diffusion_vjp, apply_fci_parallel_diffusion_field_vjp, &
        apply_fci_parallel_gradient, apply_fci_parallel_gradient_jvp, &
        apply_fci_parallel_gradient_vjp
    use fortfem_tetra_nedelec_solver_3d, only: &
        solve_tetra_nedelec_curl_mass, solve_tetra_nedelec_pml
    use fortfem_fci_field_line_tracer, only: trace_fci_field_line_rk4
    use fortfem_fci_support_geometry, only: &
        compute_fci_quadrilateral_cell_areas_2d, &
        compute_fci_curved_quadrilateral_cell_areas_2d, &
        compute_fci_polygon_cell_areas_2d, &
        compute_fci_curved_polygon_cell_areas_2d, &
        compute_fci_cubic_curved_polygon_cell_areas_2d, &
        compute_fci_quartic_curved_polygon_cell_areas_2d, &
        compute_fci_quintic_curved_polygon_cell_areas_2d, &
        compute_fci_sextic_curved_polygon_cell_areas_2d
    use fortfem_cgl_pressure_tensor, only: &
        evaluate_cgl_pressure_tensor, evaluate_cgl_pressure_tensor_jvp, &
        evaluate_cgl_pressure_tensor_vjp, evaluate_cgl_pressure_traction, &
        evaluate_cgl_pressure_traction_jvp, evaluate_cgl_pressure_traction_vjp, &
        evaluate_cgl_pressure_work, evaluate_cgl_pressure_work_jvp, &
        evaluate_cgl_pressure_work_vjp
    use fortfem_cgl_pressure_divergence, only: &
        evaluate_cgl_pressure_divergence, &
        evaluate_cgl_pressure_divergence_jvp, evaluate_cgl_pressure_divergence_vjp
    use fortfem_mixed_elasticity_residual, only: &
        assemble_mixed_elasticity_residual, &
        assemble_mixed_elasticity_residual_jvp, &
        assemble_mixed_elasticity_residual_vjp
    use fortfem_elasticity_symmetry_constraint, only: &
        assemble_elasticity_symmetry_constraint, &
        assemble_elasticity_symmetry_constraint_jvp, &
        assemble_elasticity_symmetry_constraint_vjp
    use fortfem_tree_cotree_gauge, only: &
        apply_tree_cotree_prolongation, &
        apply_tree_cotree_restriction, &
        build_tree_cotree_dof_map, &
        build_tree_cotree_gauge, &
        reduce_tree_cotree_dense_system, &
        reduce_tree_cotree_dense_system_jvp, &
        reduce_tree_cotree_dense_system_vjp, &
        tree_cotree_gauge_components, &
        tree_cotree_gauge_edges, &
        tree_cotree_gauge_t, &
        validate_tree_cotree_gauge
    use fortfem_tree_cotree_iga_parity, only: &
        diagnose_tree_cotree_iga_invariance, &
        tree_cotree_iga_parity_t
    use fortfem_beltrami_parity, only: &
        beltrami_parity_t, compare_beltrami_two_region_residual, &
        beltrami_shell_parity_t, compare_beltrami_shell_residual, &
        validate_beltrami_parity, validate_beltrami_resonance, &
        validate_beltrami_shell_parity
    use fortfem_eulerian_nonnested_residual, only: &
        assemble_eulerian_nonnested_residual, &
        assemble_eulerian_nonnested_residual_jvp
    use fortfem_force_balance_objective, only: &
        evaluate_force_balance_objective, &
        evaluate_force_balance_objective_jvp
    use fortfem_shifted_enriched_space, only: evaluate_shifted_enriched_space
    use fortfem_shifted_vector_enriched_space, only: &
        evaluate_shifted_vector_enriched_space
    use fortfem_singular_layer_matching, only: &
        assemble_singular_layer_matching, &
        assemble_singular_layer_matching_jvp
    use fortfem_api_spaces, only: &
        constant, function, function_space, vector_function_space, dirichlet_bc, &
        test_function, trial_function, vector_test_function, &
        vector_trial_function, test_function_t, trial_function_t, &
        vector_test_function_t, vector_trial_function_t
    use fortfem_api_types, only: function_t
    use fortfem_api_forms, only: &
        compile_tetra_mixed_form_csc, div, dx, form_expr_t, init_measures, &
        grad, inner, operator(*), operator(==)
    use fortfem_mixed_poisson_2d, only: solve_symbolic_mixed_poisson_rt
    use fortfem_tetra_mixed_poisson_3d, only: &
        solve_symbolic_tetra_mixed_poisson_rt
    use fortfem_tetra_discontinuous_arbitrary_order, only: &
        evaluate_tetra_discontinuous, initialize_tetra_discontinuous, &
        tetra_discontinuous_t
    use fortfem_triangle_discontinuous_dof_map, only: &
        build_triangle_discontinuous_dof_map
    use fortfem_triangle_global_dof_map, only: build_triangle_trimmed_dof_map
    use fortfem_triangle_lagrange_arbitrary_order, only: &
        evaluate_triangle_lagrange_basis, initialize_triangle_lagrange_basis, &
        triangle_lagrange_basis_t
    use fortfem_triangle_rt_arbitrary_order, only: &
        initialize_triangle_raviart_thomas, triangle_rt_basis_t
    use fortfem_triangle_vector_interpolation, only: &
        evaluate_triangle_rt_interpolant
    use fortfem_triangle_duffy_quadrature, only: triangle_duffy_quadrature
    use fortfem_api_solvers, only: &
        assemble_laplacian_system, sparse_from_dense, sparse_matrix_t, &
        sparse_direct_factor_t, sparse_direct_factor_csc, &
        sparse_direct_factor_transpose_csc, sparse_direct_free, &
        sparse_direct_solve_factored, sparse_direct_solve_factored_jvp, &
        sparse_direct_solve_factored_vjp, solve
    use fortfem_sparse_ilut, only: &
        build_sparse_ilut_row, sparse_ilut_factor_t, apply_sparse_ilut
    use fortfem_sparse_incomplete_cholesky, only: &
        build_sparse_ichol_row, sparse_incomplete_cholesky_factor_t, &
        apply_sparse_incomplete_cholesky
    use fortfem_assembly_tetra_nedelec_3d, only: &
        assemble_tetra_nedelec_pml_csc, &
        assemble_tetra_nedelec_pml_csc_jvp, &
        assemble_tetra_nedelec_pml_csc_vjp
    use fortfem_helmholtz_representation_3d, only: &
        evaluate_helmholtz_representation_torus_curved_3d
    use fortfem_laplace_representation_3d, only: &
        evaluate_laplace_representation_torus_curved_3d
    implicit none
    private

    public :: apply_tree_cotree_prolongation
    public :: apply_tree_cotree_restriction
    public :: assemble_feec_commuting_projection
    public :: assemble_feec_commuting_projection_jvp
    public :: assemble_feec_commuting_projection_vjp
    public :: assemble_feec_exact_sequence
    public :: assemble_feec_exact_sequence_jvp
    public :: assemble_feec_exact_sequence_vjp
    public :: build_tree_cotree_dof_map
    public :: build_tree_cotree_gauge
    public :: build_tetra_nedelec_dof_map
    public :: build_tetra_edge_dof_map
    public :: build_bspline_polar_feec_2d_operators
    public :: build_bspline_polar_feec_2d_extractions
    public :: build_bspline_polar_h1_extraction
    public :: diagnose_tree_cotree_iga_invariance
    public :: evaluate_bspline_basis
    public :: evaluate_field_aligned_constitutive_tensor
    public :: evaluate_field_aligned_constitutive_tensor_jvp
    public :: evaluate_field_aligned_constitutive_tensor_vjp
    public :: evaluate_periodic_bspline_basis
    public :: evaluate_tetra_nedelec_first_kind
    public :: evaluate_tetra_nedelec_first_kind_jvp
    public :: evaluate_tetra_nedelec_first_kind_vjp
    public :: evaluate_tetra_nedelec_interpolant_at_point
    public :: initialize_tetra_nedelec_first_kind
    public :: interpolate_reference_tetra_nedelec
    public :: map_tetra_nedelec_covariant
    public :: reduce_tree_cotree_dense_system
    public :: reduce_tree_cotree_dense_system_jvp
    public :: reduce_tree_cotree_dense_system_vjp
    public :: tree_cotree_gauge_components
    public :: tree_cotree_gauge_edges
    public :: tree_cotree_gauge_t
    public :: tree_cotree_iga_parity_t
    public :: tetra_nedelec_dof_count
    public :: tetra_nedelec_first_kind_t
    public :: solve_tetra_nedelec_curl_mass
    public :: solve_tetra_nedelec_pml
    public :: tetra_duffy_quadrature
    public :: assemble_bspline_polar_h1_operator_csc
    public :: assemble_bspline_polar_hcurl_operator_csc
    public :: assemble_bspline_polar_l2_mass_csc
    public :: restrict_bspline_polar_operator_csc
    public :: assemble_tensor_diffusion_matrix_3d
    public :: assemble_tensor_diffusion_matrix_3d_jvp
    public :: assemble_tensor_diffusion_matrix_3d_vjp
    public :: assemble_tensor_diffusion_matrix
    public :: assemble_tensor_diffusion_matrix_jvp
    public :: assemble_tensor_diffusion_matrix_vjp
    public :: evaluate_field_aligned_flux
    public :: evaluate_field_aligned_flux_jvp
    public :: evaluate_field_aligned_flux_vjp
    public :: apply_fci_parallel_diffusion
    public :: apply_fci_parallel_diffusion_field_vjp
    public :: apply_fci_parallel_diffusion_jvp
    public :: apply_fci_parallel_diffusion_vjp
    public :: assemble_bspline_h1_operator_csc
    public :: assemble_bspline_h1_operator_csc_jvp
    public :: assemble_bspline_h1_operator_csc_vjp
    public :: apply_fci_parallel_gradient
    public :: apply_fci_parallel_gradient_jvp
    public :: apply_fci_parallel_gradient_vjp
    public :: trace_fci_field_line_rk4
    public :: compute_fci_quadrilateral_cell_areas_2d
    public :: compute_fci_curved_quadrilateral_cell_areas_2d
    public :: compute_fci_polygon_cell_areas_2d
    public :: compute_fci_curved_polygon_cell_areas_2d
    public :: compute_fci_cubic_curved_polygon_cell_areas_2d
    public :: compute_fci_quartic_curved_polygon_cell_areas_2d
    public :: compute_fci_quintic_curved_polygon_cell_areas_2d
    public :: compute_fci_sextic_curved_polygon_cell_areas_2d
    public :: evaluate_cgl_pressure_tensor
    public :: evaluate_cgl_pressure_tensor_jvp
    public :: evaluate_cgl_pressure_tensor_vjp
    public :: evaluate_cgl_pressure_traction
    public :: evaluate_cgl_pressure_traction_jvp
    public :: evaluate_cgl_pressure_traction_vjp
    public :: evaluate_cgl_pressure_work
    public :: evaluate_cgl_pressure_work_jvp
    public :: evaluate_cgl_pressure_work_vjp
    public :: evaluate_cgl_pressure_divergence
    public :: evaluate_cgl_pressure_divergence_jvp
    public :: evaluate_cgl_pressure_divergence_vjp
    public :: assemble_mixed_elasticity_residual
    public :: assemble_mixed_elasticity_residual_jvp
    public :: assemble_mixed_elasticity_residual_vjp
    public :: assemble_elasticity_symmetry_constraint
    public :: assemble_elasticity_symmetry_constraint_jvp
    public :: assemble_elasticity_symmetry_constraint_vjp
    public :: validate_tree_cotree_gauge
    public :: beltrami_parity_t
    public :: compare_beltrami_two_region_residual
    public :: validate_beltrami_parity
    public :: validate_beltrami_resonance
    public :: beltrami_shell_parity_t
    public :: compare_beltrami_shell_residual
    public :: validate_beltrami_shell_parity
    public :: assemble_eulerian_nonnested_residual
    public :: assemble_eulerian_nonnested_residual_jvp
    public :: evaluate_force_balance_objective
    public :: evaluate_force_balance_objective_jvp
    public :: evaluate_shifted_enriched_space
    public :: evaluate_shifted_vector_enriched_space
    public :: assemble_singular_layer_matching
    public :: assemble_singular_layer_matching_jvp
    public :: function_space
    public :: constant
    public :: function
    public :: vector_function_space
    public :: dirichlet_bc
    public :: function_t
    public :: test_function
    public :: trial_function
    public :: vector_test_function
    public :: vector_trial_function
    public :: test_function_t
    public :: trial_function_t
    public :: vector_test_function_t
    public :: vector_trial_function_t
    public :: compile_tetra_mixed_form_csc
    public :: div
    public :: dx
    public :: form_expr_t
    public :: init_measures
    public :: inner
    public :: operator(*)
    public :: operator(==)
    public :: grad
    public :: solve
    public :: solve_symbolic_mixed_poisson_rt
    public :: solve_symbolic_tetra_mixed_poisson_rt
    public :: evaluate_tetra_discontinuous
    public :: initialize_tetra_discontinuous
    public :: tetra_discontinuous_t
    public :: build_triangle_discontinuous_dof_map
    public :: build_triangle_trimmed_dof_map
    public :: evaluate_triangle_lagrange_basis
    public :: initialize_triangle_lagrange_basis
    public :: triangle_lagrange_basis_t
    public :: initialize_triangle_raviart_thomas
    public :: triangle_rt_basis_t
    public :: evaluate_triangle_rt_interpolant
    public :: triangle_duffy_quadrature
    public :: assemble_laplacian_system
    public :: sparse_from_dense
    public :: sparse_matrix_t
    public :: sparse_direct_factor_t
    public :: sparse_direct_factor_csc
    public :: sparse_direct_factor_transpose_csc
    public :: sparse_direct_free
    public :: sparse_direct_solve_factored
    public :: sparse_direct_solve_factored_jvp
    public :: sparse_direct_solve_factored_vjp
    public :: build_sparse_ilut_row
    public :: sparse_ilut_factor_t
    public :: apply_sparse_ilut
    public :: build_sparse_ichol_row
    public :: sparse_incomplete_cholesky_factor_t
    public :: apply_sparse_incomplete_cholesky
    public :: assemble_tetra_nedelec_pml_csc
    public :: assemble_tetra_nedelec_pml_csc_jvp
    public :: assemble_tetra_nedelec_pml_csc_vjp
    public :: evaluate_helmholtz_representation_torus_curved_3d
    public :: evaluate_laplace_representation_torus_curved_3d

end module fortfem_feec
