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
    use fortfem_broken_feec_sequence, only: &
        assemble_broken_feec_sequence, &
        assemble_broken_feec_sequence_jvp, &
        assemble_broken_feec_sequence_vjp
    use fortfem_glued_feec_sequence, only: &
        assemble_glued_feec_sequence, &
        assemble_glued_feec_sequence_jvp, &
        assemble_glued_feec_sequence_vjp
    use fortfem_broken_skeleton_spaces, only: &
        BROKEN_SPACE_H1, BROKEN_SPACE_HCURL, BROKEN_SPACE_HDIV, &
        BROKEN_SPACE_L2, SKELETON_SPACE_SCALAR, &
        SKELETON_SPACE_TANGENTIAL, SKELETON_SPACE_NORMAL, &
        broken_space_layout_t, skeleton_space_layout_t, &
        initialize_broken_space_layout, validate_broken_space_layout, &
        broken_space_layout_maps, broken_space_layout_global_count, &
        initialize_skeleton_space_layout, validate_skeleton_space_layout, &
        skeleton_space_layout_maps, skeleton_space_layout_global_count
    use fortfem_hdg_static_condensation, only: &
        assemble_hdg_static_condensation, &
        assemble_hdg_static_condensation_jvp, &
        assemble_hdg_static_condensation_vjp
    use fortfem_hdg_global_skeleton, only: &
        assemble_hdg_global_skeleton, &
        assemble_hdg_global_skeleton_jvp, &
        assemble_hdg_global_skeleton_vjp
    use fortfem_hdg_global_skeleton_csc, only: &
        assemble_hdg_global_skeleton_csc, &
        assemble_hdg_global_skeleton_csc_jvp, &
        assemble_hdg_global_skeleton_csc_vjp
    use fortfem_block_graph_residual, only: &
        assemble_block_graph_residual, &
        assemble_block_graph_residual_jvp, &
        assemble_block_graph_residual_vjp
    use fortfem_complex_block_graph_residual, only: &
        assemble_complex_block_graph_residual, &
        assemble_complex_block_graph_residual_jvp, &
        assemble_complex_block_graph_residual_vjp
    use fortfem_complex_coupled_field_residual, only: &
        assemble_complex_coupled_field_residual, &
        assemble_complex_coupled_field_residual_jvp, &
        assemble_complex_coupled_field_residual_vjp
    use fortfem_coupled_field_residual, only: &
        assemble_coupled_field_residual, &
        assemble_coupled_field_residual_jvp, &
        assemble_coupled_field_residual_vjp
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
    use fortfem_barycentric_surface_refinement, only: &
        barycentric_refine_surface_mesh
    use fortfem_maxwell_bc_surface, only: &
        assemble_maxwell_rwg_rbc_pairing, build_maxwell_bc_transformation, &
        differentiate_maxwell_bc_transformation_jvp, &
        differentiate_maxwell_bc_transformation_vjp
    use fortfem_maxwell_rwg_surface, only: &
        assemble_maxwell_rwg_mass_matrix, build_maxwell_rwg_surface_space, &
        evaluate_maxwell_rwg_basis, map_maxwell_rwg_to_tetra_nedelec_edges
    use fortfem_maxwell_localized_rwg_surface, only: &
        evaluate_maxwell_localized_rwg_basis
    use fortfem_maxwell_surface_rt, only: &
        assemble_maxwell_surface_rt_mass_matrix, &
        build_maxwell_surface_rt_dof_map, evaluate_maxwell_surface_rt_basis, &
        evaluate_maxwell_surface_rt_global_basis
    use fortfem_maxwell_surface_rt_efie_3d, only: &
        assemble_maxwell_surface_rt_efie_3d
    use fortfem_maxwell_sphere_curved_rwg, only: &
        evaluate_maxwell_sphere_curved_localized_rwg_basis, &
        evaluate_maxwell_sphere_curved_localized_rwg_basis_jvp, &
        evaluate_maxwell_sphere_curved_localized_rwg_basis_vjp, &
        evaluate_maxwell_sphere_curved_rwg_basis, &
        evaluate_maxwell_sphere_curved_rwg_basis_jvp, &
        evaluate_maxwell_sphere_curved_rwg_basis_vjp
    use fortfem_maxwell_torus_curved_rwg, only: &
        evaluate_maxwell_torus_curved_localized_rwg_basis, &
        evaluate_maxwell_torus_curved_localized_rwg_basis_jvp, &
        evaluate_maxwell_torus_curved_localized_rwg_basis_vjp, &
        evaluate_maxwell_torus_curved_rwg_basis, &
        evaluate_maxwell_torus_curved_rwg_basis_jvp, &
        evaluate_maxwell_torus_curved_rwg_basis_vjp
    use fortfem_maxwell_efie_rwg_3d, only: assemble_maxwell_efie_rwg_3d
    use fortfem_tetra_duffy_quadrature, only: tetra_duffy_quadrature
    use fortfem_tetra_rt_arbitrary_order, only: &
        initialize_tetra_rt, tetra_rt_t
    use fortfem_tetra_rt_global_dof_map, only: &
        build_tetra_rt_basis_transform, build_tetra_rt_dof_map
    use fortfem_tetra_rt_interpolation, only: &
        interpolate_physical_tetra_rt, &
        interpolate_sampled_physical_tetra_rt, &
        interpolate_sampled_physical_tetra_rt_jvp, &
        interpolate_sampled_physical_tetra_rt_vjp, &
        tetra_rt_interpolation_points
    use fortfem_tetra_vector_samples, only: &
        tetra_vector_sample_gradients_t, tetra_vector_samples_t, &
        zero_tetra_vector_samples_like
    use fortfem_tetra_piola_maps, only: &
        map_tetra_rt_contravariant, &
        map_tetra_rt_contravariant_jvp, &
        map_tetra_rt_contravariant_vjp
    use fortfem_assembly_tetra_rt_arbitrary_order_3d, only: &
        assemble_tetra_rt_div_mass_csc, &
        assemble_tetra_rt_div_mass_csc_jvp, &
        assemble_tetra_rt_div_mass_csc_vjp, &
        assemble_tetra_rt_div_mass_element, &
        assemble_tetra_rt_div_mass_element_jvp, &
        assemble_tetra_rt_div_mass_element_vjp, &
        assemble_tetra_rt_vector_load_samples, &
        assemble_tetra_rt_vector_load_samples_jvp, &
        assemble_tetra_rt_vector_load_samples_vjp
    use fortfem_tetra_nedelec_global_dof_map, only: build_tetra_nedelec_dof_map
    use fortfem_tetra_piola_maps, only: map_tetra_nedelec_covariant
    use fortfem_assembly_bspline_polar_2d, only: &
        assemble_bspline_polar_h1_operator_csc, &
        assemble_bspline_polar_hcurl_operator_csc, &
        assemble_bspline_polar_l2_mass_csc, restrict_bspline_polar_operator_csc
    use fortfem_assembly_bspline_2d, only: &
        assemble_bspline_h1_operator_csc, &
        assemble_bspline_h1_operator_csc_jvp, &
        assemble_bspline_h1_operator_csc_vjp, &
        assemble_bspline_h1_weighted_mass_csc, &
        assemble_bspline_hcurl_operator_csc, &
        assemble_bspline_hdiv_operator_csc, &
        assemble_bspline_h1_hcurl_gradient_csc, &
        assemble_bspline_hcurl_h1_adjoint_gradient_csc, &
        assemble_bspline_hcurl_l2_curl_csc, &
        assemble_bspline_l2_hcurl_adjoint_curl_csc, &
        assemble_bspline_l2_mass_csc, &
        assemble_bspline_grad_shafranov_csc, &
        assemble_bspline_toroidal_fourier_laplacian_csc, &
        assemble_bspline_poloidal_bracket_csc, &
        apply_bspline_toroidal_poloidal_bracket, &
        apply_bspline_toroidal_poloidal_bracket_jvp, &
        apply_bspline_jorek_flux_rhs, &
        apply_bspline_jorek_flux_jvp, &
        apply_bspline_jorek_thermodynamic_rhs, &
        apply_bspline_jorek_thermodynamic_jvp, &
        apply_bspline_jorek_density_rhs, &
        apply_bspline_jorek_density_jvp, &
        project_bspline_toroidal_product, &
        build_bspline_feec_2d_operators_csc
    use fortfem_assembly_bspline_3d, only: &
        build_bspline_feec_3d_operators_csc, &
        assemble_bspline_h1_operator_3d_csc, &
        assemble_bspline_hcurl_operator_3d_csc, &
        assemble_bspline_hdiv_operator_3d_csc, &
        assemble_bspline_l2_mass_3d_csc, &
        assemble_bspline_hdiv_l2_divergence_3d_csc, &
        assemble_bspline_l2_hdiv_adjoint_divergence_3d_csc
    use fortfem_assembly_bspline_multipatch_2d, only: &
        build_bspline_feec_2d_two_patch_operators_csc
    use fortfem_assembly_bspline_multipatch_3d, only: &
        build_bspline_feec_3d_two_patch_operators_csc
    use fortfem_bspline_feec, only: &
        build_bspline_derivative_matrix, build_bspline_feec_2d_operators, &
        build_bspline_feec_3d_operators, evaluate_bspline_basis, &
        evaluate_nurbs_surface_geometry, evaluate_nurbs_surface_geometry_jvp, &
        evaluate_nurbs_surface_geometry_vjp, evaluate_nurbs_volume_geometry, &
        evaluate_nurbs_volume_geometry_jvp, evaluate_nurbs_volume_geometry_vjp, &
        map_isogeometric_h1_gradient, map_isogeometric_hcurl, &
        map_isogeometric_hdiv, map_isogeometric_l2
    use fortfem_bspline_multipatch, only: &
        BSPLINE_FACE_X_MIN, BSPLINE_FACE_X_MAX, BSPLINE_FACE_Y_MIN, &
        BSPLINE_FACE_Y_MAX, BSPLINE_FACE_Z_MIN, BSPLINE_FACE_Z_MAX, &
        build_bspline_feec_2d_interface_dofs, &
        build_bspline_feec_2d_multipatch_maps, &
        build_bspline_feec_3d_interface_dofs, &
        build_bspline_feec_3d_multipatch_maps
    use fortfem_bspline_polar, only: &
        build_bspline_polar_feec_2d_operators, &
        build_bspline_polar_feec_2d_extractions, &
        build_bspline_polar_h1_extraction, evaluate_periodic_bspline_basis
    use fortfem_field_aligned_constitutive_tensor, only: &
        evaluate_field_aligned_constitutive_tensor, &
        evaluate_field_aligned_constitutive_tensor_jvp, &
        evaluate_field_aligned_constitutive_tensor_vjp
    use fortfem_tensor_power_split, only: evaluate_tensor_power_split, &
        evaluate_tensor_power_split_jvp, evaluate_tensor_power_split_vjp
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
        assemble_fci_parallel_gradient_csc, &
        assemble_fci_parallel_support_divergence_csc, &
        apply_fci_parallel_diffusion, apply_fci_parallel_diffusion_jvp, &
        apply_fci_parallel_diffusion_vjp, apply_fci_parallel_diffusion_field_vjp, &
        apply_fci_parallel_gradient, apply_fci_parallel_gradient_jvp, &
        apply_fci_parallel_gradient_vjp, compute_fci_parallel_diffusion_diagonal
    use fortfem_fci_power_flux_ledger, only: &
        evaluate_fci_power_flux_ledger, &
        evaluate_fci_power_flux_ledger_jvp, &
        evaluate_fci_power_flux_ledger_vjp
    use fortfem_fci_field_split_preconditioner, only: &
        apply_fci_additive_field_split_preconditioner
    use fortfem_tetra_nedelec_solver_3d, only: &
        solve_tetra_nedelec_curl_mass, solve_tetra_nedelec_pml
    use fortfem_fci_field_line_tracer, only: &
        trace_fci_field_line_rk4, trace_fci_field_line_rk4_jvp
    use fortfem_fci_interpolation_map, only: &
        build_fci_bilinear_interpolation_map_2d, &
        build_fci_bilinear_interpolation_map_2d_jvp, &
        build_fci_bilinear_interpolation_map_2d_vjp, &
        build_fci_cubic_interpolation_map_1d, &
        build_fci_cubic_interpolation_map_1d_jvp, &
        build_fci_cubic_interpolation_map_1d_vjp
    use fortfem_fci_support_geometry, only: &
        compute_fci_quadrilateral_cell_areas_2d, &
        compute_fci_curved_quadrilateral_cell_areas_2d, &
        compute_fci_polygon_cell_areas_2d, &
        compute_fci_curved_polygon_cell_areas_2d, &
        compute_fci_curved_polygon_cell_areas_2d_jvp, &
        compute_fci_curved_polygon_cell_areas_2d_vjp, &
        compute_fci_cubic_curved_polygon_cell_areas_2d, &
        compute_fci_cubic_curved_polygon_cell_areas_2d_jvp, &
        compute_fci_cubic_curved_polygon_cell_areas_2d_vjp, &
        compute_fci_curved_quadrilateral_cell_areas_2d_jvp, &
        compute_fci_curved_quadrilateral_cell_areas_2d_vjp, &
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
    use fortfem_geometry_mortar_component_coupling, only: &
        assemble_geometry_mortar_component_coupling, &
        assemble_geometry_mortar_component_coupling_jvp, &
        assemble_geometry_mortar_component_coupling_vjp
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
        evaluate_force_balance_objective_jvp, &
        evaluate_force_balance_objective_vjp, &
        evaluate_force_balance_objective_hvp
    use fortfem_force_balance_residual, only: &
        assemble_force_balance_residual, &
        assemble_force_balance_residual_jvp, &
        assemble_force_balance_residual_vjp
    use fortfem_force_balance_product, only: &
        evaluate_force_balance_product, &
        evaluate_force_balance_product_jvp, &
        evaluate_force_balance_product_vjp
    use fortfem_heaviside_enrichment, only: &
        evaluate_shifted_heaviside_enrichment, &
        evaluate_shifted_heaviside_enrichment_jvp, &
        evaluate_shifted_heaviside_enrichment_vjp
    use fortfem_shifted_enriched_basis, only: &
        evaluate_shifted_enriched_basis, &
        evaluate_shifted_enriched_basis_jvp, &
        evaluate_shifted_enriched_basis_vjp
    use fortfem_shifted_enriched_space, only: &
        evaluate_shifted_enriched_space, &
        evaluate_shifted_enriched_space_jvp, &
        evaluate_shifted_enriched_space_vjp
    use fortfem_shifted_vector_enriched_basis, only: &
        evaluate_shifted_vector_enriched_basis, &
        evaluate_shifted_vector_enriched_basis_jvp, &
        evaluate_shifted_vector_enriched_basis_vjp
    use fortfem_shifted_vector_enriched_space, only: &
        evaluate_shifted_vector_enriched_space, &
        evaluate_shifted_vector_enriched_space_jvp, &
        evaluate_shifted_vector_enriched_space_vjp
    use fortfem_enriched_feec_sequence, only: &
        assemble_enriched_feec_sequence, &
        assemble_enriched_feec_sequence_jvp, &
        assemble_enriched_feec_sequence_vjp
    use fortfem_singular_layer_matching, only: &
        assemble_singular_layer_matching, &
        assemble_singular_layer_matching_jvp
    use fortfem_retained_field_split, only: &
        apply_retained_complex_field_split, &
        apply_retained_complex_field_split_jvp, &
        apply_retained_complex_field_split_vjp, &
        apply_retained_field_split, apply_retained_field_split_jvp, &
        apply_retained_field_split_vjp, factor_retained_complex_field_split, &
        factor_retained_field_split, free_retained_complex_field_split, &
        free_retained_field_split, retained_complex_field_split_t, &
        retained_field_split_t
    use fortfem_retained_coupled_schur, only: &
        assemble_retained_coupled_schur, &
        assemble_retained_coupled_schur_jvp, &
        assemble_retained_coupled_schur_vjp
    use fortfem_api_spaces, only: &
        constant, function, function_space, vector_function_space, dirichlet_bc, &
        test_function, trial_function, vector_function, vector_bc, &
        vector_test_function, &
        vector_trial_function, test_function_t, trial_function_t, &
        vector_test_function_t, vector_trial_function_t, vector_function_t, &
        vector_bc_t
    use fortfem_api_types, only: function_t
    use fortfem_api_forms, only: &
        compile_tetra_mixed_form_csc, curl, div, dx, form_expr_t, init_measures, &
        grad, inner, operator(*), operator(+), operator(==)
    use fortfem_mixed_poisson_2d, only: &
        solve_mixed_poisson_rt, solve_mixed_poisson_rt0, &
        solve_symbolic_mixed_poisson_rt
    use fortfem_mixed_rt_system, only: solve_mixed_rt_system, &
        solve_mixed_rt_system_jvp, solve_mixed_rt_system_vjp
    use fortfem_tetra_mixed_poisson_3d, only: &
        solve_symbolic_tetra_mixed_poisson_rt
    use fortfem_tetra_mixed_poisson_state_3d, only: &
        solve_tetra_mixed_poisson_state, &
        solve_tetra_mixed_poisson_state_jvp, &
        solve_tetra_mixed_poisson_state_vjp, &
        solve_tetra_mixed_poisson_sampled_state, &
        solve_tetra_mixed_poisson_sampled_state_jvp, &
        solve_tetra_mixed_poisson_sampled_state_vjp
    use fortfem_tetra_discontinuous_arbitrary_order, only: &
        evaluate_tetra_discontinuous, initialize_tetra_discontinuous, &
        tetra_discontinuous_t
    use fortfem_triangle_discontinuous_dof_map, only: &
        build_triangle_discontinuous_dof_map
    use fortfem_triangle_global_dof_map, only: build_triangle_trimmed_dof_map
    use fortfem_triangle_lagrange_arbitrary_order, only: &
        assignment(=), evaluate_triangle_lagrange_basis, &
        initialize_triangle_lagrange_basis, triangle_lagrange_basis_t, &
        triangle_lagrange_nodes
    use fortfem_triangle_nedelec_arbitrary_order, only: &
        assignment(=), evaluate_triangle_nedelec_first_kind, &
        evaluate_triangle_nedelec_first_kind_jvp, &
        evaluate_triangle_nedelec_first_kind_vjp, &
        initialize_triangle_nedelec_first_kind, triangle_nedelec_dof_count, &
        triangle_nedelec_first_kind_t
    use fortfem_triangle_feec_operators, only: &
        build_triangle_discrete_gradient
    use fortfem_triangle_rt_arbitrary_order, only: &
        evaluate_triangle_raviart_thomas, initialize_triangle_raviart_thomas, &
        triangle_rt_basis_t, triangle_rt_dof_count
    use fortfem_triangle_vector_interpolation, only: &
        evaluate_triangle_rt_interpolant, interpolate_triangle_nedelec, &
        interpolate_triangle_rt
    use fortfem_assembly_nedelec_arbitrary_order_2d, only: &
        assemble_triangle_nedelec_curl_mass_csc, &
        assemble_triangle_nedelec_curl_mass_element
    use fortfem_assembly_rt_arbitrary_order_2d, only: &
        assemble_triangle_rt_div_mass_element, &
        assemble_triangle_rt_div_mass_element_jvp, &
        assemble_triangle_rt_div_mass_element_vjp, &
        assemble_triangle_rt_div_mass_csc, &
        assemble_triangle_rt_div_mass_csc_jvp, &
        assemble_triangle_rt_div_mass_csc_vjp, &
        assemble_triangle_rt_divergence_csc, &
        assemble_triangle_rt_vector_load_samples, &
        assemble_triangle_rt_vector_load_samples_jvp, &
        assemble_triangle_rt_vector_load_samples_vjp
    use fortfem_triangle_rt_sampled_state_2d, only: &
        solve_triangle_rt_sampled_state, &
        solve_triangle_rt_sampled_state_jvp, &
        solve_triangle_rt_sampled_state_vjp
    use fortfem_triangle_duffy_quadrature, only: triangle_duffy_quadrature
    use fortfem_api_solvers, only: &
        assemble_laplacian_system, sparse_from_dense, sparse_matrix_t, &
        sparse_direct_factor_t, sparse_direct_factor_csc, &
        sparse_direct_factor_transpose_csc, sparse_direct_free, &
        sparse_direct_solve_factored, sparse_direct_solve_factored_jvp, &
        sparse_direct_solve_factored_vjp, solve
    use fortfem_sparse_direct, only: &
        sparse_direct_solve_tree_cotree, &
        sparse_direct_solve_tree_cotree_jvp, &
        sparse_direct_solve_tree_cotree_vjp
    use fortfem_sparse_ilut, only: &
        build_sparse_ilut, build_sparse_ilut_row, sparse_ilut_factor_t, &
        apply_sparse_ilut, apply_sparse_ilut_jvp, apply_sparse_ilut_vjp
    use fortfem_sparse_incomplete_cholesky, only: &
        build_sparse_ichol, build_sparse_ichol_row, &
        sparse_incomplete_cholesky_factor_t, apply_sparse_incomplete_cholesky, &
        apply_sparse_incomplete_cholesky_jvp, apply_sparse_incomplete_cholesky_vjp
    use fortfem_assembly_tetra_nedelec_3d, only: &
        assemble_tetra_nedelec_curl_mass_element, &
        assemble_tetra_nedelec_pml_csc, &
        assemble_tetra_nedelec_pml_csc_jvp, &
        assemble_tetra_nedelec_pml_csc_vjp
    implicit none
    private

    public :: apply_tree_cotree_prolongation
    public :: apply_tree_cotree_restriction
    public :: assemble_block_graph_residual
    public :: assemble_block_graph_residual_jvp
    public :: assemble_block_graph_residual_vjp
    public :: assemble_complex_block_graph_residual
    public :: assemble_complex_block_graph_residual_jvp
    public :: assemble_complex_block_graph_residual_vjp
    public :: assemble_complex_coupled_field_residual
    public :: assemble_complex_coupled_field_residual_jvp
    public :: assemble_complex_coupled_field_residual_vjp
    public :: assemble_coupled_field_residual
    public :: assemble_coupled_field_residual_jvp
    public :: assemble_coupled_field_residual_vjp
    public :: assemble_feec_commuting_projection
    public :: assemble_feec_commuting_projection_jvp
    public :: assemble_feec_commuting_projection_vjp
    public :: assemble_feec_exact_sequence
    public :: assemble_feec_exact_sequence_jvp
    public :: assemble_feec_exact_sequence_vjp
    public :: assemble_broken_feec_sequence
    public :: assemble_broken_feec_sequence_jvp
    public :: assemble_broken_feec_sequence_vjp
    public :: assemble_glued_feec_sequence
    public :: assemble_glued_feec_sequence_jvp
    public :: assemble_glued_feec_sequence_vjp
    public :: BROKEN_SPACE_H1
    public :: BROKEN_SPACE_HCURL
    public :: BROKEN_SPACE_HDIV
    public :: BROKEN_SPACE_L2
    public :: SKELETON_SPACE_SCALAR
    public :: SKELETON_SPACE_TANGENTIAL
    public :: SKELETON_SPACE_NORMAL
    public :: broken_space_layout_t
    public :: skeleton_space_layout_t
    public :: initialize_broken_space_layout
    public :: validate_broken_space_layout
    public :: broken_space_layout_maps
    public :: broken_space_layout_global_count
    public :: initialize_skeleton_space_layout
    public :: validate_skeleton_space_layout
    public :: skeleton_space_layout_maps
    public :: skeleton_space_layout_global_count
    public :: assemble_hdg_static_condensation
    public :: assemble_hdg_static_condensation_jvp
    public :: assemble_hdg_static_condensation_vjp
    public :: assemble_hdg_global_skeleton
    public :: assemble_hdg_global_skeleton_jvp
    public :: assemble_hdg_global_skeleton_vjp
    public :: assemble_hdg_global_skeleton_csc
    public :: assemble_hdg_global_skeleton_csc_jvp
    public :: assemble_hdg_global_skeleton_csc_vjp
    public :: build_tree_cotree_dof_map
    public :: build_tree_cotree_gauge
    public :: build_tetra_nedelec_dof_map
    public :: build_tetra_edge_dof_map
    public :: barycentric_refine_surface_mesh
    public :: assemble_maxwell_rwg_rbc_pairing
    public :: build_maxwell_bc_transformation
    public :: differentiate_maxwell_bc_transformation_jvp
    public :: differentiate_maxwell_bc_transformation_vjp
    public :: assemble_maxwell_rwg_mass_matrix
    public :: build_maxwell_rwg_surface_space
    public :: evaluate_maxwell_rwg_basis
    public :: map_maxwell_rwg_to_tetra_nedelec_edges
    public :: evaluate_maxwell_localized_rwg_basis
    public :: assemble_maxwell_surface_rt_mass_matrix
    public :: build_maxwell_surface_rt_dof_map
    public :: evaluate_maxwell_surface_rt_basis
    public :: evaluate_maxwell_surface_rt_global_basis
    public :: assemble_maxwell_surface_rt_efie_3d
    public :: evaluate_maxwell_sphere_curved_localized_rwg_basis
    public :: evaluate_maxwell_sphere_curved_localized_rwg_basis_jvp
    public :: evaluate_maxwell_sphere_curved_localized_rwg_basis_vjp
    public :: evaluate_maxwell_sphere_curved_rwg_basis
    public :: evaluate_maxwell_sphere_curved_rwg_basis_jvp
    public :: evaluate_maxwell_sphere_curved_rwg_basis_vjp
    public :: evaluate_maxwell_torus_curved_localized_rwg_basis
    public :: evaluate_maxwell_torus_curved_localized_rwg_basis_jvp
    public :: evaluate_maxwell_torus_curved_localized_rwg_basis_vjp
    public :: evaluate_maxwell_torus_curved_rwg_basis
    public :: evaluate_maxwell_torus_curved_rwg_basis_jvp
    public :: evaluate_maxwell_torus_curved_rwg_basis_vjp
    public :: assemble_maxwell_efie_rwg_3d
    public :: build_bspline_polar_feec_2d_operators
    public :: build_bspline_polar_feec_2d_extractions
    public :: build_bspline_polar_h1_extraction
    public :: diagnose_tree_cotree_iga_invariance
    public :: assemble_geometry_mortar_component_coupling
    public :: assemble_geometry_mortar_component_coupling_jvp
    public :: assemble_geometry_mortar_component_coupling_vjp
    public :: evaluate_bspline_basis
    public :: evaluate_field_aligned_constitutive_tensor
    public :: evaluate_field_aligned_constitutive_tensor_jvp
    public :: evaluate_field_aligned_constitutive_tensor_vjp
    public :: evaluate_tensor_power_split
    public :: evaluate_tensor_power_split_jvp
    public :: evaluate_tensor_power_split_vjp
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
    public :: initialize_tetra_rt
    public :: tetra_rt_t
    public :: build_tetra_rt_basis_transform
    public :: build_tetra_rt_dof_map
    public :: interpolate_physical_tetra_rt
    public :: interpolate_sampled_physical_tetra_rt
    public :: interpolate_sampled_physical_tetra_rt_jvp
    public :: interpolate_sampled_physical_tetra_rt_vjp
    public :: tetra_rt_interpolation_points
    public :: tetra_vector_sample_gradients_t
    public :: tetra_vector_samples_t
    public :: zero_tetra_vector_samples_like
    public :: map_tetra_rt_contravariant
    public :: map_tetra_rt_contravariant_jvp
    public :: map_tetra_rt_contravariant_vjp
    public :: assemble_tetra_rt_div_mass_csc
    public :: assemble_tetra_rt_div_mass_csc_jvp
    public :: assemble_tetra_rt_div_mass_csc_vjp
    public :: assemble_tetra_rt_div_mass_element
    public :: assemble_tetra_rt_div_mass_element_jvp
    public :: assemble_tetra_rt_div_mass_element_vjp
    public :: assemble_tetra_rt_vector_load_samples
    public :: assemble_tetra_rt_vector_load_samples_jvp
    public :: assemble_tetra_rt_vector_load_samples_vjp
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
    public :: assemble_fci_parallel_gradient_csc
    public :: assemble_fci_parallel_support_divergence_csc
    public :: assemble_bspline_h1_operator_csc
    public :: assemble_bspline_h1_operator_csc_jvp
    public :: assemble_bspline_h1_operator_csc_vjp
    public :: assemble_bspline_h1_weighted_mass_csc
    public :: assemble_bspline_hcurl_operator_csc
    public :: assemble_bspline_hdiv_operator_csc
    public :: assemble_bspline_h1_hcurl_gradient_csc
    public :: assemble_bspline_hcurl_h1_adjoint_gradient_csc
    public :: assemble_bspline_hcurl_l2_curl_csc
    public :: assemble_bspline_l2_hcurl_adjoint_curl_csc
    public :: assemble_bspline_l2_mass_csc
    public :: assemble_bspline_grad_shafranov_csc
    public :: assemble_bspline_toroidal_fourier_laplacian_csc
    public :: assemble_bspline_poloidal_bracket_csc
    public :: apply_bspline_toroidal_poloidal_bracket
    public :: apply_bspline_toroidal_poloidal_bracket_jvp
    public :: apply_bspline_jorek_flux_rhs
    public :: apply_bspline_jorek_flux_jvp
    public :: apply_bspline_jorek_thermodynamic_rhs
    public :: apply_bspline_jorek_thermodynamic_jvp
    public :: apply_bspline_jorek_density_rhs
    public :: apply_bspline_jorek_density_jvp
    public :: project_bspline_toroidal_product
    public :: build_bspline_feec_2d_operators_csc
    public :: build_bspline_derivative_matrix
    public :: build_bspline_feec_2d_operators
    public :: build_bspline_feec_3d_operators
    public :: build_bspline_feec_2d_interface_dofs
    public :: build_bspline_feec_2d_multipatch_maps
    public :: build_bspline_feec_3d_interface_dofs
    public :: build_bspline_feec_3d_multipatch_maps
    public :: build_bspline_feec_2d_two_patch_operators_csc
    public :: build_bspline_feec_3d_two_patch_operators_csc
    public :: build_bspline_feec_3d_operators_csc
    public :: assemble_bspline_h1_operator_3d_csc
    public :: assemble_bspline_hcurl_operator_3d_csc
    public :: assemble_bspline_hdiv_operator_3d_csc
    public :: assemble_bspline_l2_mass_3d_csc
    public :: assemble_bspline_hdiv_l2_divergence_3d_csc
    public :: assemble_bspline_l2_hdiv_adjoint_divergence_3d_csc
    public :: evaluate_nurbs_surface_geometry
    public :: evaluate_nurbs_surface_geometry_jvp
    public :: evaluate_nurbs_surface_geometry_vjp
    public :: evaluate_nurbs_volume_geometry
    public :: evaluate_nurbs_volume_geometry_jvp
    public :: evaluate_nurbs_volume_geometry_vjp
    public :: map_isogeometric_h1_gradient
    public :: map_isogeometric_hcurl
    public :: map_isogeometric_hdiv
    public :: map_isogeometric_l2
    public :: BSPLINE_FACE_X_MIN
    public :: BSPLINE_FACE_X_MAX
    public :: BSPLINE_FACE_Y_MIN
    public :: BSPLINE_FACE_Y_MAX
    public :: BSPLINE_FACE_Z_MIN
    public :: BSPLINE_FACE_Z_MAX
    public :: apply_fci_parallel_gradient
    public :: apply_fci_parallel_gradient_jvp
    public :: apply_fci_parallel_gradient_vjp
    public :: compute_fci_parallel_diffusion_diagonal
    public :: evaluate_fci_power_flux_ledger
    public :: evaluate_fci_power_flux_ledger_jvp
    public :: evaluate_fci_power_flux_ledger_vjp
    public :: apply_fci_additive_field_split_preconditioner
    public :: trace_fci_field_line_rk4
    public :: trace_fci_field_line_rk4_jvp
    public :: build_fci_bilinear_interpolation_map_2d
    public :: build_fci_bilinear_interpolation_map_2d_jvp
    public :: build_fci_bilinear_interpolation_map_2d_vjp
    public :: build_fci_cubic_interpolation_map_1d
    public :: build_fci_cubic_interpolation_map_1d_jvp
    public :: build_fci_cubic_interpolation_map_1d_vjp
    public :: compute_fci_quadrilateral_cell_areas_2d
    public :: compute_fci_curved_quadrilateral_cell_areas_2d
    public :: compute_fci_curved_quadrilateral_cell_areas_2d_jvp
    public :: compute_fci_curved_quadrilateral_cell_areas_2d_vjp
    public :: compute_fci_polygon_cell_areas_2d
    public :: compute_fci_curved_polygon_cell_areas_2d
    public :: compute_fci_curved_polygon_cell_areas_2d_jvp
    public :: compute_fci_curved_polygon_cell_areas_2d_vjp
    public :: compute_fci_cubic_curved_polygon_cell_areas_2d
    public :: compute_fci_cubic_curved_polygon_cell_areas_2d_jvp
    public :: compute_fci_cubic_curved_polygon_cell_areas_2d_vjp
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
    public :: evaluate_force_balance_objective_vjp
    public :: evaluate_force_balance_objective_hvp
    public :: assemble_force_balance_residual
    public :: assemble_force_balance_residual_jvp
    public :: assemble_force_balance_residual_vjp
    public :: evaluate_force_balance_product
    public :: evaluate_force_balance_product_jvp
    public :: evaluate_force_balance_product_vjp
    public :: evaluate_shifted_heaviside_enrichment
    public :: evaluate_shifted_heaviside_enrichment_jvp
    public :: evaluate_shifted_heaviside_enrichment_vjp
    public :: evaluate_shifted_enriched_basis
    public :: evaluate_shifted_enriched_basis_jvp
    public :: evaluate_shifted_enriched_basis_vjp
    public :: evaluate_shifted_enriched_space
    public :: evaluate_shifted_enriched_space_jvp
    public :: evaluate_shifted_enriched_space_vjp
    public :: evaluate_shifted_vector_enriched_basis
    public :: evaluate_shifted_vector_enriched_basis_jvp
    public :: evaluate_shifted_vector_enriched_basis_vjp
    public :: evaluate_shifted_vector_enriched_space
    public :: evaluate_shifted_vector_enriched_space_jvp
    public :: evaluate_shifted_vector_enriched_space_vjp
    public :: assemble_enriched_feec_sequence
    public :: assemble_enriched_feec_sequence_jvp
    public :: assemble_enriched_feec_sequence_vjp
    public :: assemble_singular_layer_matching
    public :: assemble_singular_layer_matching_jvp
    public :: retained_field_split_t
    public :: retained_complex_field_split_t
    public :: factor_retained_field_split
    public :: factor_retained_complex_field_split
    public :: apply_retained_field_split
    public :: apply_retained_complex_field_split
    public :: apply_retained_field_split_jvp
    public :: apply_retained_complex_field_split_jvp
    public :: apply_retained_field_split_vjp
    public :: apply_retained_complex_field_split_vjp
    public :: free_retained_field_split
    public :: free_retained_complex_field_split
    public :: assemble_retained_coupled_schur
    public :: assemble_retained_coupled_schur_jvp
    public :: assemble_retained_coupled_schur_vjp
    public :: function_space
    public :: constant
    public :: function
    public :: vector_function_space
    public :: vector_function
    public :: vector_function_t
    public :: vector_bc
    public :: vector_bc_t
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
    public :: operator(+)
    public :: operator(==)
    public :: grad
    public :: curl
    public :: solve
    public :: solve_mixed_poisson_rt0
    public :: solve_mixed_poisson_rt
    public :: solve_mixed_rt_system
    public :: solve_mixed_rt_system_jvp
    public :: solve_mixed_rt_system_vjp
    public :: solve_symbolic_mixed_poisson_rt
    public :: solve_symbolic_tetra_mixed_poisson_rt
    public :: solve_tetra_mixed_poisson_state
    public :: solve_tetra_mixed_poisson_state_jvp
    public :: solve_tetra_mixed_poisson_state_vjp
    public :: solve_tetra_mixed_poisson_sampled_state
    public :: solve_tetra_mixed_poisson_sampled_state_jvp
    public :: solve_tetra_mixed_poisson_sampled_state_vjp
    public :: evaluate_tetra_discontinuous
    public :: initialize_tetra_discontinuous
    public :: tetra_discontinuous_t
    public :: build_triangle_discontinuous_dof_map
    public :: build_triangle_trimmed_dof_map
    public :: evaluate_triangle_lagrange_basis
    public :: initialize_triangle_lagrange_basis
    public :: triangle_lagrange_basis_t
    public :: triangle_lagrange_nodes
    public :: assignment(=)
    public :: build_triangle_discrete_gradient
    public :: evaluate_triangle_nedelec_first_kind
    public :: evaluate_triangle_nedelec_first_kind_jvp
    public :: evaluate_triangle_nedelec_first_kind_vjp
    public :: initialize_triangle_nedelec_first_kind
    public :: triangle_nedelec_dof_count
    public :: triangle_nedelec_first_kind_t
    public :: interpolate_triangle_nedelec
    public :: assemble_triangle_nedelec_curl_mass_csc
    public :: assemble_triangle_nedelec_curl_mass_element
    public :: initialize_triangle_raviart_thomas
    public :: evaluate_triangle_raviart_thomas
    public :: triangle_rt_basis_t
    public :: triangle_rt_dof_count
    public :: evaluate_triangle_rt_interpolant
    public :: interpolate_triangle_rt
    public :: assemble_triangle_rt_div_mass_element
    public :: assemble_triangle_rt_div_mass_element_jvp
    public :: assemble_triangle_rt_div_mass_element_vjp
    public :: assemble_triangle_rt_div_mass_csc
    public :: assemble_triangle_rt_div_mass_csc_jvp
    public :: assemble_triangle_rt_div_mass_csc_vjp
    public :: assemble_triangle_rt_divergence_csc
    public :: assemble_triangle_rt_vector_load_samples
    public :: assemble_triangle_rt_vector_load_samples_jvp
    public :: assemble_triangle_rt_vector_load_samples_vjp
    public :: solve_triangle_rt_sampled_state
    public :: solve_triangle_rt_sampled_state_jvp
    public :: solve_triangle_rt_sampled_state_vjp
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
    public :: sparse_direct_solve_tree_cotree
    public :: sparse_direct_solve_tree_cotree_jvp
    public :: sparse_direct_solve_tree_cotree_vjp
    public :: build_sparse_ilut
    public :: build_sparse_ilut_row
    public :: sparse_ilut_factor_t
    public :: apply_sparse_ilut
    public :: apply_sparse_ilut_jvp
    public :: apply_sparse_ilut_vjp
    public :: build_sparse_ichol
    public :: build_sparse_ichol_row
    public :: sparse_incomplete_cholesky_factor_t
    public :: apply_sparse_incomplete_cholesky
    public :: apply_sparse_incomplete_cholesky_jvp
    public :: apply_sparse_incomplete_cholesky_vjp
    public :: assemble_tetra_nedelec_pml_csc
    public :: assemble_tetra_nedelec_curl_mass_element
    public :: assemble_tetra_nedelec_pml_csc_jvp
    public :: assemble_tetra_nedelec_pml_csc_vjp

end module fortfem_feec
