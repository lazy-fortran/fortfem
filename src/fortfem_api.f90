module fortfem_api
    use fortfem_interchange_samples, only: &
        compare_interchange_samples, compare_interchange_samples_jvp, &
        compare_interchange_samples_vjp, initialize_interchange_samples, &
        interchange_sample_set_t, validate_interchange_samples
    use fortfem_complex_interchange_samples, only: &
        compare_complex_interchange_samples, &
        compare_complex_interchange_samples_jvp, &
        compare_complex_interchange_samples_vjp, &
        complex_interchange_sample_set_t, &
        initialize_complex_interchange_samples, &
        validate_complex_interchange_samples
    use fortfem_wave_reflection_diagnostics, only: &
        evaluate_weighted_complex_error, &
        evaluate_weighted_complex_error_jvp, &
        evaluate_weighted_complex_error_vjp, &
        evaluate_weighted_reflection_coefficient, &
        evaluate_weighted_reflection_coefficient_jvp, &
        evaluate_weighted_reflection_coefficient_vjp
    use fortfem_singular_layer_matching, only: &
        assemble_singular_layer_matching, &
        assemble_singular_layer_matching_jvp, &
        assemble_singular_layer_matching_vjp
    use fortfem_equation_objective_metadata, only: &
        OBJECTIVE_METADATA_KIND_EQUATION, OBJECTIVE_METADATA_KIND_OBJECTIVE, &
        OBJECTIVE_METADATA_KIND_CONSTRAINT, OBJECTIVE_METADATA_UNSET_ID, &
        equation_objective_metadata_t, initialize_equation_objective_metadata, &
        validate_equation_objective_metadata, clear_equation_objective_metadata
    use fortfem_equation_objective_merit, only: &
        evaluate_equation_objective_merit, &
        evaluate_equation_objective_merit_jvp, &
        evaluate_equation_objective_merit_vjp
    use fortfem_surface_shape_objective, only: &
        evaluate_surface_shape_objective, &
        evaluate_surface_shape_objective_jvp, &
        evaluate_surface_shape_objective_vjp
    use fortfem_surface_integral_constraint, only: &
        evaluate_surface_integral_constraint, &
        evaluate_surface_integral_constraint_jvp, &
        evaluate_surface_integral_constraint_vjp
    use fortfem_linear_response_cross, only: &
        assemble_linear_response_eigen_cross_residual, &
        assemble_linear_response_eigen_cross_residual_jvp, &
        assemble_linear_response_eigen_cross_residual_vjp, &
        initialize_linear_response_cross_metadata, &
        linear_response_cross_metadata_t, &
        validate_linear_response_cross_metadata
    use fortfem_equilibrium_interchange, only: &
        equilibrium_interchange_t, equilibrium_normalization_t, &
        initialize_equilibrium_interchange, validate_equilibrium_interchange
    use fortfem_equilibrium_sample_adapter, only: &
        build_equilibrium_interchange_sample_set
    use fortfem_toroidal_diagnostic_hooks, only: &
        near_axis_diagnostic_metadata_t, &
        evaluate_boozer_like_rotational_transform, &
        evaluate_boozer_like_rotational_transform_jvp, &
        evaluate_boozer_like_rotational_transform_vjp, &
        evaluate_near_axis_diagnostic_metadata, &
        evaluate_near_axis_diagnostic_metadata_jvp, &
        evaluate_near_axis_diagnostic_metadata_vjp
    use fortfem_critical_point_metadata, only: &
        CRITICAL_CONTOUR_LIMITER, CRITICAL_CONTOUR_NONE, &
        CRITICAL_CONTOUR_SEPARATRIX, CRITICAL_POINT_DEGENERATE, &
        CRITICAL_POINT_O_POINT, CRITICAL_POINT_REGULAR, CRITICAL_POINT_X_POINT, &
        critical_point_metadata_t, evaluate_critical_point_metadata, &
        evaluate_critical_point_metadata_jvp, evaluate_critical_point_metadata_vjp, &
        validate_critical_point_metadata
    use fortfem_oracle_manifest, only: &
        oracle_manifest_schema_magic, oracle_manifest_schema_version, &
        oracle_manifest_t, oracle_normalization_t, oracle_timing_t, &
        oracle_tolerance_t, initialize_oracle_manifest, &
        validate_oracle_manifest, read_oracle_manifest, &
        write_oracle_manifest
    use fortfem_linear_response_interchange, only: &
        assemble_linear_response_operator, &
        assemble_linear_response_operator_jvp, &
        assemble_linear_response_operator_vjp, &
        evaluate_linear_response_diagnostics, &
        assemble_linear_response_residual, &
        assemble_linear_response_residual_jvp, &
        assemble_linear_response_residual_vjp, &
        initialize_linear_response_interchange, &
        linear_response_interchange_t, validate_linear_response_interchange
    use fortfem_linear_response_schema, only: &
        linear_response_schema_magic, &
        read_linear_response_interchange, write_linear_response_interchange
    use fortfem_linear_perturbation_composition, only: &
        LINEAR_PASSIVITY_HERMITIAN_NONNEGATIVE, &
        LINEAR_PASSIVITY_UNDECLARED, &
        LINEAR_RECIPROCITY_TRANSPOSE, LINEAR_RECIPROCITY_UNDECLARED, &
        assemble_linear_perturbation_operator, &
        assemble_linear_perturbation_operator_jvp, &
        assemble_linear_perturbation_operator_vjp, &
        initialize_linear_perturbation_metadata, &
        linear_perturbation_metadata_t, &
        validate_linear_perturbation_metadata
    use fortfem_nonlinear_resistive_mhd_composition, only: &
        RESISTIVE_MHD_AMPERE, RESISTIVE_MHD_ANISOTROPIC_TRANSPORT, &
        RESISTIVE_MHD_BLOCK_COUNT, RESISTIVE_MHD_FARADAY, &
        RESISTIVE_MHD_FREE_BOUNDARY, RESISTIVE_MHD_MOMENTUM, &
        RESISTIVE_MHD_PRESSURE, RESISTIVE_MHD_TENSOR, RESISTIVE_MHD_WALL, &
        assemble_nonlinear_resistive_mhd_residual, &
        assemble_nonlinear_resistive_mhd_residual_jvp, &
        assemble_nonlinear_resistive_mhd_residual_vjp, &
        nonlinear_resistive_mhd_energy_ledger_t
    use fortfem_resistive_mhd_branch_history, only: &
        compare_resistive_mhd_branch_histories, &
        evaluate_resistive_mhd_branch_diagnostics, &
        evaluate_resistive_mhd_branch_path_metric, &
        evaluate_resistive_mhd_branch_path_metric_jvp, &
        evaluate_resistive_mhd_branch_path_metric_vjp, &
        initialize_resistive_mhd_branch_history, &
        resistive_mhd_branch_diagnostics_t, &
        resistive_mhd_branch_history_t, &
        validate_resistive_mhd_branch_history
    use fortfem_pseudo_arclength_residual, only: &
        assemble_pseudo_arclength_residual, &
        assemble_pseudo_arclength_residual_jvp, &
        assemble_pseudo_arclength_residual_vjp
    use fortfem_deflated_residual, only: &
        assemble_deflated_residual, &
        assemble_deflated_residual_jvp, &
        assemble_deflated_residual_vjp
    use fortfem_pseudo_arclength_tangent, only: &
        normalize_pseudo_arclength_tangent, &
        normalize_pseudo_arclength_tangent_jvp, &
        normalize_pseudo_arclength_tangent_vjp
    use fortfem_residual_merit, only: &
        evaluate_residual_merit, evaluate_residual_merit_jvp, &
        evaluate_residual_merit_vjp
    use fortfem_continuation_event, only: &
        CONTINUATION_EVENT_NEAR_ZERO, CONTINUATION_EVENT_NONE, &
        CONTINUATION_EVENT_SIGN_CROSSING, classify_continuation_event
    use fortfem_pseudo_transient_residual, only: &
        assemble_pseudo_transient_residual, &
        assemble_pseudo_transient_residual_jvp, &
        assemble_pseudo_transient_residual_vjp
    use fortfem_eulerian_nonnested_residual, only: &
        assemble_eulerian_nonnested_residual, &
        assemble_eulerian_nonnested_residual_jvp, &
        assemble_eulerian_nonnested_residual_vjp
    use fortfem_step_reduction, only: &
        evaluate_step_reduction, evaluate_step_reduction_jvp, &
        evaluate_step_reduction_vjp
    use fortfem_generalized_eigen_residual, only: &
        assemble_generalized_eigen_residual, &
        assemble_generalized_eigen_residual_jvp, &
        assemble_generalized_eigen_residual_vjp
    use fortfem_beltrami_residual, only: &
        assemble_beltrami_residual, &
        assemble_beltrami_residual_jvp, &
        assemble_beltrami_residual_vjp
    use fortfem_beltrami_parity, only: &
        beltrami_parity_t, compare_beltrami_two_region_residual, &
        validate_beltrami_parity, validate_beltrami_resonance, &
        beltrami_shell_parity_t, compare_beltrami_shell_residual, &
        validate_beltrami_shell_parity
    use fortfem_coupled_field_residual, only: &
        assemble_coupled_field_residual, &
        assemble_coupled_field_residual_jvp, &
        assemble_coupled_field_residual_vjp
    use fortfem_block_2x2_residual, only: &
        assemble_block_2x2_residual, &
        assemble_block_2x2_residual_jvp, &
        assemble_block_2x2_residual_vjp
    use fortfem_block_graph_residual, only: &
        assemble_block_graph_residual, &
        assemble_block_graph_residual_jvp, &
        assemble_block_graph_residual_vjp
    use fortfem_block_graph_csc, only: assemble_block_graph_csc
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
    use fortfem_complex_block_graph_residual, only: &
        assemble_complex_block_graph_residual, &
        assemble_complex_block_graph_residual_jvp, &
        assemble_complex_block_graph_residual_vjp
    use fortfem_boundary_operator_contract, only: &
        BOUNDARY_OPERATOR_BACKEND_BEM, BOUNDARY_OPERATOR_BACKEND_BIEST, &
        BOUNDARY_OPERATOR_BACKEND_DTN, BOUNDARY_OPERATOR_BACKEND_FEM, &
        BOUNDARY_OPERATOR_BACKEND_NESTOR, BOUNDARY_OPERATOR_BACKEND_PML, &
        BOUNDARY_OPERATOR_BACKEND_USER, BOUNDARY_OPERATOR_BACKEND_VIRTUAL_CASING, &
        BOUNDARY_OPERATOR_TRACE_CHANNEL_MIXED, &
        BOUNDARY_OPERATOR_TRACE_CHANNEL_NORMAL, &
        BOUNDARY_OPERATOR_TRACE_CHANNEL_SCALAR, &
        BOUNDARY_OPERATOR_TRACE_CHANNEL_TANGENTIAL, &
        boundary_operator_contract_t, initialize_boundary_operator_contract, &
        initialize_boundary_operator_trace_metadata, &
        validate_boundary_operator_contract
    use fortfem_boundary_operator_parity, only: &
        boundary_operator_parity_t, compare_boundary_operator_parity, &
        compare_boundary_operator_parity_jvp, &
        compare_boundary_operator_parity_vjp, &
        validate_boundary_operator_parity
    use fortfem_surface_trace_parity_ledger, only: &
        evaluate_surface_trace_parity_ledger, &
        evaluate_surface_trace_parity_ledger_jvp, &
        evaluate_surface_trace_parity_ledger_vjp
    use fortfem_boundary_response_diagnostics, only: &
        evaluate_weighted_boundary_response_diagnostics
    use fortfem_source_trace_map, only: &
        evaluate_source_trace_map, evaluate_source_trace_map_jvp, &
        evaluate_source_trace_map_vjp, &
        evaluate_weighted_source_trace_reciprocity_defect
    use fortfem_larger_domain_parity, only: &
        larger_domain_parity_t, compare_larger_domain_solution, &
        compare_larger_domain_solution_jvp, &
        validate_larger_domain_parity
    use fortfem_boundary_trace_residual, only: &
        assemble_boundary_trace_residual, &
        assemble_boundary_trace_residual_jvp, &
        assemble_boundary_trace_residual_vjp
    use fortfem_free_boundary_port_residual, only: &
        assemble_free_boundary_port_residual, &
        assemble_free_boundary_port_residual_jvp, &
        assemble_free_boundary_port_residual_vjp
    use fortfem_free_boundary_source_response, only: &
        assemble_free_boundary_source_response, &
        assemble_free_boundary_source_response_jvp, &
        assemble_free_boundary_source_response_vjp
    use fortfem_toroidal_modal_convolution, only: &
        apply_toroidal_modal_convolution, &
        apply_toroidal_modal_convolution_jvp, &
        apply_toroidal_modal_convolution_vjp
    use fortfem_complex_coupled_field_residual, only: &
        assemble_complex_coupled_field_residual, &
        assemble_complex_coupled_field_residual_jvp, &
        assemble_complex_coupled_field_residual_vjp
    use fortfem_complex_low_rank_response, only: &
        apply_complex_low_rank_matrix, &
        apply_complex_low_rank_matrix_jvp, &
        apply_complex_low_rank_matrix_vjp, &
        compress_complex_matrix_cross, &
        complex_low_rank_matrix_t, &
        initialize_complex_low_rank_matrix, &
        materialize_complex_low_rank_matrix, &
        validate_complex_low_rank_matrix
    use fortfem_mixed_wave_wall_time, only: &
        advance_mixed_wave_wall_midpoint, &
        advance_mixed_wave_wall_midpoint_jvp, &
        advance_mixed_wave_wall_midpoint_vjp, &
        evaluate_mixed_wave_wall_energy_balance
    use fortfem_batched_vector_enrichment_differential_3d, only: &
        evaluate_batched_vector_enrichment_differential_3d, &
        evaluate_batched_vector_enrichment_differential_3d_jvp, &
        evaluate_batched_vector_enrichment_differential_3d_vjp
    use fortfem_complex_boundary_trace_residual, only: &
        assemble_complex_boundary_trace_residual, &
        assemble_complex_boundary_trace_residual_jvp, &
        assemble_complex_boundary_trace_residual_vjp
    use fortfem_mixed_elasticity_residual, only: &
        assemble_mixed_elasticity_residual, &
        assemble_mixed_elasticity_residual_jvp, &
        assemble_mixed_elasticity_residual_vjp
    use fortfem_elasticity_symmetry_constraint, only: &
        assemble_elasticity_symmetry_constraint, &
        assemble_elasticity_symmetry_constraint_jvp, &
        assemble_elasticity_symmetry_constraint_vjp
    use fortfem_symplectic_map_defect, only: &
        assemble_symplectic_map_defect, &
        assemble_symplectic_map_defect_jvp, &
        assemble_symplectic_map_defect_vjp
    use fortfem_glued_feec_sequence, only: &
        assemble_glued_feec_sequence, &
        assemble_glued_feec_sequence_jvp, &
        assemble_glued_feec_sequence_vjp
    use fortfem_glued_feec_sequence_csc, only: &
        assemble_glued_feec_sequence_csc, &
        assemble_glued_feec_sequence_csc_jvp, &
        assemble_glued_feec_sequence_csc_vjp, &
        assemble_glued_feec_sequence_csc_compositions, &
        assemble_glued_feec_sequence_csc_compositions_jvp, &
        assemble_glued_feec_sequence_csc_compositions_vjp
    use fortfem_multipatch_dof_graph, only: &
        build_multipatch_signed_dof_map
    use fortfem_cell_complex, only: cell_complex_betti_numbers, &
        cell_complex_cocycle_basis, cell_complex_cohomology_cocycle_basis, &
        cell_complex_cycle_basis, &
        cell_complex_homology_cycle_basis, &
        cell_complex_harmonic_one_forms, &
        cell_complex_euler_characteristic, cell_complex_t, &
        initialize_cell_complex, quotient_cell_complex, &
        validate_cell_complex
    use fortfem_harmonic_period_normalization, only: &
        normalize_harmonic_one_forms, normalize_harmonic_one_forms_jvp, &
        normalize_harmonic_one_forms_vjp
    use fortfem_cell_identification, only: cell_identification_classes, &
        cell_identification_t, identify_boundary_matrix, &
        initialize_cell_identification, validate_cell_identification
    use fortfem_region_interface_graph, only: &
        initialize_region_interface_graph, &
        region_interface_graph_components, region_interface_graph_incidence, &
        region_interface_graph_cycle_basis, region_interface_graph_t, &
        validate_region_interface_graph
    use fortfem_internal_manifold_graph, only: &
        initialize_internal_manifold_graph, &
        internal_manifold_graph_closed, &
        internal_manifold_graph_components, &
        internal_manifold_graph_junction_incidence, &
        internal_manifold_graph_region_incidence, &
        internal_manifold_graph_t, validate_internal_manifold_graph
    use fortfem_boundary_region_graph, only: &
        boundary_region_graph_t, initialize_boundary_region_graph, &
        validate_boundary_region_graph, boundary_region_graph_incidence, &
        boundary_region_graph_components, boundary_region_graph_cycle_basis, &
        boundary_region_graph_interface_samples, &
        boundary_region_graph_interface_metadata
    use fortfem_surface_current, only: &
        assemble_interface_surface_current, &
        assemble_interface_surface_current_jvp, &
        assemble_interface_surface_current_vjp
    use fortfem_surface_current_potential, only: &
        surface_current_potential_metadata_t, &
        initialize_surface_current_potential_metadata, &
        validate_surface_current_potential_metadata, &
        assemble_surface_current_potential, &
        assemble_surface_current_potential_jvp, &
        assemble_surface_current_potential_vjp
    use fortfem_surface_current_balance, only: &
        assemble_surface_current_junction_balance, &
        assemble_surface_current_junction_balance_jvp, &
        assemble_surface_current_junction_balance_vjp
    use fortfem_surface_current_constraints, only: &
        assemble_surface_current_loop_constraints, &
        assemble_surface_current_loop_constraints_jvp, &
        assemble_surface_current_loop_constraints_vjp
    use fortfem_period_constraints, only: &
        assemble_period_constraints, assemble_period_constraints_jvp, &
        assemble_period_constraints_vjp
    use fortfem_surface_edge_balance, only: &
        assemble_surface_edge_flux_balance, &
        assemble_surface_edge_flux_balance_jvp, &
        assemble_surface_edge_flux_balance_vjp
    use fortfem_surface_edge_flux, only: &
        assemble_surface_edge_flux, assemble_surface_edge_flux_jvp, &
        assemble_surface_edge_flux_vjp
    use fortfem_fci_terminal_boundary_ledger, only: &
        assemble_fci_terminal_boundary_ledger, &
        assemble_fci_terminal_boundary_ledger_jvp, &
        assemble_fci_terminal_boundary_ledger_vjp
    use fortfem_fci_power_flux_ledger, only: &
        evaluate_fci_power_flux_ledger, &
        evaluate_fci_power_flux_ledger_jvp, &
        evaluate_fci_power_flux_ledger_vjp
    use fortfem_volume_balance_ledger, only: &
        assemble_volume_balance_ledger, &
        assemble_volume_balance_ledger_jvp, &
        assemble_volume_balance_ledger_vjp
    use fortfem_surface_current_trace_residual, only: &
        assemble_surface_current_trace_residual, &
        assemble_surface_current_trace_residual_jvp, &
        assemble_surface_current_trace_residual_vjp
    use fortfem_regularized_surface_current_layer, only: &
        evaluate_regularized_surface_current_layer, &
        evaluate_regularized_surface_current_layer_jvp, &
        evaluate_regularized_surface_current_layer_vjp, &
        evaluate_regularized_surface_current_integral
    use fortfem_sheet_current_parity, only: compare_sheet_current_representations
    use fortfem_sheet_current_surface_parity, only: &
        compare_sheet_current_surface_representations, &
        compare_sheet_current_surface_representations_jvp
    use fortfem_broken_skeleton_spaces, only: &
        BROKEN_SPACE_H1, BROKEN_SPACE_HCURL, BROKEN_SPACE_HDIV, &
        BROKEN_SPACE_L2, SKELETON_SPACE_SCALAR, &
        SKELETON_SPACE_TANGENTIAL, SKELETON_SPACE_NORMAL, &
        broken_space_layout_t, skeleton_space_layout_t, &
        initialize_broken_space_layout, validate_broken_space_layout, &
        broken_space_layout_maps, broken_space_layout_global_count, &
        initialize_skeleton_space_layout, validate_skeleton_space_layout, &
        skeleton_space_layout_maps, skeleton_space_layout_global_count
    use fortfem_fourier_mode_registry, only: &
        evaluate_fourier_mode, evaluate_fourier_mode_jvp, &
        evaluate_fourier_mode_vjp, find_fourier_mode, &
        build_fourier_mode_closure_registry, &
        build_fourier_mode_padded_registry, build_fourier_mode_triad_map, &
        fourier_mode_conjugate_index, &
        fourier_mode_registry_t, fourier_mode_triad_closed, &
        initialize_fourier_mode_registry, &
        validate_fourier_mode_registry
    use fortfem_fourier_vector_product, only: &
        apply_fourier_bilinear_product, &
        apply_fourier_bilinear_product_jvp, &
        apply_fourier_bilinear_product_vjp, &
        assemble_fourier_vector_product, &
        assemble_fourier_vector_product_jvp, &
        assemble_fourier_vector_product_vjp
    use fortfem_fourier_mode_linear_operator, only: &
        apply_fourier_mode_linear_operator, &
        apply_fourier_mode_linear_operator_jvp, &
        apply_fourier_mode_linear_operator_vjp
    use fortfem_fourier_mode_expansion, only: &
        evaluate_fourier_mode_expansion, &
        evaluate_fourier_mode_expansion_jvp, &
        evaluate_fourier_mode_expansion_vjp, &
        evaluate_fourier_mode_expansion_hessian, &
        evaluate_fourier_mode_expansion_hvp
    use fortfem_fourier_mode_energy, only: &
        assemble_fourier_mode_energy, &
        assemble_fourier_mode_energy_jvp, &
        assemble_fourier_mode_energy_vjp, &
        fourier_coefficients_conjugate_symmetric
    use fortfem_axis_regular_fourier_modes, only: &
        AXIS_RADIAL_PARITY_EVEN, AXIS_RADIAL_PARITY_ODD, &
        axis_regular_mode_record_t, axis_regular_mode_table_t, &
        axis_regular_mode_requirements, build_axis_regular_mode_table, &
        validate_axis_regular_mode_table
    use fortfem_axis_regular_radial_basis, only: &
        evaluate_axis_regular_radial_basis, &
        evaluate_axis_regular_radial_basis_jvp, &
        evaluate_axis_regular_radial_basis_vjp
    use fortfem_fourier_zernike_basis, only: &
        FOURIER_ZERNIKE_PARITY_EVEN, FOURIER_ZERNIKE_PARITY_ODD, &
        fourier_zernike_mode_t, fourier_zernike_basis_t, &
        build_fourier_zernike_basis, validate_fourier_zernike_basis, &
        fourier_zernike_mode_requirements, evaluate_fourier_zernike_radial
    use fortfem_collocation_grid, only: &
        COLLOCATION_GRID_LINEAR, COLLOCATION_GRID_QUADRATURE, &
        COLLOCATION_GRID_CONCENTRIC, collocation_grid_t, &
        initialize_collocation_grid, validate_collocation_grid, &
        collocation_grid_metadata, collocation_grid_flat_index, &
        collocation_grid_unflatten_index, collocation_grid_chunk_bounds, &
        collocation_grid_point_count
    use fortfem_direct_fourier_transform, only: &
        direct_fourier_plan_t, initialize_direct_fourier_plan, &
        validate_direct_fourier_plan, direct_fourier_plan_metadata, &
        direct_fourier_plan_chunk_bounds, direct_fourier_forward, &
        direct_fourier_adjoint, direct_fourier_plan_sample_count, &
        direct_fourier_plan_mode_count
    use fortfem_interface_traction_balance, only: &
        assemble_normal_traction_jump, &
        assemble_normal_traction_jump_jvp, &
        assemble_normal_traction_jump_vjp, assemble_traction_jump, &
        assemble_traction_jump_jvp, assemble_traction_jump_vjp
    use fortfem_total_pressure_balance, only: &
        assemble_total_pressure_jump, assemble_total_pressure_jump_jvp, &
        assemble_total_pressure_jump_vjp
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
    use fortfem_shifted_vector_enriched_space, only: &
        evaluate_shifted_vector_enriched_space, &
        evaluate_shifted_vector_enriched_space_jvp, &
        evaluate_shifted_vector_enriched_space_vjp
    use fortfem_shifted_vector_enriched_basis, only: &
        evaluate_shifted_vector_enriched_basis, &
        evaluate_shifted_vector_enriched_basis_jvp, &
        evaluate_shifted_vector_enriched_basis_vjp
    use fortfem_vector_enrichment_differential_3d, only: &
        evaluate_vector_enrichment_differential_3d, &
        evaluate_vector_enrichment_differential_3d_jvp, &
        evaluate_vector_enrichment_differential_3d_vjp
    use fortfem_tensor_volume_work, only: &
        assemble_tensor_volume_work, assemble_tensor_volume_work_jvp, &
        assemble_tensor_volume_work_vjp
    use fortfem_force_balance_residual, only: &
        assemble_force_balance_residual, assemble_force_balance_residual_jvp, &
        assemble_force_balance_residual_vjp
    use fortfem_force_balance_objective, only: &
        evaluate_force_balance_objective, &
        evaluate_force_balance_objective_jvp, &
        evaluate_force_balance_objective_vjp, &
        evaluate_force_balance_objective_hvp
    use fortfem_force_balance_product, only: &
        evaluate_force_balance_product, evaluate_force_balance_product_jvp, &
        evaluate_force_balance_product_vjp
    use fortfem_flux_surface_average, only: &
        evaluate_flux_surface_average, evaluate_flux_surface_average_jvp, &
        evaluate_flux_surface_average_vjp
    use fortfem_equation_objective_registry, only: &
        equation_objective_block_t, equation_objective_registry_t, &
        initialize_equation_objective_registry, &
        validate_equation_objective_registry, equation_objective_registry_block, &
        equation_objective_registry_block_count, equation_objective_registry_total_rows, &
        pack_equation_objective_values, pack_equation_objective_values_jvp, &
        pack_equation_objective_values_vjp, unpack_equation_objective_values, &
        REGISTRY_KIND_EQUATION, REGISTRY_KIND_OBJECTIVE, REGISTRY_KIND_CONSTRAINT
    use fortfem_equation_objective_callbacks, only: &
        evaluate_equation_objective_callbacks, &
        evaluate_equation_objective_callbacks_jvp, &
        evaluate_equation_objective_callbacks_vjp
    use fortfem_nested_surface_geometry, only: &
        evaluate_nested_surface_geometry, &
        evaluate_nested_surface_geometry_jvp, &
        evaluate_nested_surface_geometry_vjp, &
        evaluate_nested_surface_geometry_coordinate_jvp, &
        evaluate_nested_surface_geometry_coordinate_vjp
    use fortfem_nested_geometry_differential_jet, only: &
        evaluate_nested_geometry_differential_jet, &
        evaluate_nested_geometry_differential_jet_jvp, &
        evaluate_nested_geometry_differential_jet_vjp, &
        evaluate_nested_geometry_polynomial_jet, &
        validate_nested_geometry_axis_flags
    use fortfem_tensor_diffusion_matrix, only: &
        assemble_tensor_diffusion_matrix, assemble_tensor_diffusion_matrix_jvp, &
        assemble_tensor_diffusion_matrix_vjp, assemble_tensor_diffusion_matrix_nd, &
        assemble_tensor_diffusion_matrix_nd_jvp, &
        assemble_tensor_diffusion_matrix_nd_vjp, &
        assemble_tensor_diffusion_matrix_3d, assemble_tensor_diffusion_matrix_3d_jvp, &
        assemble_tensor_diffusion_matrix_3d_vjp
    use fortfem_dissipative_cayley, only: &
        advance_dissipative_cayley, advance_dissipative_cayley_jvp, &
        advance_dissipative_cayley_vjp
    use fortfem_compatible_flux_elimination, only: &
        assemble_compatible_flux_elimination, &
        assemble_compatible_flux_elimination_jvp, &
        assemble_compatible_flux_elimination_vjp
    use fortfem_xfem_blending_correction, only: &
        evaluate_blending_corrected_enrichment, &
        evaluate_blending_corrected_enrichment_jvp, &
        evaluate_blending_corrected_enrichment_vjp, &
        evaluate_vector_blending_corrected_enrichment, &
        evaluate_vector_blending_corrected_enrichment_jvp, &
        evaluate_vector_blending_corrected_enrichment_vjp
    use fortfem_enrichment_support_diagnostics, only: &
        evaluate_enrichment_support_gram, &
        evaluate_enrichment_support_gram_jvp, &
        evaluate_enrichment_support_gram_vjp, &
        evaluate_enrichment_support_rank_condition
    use fortfem_enrichment_support_activation, only: &
        evaluate_enrichment_support_activation, &
        evaluate_enrichment_support_activation_jvp, &
        evaluate_enrichment_support_activation_vjp
    use fortfem_enrichment_support_vector_diagnostics, only: &
        evaluate_enrichment_support_vector_gram, &
        evaluate_enrichment_support_vector_gram_jvp, &
        evaluate_enrichment_support_vector_gram_vjp
    use fortfem_piola_enriched_vector, only: &
        evaluate_piola_enriched_vector_values, &
        evaluate_piola_enriched_vector_values_jvp, &
        evaluate_piola_enriched_vector_values_vjp, &
        evaluate_piola_enriched_vector_differential_3d, &
        evaluate_piola_enriched_vector_differential_3d_jvp, &
        evaluate_piola_enriched_vector_differential_3d_vjp, &
        evaluate_piola_enriched_vector_differential_2d, &
        evaluate_piola_enriched_vector_differential_2d_jvp, &
        evaluate_piola_enriched_vector_differential_2d_vjp, &
        PIOLA_COVARIANT, PIOLA_CONTRAVARIANT
    use fortfem_tree_cotree_gauge, only: &
        apply_tree_cotree_prolongation, apply_tree_cotree_restriction, &
        build_tree_cotree_dof_map, build_tree_cotree_gauge, &
        reduce_tree_cotree_dense_system, &
        reduce_tree_cotree_dense_system_jvp, reduce_tree_cotree_dense_system_vjp, &
        tree_cotree_gauge_components, &
        tree_cotree_gauge_edges, tree_cotree_gauge_t, &
        validate_tree_cotree_gauge
    use fortfem_tree_cotree_iga_parity, only: &
        diagnose_tree_cotree_iga_invariance, tree_cotree_iga_parity_t
    use fortfem_surface_triangle_geometry_3d, only: &
        evaluate_surface_triangle_geometry_3d, &
        evaluate_surface_triangle_geometry_3d_jvp, &
        evaluate_surface_triangle_geometry_3d_vjp
    use fortfem_surface_triangle_areas_3d, only: &
        assemble_surface_triangle_areas_3d, &
        assemble_surface_triangle_areas_3d_jvp, &
        assemble_surface_triangle_areas_3d_vjp
    use fortfem_surface_triangle_measures_3d, only: &
        assemble_surface_triangle_measures_3d, &
        assemble_surface_triangle_measures_3d_jvp, &
        assemble_surface_triangle_measures_3d_vjp
    use fortfem_level_set_triangle_interface_2d, only: &
        evaluate_level_set_triangle_interface_2d, &
        evaluate_level_set_triangle_interface_2d_jvp, &
        evaluate_level_set_triangle_cut_areas_2d, &
        evaluate_level_set_triangle_cut_quadrature_2d, &
        evaluate_level_set_triangle_cut_quadrature_2d_jvp, &
        evaluate_level_set_triangle_cut_moments_2d, &
        evaluate_level_set_triangle_cut_moments_2d_jvp, &
        evaluate_level_set_triangle_cut_third_moments_2d, &
        evaluate_level_set_triangle_cut_third_moments_2d_jvp, &
        evaluate_level_set_triangle_cut_fourth_moments_2d, &
        evaluate_level_set_triangle_cut_fourth_moments_2d_jvp
    use fortfem_level_set_tetra_interface_3d, only: &
        evaluate_level_set_tetra_interface_3d, &
        evaluate_level_set_tetra_interface_3d_jvp, &
        evaluate_level_set_tetra_cut_quadrature_3d, &
        evaluate_level_set_tetra_cut_quadrature_3d_jvp, &
        evaluate_level_set_tetra_cut_moments_3d, &
        evaluate_level_set_tetra_cut_moments_3d_jvp, &
        evaluate_level_set_tetra_cut_third_moments_3d, &
        evaluate_level_set_tetra_cut_third_moments_3d_jvp, &
        evaluate_level_set_tetra_cut_fourth_moments_3d, &
        evaluate_level_set_tetra_cut_fourth_moments_3d_jvp
    use fortfem_fci_terminal_segment_2d, only: &
        find_fci_first_hit_segment_2d, find_fci_first_hit_segment_2d_jvp
    use fortfem_fci_terminal_triangle_3d, only: &
        find_fci_first_hit_triangle_3d, find_fci_first_hit_triangle_3d_jvp
    use fortfem_fci_terminal_polyline_3d, only: &
        find_fci_first_hit_polyline_3d, find_fci_first_hit_polyline_3d_jvp
    use fortfem_fci_terminal_boundary_flux, only: &
        assemble_fci_terminal_boundary_flux, &
        assemble_fci_terminal_boundary_flux_jvp, &
        assemble_fci_terminal_boundary_flux_vjp
    use fortfem_nonlinear_surface_flux, only: &
        assemble_nonlinear_surface_flux, &
        assemble_nonlinear_surface_flux_jvp, &
        assemble_nonlinear_surface_flux_vjp
    use fortfem_tetra_affine_map, only: &
        invert_tetra_affine_map, invert_tetra_affine_map_jvp, &
        invert_tetra_affine_map_vjp
    use fortfem_tetra_vector_evaluation, only: &
        evaluate_tetra_nedelec_interpolant, &
        evaluate_tetra_nedelec_interpolant_jvp, &
        evaluate_tetra_nedelec_interpolant_vjp, &
        evaluate_tetra_nedelec_interpolant_at_point, &
        evaluate_tetra_nedelec_interpolant_at_point_jvp, &
        evaluate_tetra_nedelec_interpolant_at_point_vjp, &
        evaluate_tetra_rt_interpolant, evaluate_tetra_rt_interpolant_jvp, &
        evaluate_tetra_rt_interpolant_vjp, &
        evaluate_tetra_rt_interpolant_at_point, &
        evaluate_tetra_rt_interpolant_at_point_jvp, &
        evaluate_tetra_rt_interpolant_at_point_vjp
    use fortfem_triangle_affine_map, only: &
        invert_triangle_affine_map, invert_triangle_affine_map_jvp, &
        invert_triangle_affine_map_vjp
    use fortfem_triangle_full_vector_sampled_state_2d, only: &
        solve_triangle_bdm_sampled_state, &
        solve_triangle_bdm_sampled_state_jvp, &
        solve_triangle_bdm_sampled_state_vjp, &
        solve_triangle_nedelec_second_sampled_state, &
        solve_triangle_nedelec_second_sampled_state_jvp, &
        solve_triangle_nedelec_second_sampled_state_vjp
    use fortfem_triangle_rt_sampled_state_2d, only: &
        solve_triangle_rt_sampled_state, &
        solve_triangle_rt_sampled_state_jvp, &
        solve_triangle_rt_sampled_state_vjp
    use fortfem_triangle_nedelec_sampled_state_2d, only: &
        solve_triangle_nedelec_sampled_state, &
        solve_triangle_nedelec_sampled_state_jvp, &
        solve_triangle_nedelec_sampled_state_vjp
    use fortfem_tetra_nedelec_sampled_state_3d, only: &
        solve_tetra_nedelec_sampled_state, &
        solve_tetra_nedelec_sampled_state_jvp, &
        solve_tetra_nedelec_sampled_state_vjp
    use fortfem_tetra_lagrange_state_3d, only: &
        solve_tetra_lagrange_sampled_state, &
        solve_tetra_lagrange_sampled_state_jvp, &
        solve_tetra_lagrange_sampled_state_vjp, &
        solve_tetra_lagrange_state, solve_tetra_lagrange_state_jvp, &
        solve_tetra_lagrange_state_vjp
    use fortfem_tetra_lagrange_curvilinear_pml_state_3d, only: &
        solve_tetra_lagrange_curvilinear_pml, &
        solve_tetra_lagrange_curvilinear_pml_jvp, &
        solve_tetra_lagrange_curvilinear_pml_vjp
    use fortfem_tetra_mixed_poisson_state_3d, only: &
        solve_tetra_mixed_poisson_sampled_state, &
        solve_tetra_mixed_poisson_sampled_state_jvp, &
        solve_tetra_mixed_poisson_sampled_state_vjp, &
        solve_tetra_mixed_poisson_state, &
        solve_tetra_mixed_poisson_state_jvp, &
        solve_tetra_mixed_poisson_state_vjp
    use fortfem_mixed_rt_system, only: solve_mixed_rt_system, &
        solve_mixed_rt_system_jvp, solve_mixed_rt_system_vjp
    use fortfem_assembly_bspline_polar_2d, only: &
        assemble_bspline_polar_hcurl_operator_csc, &
        assemble_bspline_polar_h1_operator_csc, &
        assemble_bspline_polar_l2_mass_csc, &
        restrict_bspline_polar_operator_csc
    use fortfem_magnetic_curvilinear_coefficients_2d, only: &
        scalar_reluctivity_curvilinear_fourier_coefficients
    use fortfem_cartesian_pml_geometry, only: &
        build_cartesian_pml_element_stretch, &
        build_cartesian_pml_element_stretch_jvp, &
        build_cartesian_pml_element_stretch_vjp
    use fortfem_curvilinear_pml_geometry, only: &
        build_curvilinear_normal_pml_element_stretch, &
        build_curvilinear_normal_pml_element_stretch_jvp, &
        build_curvilinear_normal_pml_element_stretch_vjp
    use fortfem_assembly_full_vector_arbitrary_order_2d, only: &
        assemble_triangle_bdm_div_mass_csc, &
        assemble_triangle_bdm_div_mass_csc_jvp, &
        assemble_triangle_bdm_div_mass_csc_vjp, &
        assemble_triangle_bdm_div_mass_element, &
        assemble_triangle_bdm_div_mass_element_jvp, &
        assemble_triangle_bdm_div_mass_element_vjp, &
        assemble_triangle_bdm_vector_load_samples, &
        assemble_triangle_bdm_vector_load_samples_jvp, &
        assemble_triangle_bdm_vector_load_samples_vjp, &
        assemble_triangle_nedelec_second_curl_mass_csc, &
        assemble_triangle_nedelec_second_curl_mass_csc_jvp, &
        assemble_triangle_nedelec_second_curl_mass_csc_vjp, &
        assemble_triangle_nedelec_second_curl_mass_element, &
        assemble_triangle_nedelec_second_curl_mass_element_jvp, &
        assemble_triangle_nedelec_second_curl_mass_element_vjp, &
        assemble_triangle_nedelec_second_vector_load_samples, &
        assemble_triangle_nedelec_second_vector_load_samples_jvp, &
        assemble_triangle_nedelec_second_vector_load_samples_vjp
    use fortfem_assembly_nedelec_arbitrary_order_2d, only: &
        assemble_triangle_nedelec_curl_csc, &
        assemble_triangle_nedelec_curl_mass_csc, &
        assemble_triangle_nedelec_curl_mass_csc_jvp, &
        assemble_triangle_nedelec_curl_mass_csc_vjp, &
        assemble_triangle_nedelec_curl_mass_element, &
        assemble_triangle_nedelec_curl_mass_element_jvp, &
        assemble_triangle_nedelec_curl_mass_element_vjp, &
        assemble_triangle_nedelec_vector_load_samples, &
        assemble_triangle_nedelec_vector_load_samples_jvp, &
        assemble_triangle_nedelec_vector_load_samples_vjp
    use fortfem_assembly_rt_arbitrary_order_2d, only: &
        assemble_triangle_rt_div_mass_csc, &
        assemble_triangle_rt_div_mass_csc_jvp, &
        assemble_triangle_rt_div_mass_csc_vjp, &
        assemble_triangle_rt_div_mass_element, &
        assemble_triangle_rt_div_mass_element_jvp, &
        assemble_triangle_rt_div_mass_element_vjp, &
        assemble_triangle_rt_divergence_csc, &
        assemble_triangle_rt_vector_load_samples, &
        assemble_triangle_rt_vector_load_samples_jvp, &
        assemble_triangle_rt_vector_load_samples_vjp
    use fortfem_assembly_tetra_rt_arbitrary_order_3d, only: &
        assemble_tetra_rt_div_mass_csc, assemble_tetra_rt_div_mass_csc_jvp, &
        assemble_tetra_rt_div_mass_csc_vjp, &
        assemble_tetra_rt_div_mass_element, &
        assemble_tetra_rt_div_mass_element_jvp, &
        assemble_tetra_rt_div_mass_element_vjp, &
        assemble_tetra_rt_divergence_csc, &
        assemble_tetra_rt_vector_load_samples, &
        assemble_tetra_rt_vector_load_samples_jvp, &
        assemble_tetra_rt_vector_load_samples_vjp
    use fortfem_assembly_tetra_lagrange_arbitrary_order_3d, only: &
        assemble_tetra_lagrange_scalar_load, &
        assemble_tetra_lagrange_scalar_load_samples, &
        assemble_tetra_lagrange_scalar_load_samples_jvp, &
        assemble_tetra_lagrange_scalar_load_samples_vjp, &
        assemble_tetra_lagrange_stiffness_csc, &
        assemble_tetra_lagrange_stiffness_csc_jvp, &
        assemble_tetra_lagrange_stiffness_csc_vjp, &
        assemble_tetra_lagrange_stiffness_element, &
        assemble_tetra_lagrange_stiffness_element_jvp, &
        assemble_tetra_lagrange_stiffness_element_vjp
    use fortfem_assembly_tetra_lagrange_curvilinear_pml_3d, only: &
        assemble_tetra_lagrange_curvilinear_pml_element, &
        assemble_tetra_lagrange_curvilinear_pml_element_jvp, &
        assemble_tetra_lagrange_curvilinear_pml_element_vjp, &
        assemble_tetra_lagrange_curvilinear_pml_csc, &
        assemble_tetra_lagrange_curvilinear_pml_csc_jvp, &
        assemble_tetra_lagrange_curvilinear_pml_csc_vjp
    use fortfem_assembly_tetra_lagrange_geometry_pml_3d, only: &
        assemble_tetra_lagrange_geometry_pml_csc, &
        assemble_tetra_lagrange_geometry_pml_csc_jvp, &
        assemble_tetra_lagrange_geometry_pml_csc_vjp
    use fortfem_assembly_tetra_nedelec_3d, only: &
        assemble_tetra_nedelec_curl_mass_element, &
        assemble_tetra_nedelec_curl_mass_element_jvp, &
        assemble_tetra_nedelec_curl_mass_element_vjp, &
        assemble_tetra_nedelec_curl_mass_csc, &
        assemble_tetra_nedelec_curl_mass_csc_jvp, &
        assemble_tetra_nedelec_curl_mass_csc_vjp, &
        assemble_tetra_nedelec_pml_element, &
        assemble_tetra_nedelec_pml_element_jvp, &
        assemble_tetra_nedelec_pml_element_vjp, &
        assemble_tetra_nedelec_curvilinear_pml_element, &
        assemble_tetra_nedelec_curvilinear_pml_element_jvp, &
        assemble_tetra_nedelec_curvilinear_pml_element_vjp, &
        assemble_tetra_nedelec_curvilinear_pml_csc, &
        assemble_tetra_nedelec_curvilinear_pml_csc_jvp, &
        assemble_tetra_nedelec_curvilinear_pml_csc_vjp, &
        assemble_tetra_nedelec_pml_csc, &
        assemble_tetra_nedelec_pml_csc_jvp, &
        assemble_tetra_nedelec_pml_csc_vjp, &
        assemble_tetra_nedelec_vector_load, &
        assemble_tetra_nedelec_vector_load_samples, &
        assemble_tetra_nedelec_vector_load_samples_jvp, &
        assemble_tetra_nedelec_vector_load_samples_vjp, &
        assemble_tetra_nedelec_weighted_csc
    use fortfem_assembly_tetra_nedelec_geometry_pml_3d, only: &
        assemble_tetra_nedelec_geometry_pml_csc, &
        assemble_tetra_nedelec_geometry_pml_csc_jvp, &
        assemble_tetra_nedelec_geometry_pml_csc_vjp
    use fortfem_kinds
    use fortfem_bspline_polar, only: &
        build_bspline_polar_feec_2d_extractions, &
        build_bspline_polar_feec_2d_operators, &
        build_bspline_polar_h1_extraction, evaluate_periodic_bspline_basis
    use fortfem_bspline_multipatch, only: &
        BSPLINE_FACE_X_MAX, BSPLINE_FACE_X_MIN, BSPLINE_FACE_Y_MAX, &
        BSPLINE_FACE_Y_MIN, BSPLINE_FACE_Z_MAX, BSPLINE_FACE_Z_MIN, &
        build_bspline_feec_2d_interface_dofs, &
        build_bspline_feec_2d_multipatch_maps, &
        build_bspline_feec_2d_two_patch_maps, &
        build_bspline_feec_3d_interface_dofs, &
        build_bspline_feec_3d_multipatch_maps, &
        build_bspline_feec_3d_two_patch_maps
    use fortfem_assembly_bspline_multipatch_2d, only: &
        build_bspline_feec_2d_two_patch_operators_csc
    use fortfem_assembly_bspline_multipatch_3d, only: &
        build_bspline_feec_3d_two_patch_operators_csc
    use fortfem_boundary, only: boundary_t
    use fortfem_api_types
    use fortfem_api_mesh
    use fortfem_api_spaces
    use fortfem_api_forms
    use fortfem_api_solvers
    use fortfem_fci_parallel_operator, only: &
        assemble_fci_parallel_gradient_csc, &
        assemble_fci_parallel_support_divergence_csc, &
        apply_fci_parallel_gradient, apply_fci_parallel_diffusion, &
        apply_fci_parallel_diffusion_jvp, apply_fci_parallel_diffusion_vjp, &
        apply_fci_parallel_diffusion_field_vjp, &
        compute_fci_parallel_diffusion_diagonal, &
        compute_fci_anisotropic_diffusion_diagonal, &
        apply_fci_parallel_jacobi_preconditioner, &
        apply_fci_anisotropic_jacobi_preconditioner, &
        apply_fci_anisotropic_diffusion, &
        apply_fci_anisotropic_diffusion_field_vjp, &
        apply_fci_parallel_gradient_jvp, &
        apply_fci_parallel_gradient_vjp
    use fortfem_fci_plane_multigrid, only: &
        apply_fci_plane_two_level_vcycle, &
        factor_fci_plane_coarse_operator, &
        apply_fci_plane_two_level_vcycle_factored, &
        apply_fci_plane_multilevel_vcycle, &
        apply_fci_plane_multilevel_wcycle, &
        apply_fci_plane_two_level_vcycles, apply_fci_plane_two_level_vcycles_ragged
    use fortfem_fci_field_split_preconditioner, only: &
        apply_fci_additive_field_split_preconditioner
    use fortfem_fci_field_line_tracer, only: trace_fci_field_line_rk4
    use fortfem_fci_field_line_tracer, only: trace_fci_field_line_rk4_jvp
    use fortfem_fci_interpolation_map, only: &
        build_fci_linear_interpolation_map_1d, &
        build_fci_linear_interpolation_map_1d_jvp, &
        build_fci_linear_interpolation_map_1d_vjp, &
        build_fci_quadratic_interpolation_map_1d, &
        build_fci_quadratic_interpolation_map_1d_jvp, &
        build_fci_quadratic_interpolation_map_1d_vjp, &
        build_fci_cubic_interpolation_map_1d, &
        build_fci_cubic_interpolation_map_1d_jvp, &
        build_fci_cubic_interpolation_map_1d_vjp, &
        build_fci_quartic_interpolation_map_1d, &
        build_fci_quartic_interpolation_map_1d_jvp, &
        build_fci_quartic_interpolation_map_1d_vjp, &
        build_fci_quintic_interpolation_map_1d, &
        build_fci_quintic_interpolation_map_1d_jvp, &
        build_fci_quintic_interpolation_map_1d_vjp, &
        build_fci_sextic_interpolation_map_1d, &
        build_fci_sextic_interpolation_map_1d_jvp, &
        build_fci_sextic_interpolation_map_1d_vjp, &
        build_fci_quadratic_interpolation_maps_1d, &
        build_fci_quadratic_interpolation_maps_1d_jvp, &
        build_fci_quadratic_interpolation_maps_1d_vjp, &
        build_fci_bilinear_interpolation_map_2d, &
        build_fci_bilinear_interpolation_map_2d_jvp, &
        build_fci_bilinear_interpolation_map_2d_vjp, &
        build_fci_bilinear_interpolation_maps_2d, &
        build_fci_bilinear_interpolation_maps_2d_jvp, &
        build_fci_bilinear_interpolation_maps_2d_vjp, &
        build_fci_triangle_interpolation_map_2d, &
        build_fci_triangle_interpolation_map_2d_jvp, &
        build_fci_triangle_interpolation_map_2d_vjp, &
        build_fci_triangle_interpolation_maps_2d, &
        build_fci_triangle_interpolation_maps_2d_jvp, &
        build_fci_triangle_interpolation_maps_2d_vjp
    use fortfem_fci_support_geometry, only: &
        compute_fci_staggered_flux_box_volumes, &
        compute_fci_staggered_flux_box_volumes_jvp, &
        compute_fci_staggered_flux_box_volumes_vjp, &
        compute_fci_quadrilateral_cell_areas_2d, &
        compute_fci_quadrilateral_cell_areas_2d_jvp, &
        compute_fci_quadrilateral_cell_areas_2d_vjp, &
        compute_fci_polygon_cell_areas_2d, &
        compute_fci_polygon_cell_areas_2d_jvp, &
        compute_fci_polygon_cell_areas_2d_vjp, &
        compute_fci_curved_polygon_cell_areas_2d, &
        compute_fci_curved_polygon_cell_areas_2d_jvp, &
        compute_fci_curved_polygon_cell_areas_2d_vjp, &
        compute_fci_cubic_curved_polygon_cell_areas_2d, &
        compute_fci_cubic_curved_polygon_cell_areas_2d_jvp, &
        compute_fci_cubic_curved_polygon_cell_areas_2d_vjp, &
        compute_fci_quartic_curved_polygon_cell_areas_2d, &
        compute_fci_quartic_curved_polygon_cell_areas_2d_jvp, &
        compute_fci_quartic_curved_polygon_cell_areas_2d_vjp, &
        compute_fci_quintic_curved_polygon_cell_areas_2d, &
        compute_fci_quintic_curved_polygon_cell_areas_2d_jvp, &
        compute_fci_quintic_curved_polygon_cell_areas_2d_vjp, &
        compute_fci_sextic_curved_polygon_cell_areas_2d, &
        compute_fci_sextic_curved_polygon_cell_areas_2d_jvp, &
        compute_fci_sextic_curved_polygon_cell_areas_2d_vjp, &
        compute_fci_curved_quadrilateral_cell_areas_2d, &
        compute_fci_curved_quadrilateral_cell_areas_2d_jvp, &
        compute_fci_curved_quadrilateral_cell_areas_2d_vjp
    use fortfem_mixed_wave_time, only: advance_mixed_wave_midpoint, &
        assemble_mixed_wave_midpoint_map, &
        advance_mixed_wave_midpoint_jvp, advance_mixed_wave_midpoint_vjp, &
        advance_mixed_wave_symplectic_euler, &
        advance_mixed_wave_symplectic_euler_jvp, &
        advance_mixed_wave_symplectic_euler_vjp, advance_mixed_wave_strang, &
        advance_mixed_wave_strang_jvp, advance_mixed_wave_strang_vjp
    use fortfem_quadratic_average_vector_field, only: &
        advance_quadratic_avf, advance_quadratic_avf_jvp, &
        advance_quadratic_avf_vjp
    use fortfem_cgl_pressure_tensor, only: &
        evaluate_cgl_pressure_tensor, evaluate_cgl_pressure_tensor_jvp, &
        evaluate_cgl_pressure_tensor_vjp, evaluate_cgl_pressure_traction, &
        evaluate_cgl_pressure_traction_jvp, evaluate_cgl_pressure_traction_vjp, &
        evaluate_cgl_pressure_work, evaluate_cgl_pressure_work_jvp, &
        evaluate_cgl_pressure_work_vjp
    use fortfem_field_aligned_flux, only: &
        evaluate_field_aligned_flux, evaluate_field_aligned_flux_jvp, &
        evaluate_field_aligned_flux_vjp
    use fortfem_field_aligned_constitutive_tensor, only: &
        evaluate_field_aligned_constitutive_tensor, &
        evaluate_field_aligned_constitutive_tensor_jvp, &
        evaluate_field_aligned_constitutive_tensor_vjp
    use fortfem_field_aligned_tensor_pullback, only: &
        pullback_field_aligned_tensor, pullback_field_aligned_tensor_jvp, &
        pullback_field_aligned_tensor_vjp
    use fortfem_tensor_power_split, only: evaluate_tensor_power_split, &
        evaluate_tensor_power_split_jvp, evaluate_tensor_power_split_vjp
    use fortfem_spherical_harmonics, only: spherical_harmonic, &
        spherical_harmonic_product_coefficient, &
        spherical_harmonic_theta_derivative, spherical_harmonic_phi_derivative
    use fortfem_toroidal_harmonics, only: toroidal_p, toroidal_q, &
        toroidal_p_derivative, toroidal_q_derivative, &
        toroidal_p_second_derivative, toroidal_q_second_derivative
    use fortfem_nestor_fourier_response, only: &
        apply_nestor_fourier_response_map, &
        apply_nestor_fourier_response_map_jvp, &
        apply_nestor_fourier_response_map_vjp, &
        evaluate_nestor_fourier_reciprocity_defect
    use fortfem_interface_traces, only: &
        compute_interface_scalar_jump_average, compute_interface_vector_traces
    use fortfem_surface_delta_load, only: &
        assemble_surface_delta_load, assemble_surface_delta_load_jvp, &
        assemble_surface_delta_load_vjp, assemble_surface_vector_delta_load, &
        assemble_surface_vector_delta_load_jvp, &
        assemble_surface_vector_delta_load_vjp
    use fortfem_interface_jump_penalty, only: assemble_interface_jump_penalty, &
        assemble_interface_jump_penalty_jvp, assemble_interface_jump_penalty_vjp
    use fortfem_nitsche_interface, only: assemble_symmetric_nitsche_interface, &
        assemble_symmetric_nitsche_interface_jvp, &
        assemble_symmetric_nitsche_interface_vjp
    use fortfem_scalar_sipg_interface, only: assemble_scalar_sipg_interface, &
        assemble_scalar_sipg_interface_jvp, assemble_scalar_sipg_interface_vjp
    use fortfem_vector_jump_penalty, only: assemble_vector_jump_penalty, &
        assemble_vector_jump_penalty_jvp, assemble_vector_jump_penalty_vjp, &
        assemble_vector_sipg_interface, assemble_vector_sipg_interface_jvp, &
        assemble_vector_sipg_interface_vjp
    use fortfem_hdg_static_condensation, only: assemble_hdg_static_condensation, &
        assemble_hdg_static_condensation_jvp, assemble_hdg_static_condensation_vjp
    use fortfem_hdg_global_skeleton, only: assemble_hdg_global_skeleton, &
        assemble_hdg_global_skeleton_jvp, assemble_hdg_global_skeleton_vjp
    use fortfem_hdg_global_skeleton_csc, only: &
        assemble_hdg_global_skeleton_csc, &
        assemble_hdg_global_skeleton_csc_jvp, &
        assemble_hdg_global_skeleton_csc_vjp
    use fortfem_feec_exact_sequence, only: assemble_feec_exact_sequence, &
        assemble_feec_exact_sequence_jvp, assemble_feec_exact_sequence_vjp
    use fortfem_enriched_feec_sequence, only: &
        assemble_enriched_feec_sequence, assemble_enriched_feec_sequence_jvp, &
        assemble_enriched_feec_sequence_vjp
    use fortfem_broken_feec_sequence, only: &
        assemble_broken_feec_sequence, assemble_broken_feec_sequence_jvp, &
        assemble_broken_feec_sequence_vjp
    use fortfem_fitted_trace_constraint, only: &
        assemble_fitted_trace_constraint, &
        assemble_fitted_trace_constraint_jvp, &
        assemble_fitted_trace_constraint_vjp
    use fortfem_partition_layout, only: &
        partition_layout_t, initialize_partition_layout, &
        validate_partition_layout, partition_layout_maps, &
        partition_layout_global_count, partition_layout_owned_count, &
        partition_layout_ghost_count, assemble_partitioned_sum, &
        assemble_partitioned_sum_jvp, assemble_partitioned_sum_vjp
    use fortfem_distributed_trace_ownership, only: &
        distributed_trace_layout_t, initialize_distributed_trace_layout, &
        validate_distributed_trace_layout, &
        distributed_trace_layout_partition_count, &
        distributed_trace_layout_global_count, &
        distributed_trace_layout_local_count, &
        distributed_trace_layout_component_count, &
        assemble_distributed_trace_reduction, &
        assemble_distributed_trace_reduction_jvp, &
        assemble_distributed_trace_reduction_vjp
    use fortfem_physical_trace_ownership, only: &
        physical_trace_ownership_t, initialize_physical_trace_ownership, &
        validate_physical_trace_ownership, physical_trace_ownership_maps, &
        physical_trace_ownership_dimension, physical_trace_ownership_point_count, &
        physical_trace_ownership_rank, compare_physical_trace_coordinates
    use fortfem_physical_trace_reconciliation, only: &
        physical_trace_reconciliation_t, initialize_physical_trace_reconciliation, &
        validate_physical_trace_reconciliation, physical_trace_reconciliation_maps, &
        reconcile_physical_trace_values, reconcile_physical_trace_values_jvp, &
        reconcile_physical_trace_values_vjp
    use fortfem_physical_trace_owner_selection, only: &
        physical_trace_owner_selection_t, initialize_physical_trace_owner_selection, &
        validate_physical_trace_owner_selection, physical_trace_owner_selection_maps, &
        gather_physical_trace_values, gather_physical_trace_values_jvp, &
        gather_physical_trace_values_vjp
    use fortfem_feec_commuting_projection, only: &
        assemble_feec_commuting_projection, &
        assemble_feec_commuting_projection_jvp, &
        assemble_feec_commuting_projection_vjp
    use fortfem_scalar_numerical_flux, only: assemble_scalar_numerical_flux, &
        assemble_scalar_numerical_flux_jvp, assemble_scalar_numerical_flux_vjp, &
        NUMERICAL_FLUX_CENTRAL, NUMERICAL_FLUX_UPWIND, &
        NUMERICAL_FLUX_LAX_FRIEDRICHS
    use fortfem_vector_numerical_flux, only: assemble_vector_numerical_flux, &
        assemble_vector_numerical_flux_jvp, assemble_vector_numerical_flux_vjp, &
        assemble_vector_entropy_stable_flux, &
        assemble_vector_entropy_stable_flux_jvp, &
        assemble_vector_entropy_stable_flux_vjp
    use fortfem_mortar_trace_coupling, only: assemble_mortar_trace_coupling, &
        assemble_mortar_trace_coupling_jvp, assemble_mortar_trace_coupling_vjp
    use fortfem_geometry_mortar_trace_coupling, only: &
        assemble_geometry_mortar_trace_coupling, &
        assemble_geometry_mortar_trace_coupling_jvp, &
        assemble_geometry_mortar_trace_coupling_vjp
    use fortfem_geometry_mortar_component_coupling, only: &
        assemble_geometry_mortar_component_coupling, &
        assemble_geometry_mortar_component_coupling_jvp, &
        assemble_geometry_mortar_component_coupling_vjp
    use fortfem_physical_surface_geometry, only: &
        sample_physical_surface_geometry, &
        sample_physical_surface_geometry_jvp, &
        sample_physical_surface_geometry_vjp
    use fortfem_surface_vector_trace, only: &
        evaluate_surface_vector_trace, &
        evaluate_surface_vector_trace_jvp, &
        evaluate_surface_vector_trace_vjp
    use fortfem_fci_boundary_patch_mortar, only: &
        assemble_fci_boundary_patch_mortar, &
        assemble_fci_boundary_patch_mortar_jvp, &
        assemble_fci_boundary_patch_mortar_vjp
    use fortfem_cgl_pressure_divergence, only: &
        evaluate_cgl_pressure_divergence, evaluate_cgl_pressure_divergence_jvp, &
        evaluate_cgl_pressure_divergence_vjp
    use fortfem_api_plot
    use fortfem_torus_surface_mesh, only: generate_torus_surface_mesh
    use fortfem_solid_torus_tetra_mesh, only: generate_solid_torus_tetra_mesh
    use fortfem_structured_tetra_box_mesh, only: &
        generate_structured_tetra_box_mesh
    use fortfem_planar_nedelec_maxwell_dtn, only: &
        assemble_planar_nedelec_maxwell_dtn_form, &
        build_planar_nedelec_trace_sampling, &
        pullback_planar_maxwell_dtn_form, &
        pullback_planar_maxwell_dtn_form_jvp, &
        pullback_planar_maxwell_dtn_form_vjp
    use fortfem_torus_curved_panel, only: &
        evaluate_torus_curved_panel, evaluate_torus_curved_panel_jvp, &
        evaluate_torus_curved_panel_vjp
    use fortfem_laplace_torus_panel_pair_ad_3d, only: &
        integrate_laplace_torus_panel_p0_3d, &
        integrate_laplace_torus_panel_p0_3d_jvp, &
        integrate_laplace_torus_panel_p0_3d_vjp
    use fortfem_laplace_sphere_panel_pair_ad_3d, only: &
        integrate_laplace_sphere_panel_p0_3d, &
        integrate_laplace_sphere_panel_p0_3d_jvp, &
        integrate_laplace_sphere_panel_p0_3d_vjp
    use fortfem_helmholtz_torus_panel_pair_ad_3d, only: &
        integrate_helmholtz_torus_panel_p0_3d, &
        integrate_helmholtz_torus_panel_p0_3d_jvp, &
        integrate_helmholtz_torus_panel_p0_3d_vjp
    use fortfem_helmholtz_sphere_panel_pair_ad_3d, only: &
        integrate_helmholtz_sphere_panel_p0_3d, &
        integrate_helmholtz_sphere_panel_p0_3d_jvp, &
        integrate_helmholtz_sphere_panel_p0_3d_vjp
    use fortfem_sphere_surface_mesh, only: generate_sphere_surface_mesh
    use fortfem_sphere_curved_panel, only: &
        evaluate_sphere_curved_panel, evaluate_sphere_curved_panel_jvp, &
        evaluate_sphere_curved_panel_vjp, invert_sphere_curved_panel
    use fortfem_maxwell_sphere_curved_rwg, only: &
        assemble_maxwell_sphere_curved_efie_rwg_3d, &
        assemble_maxwell_sphere_curved_efie_imaginary_rwg_3d, &
        assemble_maxwell_sphere_efie_propagating_impedance_jvp, &
        assemble_maxwell_sphere_efie_propagating_impedance_vjp, &
        assemble_maxwell_sphere_efie_imaginary_impedance_jvp, &
        assemble_maxwell_sphere_efie_imaginary_impedance_vjp, &
        assemble_maxwell_sphere_efie_wave_number_jvp, &
        assemble_maxwell_sphere_efie_wave_number_vjp, &
        assemble_maxwell_sphere_efie_imaginary_decay_jvp, &
        assemble_maxwell_sphere_efie_imaginary_decay_vjp, &
        assemble_maxwell_sphere_curved_efie_bc_imaginary_3d, &
        assemble_maxwell_sphere_curved_regularized_cfie_rwg_3d, &
        assemble_maxwell_sphere_curved_regularized_cfie_rhs_rwg_3d, &
        assemble_maxwell_sphere_curved_mfie_exterior_trace_rwg_rbc_3d, &
        assemble_maxwell_sphere_curved_mfie_rwg_rbc_3d, &
        assemble_maxwell_sphere_curved_plane_wave_rhs_rwg_3d, &
        assemble_maxwell_sphere_curved_plane_wave_rhs_bc_3d, &
        assemble_maxwell_sphere_curved_plane_wave_rhs_bc_3d_jvp, &
        assemble_maxwell_sphere_curved_plane_wave_rhs_bc_3d_vjp, &
        assemble_maxwell_sphere_curved_potential_operators_rwg_3d, &
        assemble_maxwell_sphere_curved_rwg_mass_matrix, &
        assemble_maxwell_sphere_curved_rwg_mass_matrix_jvp, &
        assemble_maxwell_sphere_curved_rwg_mass_matrix_vjp, &
        assemble_maxwell_sphere_curved_rwg_rbc_pairing, &
        assemble_maxwell_sphere_curved_rwg_rbc_pairing_jvp, &
        assemble_maxwell_sphere_curved_rwg_rbc_pairing_vjp, &
        assemble_maxwell_sphere_curved_vector_potential_rwg_3d, &
        evaluate_maxwell_sphere_curved_far_field_rwg_3d, &
        evaluate_maxwell_sphere_curved_far_field_rwg_3d_jvp, &
        evaluate_maxwell_sphere_curved_far_field_rwg_3d_vjp, &
        evaluate_maxwell_sphere_curved_localized_rwg_basis, &
        evaluate_maxwell_sphere_curved_localized_rwg_basis_jvp, &
        evaluate_maxwell_sphere_curved_localized_rwg_basis_vjp, &
        evaluate_maxwell_sphere_curved_magnetic_field_rwg_3d, &
        evaluate_maxwell_sphere_curved_magnetic_field_rwg_3d_jvp, &
        evaluate_maxwell_sphere_curved_magnetic_field_rwg_3d_vjp, &
        evaluate_maxwell_sphere_curved_rwg_basis, &
        evaluate_maxwell_sphere_curved_rwg_basis_jvp, &
        evaluate_maxwell_sphere_curved_rwg_basis_vjp, &
        integrate_maxwell_sphere_curved_adjacent_rwg_pair_3d, &
        integrate_maxwell_sphere_curved_coincident_rwg_pair_3d, &
        solve_maxwell_pec_sphere_curved_efie_rwg_3d, &
        solve_maxwell_pec_sphere_curved_regularized_cfie_rwg_3d
    use fortfem_barycentric_surface_refinement, only: &
        barycentric_refine_surface_mesh, &
        barycentric_refine_sphere_surface_mesh_jvp, &
        barycentric_refine_sphere_surface_mesh_vjp, &
        barycentric_refine_torus_surface_mesh, &
        barycentric_refine_torus_surface_mesh_jvp, &
        barycentric_refine_torus_surface_mesh_vjp
    use fortfem_maxwell_localized_rwg_surface, only: &
        evaluate_maxwell_localized_rwg_basis
    use fortfem_maxwell_bc_surface, only: &
        assemble_maxwell_rwg_rbc_pairing, build_maxwell_bc_transformation, &
        differentiate_maxwell_bc_transformation_jvp, &
        differentiate_maxwell_bc_transformation_vjp
    use fortfem_maxwell_magnetic_rwg_3d, only: &
        evaluate_maxwell_magnetic_field_rwg_3d, &
        evaluate_maxwell_magnetic_field_rwg_3d_jvp, &
        evaluate_maxwell_magnetic_field_rwg_3d_vjp, &
        evaluate_maxwell_magnetic_field_rwg_3d_targets, &
        evaluate_maxwell_magnetic_field_rwg_3d_targets_jvp, &
        evaluate_maxwell_magnetic_field_rwg_3d_targets_vjp
    use fortfem_maxwell_virtual_casing, only: &
        evaluate_maxwell_virtual_casing_rwg_3d, &
        evaluate_maxwell_virtual_casing_rwg_3d_jvp, &
        evaluate_maxwell_virtual_casing_rwg_3d_vjp
    use fortfem_resistive_wall_midpoint, only: &
        advance_resistive_wall_midpoint, advance_resistive_wall_midpoint_jvp, &
        advance_resistive_wall_midpoint_vjp, evaluate_resistive_wall_energy_balance
    use fortfem_toroidal_spectral_trace, only: &
        evaluate_toroidal_spectral_trace, evaluate_toroidal_spectral_trace_jvp, &
        evaluate_toroidal_spectral_trace_vjp, &
        evaluate_toroidal_spectral_trace_grid, &
        evaluate_toroidal_spectral_trace_grid_jvp, &
        evaluate_toroidal_spectral_trace_grid_vjp, &
        solve_toroidal_spectral_neumann, &
        solve_toroidal_spectral_neumann_jvp, &
        solve_toroidal_spectral_neumann_vjp
    use fortfem_toroidal_spectral_metadata, only: &
        analyze_toroidal_spectral_modes, &
        analyze_toroidal_spectral_modes_jvp, &
        analyze_toroidal_spectral_modes_vjp
    use fortfem_wall_response_condensation, only: &
        condense_wall_response_blocks, condense_wall_response_blocks_jvp, &
        condense_wall_response_blocks_vjp
    use fortfem_generalized_debye_source, only: &
        assemble_generalized_debye_source_residual, &
        assemble_generalized_debye_source_residual_jvp, &
        assemble_generalized_debye_source_residual_vjp, &
        assemble_generalized_debye_source_second_kind, &
        assemble_generalized_debye_source_second_kind_jvp, &
        assemble_generalized_debye_source_second_kind_vjp
    use fortfem_maxwell_mfie_rwg_rbc_3d, only: &
        assemble_maxwell_mfie_rwg_rbc_3d
    use fortfem_maxwell_efie_bc_3d, only: &
        assemble_maxwell_bc_potential_operators_3d, &
        assemble_maxwell_bc_scalar_potential_3d, &
        assemble_maxwell_efie_bc_3d, &
        assemble_maxwell_efie_bc_imaginary_3d, &
        build_maxwell_bc_panel_divergence, &
        build_maxwell_bc_to_refined_rwg
    use fortfem_maxwell_cfie_regularized_3d, only: &
        assemble_maxwell_plane_wave_rhs_bc_3d, &
        assemble_maxwell_regularized_cfie_rwg_3d, &
        solve_maxwell_pec_regularized_cfie_rwg_3d, &
        solve_maxwell_pec_regularized_cfie_rwg_multiple_3d
    ! The offset-trace derivative actions optionally expose wave-number
    ! directions and cotangents in addition to the offset parameter.
    use fortfem_maxwell_torus_curved_rwg, only: &
        assemble_maxwell_torus_curved_efie_rwg_3d, &
        assemble_maxwell_torus_curved_efie_imaginary_rwg_3d, &
        assemble_maxwell_torus_efie_imaginary_impedance_jvp, &
        assemble_maxwell_torus_efie_imaginary_impedance_vjp, &
        assemble_maxwell_torus_efie_propagating_impedance_jvp, &
        assemble_maxwell_torus_efie_propagating_impedance_vjp, &
        assemble_maxwell_torus_efie_imaginary_decay_jvp, &
        assemble_maxwell_torus_efie_imaginary_decay_vjp, &
        assemble_maxwell_torus_efie_wave_number_jvp, &
        assemble_maxwell_torus_efie_wave_number_vjp, &
        assemble_maxwell_torus_curved_efie_bc_imaginary_3d, &
        assemble_maxwell_torus_curved_mfie_offset_trace_rwg_rbc_3d, &
        assemble_maxwell_torus_mfie_offset_jvp, &
        assemble_maxwell_torus_mfie_offset_vjp, &
        assemble_maxwell_torus_mfie_offset_geometry_jvp, &
        assemble_maxwell_torus_mfie_offset_geometry_vjp, &
        assemble_maxwell_torus_curved_mfie_rwg_rbc_3d, &
        assemble_maxwell_torus_curved_plane_wave_rhs_rwg_3d, &
        assemble_maxwell_torus_curved_plane_wave_rhs_rwg_3d_jvp, &
        assemble_maxwell_torus_curved_plane_wave_rhs_rwg_3d_vjp, &
        assemble_maxwell_torus_curved_plane_wave_rhs_bc_3d, &
        assemble_maxwell_torus_curved_potential_operators_rwg_3d, &
        assemble_maxwell_torus_curved_regularized_cfie_rwg_3d, &
        assemble_maxwell_torus_curved_rwg_mass_matrix, &
        assemble_maxwell_torus_curved_rwg_mass_matrix_jvp, &
        assemble_maxwell_torus_curved_rwg_mass_matrix_vjp, &
        assemble_maxwell_torus_curved_rwg_rbc_pairing, &
        assemble_maxwell_torus_curved_rwg_rbc_pairing_jvp, &
        assemble_maxwell_torus_curved_rwg_rbc_pairing_vjp, &
        evaluate_maxwell_torus_curved_far_field_rwg_3d, &
        evaluate_maxwell_torus_curved_far_field_rwg_3d_jvp, &
        evaluate_maxwell_torus_curved_far_field_rwg_3d_vjp, &
        evaluate_maxwell_torus_curved_magnetic_field_rwg_3d, &
        evaluate_maxwell_torus_curved_magnetic_field_rwg_3d_jvp, &
        evaluate_maxwell_torus_magnetic_geometry_jvp, &
        evaluate_maxwell_torus_magnetic_geometry_vjp, &
        evaluate_maxwell_torus_curved_magnetic_field_rwg_3d_vjp, &
        evaluate_maxwell_torus_curved_localized_rwg_basis, &
        evaluate_maxwell_torus_curved_localized_rwg_basis_jvp, &
        evaluate_maxwell_torus_curved_localized_rwg_basis_vjp, &
        evaluate_maxwell_torus_curved_rwg_basis, &
        evaluate_maxwell_torus_curved_rwg_basis_jvp, &
        evaluate_maxwell_torus_curved_rwg_basis_vjp, &
        assemble_maxwell_torus_curved_plane_wave_rhs_bc_3d_jvp, &
        assemble_maxwell_torus_curved_plane_wave_rhs_bc_3d_vjp, &
        integrate_maxwell_torus_curved_adjacent_rwg_pair_3d, &
        integrate_maxwell_torus_curved_coincident_rwg_pair_3d, &
        solve_maxwell_pec_torus_curved_regularized_cfie_rwg_3d, &
        solve_maxwell_pec_torus_curved_regularized_cfie_rwg_multiple_3d
    use fortfem_maxwell_curved_dtn, only: &
        apply_maxwell_trace_to_flux, apply_maxwell_trace_to_flux_jvp, &
        apply_maxwell_trace_to_flux_map, apply_maxwell_trace_to_flux_vjp, &
        apply_maxwell_weak_trace_reconstruction, &
        apply_maxwell_weak_trace_reconstruction_jvp, &
        apply_maxwell_weak_trace_reconstruction_vjp, &
        assemble_maxwell_trace_to_flux_map, &
        assemble_maxwell_trace_to_flux_map_jvp, &
        assemble_maxwell_trace_to_flux_map_vjp, &
        assemble_maxwell_weak_trace_reconstruction, &
        assemble_maxwell_torus_curved_dtn_rwg_3d
    use fortfem_toroidal_coordinates, only: cartesian_to_toroidal, &
        cartesian_to_toroidal_jvp, cartesian_to_toroidal_vjp, &
        toroidal_point_to_cartesian, toroidal_point_to_cartesian_jvp, &
        toroidal_point_to_cartesian_vjp, toroidal_vector_to_cartesian, &
        toroidal_vector_to_cartesian_jvp, toroidal_vector_to_cartesian_vjp
    use fortfem_laplace_representation_3d, only: &
        evaluate_laplace_representation_triangles_3d, &
        evaluate_laplace_representation_torus_curved_3d
    use fortfem_laplace_representation_ad_3d, only: &
        evaluate_laplace_representation_torus_curved_3d_jvp, &
        evaluate_laplace_representation_torus_curved_3d_vjp, &
        evaluate_laplace_representation_torus_curved_3d_geometry_jvp, &
        evaluate_laplace_representation_torus_curved_3d_geometry_vjp
    use fortfem_adaptive_surface_bem, only: &
        estimate_helmholtz_p0_two_level_residual_3d, &
        estimate_laplace_p0_two_level_residual_3d, mark_bem_dorfler, &
        refine_surface_mesh_marked
    use fortfem_laplace_torus_curved_bem_3d, only: &
        assemble_laplace_torus_curved_calderon_3d, &
        assemble_laplace_torus_curved_dtn_3d, &
        solve_laplace_bem_dtn_torus_curved_3d
    use fortfem_laplace_torus_curved_bem_ad_3d, only: &
        assemble_laplace_torus_curved_dtn_3d_geometry_jvp, &
        assemble_laplace_torus_curved_dtn_3d_geometry_vjp
    use fortfem_laplace_torus_curved_fem_bem_coupling_3d, only: &
        assemble_laplace_fem_bem_costabel_torus_curved_3d, &
        solve_laplace_fem_bem_costabel_torus_curved_3d
    use fortfem_laplace_galerkin_3d, only: &
        assemble_laplace_calderon_p1_p0_3d, &
        assemble_laplace_single_layer_p0_3d, &
        assemble_laplace_single_layer_p0_adaptive_3d, &
        solve_laplace_dirichlet_p0_3d
    use fortfem_laplace_bem_state_ad_3d, only: &
        solve_laplace_dirichlet_p0_3d_jvp, &
        solve_laplace_dirichlet_p0_3d_vjp
    use fortfem_laplace_panel_pair_3d, only: &
        assemble_laplace_single_layer_p0_3d_jvp, &
        assemble_laplace_single_layer_p0_3d_vjp, &
        integrate_laplace_single_layer_regular_panel_pair_p0_3d, &
        integrate_laplace_single_layer_regular_panel_pair_p0_3d_jvp, &
        integrate_laplace_single_layer_regular_panel_pair_p0_3d_vjp
    use fortfem_laplace_fem_bem_coupling_3d, only: &
        assemble_laplace_fem_bem_costabel_3d, &
        solve_laplace_fem_bem_costabel_3d, &
        solve_laplace_fem_bem_johnson_nedelec_3d
    use fortfem_laplace_hierarchical_3d, only: &
        apply_laplace_single_layer_p0_hierarchical_3d
    use fortfem_helmholtz_hierarchical_3d, only: &
        apply_helmholtz_cfie_p0_hierarchical_3d, &
        apply_helmholtz_single_layer_p0_hierarchical_3d, &
        solve_helmholtz_cfie_p0_hierarchical_3d, &
        solve_helmholtz_dirichlet_p0_hierarchical_3d
    use fortfem_helmholtz_representation_3d, only: &
        evaluate_helmholtz_representation_triangles_3d, &
        evaluate_helmholtz_representation_torus_curved_3d
    use fortfem_helmholtz_representation_ad_3d, only: &
        evaluate_helmholtz_representation_torus_curved_3d_geometry_jvp, &
        evaluate_helmholtz_representation_torus_curved_3d_geometry_vjp
    use fortfem_helmholtz_torus_curved_bem_3d, only: &
        assemble_helmholtz_torus_curved_calderon_3d, &
        assemble_helmholtz_torus_curved_dtn_3d, &
        solve_helmholtz_bem_dtn_torus_curved_3d
    use fortfem_helmholtz_torus_curved_bem_ad_3d, only: &
        assemble_helmholtz_torus_curved_dtn_3d_geometry_jvp, &
        assemble_helmholtz_torus_curved_dtn_3d_geometry_vjp
    use fortfem_helmholtz_torus_curved_fem_bem_coupling_3d, only: &
        assemble_helmholtz_fem_bem_costabel_torus_curved_3d, &
        solve_helmholtz_fem_bem_costabel_torus_curved_3d
    use fortfem_helmholtz_galerkin_3d, only: &
        assemble_helmholtz_calderon_p1_p0_3d, &
        assemble_helmholtz_double_layer_p0_3d, &
        assemble_helmholtz_single_layer_p0_adaptive_3d, &
        assemble_helmholtz_single_layer_p0_3d, &
        evaluate_helmholtz_cfie_p0_3d, solve_helmholtz_cfie_p0_3d, &
        solve_helmholtz_dirichlet_p0_3d
    use fortfem_helmholtz_bem_state_ad_3d, only: &
        solve_helmholtz_dirichlet_p0_3d_jvp, &
        solve_helmholtz_dirichlet_p0_3d_vjp
    use fortfem_helmholtz_panel_pair_3d, only: &
        assemble_helmholtz_single_layer_p0_3d_jvp, &
        assemble_helmholtz_single_layer_p0_3d_vjp, &
        integrate_helmholtz_single_layer_regular_panel_pair_p0_3d, &
        integrate_helmholtz_single_layer_regular_panel_pair_p0_3d_jvp, &
        integrate_helmholtz_single_layer_regular_panel_pair_p0_3d_vjp
    use fortfem_helmholtz_fem_bem_coupling_3d, only: &
        assemble_helmholtz_fem_bem_costabel_3d, &
        solve_helmholtz_fem_bem_costabel_3d
    use fortfem_maxwell_rwg_surface, only: &
        assemble_maxwell_rwg_mass_matrix, &
        build_maxwell_rwg_surface_space, evaluate_maxwell_rwg_basis, &
        map_maxwell_rwg_to_tetra_nedelec_edges
    use fortfem_maxwell_surface_rt, only: &
        assemble_maxwell_surface_rt_mass_matrix, &
        build_maxwell_surface_rt_dof_map, evaluate_maxwell_surface_rt_basis, &
        evaluate_maxwell_surface_rt_global_basis
    use fortfem_maxwell_surface_rt_efie_3d, only: &
        assemble_maxwell_surface_rt_efie_3d
    use fortfem_bspline_feec, only: &
        build_bspline_derivative_matrix, build_bspline_feec_2d_operators, &
        build_bspline_feec_3d_operators, evaluate_bspline_basis, &
        evaluate_nurbs_surface_geometry, evaluate_nurbs_surface_geometry_jvp, &
        evaluate_nurbs_surface_geometry_vjp, map_isogeometric_h1_gradient, &
        evaluate_nurbs_volume_geometry, evaluate_nurbs_volume_geometry_jvp, &
        evaluate_nurbs_volume_geometry_vjp, &
        map_isogeometric_hcurl, map_isogeometric_hdiv, map_isogeometric_l2
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
        apply_bspline_jorek_flux_rhs, &
        apply_bspline_jorek_flux_jvp, &
        apply_bspline_jorek_thermodynamic_rhs, &
        apply_bspline_jorek_thermodynamic_jvp, &
        apply_bspline_jorek_density_rhs, &
        apply_bspline_jorek_density_jvp, &
        project_bspline_toroidal_product, &
        advance_bspline_jorek_poloidal_flux_midpoint, &
        advance_bspline_jorek_poloidal_flux_midpoint_steps, &
        apply_bspline_toroidal_poloidal_bracket, &
        apply_bspline_toroidal_poloidal_bracket_jvp, &
        apply_toroidal_fourier_derivative, &
        build_bspline_feec_2d_operators_csc, scalar_weight_2d, tensor_weight_2d
    use fortfem_assembly_bspline_3d, only: &
        assemble_bspline_h1_operator_3d_csc, &
        assemble_bspline_hcurl_operator_3d_csc, &
        assemble_bspline_hdiv_operator_3d_csc, &
        assemble_bspline_hdiv_l2_divergence_3d_csc, &
        assemble_bspline_l2_hdiv_adjoint_divergence_3d_csc, &
        assemble_bspline_l2_mass_3d_csc, &
        build_bspline_feec_3d_operators_csc
    use fortfem_maxwell_efie_rwg_3d, only: &
        assemble_maxwell_efie_rwg_3d, &
        assemble_maxwell_plane_wave_rhs_rwg_3d, &
        assemble_maxwell_rwg_potential_operators_3d, &
        evaluate_maxwell_efie_far_field_rwg_3d, &
        evaluate_maxwell_efie_field_rwg_3d, solve_maxwell_pec_efie_rwg_3d
    use fortfem_maxwell_fem_bem_coupling_3d, only: &
        assemble_maxwell_fem_bem_boundary_matrix_3d, &
        assemble_maxwell_fem_bem_system_3d, &
        assemble_maxwell_fem_bem_torus_curved_system_3d, &
        assemble_maxwell_rwg_nedelec_coupling_3d, &
        assemble_maxwell_torus_curved_rwg_nedelec_coupling_3d, &
        solve_maxwell_fem_bem_system_3d, &
        solve_maxwell_fem_bem_torus_curved_system_3d
    use fortfem_maxwell_fem_bem_state_ad, only: &
        solve_maxwell_fem_bem_linear_state, &
        solve_maxwell_fem_bem_linear_state_jvp, &
        solve_maxwell_fem_bem_linear_state_vjp
    use fortfem_planar_helmholtz_dtn, only: &
        apply_planar_helmholtz_dtn, apply_planar_helmholtz_dtn_jvp, &
        apply_planar_helmholtz_dtn_vjp, assemble_planar_helmholtz_dtn_form, &
        assemble_planar_helmholtz_dtn_form_jvp, &
        assemble_planar_helmholtz_dtn_form_vjp
    use fortfem_planar_maxwell_dtn, only: apply_planar_maxwell_dtn, &
        apply_planar_maxwell_dtn_jvp, apply_planar_maxwell_dtn_vjp, &
        assemble_planar_maxwell_dtn_form, &
        assemble_planar_maxwell_dtn_form_jvp, &
        assemble_planar_maxwell_dtn_form_vjp
    use fortfem_planar_maxwell_dtn_system, only: &
        solve_planar_maxwell_dtn_system, &
        solve_planar_maxwell_dtn_system_jvp, &
        solve_planar_maxwell_dtn_system_vjp
    use fortfem_planar_acoustic_displacement_dtn, only: &
        apply_planar_acoustic_displacement_dtn, &
        assemble_planar_acoustic_displacement_dtn_form
    use fortfem_elasticity_planar_acoustic_dtn_2d, only: &
        solve_elasticity_planar_acoustic_dtn_p1
    use fortfem_elasticity_curved_acoustic_ntd_2d, only: &
        solve_elasticity_curved_acoustic_ntd_p1
    use fortfem_scalar_helmholtz_planar_dtn_2d, only: &
        solve_scalar_helmholtz_planar_dtn_p1, &
        solve_scalar_helmholtz_planar_dtn_p1_jvp, &
        solve_scalar_helmholtz_planar_dtn_p1_vjp
    use fortfem_circular_dtn_2d, only: apply_circular_helmholtz_dtn, &
        circular_helmholtz_dtn_eigenvalue
    use fortfem_circular_dtn_2d_ad, only: &
        apply_circular_helmholtz_dtn_jvp, apply_circular_helmholtz_dtn_vjp, &
        circular_helmholtz_dtn_eigenvalue_jvp, &
        circular_helmholtz_dtn_eigenvalue_vjp
    use fortfem_spherical_helmholtz_dtn, only: &
        apply_spherical_helmholtz_dtn, apply_spherical_helmholtz_dtn_jvp, &
        apply_spherical_helmholtz_dtn_vjp, &
        spherical_helmholtz_dtn_eigenvalue, &
        spherical_helmholtz_dtn_eigenvalue_jvp, &
        spherical_helmholtz_dtn_eigenvalue_vjp
    use fortfem_spherical_maxwell_dtn, only: &
        apply_spherical_maxwell_dtn, apply_spherical_maxwell_dtn_jvp, &
        apply_spherical_maxwell_dtn_vjp, spherical_maxwell_dtn_eigenvalues, &
        spherical_maxwell_dtn_eigenvalues_jvp, &
        spherical_maxwell_dtn_eigenvalues_vjp
    use fortfem_cartesian_helmholtz_pml, only: &
        cartesian_curl_curl_pml_coefficients, &
        cartesian_curl_curl_pml_coefficients_jvp, &
        cartesian_curl_curl_pml_coefficients_vjp, &
        cartesian_scalar_helmholtz_pml_coefficients, &
        cartesian_scalar_helmholtz_pml_coefficients_jvp, &
        cartesian_scalar_helmholtz_pml_coefficients_vjp
    use fortfem_curvilinear_helmholtz_pml, only: &
        curvilinear_curl_curl_pml_coefficients, &
        curvilinear_curl_curl_pml_coefficients_jvp, &
        curvilinear_curl_curl_pml_coefficients_vjp, &
        curvilinear_scalar_helmholtz_pml_coefficients, &
        curvilinear_scalar_helmholtz_pml_coefficients_jvp, &
        curvilinear_scalar_helmholtz_pml_coefficients_vjp
    use fortfem_scalar_helmholtz_pml_slab_1d, only: &
        solve_scalar_helmholtz_pml_slab_1d
    use fortfem_scalar_helmholtz_pml_slab_1d_ad, only: &
        solve_scalar_helmholtz_pml_slab_1d_jvp, &
        solve_scalar_helmholtz_pml_slab_1d_vjp
    use fortfem_scalar_helmholtz_pml_2d, only: &
        solve_scalar_helmholtz_pml_p1_2d
    use fortfem_scalar_helmholtz_pml_2d_ad, only: &
        solve_scalar_helmholtz_pml_p1_2d_jvp, &
        solve_scalar_helmholtz_pml_p1_2d_vjp
    use fortfem_scalar_helmholtz_pml_3d, only: &
        solve_scalar_helmholtz_pml_p1_3d
    use fortfem_scalar_helmholtz_pml_3d_ad, only: &
        solve_scalar_helmholtz_pml_p1_3d_jvp, &
        solve_scalar_helmholtz_pml_p1_3d_vjp
    use fortfem_toroidal_poisson_dtn, only: &
        evaluate_toroidal_harmonic_p, evaluate_toroidal_ampere_field_p, &
        toroidal_poisson_exterior_dtn_p
    use fortfem_toroidal_poisson_dtn_ad, only: &
        evaluate_toroidal_harmonic_p_jvp, evaluate_toroidal_harmonic_p_vjp, &
        evaluate_toroidal_ampere_field_p_jvp, &
        evaluate_toroidal_ampere_field_p_vjp, &
        toroidal_poisson_exterior_dtn_p_jvp, &
        toroidal_poisson_exterior_dtn_p_vjp
    use fortfem_laplace_boundary_operators_2d, only: &
        assemble_laplace_adjoint_double_layer_constant, &
        assemble_laplace_double_layer_constant, &
        assemble_laplace_double_layer_mixed_linear, &
        assemble_laplace_hypersingular_linear, &
        assemble_laplace_single_layer_constant
    use fortfem_laplace_symmetric_coupling_2d, only: &
        assemble_helmholtz_symmetric_coupling_p1_p0, &
        assemble_laplace_symmetric_coupling_p1_p0, &
        solve_helmholtz_symmetric_coupling_p1_p0, &
        solve_laplace_symmetric_coupling_p1_p0
    use fortfem_helmholtz_boundary_operators_2d, only: &
        assemble_helmholtz_adjoint_double_layer_constant, &
        assemble_helmholtz_double_layer_constant, &
        assemble_helmholtz_double_layer_mixed_linear, &
        assemble_helmholtz_hypersingular_linear, &
        assemble_helmholtz_single_layer_constant, &
        assemble_helmholtz_single_layer_linear
    use fortfem_curved_acoustic_displacement_ntd_2d, only: &
        apply_curved_acoustic_displacement_ntd_2d, &
        assemble_curved_acoustic_displacement_ntd_form_2d
    use fortfem_helmholtz_exterior_2d, only: &
        evaluate_helmholtz_combined_potential_adaptive_constant, &
        evaluate_helmholtz_combined_potential_constant, &
        solve_helmholtz_cfie_constant
    use fortfem_triangle_duffy_quadrature, only: triangle_duffy_quadrature
    use fortfem_tetra_duffy_quadrature, only: tetra_duffy_quadrature
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
    use fortfem_triangle_nedelec_second_kind, only: &
        assignment(=), evaluate_triangle_nedelec_second_kind, &
        evaluate_triangle_nedelec_second_kind_jvp, &
        evaluate_triangle_nedelec_second_kind_vjp, &
        initialize_triangle_nedelec_second_kind, &
        triangle_nedelec_second_kind_dof_count, &
        triangle_nedelec_second_kind_t
    use fortfem_tetra_nedelec_first_order, only: &
        evaluate_tetra_nedelec_first_order
    use fortfem_tetra_nedelec_arbitrary_order, only: &
        assignment(=), evaluate_tetra_nedelec_first_kind, &
        evaluate_tetra_nedelec_first_kind_jvp, &
        evaluate_tetra_nedelec_first_kind_vjp, &
        initialize_tetra_nedelec_first_kind, tetra_nedelec_dof_count, &
        tetra_nedelec_first_kind_t
    use fortfem_tetra_discontinuous_arbitrary_order, only: &
        assignment(=), evaluate_tetra_discontinuous, &
        evaluate_tetra_discontinuous_jvp, evaluate_tetra_discontinuous_vjp, &
        initialize_tetra_discontinuous, tetra_discontinuous_dof_count, &
        tetra_discontinuous_t
    use fortfem_tetra_discontinuous_projection, only: &
        project_physical_tetra_discontinuous, &
        project_sampled_physical_tetra_discontinuous, &
        project_sampled_physical_tetra_discontinuous_jvp, &
        project_sampled_physical_tetra_discontinuous_vjp
    use fortfem_tetra_lagrange_arbitrary_order, only: &
        assignment(=), evaluate_tetra_lagrange, evaluate_tetra_lagrange_jvp, &
        evaluate_tetra_lagrange_vjp, initialize_tetra_lagrange, &
        tetra_lagrange_barycentric_indices, tetra_lagrange_dof_count, &
        tetra_lagrange_nodes, tetra_lagrange_t
    use fortfem_tetra_lagrange_global_dof_map, only: &
        build_tetra_lagrange_dof_map
    use fortfem_tetra_feec_operators, only: &
        build_tetra_discrete_curl, build_tetra_discrete_gradient
    use fortfem_tetra_rt_arbitrary_order, only: &
        assignment(=), evaluate_tetra_rt, evaluate_tetra_rt_jvp, &
        evaluate_tetra_rt_vjp, initialize_tetra_rt, &
        tetra_rt_dof_count, tetra_rt_t
    use fortfem_tetra_rt_global_dof_map, only: &
        build_tetra_rt_basis_transform, build_tetra_rt_dof_map
    use fortfem_tetra_piola_maps, only: map_tetra_nedelec_covariant, &
        map_tetra_nedelec_covariant_jvp, map_tetra_nedelec_covariant_vjp, &
        map_tetra_rt_contravariant, map_tetra_rt_contravariant_jvp, &
        map_tetra_rt_contravariant_vjp
    use fortfem_tetra_edge_dof_map, only: build_tetra_edge_dof_map
    use fortfem_tetra_nedelec_global_dof_map, only: &
        build_tetra_nedelec_basis_transform, build_tetra_nedelec_dof_map
    use fortfem_tetra_nedelec_interpolation, only: &
        interpolate_physical_tetra_nedelec, &
        interpolate_reference_tetra_nedelec, &
        interpolate_sampled_physical_tetra_nedelec, &
        interpolate_sampled_physical_tetra_nedelec_jvp, &
        interpolate_sampled_physical_tetra_nedelec_vjp, &
        tetra_nedelec_interpolation_points
    use fortfem_tetra_vector_samples, only: &
        assignment(=), tetra_vector_sample_gradients_t, &
        tetra_vector_samples_t, zero_tetra_vector_samples_like
    use fortfem_tetra_rt_interpolation, only: &
        interpolate_physical_tetra_rt, interpolate_reference_tetra_rt, &
        interpolate_sampled_physical_tetra_rt, &
        interpolate_sampled_physical_tetra_rt_jvp, &
        interpolate_sampled_physical_tetra_rt_vjp, &
        tetra_rt_interpolation_points
    use fortfem_triangle_rt_arbitrary_order, only: &
        assignment(=), evaluate_triangle_raviart_thomas, &
        evaluate_triangle_raviart_thomas_jvp, &
        evaluate_triangle_raviart_thomas_vjp, &
        initialize_triangle_raviart_thomas, triangle_rt_basis_t, &
        triangle_rt_dof_count
    use fortfem_triangle_bdm_arbitrary_order, only: &
        assignment(=), evaluate_triangle_bdm, evaluate_triangle_bdm_jvp, &
        evaluate_triangle_bdm_vjp, initialize_triangle_bdm, &
        triangle_bdm_basis_t, triangle_bdm_dof_count
    use fortfem_edge_moment_orientation, only: apply_edge_moment_orientation
    use fortfem_triangle_piola_maps, only: &
        map_triangle_nedelec_covariant, &
        map_triangle_nedelec_covariant_jvp, &
        map_triangle_nedelec_covariant_vjp, &
        map_triangle_rt_contravariant, &
        map_triangle_rt_contravariant_jvp, &
        map_triangle_rt_contravariant_vjp
    use fortfem_triangle_feec_operators, only: &
        build_triangle_discrete_gradient
    use fortfem_triangle_global_dof_map, only: &
        build_triangle_full_vector_dof_map, build_triangle_trimmed_dof_map
    use fortfem_triangle_discontinuous_dof_map, only: &
        build_triangle_discontinuous_dof_map
    use fortfem_triangle_vector_interpolation, only: &
        evaluate_triangle_bdm_interpolant, &
        evaluate_triangle_bdm_interpolant_at_point, &
        evaluate_triangle_bdm_interpolant_at_point_jvp, &
        evaluate_triangle_bdm_interpolant_at_point_vjp, &
        evaluate_triangle_bdm_interpolant_jvp, &
        evaluate_triangle_bdm_interpolant_vjp, &
        evaluate_triangle_nedelec_interpolant, &
        evaluate_triangle_nedelec_interpolant_at_point, &
        evaluate_triangle_nedelec_interpolant_at_point_jvp, &
        evaluate_triangle_nedelec_interpolant_at_point_vjp, &
        evaluate_triangle_nedelec_interpolant_jvp, &
        evaluate_triangle_nedelec_interpolant_vjp, &
        evaluate_triangle_nedelec_second_kind_interpolant, &
        evaluate_triangle_nedelec_second_kind_interpolant_jvp, &
        evaluate_triangle_nedelec_second_kind_interpolant_vjp, &
        evaluate_triangle_nedelec_second_interpolant_at_point, &
        evaluate_triangle_nedelec_second_interpolant_at_point_jvp, &
        evaluate_triangle_nedelec_second_interpolant_at_point_vjp, &
        evaluate_triangle_rt_interpolant, &
        evaluate_triangle_rt_interpolant_at_point, &
        evaluate_triangle_rt_interpolant_at_point_jvp, &
        evaluate_triangle_rt_interpolant_at_point_vjp, &
        evaluate_triangle_rt_interpolant_jvp, &
        evaluate_triangle_rt_interpolant_vjp, interpolate_triangle_bdm, &
        interpolate_triangle_nedelec, &
        interpolate_triangle_nedelec_second_kind, interpolate_triangle_rt
    use fortfem_edge_interpolation_2d, only: &
        interpolate_axisymmetric_rt_edge_dofs, &
        interpolate_nedelec_edge_dofs, interpolate_rt_edge_dofs
    use fortfem_rt_field_2d, only: evaluate_rt_field_2d, &
        reconstruct_axisymmetric_fourier_toroidal, rt_l2_norm
    use fortfem_sparse_direct, only: sparse_direct_factor_t, &
        sparse_direct_factor_adjoint_csc, sparse_direct_factor_csc, &
        sparse_direct_factor_transpose_csc, sparse_direct_free, &
        sparse_direct_solve_csc, sparse_direct_solve_constrained, &
        sparse_direct_solve_constrained_jvp, &
        sparse_direct_solve_constrained_vjp, sparse_direct_solve_factored, &
        sparse_direct_solve_factored_jvp, sparse_direct_solve_factored_vjp, &
        sparse_direct_solve_tree_cotree, &
        sparse_direct_solve_tree_cotree_jvp, &
        sparse_direct_solve_tree_cotree_vjp
    use fortfem_mixed_poisson_2d, only: &
        solve_mixed_poisson_rt, solve_mixed_poisson_rt0, &
        solve_symbolic_mixed_poisson_rt
    use fortfem_magnetic_box_3d, only: solve_magnetic_box_3d
    ! Public arbitrary-order H(curl) solve, including optional homogeneous PEC.
    use fortfem_tetra_nedelec_solver_3d, only: &
        solve_tetra_nedelec_curl_mass, solve_tetra_nedelec_pml, &
        solve_tetra_nedelec_curvilinear_pml, &
        solve_tetra_nedelec_weighted_curl_mass
    use fortfem_tetra_nedelec_pml_state_3d, only: &
        solve_tetra_nedelec_pml_jvp, solve_tetra_nedelec_pml_vjp, &
        solve_tetra_nedelec_curvilinear_pml_jvp, &
        solve_tetra_nedelec_curvilinear_pml_vjp
    ! Public arbitrary-order H(div) solve, including optional zero normal trace.
    use fortfem_tetra_rt_solver_3d, only: solve_tetra_rt_div_mass
    use fortfem_tetra_mixed_poisson_3d, only: &
        assemble_tetra_dg_source_load_samples, &
        assemble_tetra_dg_source_load_samples_jvp, &
        assemble_tetra_dg_source_load_samples_vjp, &
        solve_symbolic_tetra_mixed_poisson_rt
    use fortfem_tetra_lagrange_solver_3d, only: &
        assignment(=), evaluate_tetra_lagrange_solution, &
        evaluate_tetra_lagrange_solution_at_point, &
        evaluate_tetra_lagrange_solution_at_point_jvp, &
        evaluate_tetra_lagrange_solution_at_point_vjp, &
        evaluate_tetra_lagrange_solution_jvp, &
        evaluate_tetra_lagrange_solution_prepared, &
        evaluate_tetra_lagrange_solution_vjp, &
        initialize_tetra_lagrange_solution_evaluator, &
        solve_tetra_lagrange_diffusion_reaction, &
        solve_tetra_lagrange_poisson, tetra_lagrange_solution_evaluator_t
    use fortfem_advanced_solvers, only: solver_options_t, solver_stats_t,   &
        solver_options, cg_solve, pcg_solve, bicgstab_solve, gmres_solve,   &
        cg_solve_jvp, cg_solve_vjp, pcg_solve_jvp, pcg_solve_vjp,           &
        bicgstab_solve_jvp, bicgstab_solve_vjp, gmres_solve_jvp,            &
        gmres_solve_vjp, jacobi_preconditioner, ilu_preconditioner, &
        ichol_preconditioner, ichol_controlled_preconditioner
    use fortfem_solver_resource_budget, only: &
        solver_resource_budget_t, initialize_solver_resource_budget, &
        validate_solver_resource_budget, evaluate_solver_resource_usage
    use fortfem_incomplete_cholesky, only: apply_incomplete_cholesky, &
        build_incomplete_cholesky, incomplete_cholesky_factor_t
    use fortfem_sparse_incomplete_cholesky, only: &
        apply_sparse_incomplete_cholesky, &
        apply_sparse_incomplete_cholesky_jvp, &
        apply_sparse_incomplete_cholesky_vjp, &
        build_sparse_ichol, build_sparse_ichol_row, &
        build_sparse_incomplete_cholesky, &
        sparse_incomplete_cholesky_factor_t
    use fortfem_sparse_incomplete_lu, only: apply_sparse_incomplete_lu, &
        apply_sparse_incomplete_lu_jvp, apply_sparse_incomplete_lu_vjp, &
        build_sparse_incomplete_lu, sparse_incomplete_lu_factor_t
    use fortfem_sparse_ilut, only: apply_sparse_ilut, apply_sparse_ilut_jvp, &
        apply_sparse_ilut_vjp, build_sparse_ilut, build_sparse_ilut_row, &
        sparse_ilut_factor_t
    implicit none

    private

    ! Public types
    public :: cell_complex_t
    public :: initialize_cell_complex
    public :: validate_cell_complex
    public :: cell_complex_euler_characteristic
    public :: cell_complex_betti_numbers
    public :: cell_complex_cycle_basis
    public :: cell_complex_homology_cycle_basis
    public :: cell_complex_cocycle_basis
    public :: cell_complex_cohomology_cocycle_basis
    public :: cell_complex_harmonic_one_forms
    public :: quotient_cell_complex
    public :: partition_layout_t
    public :: initialize_partition_layout
    public :: validate_partition_layout
    public :: partition_layout_maps
    public :: partition_layout_global_count
    public :: partition_layout_owned_count
    public :: partition_layout_ghost_count
    public :: assemble_partitioned_sum
    public :: assemble_partitioned_sum_jvp
    public :: assemble_partitioned_sum_vjp
    public :: distributed_trace_layout_t
    public :: initialize_distributed_trace_layout
    public :: validate_distributed_trace_layout
    public :: distributed_trace_layout_partition_count
    public :: distributed_trace_layout_global_count
    public :: distributed_trace_layout_local_count
    public :: distributed_trace_layout_component_count
    public :: assemble_distributed_trace_reduction
    public :: assemble_distributed_trace_reduction_jvp
    public :: assemble_distributed_trace_reduction_vjp
    public :: physical_trace_ownership_t
    public :: initialize_physical_trace_ownership
    public :: validate_physical_trace_ownership
    public :: physical_trace_ownership_maps
    public :: physical_trace_ownership_dimension
    public :: physical_trace_ownership_point_count
    public :: physical_trace_ownership_rank
    public :: compare_physical_trace_coordinates
    public :: physical_trace_reconciliation_t
    public :: initialize_physical_trace_reconciliation
    public :: validate_physical_trace_reconciliation
    public :: physical_trace_reconciliation_maps
    public :: reconcile_physical_trace_values
    public :: reconcile_physical_trace_values_jvp
    public :: reconcile_physical_trace_values_vjp
    public :: physical_trace_owner_selection_t
    public :: initialize_physical_trace_owner_selection
    public :: validate_physical_trace_owner_selection
    public :: physical_trace_owner_selection_maps
    public :: gather_physical_trace_values
    public :: gather_physical_trace_values_jvp
    public :: gather_physical_trace_values_vjp
    public :: assemble_pseudo_arclength_residual
    public :: assemble_pseudo_arclength_residual_jvp
    public :: assemble_pseudo_arclength_residual_vjp
    public :: assemble_deflated_residual
    public :: assemble_deflated_residual_jvp
    public :: assemble_deflated_residual_vjp
    public :: normalize_pseudo_arclength_tangent
    public :: normalize_pseudo_arclength_tangent_jvp
    public :: normalize_pseudo_arclength_tangent_vjp
    public :: evaluate_residual_merit
    public :: evaluate_residual_merit_jvp
    public :: evaluate_residual_merit_vjp
    public :: CONTINUATION_EVENT_NONE
    public :: CONTINUATION_EVENT_SIGN_CROSSING
    public :: CONTINUATION_EVENT_NEAR_ZERO
    public :: classify_continuation_event
    public :: assemble_pseudo_transient_residual
    public :: assemble_pseudo_transient_residual_jvp
    public :: assemble_pseudo_transient_residual_vjp
    public :: assemble_eulerian_nonnested_residual
    public :: assemble_eulerian_nonnested_residual_jvp
    public :: assemble_eulerian_nonnested_residual_vjp
    public :: evaluate_step_reduction
    public :: evaluate_step_reduction_jvp
    public :: evaluate_step_reduction_vjp
    public :: normalize_harmonic_one_forms
    public :: normalize_harmonic_one_forms_jvp
    public :: normalize_harmonic_one_forms_vjp
    public :: cell_identification_t
    public :: initialize_cell_identification
    public :: validate_cell_identification
    public :: cell_identification_classes
    public :: identify_boundary_matrix
    public :: region_interface_graph_t
    public :: initialize_region_interface_graph
    public :: validate_region_interface_graph
    public :: region_interface_graph_incidence
    public :: region_interface_graph_components
    public :: region_interface_graph_cycle_basis
    public :: internal_manifold_graph_t
    public :: initialize_internal_manifold_graph
    public :: validate_internal_manifold_graph
    public :: internal_manifold_graph_region_incidence
    public :: internal_manifold_graph_junction_incidence
    public :: internal_manifold_graph_closed
    public :: internal_manifold_graph_components
    public :: boundary_region_graph_t
    public :: initialize_boundary_region_graph
    public :: validate_boundary_region_graph
    public :: boundary_region_graph_incidence
    public :: boundary_region_graph_components
    public :: boundary_region_graph_cycle_basis
    public :: boundary_region_graph_interface_samples
    public :: boundary_region_graph_interface_metadata
    public :: assemble_interface_surface_current
    public :: assemble_interface_surface_current_jvp
    public :: assemble_interface_surface_current_vjp
    public :: surface_current_potential_metadata_t
    public :: initialize_surface_current_potential_metadata
    public :: validate_surface_current_potential_metadata
    public :: assemble_surface_current_potential
    public :: assemble_surface_current_potential_jvp
    public :: assemble_surface_current_potential_vjp
    public :: assemble_surface_current_junction_balance
    public :: assemble_surface_current_junction_balance_jvp
    public :: assemble_surface_current_junction_balance_vjp
    public :: assemble_surface_current_loop_constraints
    public :: assemble_surface_current_loop_constraints_jvp
    public :: assemble_surface_current_loop_constraints_vjp
    public :: assemble_period_constraints
    public :: assemble_period_constraints_jvp
    public :: assemble_period_constraints_vjp
    public :: assemble_surface_edge_flux_balance
    public :: assemble_surface_edge_flux_balance_jvp
    public :: assemble_surface_edge_flux_balance_vjp
    public :: assemble_surface_edge_flux
    public :: assemble_surface_edge_flux_jvp
    public :: assemble_surface_edge_flux_vjp
    public :: assemble_surface_current_trace_residual
    public :: assemble_surface_current_trace_residual_jvp
    public :: assemble_surface_current_trace_residual_vjp
    public :: evaluate_regularized_surface_current_layer
    public :: evaluate_regularized_surface_current_layer_jvp
    public :: evaluate_regularized_surface_current_layer_vjp
    public :: evaluate_regularized_surface_current_integral
    public :: compare_sheet_current_representations
    public :: compare_sheet_current_surface_representations
    public :: compare_sheet_current_surface_representations_jvp
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
    public :: fourier_mode_registry_t
    public :: initialize_fourier_mode_registry
    public :: validate_fourier_mode_registry
    public :: find_fourier_mode
    public :: fourier_mode_conjugate_index
    public :: fourier_mode_triad_closed
    public :: build_fourier_mode_triad_map
    public :: build_fourier_mode_padded_registry
    public :: build_fourier_mode_closure_registry
    public :: evaluate_fourier_mode
    public :: evaluate_fourier_mode_jvp
    public :: evaluate_fourier_mode_vjp
    public :: evaluate_nested_surface_geometry
    public :: evaluate_nested_surface_geometry_jvp
    public :: evaluate_nested_surface_geometry_vjp
    public :: evaluate_nested_surface_geometry_coordinate_jvp
    public :: evaluate_nested_surface_geometry_coordinate_vjp
    public :: evaluate_nested_geometry_differential_jet
    public :: evaluate_nested_geometry_differential_jet_jvp
    public :: evaluate_nested_geometry_differential_jet_vjp
    public :: evaluate_nested_geometry_polynomial_jet
    public :: validate_nested_geometry_axis_flags
    public :: assemble_fourier_vector_product
    public :: assemble_fourier_vector_product_jvp
    public :: assemble_fourier_vector_product_vjp
    public :: apply_fourier_bilinear_product
    public :: apply_fourier_bilinear_product_jvp
    public :: apply_fourier_bilinear_product_vjp
    public :: apply_fourier_mode_linear_operator
    public :: apply_fourier_mode_linear_operator_jvp
    public :: apply_fourier_mode_linear_operator_vjp
    public :: evaluate_fourier_mode_expansion
    public :: evaluate_fourier_mode_expansion_jvp
    public :: evaluate_fourier_mode_expansion_vjp
    public :: evaluate_fourier_mode_expansion_hessian
    public :: evaluate_fourier_mode_expansion_hvp
    public :: assemble_fourier_mode_energy
    public :: assemble_fourier_mode_energy_jvp
    public :: assemble_fourier_mode_energy_vjp
    public :: fourier_coefficients_conjugate_symmetric
    public :: AXIS_RADIAL_PARITY_EVEN
    public :: AXIS_RADIAL_PARITY_ODD
    public :: axis_regular_mode_record_t
    public :: axis_regular_mode_table_t
    public :: axis_regular_mode_requirements
    public :: build_axis_regular_mode_table
    public :: validate_axis_regular_mode_table
    public :: evaluate_axis_regular_radial_basis
    public :: evaluate_axis_regular_radial_basis_jvp
    public :: evaluate_axis_regular_radial_basis_vjp
    public :: FOURIER_ZERNIKE_PARITY_EVEN
    public :: FOURIER_ZERNIKE_PARITY_ODD
    public :: fourier_zernike_mode_t
    public :: fourier_zernike_basis_t
    public :: build_fourier_zernike_basis
    public :: validate_fourier_zernike_basis
    public :: fourier_zernike_mode_requirements
    public :: evaluate_fourier_zernike_radial
    public :: COLLOCATION_GRID_LINEAR
    public :: COLLOCATION_GRID_QUADRATURE
    public :: COLLOCATION_GRID_CONCENTRIC
    public :: collocation_grid_t
    public :: initialize_collocation_grid
    public :: validate_collocation_grid
    public :: collocation_grid_metadata
    public :: collocation_grid_flat_index
    public :: collocation_grid_unflatten_index
    public :: collocation_grid_chunk_bounds
    public :: collocation_grid_point_count
    public :: direct_fourier_plan_t
    public :: initialize_direct_fourier_plan
    public :: validate_direct_fourier_plan
    public :: direct_fourier_plan_metadata
    public :: direct_fourier_plan_chunk_bounds
    public :: direct_fourier_forward
    public :: direct_fourier_adjoint
    public :: direct_fourier_plan_sample_count
    public :: direct_fourier_plan_mode_count
    public :: assemble_normal_traction_jump
    public :: assemble_normal_traction_jump_jvp
    public :: assemble_normal_traction_jump_vjp
    public :: assemble_traction_jump
    public :: assemble_traction_jump_jvp
    public :: assemble_traction_jump_vjp
    public :: assemble_total_pressure_jump
    public :: assemble_total_pressure_jump_jvp
    public :: assemble_total_pressure_jump_vjp
    public :: evaluate_shifted_heaviside_enrichment
    public :: evaluate_shifted_heaviside_enrichment_jvp
    public :: evaluate_shifted_heaviside_enrichment_vjp
    public :: evaluate_shifted_enriched_basis
    public :: evaluate_shifted_enriched_basis_jvp
    public :: evaluate_shifted_enriched_basis_vjp
    public :: evaluate_shifted_enriched_space
    public :: evaluate_shifted_enriched_space_jvp
    public :: evaluate_shifted_enriched_space_vjp
    public :: evaluate_shifted_vector_enriched_space
    public :: evaluate_shifted_vector_enriched_space_jvp
    public :: evaluate_shifted_vector_enriched_space_vjp
    public :: evaluate_shifted_vector_enriched_basis
    public :: evaluate_shifted_vector_enriched_basis_jvp
    public :: evaluate_shifted_vector_enriched_basis_vjp
    public :: evaluate_vector_enrichment_differential_3d
    public :: evaluate_vector_enrichment_differential_3d_jvp
    public :: evaluate_vector_enrichment_differential_3d_vjp
    public :: assemble_tensor_volume_work
    public :: assemble_tensor_volume_work_jvp
    public :: assemble_tensor_volume_work_vjp
    public :: assemble_force_balance_residual
    public :: assemble_force_balance_residual_jvp
    public :: assemble_force_balance_residual_vjp
    public :: evaluate_force_balance_objective
    public :: evaluate_force_balance_objective_jvp
    public :: evaluate_force_balance_objective_vjp
    public :: evaluate_force_balance_objective_hvp
    public :: evaluate_force_balance_product
    public :: evaluate_force_balance_product_jvp
    public :: evaluate_force_balance_product_vjp
    public :: evaluate_flux_surface_average
    public :: evaluate_flux_surface_average_jvp
    public :: evaluate_flux_surface_average_vjp
    public :: equation_objective_block_t
    public :: equation_objective_registry_t
    public :: initialize_equation_objective_registry
    public :: validate_equation_objective_registry
    public :: equation_objective_registry_block
    public :: equation_objective_registry_block_count
    public :: equation_objective_registry_total_rows
    public :: pack_equation_objective_values
    public :: pack_equation_objective_values_jvp
    public :: pack_equation_objective_values_vjp
    public :: unpack_equation_objective_values
    public :: evaluate_equation_objective_callbacks
    public :: evaluate_equation_objective_callbacks_jvp
    public :: evaluate_equation_objective_callbacks_vjp
    public :: REGISTRY_KIND_EQUATION
    public :: REGISTRY_KIND_OBJECTIVE
    public :: REGISTRY_KIND_CONSTRAINT
    public :: assemble_tensor_diffusion_matrix
    public :: assemble_tensor_diffusion_matrix_jvp
    public :: assemble_tensor_diffusion_matrix_vjp
    public :: assemble_tensor_diffusion_matrix_nd
    public :: assemble_tensor_diffusion_matrix_nd_jvp
    public :: assemble_tensor_diffusion_matrix_nd_vjp
    public :: assemble_tensor_diffusion_matrix_3d
    public :: assemble_tensor_diffusion_matrix_3d_jvp
    public :: assemble_tensor_diffusion_matrix_3d_vjp
    public :: advance_dissipative_cayley
    public :: advance_dissipative_cayley_jvp
    public :: advance_dissipative_cayley_vjp
    public :: assemble_compatible_flux_elimination
    public :: assemble_compatible_flux_elimination_jvp
    public :: assemble_compatible_flux_elimination_vjp
    public :: evaluate_blending_corrected_enrichment
    public :: evaluate_blending_corrected_enrichment_jvp
    public :: evaluate_blending_corrected_enrichment_vjp
    public :: evaluate_vector_blending_corrected_enrichment
    public :: evaluate_vector_blending_corrected_enrichment_jvp
    public :: evaluate_vector_blending_corrected_enrichment_vjp
    public :: evaluate_enrichment_support_gram
    public :: evaluate_enrichment_support_gram_jvp
    public :: evaluate_enrichment_support_gram_vjp
    public :: evaluate_enrichment_support_rank_condition
    public :: evaluate_enrichment_support_activation
    public :: evaluate_enrichment_support_activation_jvp
    public :: evaluate_enrichment_support_activation_vjp
    public :: evaluate_enrichment_support_vector_gram
    public :: evaluate_enrichment_support_vector_gram_jvp
    public :: evaluate_enrichment_support_vector_gram_vjp
    public :: evaluate_piola_enriched_vector_values
    public :: evaluate_piola_enriched_vector_values_jvp
    public :: evaluate_piola_enriched_vector_values_vjp
    public :: evaluate_piola_enriched_vector_differential_3d
    public :: evaluate_piola_enriched_vector_differential_3d_jvp
    public :: evaluate_piola_enriched_vector_differential_3d_vjp
    public :: evaluate_piola_enriched_vector_differential_2d
    public :: evaluate_piola_enriched_vector_differential_2d_jvp
    public :: evaluate_piola_enriched_vector_differential_2d_vjp
    public :: PIOLA_COVARIANT, PIOLA_CONTRAVARIANT
    public :: tree_cotree_gauge_t
    public :: build_tree_cotree_gauge
    public :: build_tree_cotree_dof_map
    public :: validate_tree_cotree_gauge
    public :: tree_cotree_gauge_components
    public :: tree_cotree_gauge_edges
    public :: apply_tree_cotree_restriction
    public :: apply_tree_cotree_prolongation
    public :: reduce_tree_cotree_dense_system
    public :: reduce_tree_cotree_dense_system_jvp
    public :: reduce_tree_cotree_dense_system_vjp
    public :: tree_cotree_iga_parity_t
    public :: diagnose_tree_cotree_iga_invariance
    public :: incomplete_cholesky_factor_t
    public :: build_incomplete_cholesky
    public :: apply_incomplete_cholesky
    public :: sparse_incomplete_cholesky_factor_t
    public :: build_sparse_ichol
    public :: build_sparse_ichol_row
    public :: build_sparse_incomplete_cholesky
    public :: apply_sparse_incomplete_cholesky
    public :: apply_sparse_incomplete_cholesky_jvp
    public :: apply_sparse_incomplete_cholesky_vjp
    public :: sparse_incomplete_lu_factor_t
    public :: build_sparse_incomplete_lu
    public :: apply_sparse_incomplete_lu
    public :: apply_sparse_incomplete_lu_jvp
    public :: apply_sparse_incomplete_lu_vjp
    public :: sparse_ilut_factor_t
    public :: build_sparse_ilut
    public :: build_sparse_ilut_row
    public :: apply_sparse_ilut
    public :: apply_sparse_ilut_jvp
    public :: apply_sparse_ilut_vjp
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
    public :: vector_neumann_bc_t
    public :: neumann_bc_t
    public :: boundary_t
    public :: simple_expression_t
    public :: cell_coefficient_t
    public :: cell_tensor_coefficient_t
    public :: cell_vector_source_t
    public :: form_expr_t
    public :: form_equation_t
    public :: scalar_reluctivity_curvilinear_fourier_coefficients

    ! Public mesh constructors and refinement
    public :: unit_square_mesh
    public :: rectangle_mesh
    public :: unit_disk_mesh
    public :: generate_torus_surface_mesh
    public :: generate_solid_torus_tetra_mesh
    public :: generate_structured_tetra_box_mesh
    public :: evaluate_torus_curved_panel
    public :: evaluate_torus_curved_panel_jvp
    public :: evaluate_torus_curved_panel_vjp
    public :: integrate_laplace_torus_panel_p0_3d
    public :: integrate_laplace_torus_panel_p0_3d_jvp
    public :: integrate_laplace_torus_panel_p0_3d_vjp
    public :: integrate_laplace_sphere_panel_p0_3d
    public :: integrate_laplace_sphere_panel_p0_3d_jvp
    public :: integrate_laplace_sphere_panel_p0_3d_vjp
    public :: integrate_helmholtz_torus_panel_p0_3d
    public :: integrate_helmholtz_torus_panel_p0_3d_jvp
    public :: integrate_helmholtz_torus_panel_p0_3d_vjp
    public :: integrate_helmholtz_sphere_panel_p0_3d
    public :: integrate_helmholtz_sphere_panel_p0_3d_jvp
    public :: integrate_helmholtz_sphere_panel_p0_3d_vjp
    public :: generate_sphere_surface_mesh
    public :: evaluate_sphere_curved_panel
    public :: evaluate_sphere_curved_panel_jvp
    public :: evaluate_sphere_curved_panel_vjp
    public :: invert_sphere_curved_panel
    public :: evaluate_maxwell_sphere_curved_rwg_basis
    public :: assemble_maxwell_sphere_curved_efie_rwg_3d
    public :: assemble_maxwell_sphere_curved_efie_imaginary_rwg_3d
    public :: assemble_maxwell_sphere_efie_propagating_impedance_jvp
    public :: assemble_maxwell_sphere_efie_propagating_impedance_vjp
    public :: assemble_maxwell_sphere_efie_imaginary_impedance_jvp
    public :: assemble_maxwell_sphere_efie_imaginary_impedance_vjp
    public :: assemble_maxwell_sphere_efie_wave_number_jvp
    public :: assemble_maxwell_sphere_efie_wave_number_vjp
    public :: assemble_maxwell_sphere_efie_imaginary_decay_jvp
    public :: assemble_maxwell_sphere_efie_imaginary_decay_vjp
    public :: assemble_maxwell_sphere_curved_efie_bc_imaginary_3d
    public :: assemble_maxwell_sphere_curved_regularized_cfie_rwg_3d
    public :: assemble_maxwell_sphere_curved_regularized_cfie_rhs_rwg_3d
    public :: assemble_maxwell_sphere_curved_mfie_exterior_trace_rwg_rbc_3d
    public :: assemble_maxwell_sphere_curved_mfie_rwg_rbc_3d
    public :: solve_maxwell_pec_sphere_curved_efie_rwg_3d
    public :: solve_maxwell_pec_sphere_curved_regularized_cfie_rwg_3d
    public :: assemble_maxwell_sphere_curved_rwg_mass_matrix
    public :: assemble_maxwell_sphere_curved_rwg_mass_matrix_jvp
    public :: assemble_maxwell_sphere_curved_rwg_mass_matrix_vjp
    public :: assemble_maxwell_sphere_curved_rwg_rbc_pairing
    public :: assemble_maxwell_sphere_curved_rwg_rbc_pairing_jvp
    public :: assemble_maxwell_sphere_curved_rwg_rbc_pairing_vjp
    public :: assemble_maxwell_sphere_curved_plane_wave_rhs_rwg_3d
    public :: assemble_maxwell_sphere_curved_plane_wave_rhs_bc_3d
    public :: assemble_maxwell_sphere_curved_plane_wave_rhs_bc_3d_jvp
    public :: assemble_maxwell_sphere_curved_plane_wave_rhs_bc_3d_vjp
    public :: assemble_maxwell_sphere_curved_potential_operators_rwg_3d
    public :: assemble_maxwell_sphere_curved_vector_potential_rwg_3d
    public :: evaluate_maxwell_sphere_curved_far_field_rwg_3d
    public :: evaluate_maxwell_sphere_curved_far_field_rwg_3d_jvp
    public :: evaluate_maxwell_sphere_curved_far_field_rwg_3d_vjp
    public :: evaluate_maxwell_sphere_curved_localized_rwg_basis
    public :: evaluate_maxwell_sphere_curved_localized_rwg_basis_jvp
    public :: evaluate_maxwell_sphere_curved_localized_rwg_basis_vjp
    public :: evaluate_maxwell_sphere_curved_rwg_basis_jvp
    public :: evaluate_maxwell_sphere_curved_rwg_basis_vjp
    public :: evaluate_maxwell_sphere_curved_magnetic_field_rwg_3d
    public :: evaluate_maxwell_sphere_curved_magnetic_field_rwg_3d_jvp
    public :: evaluate_maxwell_sphere_curved_magnetic_field_rwg_3d_vjp
    public :: integrate_maxwell_sphere_curved_coincident_rwg_pair_3d
    public :: integrate_maxwell_sphere_curved_adjacent_rwg_pair_3d
    public :: barycentric_refine_surface_mesh
    public :: barycentric_refine_sphere_surface_mesh_jvp
    public :: barycentric_refine_sphere_surface_mesh_vjp
    public :: barycentric_refine_torus_surface_mesh
    public :: barycentric_refine_torus_surface_mesh_jvp
    public :: barycentric_refine_torus_surface_mesh_vjp
    public :: evaluate_maxwell_localized_rwg_basis
    public :: build_maxwell_bc_transformation
    public :: differentiate_maxwell_bc_transformation_jvp
    public :: differentiate_maxwell_bc_transformation_vjp
    public :: assemble_maxwell_rwg_rbc_pairing
    public :: evaluate_maxwell_magnetic_field_rwg_3d
    public :: evaluate_maxwell_magnetic_field_rwg_3d_jvp
    public :: evaluate_maxwell_magnetic_field_rwg_3d_vjp
    public :: evaluate_maxwell_magnetic_field_rwg_3d_targets
    public :: evaluate_maxwell_magnetic_field_rwg_3d_targets_jvp
    public :: evaluate_maxwell_magnetic_field_rwg_3d_targets_vjp
    public :: evaluate_maxwell_virtual_casing_rwg_3d
    public :: evaluate_maxwell_virtual_casing_rwg_3d_jvp
    public :: evaluate_maxwell_virtual_casing_rwg_3d_vjp
    public :: advance_resistive_wall_midpoint
    public :: advance_resistive_wall_midpoint_jvp
    public :: advance_resistive_wall_midpoint_vjp
    public :: evaluate_resistive_wall_energy_balance
    public :: evaluate_toroidal_spectral_trace
    public :: evaluate_toroidal_spectral_trace_jvp
    public :: evaluate_toroidal_spectral_trace_vjp
    public :: evaluate_toroidal_spectral_trace_grid
    public :: evaluate_toroidal_spectral_trace_grid_jvp
    public :: evaluate_toroidal_spectral_trace_grid_vjp
    public :: solve_toroidal_spectral_neumann
    public :: solve_toroidal_spectral_neumann_jvp
    public :: solve_toroidal_spectral_neumann_vjp
    public :: analyze_toroidal_spectral_modes
    public :: analyze_toroidal_spectral_modes_jvp
    public :: analyze_toroidal_spectral_modes_vjp
    public :: condense_wall_response_blocks
    public :: condense_wall_response_blocks_jvp
    public :: condense_wall_response_blocks_vjp
    public :: assemble_generalized_debye_source_residual
    public :: assemble_generalized_debye_source_residual_jvp
    public :: assemble_generalized_debye_source_residual_vjp
    public :: assemble_generalized_debye_source_second_kind
    public :: assemble_generalized_debye_source_second_kind_jvp
    public :: assemble_generalized_debye_source_second_kind_vjp
    public :: assemble_maxwell_mfie_rwg_rbc_3d
    public :: assemble_maxwell_bc_scalar_potential_3d
    public :: assemble_maxwell_bc_potential_operators_3d
    public :: assemble_maxwell_efie_bc_3d
    public :: assemble_maxwell_efie_bc_imaginary_3d
    public :: build_maxwell_bc_panel_divergence
    public :: build_maxwell_bc_to_refined_rwg
    public :: assemble_maxwell_regularized_cfie_rwg_3d
    public :: assemble_maxwell_plane_wave_rhs_bc_3d
    public :: solve_maxwell_pec_regularized_cfie_rwg_3d
    public :: solve_maxwell_pec_regularized_cfie_rwg_multiple_3d
    public :: assemble_maxwell_torus_curved_rwg_mass_matrix
    public :: assemble_maxwell_torus_curved_rwg_mass_matrix_jvp
    public :: assemble_maxwell_torus_curved_rwg_mass_matrix_vjp
    public :: assemble_maxwell_torus_curved_rwg_rbc_pairing
    public :: assemble_maxwell_torus_curved_rwg_rbc_pairing_jvp
    public :: assemble_maxwell_torus_curved_rwg_rbc_pairing_vjp
    public :: assemble_maxwell_torus_curved_efie_rwg_3d
    public :: assemble_maxwell_torus_curved_efie_imaginary_rwg_3d
    public :: assemble_maxwell_torus_efie_imaginary_impedance_jvp
    public :: assemble_maxwell_torus_efie_imaginary_impedance_vjp
    public :: assemble_maxwell_torus_efie_propagating_impedance_jvp
    public :: assemble_maxwell_torus_efie_propagating_impedance_vjp
    public :: assemble_maxwell_torus_efie_imaginary_decay_jvp
    public :: assemble_maxwell_torus_efie_imaginary_decay_vjp
    public :: assemble_maxwell_torus_efie_wave_number_jvp
    public :: assemble_maxwell_torus_efie_wave_number_vjp
    public :: assemble_maxwell_torus_curved_efie_bc_imaginary_3d
    public :: assemble_maxwell_torus_curved_mfie_offset_trace_rwg_rbc_3d
    public :: assemble_maxwell_torus_mfie_offset_jvp
    public :: assemble_maxwell_torus_mfie_offset_vjp
    public :: assemble_maxwell_torus_mfie_offset_geometry_jvp
    public :: assemble_maxwell_torus_mfie_offset_geometry_vjp
    public :: assemble_maxwell_torus_curved_mfie_rwg_rbc_3d
    public :: assemble_maxwell_torus_curved_plane_wave_rhs_rwg_3d
    public :: assemble_maxwell_torus_curved_plane_wave_rhs_rwg_3d_jvp
    public :: assemble_maxwell_torus_curved_plane_wave_rhs_rwg_3d_vjp
    public :: assemble_maxwell_torus_curved_plane_wave_rhs_bc_3d
    public :: assemble_maxwell_torus_curved_plane_wave_rhs_bc_3d_jvp
    public :: assemble_maxwell_torus_curved_plane_wave_rhs_bc_3d_vjp
    public :: assemble_maxwell_torus_curved_potential_operators_rwg_3d
    public :: assemble_maxwell_torus_curved_regularized_cfie_rwg_3d
    public :: evaluate_maxwell_torus_curved_far_field_rwg_3d
    public :: evaluate_maxwell_torus_curved_far_field_rwg_3d_jvp
    public :: evaluate_maxwell_torus_curved_far_field_rwg_3d_vjp
    public :: evaluate_maxwell_torus_curved_magnetic_field_rwg_3d
    public :: evaluate_maxwell_torus_curved_magnetic_field_rwg_3d_jvp
    public :: evaluate_maxwell_torus_magnetic_geometry_jvp
    public :: evaluate_maxwell_torus_magnetic_geometry_vjp
    public :: evaluate_maxwell_torus_curved_magnetic_field_rwg_3d_vjp
    public :: evaluate_maxwell_torus_curved_localized_rwg_basis
    public :: evaluate_maxwell_torus_curved_localized_rwg_basis_jvp
    public :: evaluate_maxwell_torus_curved_localized_rwg_basis_vjp
    public :: evaluate_maxwell_torus_curved_rwg_basis
    public :: evaluate_maxwell_torus_curved_rwg_basis_jvp
    public :: evaluate_maxwell_torus_curved_rwg_basis_vjp
    public :: integrate_maxwell_torus_curved_adjacent_rwg_pair_3d
    public :: integrate_maxwell_torus_curved_coincident_rwg_pair_3d
    public :: solve_maxwell_pec_torus_curved_regularized_cfie_rwg_3d
    public :: solve_maxwell_pec_torus_curved_regularized_cfie_rwg_multiple_3d
    public :: apply_maxwell_trace_to_flux
    public :: apply_maxwell_trace_to_flux_jvp
    public :: apply_maxwell_trace_to_flux_map
    public :: apply_maxwell_trace_to_flux_vjp
    public :: apply_maxwell_weak_trace_reconstruction
    public :: apply_maxwell_weak_trace_reconstruction_jvp
    public :: apply_maxwell_weak_trace_reconstruction_vjp
    public :: assemble_maxwell_trace_to_flux_map
    public :: assemble_maxwell_trace_to_flux_map_jvp
    public :: assemble_maxwell_trace_to_flux_map_vjp
    public :: assemble_maxwell_weak_trace_reconstruction
    public :: assemble_maxwell_torus_curved_dtn_rwg_3d
    public :: cartesian_to_toroidal
    public :: cartesian_to_toroidal_jvp
    public :: cartesian_to_toroidal_vjp
    public :: toroidal_point_to_cartesian
    public :: toroidal_point_to_cartesian_jvp
    public :: toroidal_point_to_cartesian_vjp
    public :: toroidal_vector_to_cartesian
    public :: toroidal_vector_to_cartesian_jvp
    public :: toroidal_vector_to_cartesian_vjp
    public :: evaluate_laplace_representation_triangles_3d
    public :: evaluate_laplace_representation_torus_curved_3d
    public :: evaluate_laplace_representation_torus_curved_3d_jvp
    public :: evaluate_laplace_representation_torus_curved_3d_vjp
    public :: evaluate_laplace_representation_torus_curved_3d_geometry_jvp
    public :: evaluate_laplace_representation_torus_curved_3d_geometry_vjp
    public :: assemble_laplace_torus_curved_calderon_3d
    public :: assemble_laplace_torus_curved_dtn_3d
    public :: assemble_laplace_torus_curved_dtn_3d_geometry_jvp
    public :: assemble_laplace_torus_curved_dtn_3d_geometry_vjp
    public :: solve_laplace_bem_dtn_torus_curved_3d
    public :: assemble_laplace_fem_bem_costabel_torus_curved_3d
    public :: solve_laplace_fem_bem_costabel_torus_curved_3d
    public :: assemble_laplace_single_layer_p0_3d
    public :: assemble_laplace_single_layer_p0_3d_jvp
    public :: assemble_laplace_single_layer_p0_3d_vjp
    public :: integrate_laplace_single_layer_regular_panel_pair_p0_3d
    public :: integrate_laplace_single_layer_regular_panel_pair_p0_3d_jvp
    public :: integrate_laplace_single_layer_regular_panel_pair_p0_3d_vjp
    public :: assemble_laplace_single_layer_p0_adaptive_3d
    public :: assemble_laplace_calderon_p1_p0_3d
    public :: solve_laplace_dirichlet_p0_3d
    public :: solve_laplace_dirichlet_p0_3d_jvp
    public :: solve_laplace_dirichlet_p0_3d_vjp
    public :: estimate_helmholtz_p0_two_level_residual_3d
    public :: estimate_laplace_p0_two_level_residual_3d
    public :: mark_bem_dorfler
    public :: refine_surface_mesh_marked
    public :: assemble_laplace_fem_bem_costabel_3d
    public :: solve_laplace_fem_bem_costabel_3d
    public :: solve_laplace_fem_bem_johnson_nedelec_3d
    public :: apply_laplace_single_layer_p0_hierarchical_3d
    public :: apply_helmholtz_single_layer_p0_hierarchical_3d
    public :: apply_helmholtz_cfie_p0_hierarchical_3d
    public :: solve_helmholtz_dirichlet_p0_hierarchical_3d
    public :: solve_helmholtz_cfie_p0_hierarchical_3d
    public :: evaluate_helmholtz_representation_triangles_3d
    public :: evaluate_helmholtz_representation_torus_curved_3d
    public :: evaluate_helmholtz_representation_torus_curved_3d_geometry_jvp
    public :: evaluate_helmholtz_representation_torus_curved_3d_geometry_vjp
    public :: assemble_helmholtz_torus_curved_calderon_3d
    public :: assemble_helmholtz_torus_curved_dtn_3d
    public :: assemble_helmholtz_torus_curved_dtn_3d_geometry_jvp
    public :: assemble_helmholtz_torus_curved_dtn_3d_geometry_vjp
    public :: solve_helmholtz_bem_dtn_torus_curved_3d
    public :: assemble_helmholtz_fem_bem_costabel_torus_curved_3d
    public :: solve_helmholtz_fem_bem_costabel_torus_curved_3d
    public :: assemble_helmholtz_single_layer_p0_3d
    public :: assemble_helmholtz_single_layer_p0_3d_jvp
    public :: assemble_helmholtz_single_layer_p0_3d_vjp
    public :: integrate_helmholtz_single_layer_regular_panel_pair_p0_3d
    public :: integrate_helmholtz_single_layer_regular_panel_pair_p0_3d_jvp
    public :: integrate_helmholtz_single_layer_regular_panel_pair_p0_3d_vjp
    public :: assemble_helmholtz_calderon_p1_p0_3d
    public :: assemble_helmholtz_fem_bem_costabel_3d
    public :: assemble_helmholtz_single_layer_p0_adaptive_3d
    public :: assemble_helmholtz_double_layer_p0_3d
    public :: evaluate_helmholtz_cfie_p0_3d
    public :: solve_helmholtz_cfie_p0_3d
    public :: solve_helmholtz_dirichlet_p0_3d
    public :: solve_helmholtz_dirichlet_p0_3d_jvp
    public :: solve_helmholtz_dirichlet_p0_3d_vjp
    public :: solve_helmholtz_fem_bem_costabel_3d
    public :: build_maxwell_rwg_surface_space
    public :: assemble_maxwell_rwg_mass_matrix
    public :: evaluate_maxwell_rwg_basis
    public :: build_maxwell_surface_rt_dof_map
    public :: assemble_maxwell_surface_rt_mass_matrix
    public :: assemble_maxwell_surface_rt_efie_3d
    public :: build_bspline_derivative_matrix
    public :: build_bspline_feec_2d_operators
    public :: build_bspline_feec_3d_operators
    public :: evaluate_bspline_basis
    public :: evaluate_nurbs_surface_geometry
    public :: evaluate_nurbs_surface_geometry_jvp
    public :: evaluate_nurbs_surface_geometry_vjp
    public :: evaluate_surface_triangle_geometry_3d
    public :: evaluate_surface_triangle_geometry_3d_jvp
    public :: evaluate_surface_triangle_geometry_3d_vjp
    public :: assemble_surface_triangle_areas_3d
    public :: assemble_surface_triangle_areas_3d_jvp
    public :: assemble_surface_triangle_areas_3d_vjp
    public :: assemble_surface_triangle_measures_3d
    public :: assemble_surface_triangle_measures_3d_jvp
    public :: assemble_surface_triangle_measures_3d_vjp
    public :: evaluate_level_set_triangle_interface_2d
    public :: evaluate_level_set_triangle_interface_2d_jvp
    public :: evaluate_level_set_triangle_cut_areas_2d
    public :: evaluate_level_set_triangle_cut_quadrature_2d
    public :: evaluate_level_set_triangle_cut_quadrature_2d_jvp
    public :: evaluate_level_set_triangle_cut_moments_2d
    public :: evaluate_level_set_triangle_cut_moments_2d_jvp
    public :: evaluate_level_set_triangle_cut_third_moments_2d
    public :: evaluate_level_set_triangle_cut_third_moments_2d_jvp
    public :: evaluate_level_set_triangle_cut_fourth_moments_2d
    public :: evaluate_level_set_triangle_cut_fourth_moments_2d_jvp
    public :: evaluate_level_set_tetra_interface_3d
    public :: evaluate_level_set_tetra_interface_3d_jvp
    public :: evaluate_level_set_tetra_cut_quadrature_3d
    public :: evaluate_level_set_tetra_cut_quadrature_3d_jvp
    public :: evaluate_level_set_tetra_cut_moments_3d
    public :: evaluate_level_set_tetra_cut_moments_3d_jvp
    public :: evaluate_level_set_tetra_cut_third_moments_3d
    public :: evaluate_level_set_tetra_cut_third_moments_3d_jvp
    public :: evaluate_level_set_tetra_cut_fourth_moments_3d
    public :: evaluate_level_set_tetra_cut_fourth_moments_3d_jvp
    public :: find_fci_first_hit_segment_2d
    public :: find_fci_first_hit_segment_2d_jvp
    public :: find_fci_first_hit_triangle_3d
    public :: find_fci_first_hit_triangle_3d_jvp
    public :: find_fci_first_hit_polyline_3d
    public :: find_fci_first_hit_polyline_3d_jvp
    public :: assemble_fci_terminal_boundary_flux
    public :: assemble_fci_terminal_boundary_flux_jvp
    public :: assemble_fci_terminal_boundary_flux_vjp
    public :: assemble_fci_terminal_boundary_ledger
    public :: assemble_fci_terminal_boundary_ledger_jvp
    public :: assemble_fci_terminal_boundary_ledger_vjp
    public :: evaluate_fci_power_flux_ledger
    public :: evaluate_fci_power_flux_ledger_jvp
    public :: evaluate_fci_power_flux_ledger_vjp
    public :: assemble_volume_balance_ledger
    public :: assemble_volume_balance_ledger_jvp
    public :: assemble_volume_balance_ledger_vjp
    public :: assemble_nonlinear_surface_flux
    public :: assemble_nonlinear_surface_flux_jvp
    public :: assemble_nonlinear_surface_flux_vjp
    public :: evaluate_nurbs_volume_geometry
    public :: evaluate_nurbs_volume_geometry_jvp
    public :: evaluate_nurbs_volume_geometry_vjp
    public :: map_isogeometric_h1_gradient
    public :: map_isogeometric_hcurl
    public :: map_isogeometric_hdiv
    public :: map_isogeometric_l2
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
    public :: apply_bspline_jorek_flux_rhs
    public :: apply_bspline_jorek_flux_jvp
    public :: apply_bspline_jorek_thermodynamic_rhs
    public :: apply_bspline_jorek_thermodynamic_jvp
    public :: apply_bspline_jorek_density_rhs
    public :: apply_bspline_jorek_density_jvp
    public :: project_bspline_toroidal_product
    public :: advance_bspline_jorek_poloidal_flux_midpoint
    public :: advance_bspline_jorek_poloidal_flux_midpoint_steps
    public :: build_bspline_feec_2d_operators_csc
    public :: build_bspline_feec_3d_operators_csc
    public :: assemble_bspline_h1_operator_3d_csc
    public :: assemble_bspline_hcurl_operator_3d_csc
    public :: assemble_bspline_hdiv_operator_3d_csc
    public :: assemble_bspline_hdiv_l2_divergence_3d_csc
    public :: assemble_bspline_l2_hdiv_adjoint_divergence_3d_csc
    public :: assemble_bspline_l2_mass_3d_csc
    public :: apply_bspline_toroidal_poloidal_bracket
    public :: apply_bspline_toroidal_poloidal_bracket_jvp
    public :: apply_toroidal_fourier_derivative
    public :: BSPLINE_FACE_X_MAX
    public :: BSPLINE_FACE_X_MIN
    public :: BSPLINE_FACE_Y_MAX
    public :: BSPLINE_FACE_Y_MIN
    public :: BSPLINE_FACE_Z_MAX
    public :: BSPLINE_FACE_Z_MIN
    public :: build_bspline_feec_2d_interface_dofs
    public :: build_bspline_feec_2d_multipatch_maps
    public :: build_bspline_feec_2d_two_patch_maps
    public :: build_bspline_feec_2d_two_patch_operators_csc
    public :: build_bspline_feec_3d_interface_dofs
    public :: build_bspline_feec_3d_multipatch_maps
    public :: build_bspline_feec_3d_two_patch_maps
    public :: build_bspline_feec_3d_two_patch_operators_csc
    public :: scalar_weight_2d
    public :: tensor_weight_2d
    public :: evaluate_maxwell_surface_rt_basis
    public :: evaluate_maxwell_surface_rt_global_basis
    public :: map_maxwell_rwg_to_tetra_nedelec_edges
    public :: assemble_maxwell_efie_rwg_3d
    public :: assemble_maxwell_plane_wave_rhs_rwg_3d
    public :: assemble_maxwell_rwg_potential_operators_3d
    public :: evaluate_maxwell_efie_field_rwg_3d
    public :: evaluate_maxwell_efie_far_field_rwg_3d
    public :: solve_maxwell_pec_efie_rwg_3d
    public :: assemble_maxwell_fem_bem_boundary_matrix_3d
    public :: assemble_maxwell_fem_bem_system_3d
    public :: assemble_maxwell_fem_bem_torus_curved_system_3d
    public :: assemble_maxwell_rwg_nedelec_coupling_3d
    public :: assemble_maxwell_torus_curved_rwg_nedelec_coupling_3d
    public :: solve_maxwell_fem_bem_system_3d
    public :: solve_maxwell_fem_bem_torus_curved_system_3d
    public :: solve_maxwell_fem_bem_linear_state
    public :: solve_maxwell_fem_bem_linear_state_jvp
    public :: solve_maxwell_fem_bem_linear_state_vjp
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
    public :: cell_coefficient
    public :: cell_tensor_coefficient
    public :: cell_vector_source
    public :: dirichlet_bc
    public :: dirichlet_bc_on_boundary
    public :: vector_bc
    public :: vector_bc_edge_moments
    public :: vector_bc_on_edges
    public :: vector_neumann_bc_on_edges
    public :: neumann_bc_constant
    public :: neumann_bc_on_boundary

    ! Form operations
    public :: inner
    public :: grad
    public :: curl
    public :: div
    public :: dx
    public :: compile_form
    public :: compile_form_matrix
    public :: compile_form_vector
    public :: compile_vector_form_element
    public :: compile_vector_form_csc
    public :: compile_mixed_form_csc
    public :: compile_tetra_mixed_form_csc
    public :: compile_tetra_rt_form_csc
    public :: compile_vector_form_rhs
    public :: init_measures
    public :: operator(*)
    public :: operator(+)
    public :: operator(==)

    ! Solver interface
    public :: assemble_laplacian_system
    public :: solve
    public :: solve_mixed_poisson_rt0
    public :: solve_mixed_poisson_rt
    public :: solve_mixed_rt_system
    public :: solve_mixed_rt_system_jvp
    public :: solve_mixed_rt_system_vjp
    public :: solve_tetra_mixed_poisson_state
    public :: solve_tetra_mixed_poisson_state_jvp
    public :: solve_tetra_mixed_poisson_state_vjp
    public :: solve_tetra_mixed_poisson_sampled_state
    public :: solve_tetra_mixed_poisson_sampled_state_jvp
    public :: solve_tetra_mixed_poisson_sampled_state_vjp
    public :: solve_symbolic_mixed_poisson_rt
    public :: solve_magnetic_box_3d
    public :: solve_tetra_nedelec_curl_mass
    public :: solve_tetra_nedelec_sampled_state
    public :: solve_tetra_nedelec_sampled_state_jvp
    public :: solve_tetra_nedelec_sampled_state_vjp
    public :: solve_triangle_nedelec_sampled_state
    public :: solve_triangle_nedelec_sampled_state_jvp
    public :: solve_triangle_nedelec_sampled_state_vjp
    public :: solve_triangle_nedelec_second_sampled_state
    public :: solve_triangle_nedelec_second_sampled_state_jvp
    public :: solve_triangle_nedelec_second_sampled_state_vjp
    public :: solve_triangle_bdm_sampled_state
    public :: solve_triangle_bdm_sampled_state_jvp
    public :: solve_triangle_bdm_sampled_state_vjp
    public :: invert_triangle_affine_map
    public :: invert_triangle_affine_map_jvp
    public :: invert_triangle_affine_map_vjp
    public :: evaluate_tetra_nedelec_interpolant
    public :: evaluate_tetra_nedelec_interpolant_jvp
    public :: evaluate_tetra_nedelec_interpolant_vjp
    public :: evaluate_tetra_nedelec_interpolant_at_point
    public :: evaluate_tetra_nedelec_interpolant_at_point_jvp
    public :: evaluate_tetra_nedelec_interpolant_at_point_vjp
    public :: evaluate_tetra_rt_interpolant
    public :: evaluate_tetra_rt_interpolant_jvp
    public :: evaluate_tetra_rt_interpolant_vjp
    public :: evaluate_tetra_rt_interpolant_at_point
    public :: evaluate_tetra_rt_interpolant_at_point_jvp
    public :: evaluate_tetra_rt_interpolant_at_point_vjp
    public :: invert_tetra_affine_map
    public :: invert_tetra_affine_map_jvp
    public :: invert_tetra_affine_map_vjp
    public :: solve_triangle_rt_sampled_state
    public :: solve_triangle_rt_sampled_state_jvp
    public :: solve_triangle_rt_sampled_state_vjp
    public :: solve_tetra_nedelec_pml
    public :: solve_tetra_nedelec_curvilinear_pml
    public :: solve_tetra_nedelec_pml_jvp
    public :: solve_tetra_nedelec_pml_vjp
    public :: solve_tetra_nedelec_curvilinear_pml_jvp
    public :: solve_tetra_nedelec_curvilinear_pml_vjp
    public :: solve_tetra_nedelec_weighted_curl_mass
    public :: solve_tetra_rt_div_mass
    public :: solve_symbolic_tetra_mixed_poisson_rt
    public :: assemble_tetra_dg_source_load_samples
    public :: assemble_tetra_dg_source_load_samples_jvp
    public :: assemble_tetra_dg_source_load_samples_vjp
    public :: evaluate_tetra_lagrange_solution
    public :: evaluate_tetra_lagrange_solution_at_point
    public :: evaluate_tetra_lagrange_solution_at_point_jvp
    public :: evaluate_tetra_lagrange_solution_at_point_vjp
    public :: evaluate_tetra_lagrange_solution_jvp
    public :: evaluate_tetra_lagrange_solution_vjp
    public :: evaluate_tetra_lagrange_solution_prepared
    public :: initialize_tetra_lagrange_solution_evaluator
    public :: solve_tetra_lagrange_diffusion_reaction
    public :: solve_tetra_lagrange_poisson
    public :: solve_tetra_lagrange_state
    public :: solve_tetra_lagrange_state_jvp
    public :: solve_tetra_lagrange_state_vjp
    public :: solve_tetra_lagrange_curvilinear_pml
    public :: solve_tetra_lagrange_curvilinear_pml_jvp
    public :: solve_tetra_lagrange_curvilinear_pml_vjp
    public :: solve_tetra_lagrange_sampled_state
    public :: solve_tetra_lagrange_sampled_state_jvp
    public :: solve_tetra_lagrange_sampled_state_vjp
    public :: tetra_lagrange_solution_evaluator_t
    public :: solve_mixed_bc
    public :: solve_neumann
    public :: compute_boundary_integral
    public :: apply_planar_helmholtz_dtn
    public :: apply_planar_helmholtz_dtn_jvp
    public :: apply_planar_helmholtz_dtn_vjp
    public :: apply_planar_maxwell_dtn
    public :: apply_planar_maxwell_dtn_jvp
    public :: apply_planar_maxwell_dtn_vjp
    public :: assemble_planar_maxwell_dtn_form
    public :: assemble_planar_maxwell_dtn_form_jvp
    public :: assemble_planar_maxwell_dtn_form_vjp
    public :: solve_planar_maxwell_dtn_system
    public :: solve_planar_maxwell_dtn_system_jvp
    public :: solve_planar_maxwell_dtn_system_vjp
    public :: assemble_planar_nedelec_maxwell_dtn_form
    public :: build_planar_nedelec_trace_sampling
    public :: pullback_planar_maxwell_dtn_form
    public :: pullback_planar_maxwell_dtn_form_jvp
    public :: pullback_planar_maxwell_dtn_form_vjp
    public :: apply_planar_acoustic_displacement_dtn
    public :: assemble_planar_acoustic_displacement_dtn_form
    public :: solve_elasticity_planar_acoustic_dtn_p1
    public :: solve_elasticity_curved_acoustic_ntd_p1
    public :: assemble_planar_helmholtz_dtn_form
    public :: assemble_planar_helmholtz_dtn_form_jvp
    public :: assemble_planar_helmholtz_dtn_form_vjp
    public :: solve_scalar_helmholtz_planar_dtn_p1
    public :: solve_scalar_helmholtz_planar_dtn_p1_jvp
    public :: solve_scalar_helmholtz_planar_dtn_p1_vjp
    public :: apply_circular_helmholtz_dtn
    public :: apply_circular_helmholtz_dtn_jvp
    public :: apply_circular_helmholtz_dtn_vjp
    public :: circular_helmholtz_dtn_eigenvalue
    public :: circular_helmholtz_dtn_eigenvalue_jvp
    public :: circular_helmholtz_dtn_eigenvalue_vjp
    public :: apply_spherical_helmholtz_dtn
    public :: apply_spherical_helmholtz_dtn_jvp
    public :: apply_spherical_helmholtz_dtn_vjp
    public :: spherical_helmholtz_dtn_eigenvalue
    public :: spherical_helmholtz_dtn_eigenvalue_jvp
    public :: spherical_helmholtz_dtn_eigenvalue_vjp
    public :: apply_spherical_maxwell_dtn
    public :: apply_spherical_maxwell_dtn_jvp
    public :: apply_spherical_maxwell_dtn_vjp
    public :: spherical_maxwell_dtn_eigenvalues
    public :: spherical_maxwell_dtn_eigenvalues_jvp
    public :: spherical_maxwell_dtn_eigenvalues_vjp
    public :: cartesian_curl_curl_pml_coefficients
    public :: cartesian_curl_curl_pml_coefficients_jvp
    public :: cartesian_curl_curl_pml_coefficients_vjp
    public :: cartesian_scalar_helmholtz_pml_coefficients
    public :: cartesian_scalar_helmholtz_pml_coefficients_jvp
    public :: cartesian_scalar_helmholtz_pml_coefficients_vjp
    public :: curvilinear_curl_curl_pml_coefficients
    public :: curvilinear_curl_curl_pml_coefficients_jvp
    public :: curvilinear_curl_curl_pml_coefficients_vjp
    public :: curvilinear_scalar_helmholtz_pml_coefficients
    public :: curvilinear_scalar_helmholtz_pml_coefficients_jvp
    public :: curvilinear_scalar_helmholtz_pml_coefficients_vjp
    public :: build_cartesian_pml_element_stretch
    public :: build_cartesian_pml_element_stretch_jvp
    public :: build_cartesian_pml_element_stretch_vjp
    public :: build_curvilinear_normal_pml_element_stretch
    public :: build_curvilinear_normal_pml_element_stretch_jvp
    public :: build_curvilinear_normal_pml_element_stretch_vjp
    public :: solve_scalar_helmholtz_pml_slab_1d
    public :: solve_scalar_helmholtz_pml_slab_1d_jvp
    public :: solve_scalar_helmholtz_pml_slab_1d_vjp
    public :: solve_scalar_helmholtz_pml_p1_2d
    public :: solve_scalar_helmholtz_pml_p1_2d_jvp
    public :: solve_scalar_helmholtz_pml_p1_2d_vjp
    public :: solve_scalar_helmholtz_pml_p1_3d
    public :: solve_scalar_helmholtz_pml_p1_3d_jvp
    public :: solve_scalar_helmholtz_pml_p1_3d_vjp
    public :: evaluate_toroidal_harmonic_p
    public :: evaluate_toroidal_harmonic_p_jvp
    public :: evaluate_toroidal_harmonic_p_vjp
    public :: evaluate_toroidal_ampere_field_p
    public :: evaluate_toroidal_ampere_field_p_jvp
    public :: evaluate_toroidal_ampere_field_p_vjp
    public :: toroidal_poisson_exterior_dtn_p
    public :: toroidal_poisson_exterior_dtn_p_jvp
    public :: toroidal_poisson_exterior_dtn_p_vjp
    public :: evaluate_helmholtz_combined_potential_adaptive_constant
    public :: evaluate_helmholtz_combined_potential_constant
    public :: solve_helmholtz_cfie_constant
    public :: triangle_duffy_quadrature
    public :: tetra_duffy_quadrature
    public :: assignment(=)
    public :: evaluate_triangle_lagrange_basis
    public :: initialize_triangle_lagrange_basis
    public :: triangle_lagrange_basis_t
    public :: triangle_lagrange_nodes
    public :: evaluate_triangle_nedelec_first_kind
    public :: evaluate_triangle_nedelec_first_kind_jvp
    public :: evaluate_triangle_nedelec_first_kind_vjp
    public :: initialize_triangle_nedelec_first_kind
    public :: triangle_nedelec_dof_count
    public :: triangle_nedelec_first_kind_t
    public :: evaluate_triangle_nedelec_second_kind
    public :: evaluate_triangle_nedelec_second_kind_jvp
    public :: evaluate_triangle_nedelec_second_kind_vjp
    public :: initialize_triangle_nedelec_second_kind
    public :: triangle_nedelec_second_kind_dof_count
    public :: triangle_nedelec_second_kind_t
    public :: evaluate_tetra_nedelec_first_order
    public :: evaluate_tetra_nedelec_first_kind
    public :: evaluate_tetra_nedelec_first_kind_jvp
    public :: evaluate_tetra_nedelec_first_kind_vjp
    public :: initialize_tetra_nedelec_first_kind
    public :: tetra_nedelec_dof_count
    public :: tetra_nedelec_first_kind_t
    public :: evaluate_tetra_discontinuous
    public :: evaluate_tetra_discontinuous_jvp
    public :: evaluate_tetra_discontinuous_vjp
    public :: initialize_tetra_discontinuous
    public :: tetra_discontinuous_dof_count
    public :: tetra_discontinuous_t
    public :: project_physical_tetra_discontinuous
    public :: project_sampled_physical_tetra_discontinuous
    public :: project_sampled_physical_tetra_discontinuous_jvp
    public :: project_sampled_physical_tetra_discontinuous_vjp
    public :: evaluate_tetra_lagrange
    public :: evaluate_tetra_lagrange_jvp
    public :: evaluate_tetra_lagrange_vjp
    public :: initialize_tetra_lagrange
    public :: tetra_lagrange_dof_count
    public :: tetra_lagrange_barycentric_indices
    public :: tetra_lagrange_nodes
    public :: tetra_lagrange_t
    public :: build_tetra_lagrange_dof_map
    public :: build_tetra_discrete_gradient
    public :: build_tetra_discrete_curl
    public :: interpolate_reference_tetra_rt
    public :: interpolate_physical_tetra_rt
    public :: interpolate_sampled_physical_tetra_rt
    public :: interpolate_sampled_physical_tetra_rt_jvp
    public :: interpolate_sampled_physical_tetra_rt_vjp
    public :: tetra_rt_interpolation_points
    public :: evaluate_tetra_rt
    public :: evaluate_tetra_rt_jvp
    public :: evaluate_tetra_rt_vjp
    public :: initialize_tetra_rt
    public :: tetra_rt_dof_count
    public :: tetra_rt_t
    public :: build_tetra_rt_basis_transform
    public :: build_tetra_rt_dof_map
    public :: map_tetra_nedelec_covariant
    public :: map_tetra_nedelec_covariant_jvp
    public :: map_tetra_nedelec_covariant_vjp
    public :: map_tetra_rt_contravariant
    public :: map_tetra_rt_contravariant_jvp
    public :: map_tetra_rt_contravariant_vjp
    public :: build_tetra_edge_dof_map
    public :: build_tetra_nedelec_basis_transform
    public :: build_tetra_nedelec_dof_map
    public :: interpolate_reference_tetra_nedelec
    public :: interpolate_physical_tetra_nedelec
    public :: interpolate_sampled_physical_tetra_nedelec
    public :: interpolate_sampled_physical_tetra_nedelec_jvp
    public :: interpolate_sampled_physical_tetra_nedelec_vjp
    public :: tetra_nedelec_interpolation_points
    public :: tetra_vector_sample_gradients_t
    public :: tetra_vector_samples_t
    public :: zero_tetra_vector_samples_like
    public :: evaluate_triangle_bdm
    public :: evaluate_triangle_bdm_jvp
    public :: evaluate_triangle_bdm_vjp
    public :: initialize_triangle_bdm
    public :: triangle_bdm_basis_t
    public :: triangle_bdm_dof_count
    public :: evaluate_triangle_raviart_thomas
    public :: evaluate_triangle_raviart_thomas_jvp
    public :: evaluate_triangle_raviart_thomas_vjp
    public :: initialize_triangle_raviart_thomas
    public :: triangle_rt_basis_t
    public :: triangle_rt_dof_count
    public :: apply_edge_moment_orientation
    public :: map_triangle_nedelec_covariant
    public :: map_triangle_nedelec_covariant_jvp
    public :: map_triangle_nedelec_covariant_vjp
    public :: map_triangle_rt_contravariant
    public :: map_triangle_rt_contravariant_jvp
    public :: map_triangle_rt_contravariant_vjp
    public :: build_triangle_discrete_gradient
    public :: build_triangle_trimmed_dof_map
    public :: build_triangle_full_vector_dof_map
    public :: build_triangle_discontinuous_dof_map
    public :: interpolate_triangle_nedelec
    public :: evaluate_triangle_nedelec_interpolant
    public :: evaluate_triangle_nedelec_interpolant_at_point
    public :: evaluate_triangle_nedelec_interpolant_at_point_jvp
    public :: evaluate_triangle_nedelec_interpolant_at_point_vjp
    public :: evaluate_triangle_nedelec_interpolant_jvp
    public :: evaluate_triangle_nedelec_interpolant_vjp
    public :: evaluate_triangle_rt_interpolant
    public :: evaluate_triangle_rt_interpolant_at_point
    public :: evaluate_triangle_rt_interpolant_at_point_jvp
    public :: evaluate_triangle_rt_interpolant_at_point_vjp
    public :: evaluate_triangle_rt_interpolant_jvp
    public :: evaluate_triangle_rt_interpolant_vjp
    public :: interpolate_triangle_rt
    public :: evaluate_triangle_bdm_interpolant
    public :: evaluate_triangle_bdm_interpolant_at_point
    public :: evaluate_triangle_bdm_interpolant_at_point_jvp
    public :: evaluate_triangle_bdm_interpolant_at_point_vjp
    public :: evaluate_triangle_bdm_interpolant_jvp
    public :: evaluate_triangle_bdm_interpolant_vjp
    public :: evaluate_triangle_nedelec_second_kind_interpolant
    public :: evaluate_triangle_nedelec_second_kind_interpolant_jvp
    public :: evaluate_triangle_nedelec_second_kind_interpolant_vjp
    public :: evaluate_triangle_nedelec_second_interpolant_at_point
    public :: evaluate_triangle_nedelec_second_interpolant_at_point_jvp
    public :: evaluate_triangle_nedelec_second_interpolant_at_point_vjp
    public :: interpolate_triangle_bdm
    public :: interpolate_triangle_nedelec_second_kind
    public :: assemble_helmholtz_adjoint_double_layer_constant
    public :: assemble_helmholtz_double_layer_constant
    public :: assemble_helmholtz_double_layer_mixed_linear
    public :: assemble_helmholtz_hypersingular_linear
    public :: assemble_helmholtz_single_layer_constant
    public :: assemble_triangle_bdm_div_mass_csc
    public :: assemble_triangle_bdm_div_mass_csc_jvp
    public :: assemble_triangle_bdm_div_mass_csc_vjp
    public :: assemble_triangle_bdm_div_mass_element
    public :: assemble_triangle_bdm_div_mass_element_jvp
    public :: assemble_triangle_bdm_div_mass_element_vjp
    public :: assemble_triangle_bdm_vector_load_samples
    public :: assemble_triangle_bdm_vector_load_samples_jvp
    public :: assemble_triangle_bdm_vector_load_samples_vjp
    public :: assemble_triangle_nedelec_second_curl_mass_csc
    public :: assemble_triangle_nedelec_second_curl_mass_csc_jvp
    public :: assemble_triangle_nedelec_second_curl_mass_csc_vjp
    public :: assemble_triangle_nedelec_second_curl_mass_element
    public :: assemble_triangle_nedelec_second_curl_mass_element_jvp
    public :: assemble_triangle_nedelec_second_curl_mass_element_vjp
    public :: assemble_triangle_nedelec_second_vector_load_samples
    public :: assemble_triangle_nedelec_second_vector_load_samples_jvp
    public :: assemble_triangle_nedelec_second_vector_load_samples_vjp
    public :: assemble_triangle_nedelec_curl_mass_csc
    public :: assemble_triangle_nedelec_curl_mass_csc_jvp
    public :: assemble_triangle_nedelec_curl_mass_csc_vjp
    public :: assemble_triangle_nedelec_curl_mass_element
    public :: assemble_triangle_nedelec_curl_mass_element_jvp
    public :: assemble_triangle_nedelec_curl_mass_element_vjp
    public :: assemble_triangle_nedelec_curl_csc
    public :: assemble_triangle_nedelec_vector_load_samples
    public :: assemble_triangle_nedelec_vector_load_samples_jvp
    public :: assemble_triangle_nedelec_vector_load_samples_vjp
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
    public :: assemble_tetra_rt_div_mass_csc
    public :: assemble_tetra_rt_div_mass_csc_jvp
    public :: assemble_tetra_rt_div_mass_csc_vjp
    public :: assemble_tetra_rt_div_mass_element
    public :: assemble_tetra_rt_div_mass_element_jvp
    public :: assemble_tetra_rt_div_mass_element_vjp
    public :: assemble_tetra_rt_divergence_csc
    public :: assemble_tetra_rt_vector_load_samples
    public :: assemble_tetra_rt_vector_load_samples_jvp
    public :: assemble_tetra_rt_vector_load_samples_vjp
    public :: assemble_tetra_lagrange_stiffness_csc
    public :: assemble_tetra_lagrange_stiffness_csc_jvp
    public :: assemble_tetra_lagrange_stiffness_csc_vjp
    public :: assemble_tetra_lagrange_stiffness_element
    public :: assemble_tetra_lagrange_stiffness_element_jvp
    public :: assemble_tetra_lagrange_stiffness_element_vjp
    public :: assemble_tetra_lagrange_curvilinear_pml_element
    public :: assemble_tetra_lagrange_curvilinear_pml_element_jvp
    public :: assemble_tetra_lagrange_curvilinear_pml_element_vjp
    public :: assemble_tetra_lagrange_curvilinear_pml_csc
    public :: assemble_tetra_lagrange_curvilinear_pml_csc_jvp
    public :: assemble_tetra_lagrange_curvilinear_pml_csc_vjp
    public :: assemble_tetra_lagrange_geometry_pml_csc
    public :: assemble_tetra_lagrange_geometry_pml_csc_jvp
    public :: assemble_tetra_lagrange_geometry_pml_csc_vjp
    public :: assemble_tetra_lagrange_scalar_load
    public :: assemble_tetra_lagrange_scalar_load_samples
    public :: assemble_tetra_lagrange_scalar_load_samples_jvp
    public :: assemble_tetra_lagrange_scalar_load_samples_vjp
    public :: assemble_tetra_nedelec_curl_mass_csc
    public :: assemble_tetra_nedelec_curl_mass_csc_jvp
    public :: assemble_tetra_nedelec_curl_mass_csc_vjp
    public :: assemble_tetra_nedelec_curl_mass_element
    public :: assemble_tetra_nedelec_curl_mass_element_jvp
    public :: assemble_tetra_nedelec_curl_mass_element_vjp
    public :: assemble_tetra_nedelec_pml_element
    public :: assemble_tetra_nedelec_pml_element_jvp
    public :: assemble_tetra_nedelec_pml_element_vjp
    public :: assemble_tetra_nedelec_curvilinear_pml_element
    public :: assemble_tetra_nedelec_curvilinear_pml_element_jvp
    public :: assemble_tetra_nedelec_curvilinear_pml_element_vjp
    public :: assemble_tetra_nedelec_curvilinear_pml_csc
    public :: assemble_tetra_nedelec_curvilinear_pml_csc_jvp
    public :: assemble_tetra_nedelec_curvilinear_pml_csc_vjp
    public :: assemble_tetra_nedelec_geometry_pml_csc
    public :: assemble_tetra_nedelec_geometry_pml_csc_jvp
    public :: assemble_tetra_nedelec_geometry_pml_csc_vjp
    public :: assemble_tetra_nedelec_pml_csc
    public :: assemble_tetra_nedelec_pml_csc_jvp
    public :: assemble_tetra_nedelec_pml_csc_vjp
    public :: assemble_tetra_nedelec_weighted_csc
    public :: assemble_tetra_nedelec_vector_load
    public :: assemble_tetra_nedelec_vector_load_samples
    public :: assemble_tetra_nedelec_vector_load_samples_jvp
    public :: assemble_tetra_nedelec_vector_load_samples_vjp
    public :: assemble_helmholtz_single_layer_linear
    public :: assemble_laplace_adjoint_double_layer_constant
    public :: assemble_laplace_double_layer_constant
    public :: assemble_laplace_double_layer_mixed_linear
    public :: assemble_laplace_hypersingular_linear
    public :: assemble_laplace_single_layer_constant
    public :: assemble_laplace_symmetric_coupling_p1_p0
    public :: assemble_helmholtz_symmetric_coupling_p1_p0
    public :: apply_curved_acoustic_displacement_ntd_2d
    public :: assemble_curved_acoustic_displacement_ntd_form_2d
    public :: solve_helmholtz_symmetric_coupling_p1_p0
    public :: solve_laplace_symmetric_coupling_p1_p0
    public :: interpolate_nedelec_edge_dofs, interpolate_rt_edge_dofs
    public :: interpolate_axisymmetric_rt_edge_dofs
    public :: evaluate_rt_field_2d
    public :: reconstruct_axisymmetric_fourier_toroidal, rt_l2_norm
    public :: build_bspline_polar_h1_extraction
    public :: build_bspline_polar_feec_2d_operators
    public :: build_bspline_polar_feec_2d_extractions
    public :: evaluate_periodic_bspline_basis
    public :: restrict_bspline_polar_operator_csc
    public :: assemble_bspline_polar_hcurl_operator_csc
    public :: assemble_bspline_polar_h1_operator_csc
    public :: assemble_bspline_polar_l2_mass_csc
    public :: sparse_direct_factor_t, sparse_direct_factor_csc
    public :: sparse_direct_factor_adjoint_csc
    public :: sparse_direct_factor_transpose_csc, sparse_direct_solve_csc
    public :: sparse_direct_solve_constrained
    public :: sparse_direct_solve_constrained_jvp
    public :: sparse_direct_solve_constrained_vjp
    public :: sparse_direct_solve_tree_cotree
    public :: sparse_direct_solve_tree_cotree_jvp
    public :: sparse_direct_solve_tree_cotree_vjp
    public :: sparse_direct_solve_factored
    public :: sparse_direct_solve_factored_jvp
    public :: sparse_direct_solve_factored_vjp, sparse_direct_free
    public :: sparse_matrix_t, sparse_from_dense, spmv, spmv_jvp, spmv_vjp

    ! Advanced solver types and functions
    public :: solver_options_t, solver_stats_t
    public :: solver_resource_budget_t
    public :: initialize_solver_resource_budget
    public :: validate_solver_resource_budget
    public :: evaluate_solver_resource_usage
    public :: solver_options
    public :: solve_sparse
    public :: cg_solve, pcg_solve, bicgstab_solve, gmres_solve
    public :: cg_solve_jvp, cg_solve_vjp
    public :: pcg_solve_jvp, pcg_solve_vjp
    public :: bicgstab_solve_jvp, bicgstab_solve_vjp
    public :: gmres_solve_jvp, gmres_solve_vjp
    public :: jacobi_preconditioner, ilu_preconditioner, ichol_preconditioner
    public :: ichol_controlled_preconditioner

    ! Field-coordinate-independent (FCI) support operators
    public :: assemble_fci_parallel_gradient_csc
    public :: assemble_fci_parallel_support_divergence_csc
    public :: apply_fci_parallel_gradient
    public :: apply_fci_parallel_diffusion
    public :: apply_fci_parallel_diffusion_jvp
    public :: apply_fci_parallel_diffusion_vjp
    public :: apply_fci_parallel_diffusion_field_vjp
    public :: compute_fci_parallel_diffusion_diagonal
    public :: compute_fci_anisotropic_diffusion_diagonal
    public :: apply_fci_parallel_jacobi_preconditioner
    public :: apply_fci_anisotropic_jacobi_preconditioner
    public :: apply_fci_anisotropic_diffusion
    public :: apply_fci_anisotropic_diffusion_field_vjp
    public :: apply_fci_parallel_gradient_jvp
    public :: apply_fci_parallel_gradient_vjp
    public :: apply_fci_plane_two_level_vcycle
    public :: factor_fci_plane_coarse_operator
    public :: apply_fci_plane_two_level_vcycle_factored
    public :: apply_fci_plane_multilevel_vcycle
    public :: apply_fci_plane_multilevel_wcycle
    public :: apply_fci_plane_two_level_vcycles
    public :: apply_fci_plane_two_level_vcycles_ragged
    public :: apply_fci_additive_field_split_preconditioner
    public :: trace_fci_field_line_rk4
    public :: trace_fci_field_line_rk4_jvp
    public :: build_fci_linear_interpolation_map_1d
    public :: build_fci_linear_interpolation_map_1d_jvp
    public :: build_fci_linear_interpolation_map_1d_vjp
    public :: build_fci_quadratic_interpolation_map_1d
    public :: build_fci_quadratic_interpolation_map_1d_jvp
    public :: build_fci_quadratic_interpolation_map_1d_vjp
    public :: build_fci_cubic_interpolation_map_1d
    public :: build_fci_cubic_interpolation_map_1d_jvp
    public :: build_fci_cubic_interpolation_map_1d_vjp
    public :: build_fci_quartic_interpolation_map_1d
    public :: build_fci_quartic_interpolation_map_1d_jvp
    public :: build_fci_quartic_interpolation_map_1d_vjp
    public :: build_fci_quintic_interpolation_map_1d
    public :: build_fci_quintic_interpolation_map_1d_jvp
    public :: build_fci_quintic_interpolation_map_1d_vjp
    public :: build_fci_sextic_interpolation_map_1d
    public :: build_fci_sextic_interpolation_map_1d_jvp
    public :: build_fci_sextic_interpolation_map_1d_vjp
    public :: build_fci_quadratic_interpolation_maps_1d
    public :: build_fci_quadratic_interpolation_maps_1d_jvp
    public :: build_fci_quadratic_interpolation_maps_1d_vjp
    public :: build_fci_bilinear_interpolation_map_2d
    public :: build_fci_bilinear_interpolation_map_2d_jvp
    public :: build_fci_bilinear_interpolation_map_2d_vjp
    public :: build_fci_bilinear_interpolation_maps_2d
    public :: build_fci_bilinear_interpolation_maps_2d_jvp
    public :: build_fci_bilinear_interpolation_maps_2d_vjp
    public :: build_fci_triangle_interpolation_map_2d
    public :: build_fci_triangle_interpolation_map_2d_jvp
    public :: build_fci_triangle_interpolation_map_2d_vjp
    public :: build_fci_triangle_interpolation_maps_2d
    public :: build_fci_triangle_interpolation_maps_2d_jvp
    public :: build_fci_triangle_interpolation_maps_2d_vjp
    public :: compute_fci_staggered_flux_box_volumes
    public :: compute_fci_staggered_flux_box_volumes_jvp
    public :: compute_fci_staggered_flux_box_volumes_vjp
    public :: compute_fci_quadrilateral_cell_areas_2d
    public :: compute_fci_quadrilateral_cell_areas_2d_jvp
    public :: compute_fci_quadrilateral_cell_areas_2d_vjp
    public :: compute_fci_polygon_cell_areas_2d
    public :: compute_fci_polygon_cell_areas_2d_jvp
    public :: compute_fci_polygon_cell_areas_2d_vjp
    public :: compute_fci_curved_polygon_cell_areas_2d
    public :: compute_fci_curved_polygon_cell_areas_2d_jvp
    public :: compute_fci_curved_polygon_cell_areas_2d_vjp
    public :: compute_fci_cubic_curved_polygon_cell_areas_2d
    public :: compute_fci_cubic_curved_polygon_cell_areas_2d_jvp
    public :: compute_fci_cubic_curved_polygon_cell_areas_2d_vjp
    public :: compute_fci_quartic_curved_polygon_cell_areas_2d
    public :: compute_fci_quartic_curved_polygon_cell_areas_2d_jvp
    public :: compute_fci_quartic_curved_polygon_cell_areas_2d_vjp
    public :: compute_fci_quintic_curved_polygon_cell_areas_2d
    public :: compute_fci_quintic_curved_polygon_cell_areas_2d_jvp
    public :: compute_fci_quintic_curved_polygon_cell_areas_2d_vjp
    public :: compute_fci_sextic_curved_polygon_cell_areas_2d
    public :: compute_fci_sextic_curved_polygon_cell_areas_2d_jvp
    public :: compute_fci_sextic_curved_polygon_cell_areas_2d_vjp
    public :: compute_fci_curved_quadrilateral_cell_areas_2d
    public :: compute_fci_curved_quadrilateral_cell_areas_2d_jvp
    public :: compute_fci_curved_quadrilateral_cell_areas_2d_vjp
    public :: advance_mixed_wave_midpoint
    public :: assemble_mixed_wave_midpoint_map
    public :: advance_mixed_wave_midpoint_jvp
    public :: advance_mixed_wave_midpoint_vjp
    public :: advance_mixed_wave_symplectic_euler
    public :: advance_mixed_wave_symplectic_euler_jvp
    public :: advance_mixed_wave_symplectic_euler_vjp
    public :: advance_mixed_wave_strang
    public :: advance_mixed_wave_strang_jvp
    public :: advance_mixed_wave_strang_vjp
    public :: advance_quadratic_avf
    public :: advance_quadratic_avf_jvp
    public :: advance_quadratic_avf_vjp
    public :: evaluate_cgl_pressure_tensor
    public :: evaluate_cgl_pressure_tensor_jvp
    public :: evaluate_cgl_pressure_tensor_vjp
    public :: evaluate_cgl_pressure_traction
    public :: evaluate_cgl_pressure_traction_jvp
    public :: evaluate_cgl_pressure_traction_vjp
    public :: evaluate_cgl_pressure_work
    public :: evaluate_cgl_pressure_work_jvp
    public :: evaluate_cgl_pressure_work_vjp
    public :: evaluate_field_aligned_flux
    public :: evaluate_field_aligned_flux_jvp
    public :: evaluate_field_aligned_flux_vjp
    public :: evaluate_field_aligned_constitutive_tensor
    public :: evaluate_field_aligned_constitutive_tensor_jvp
    public :: evaluate_field_aligned_constitutive_tensor_vjp
    public :: pullback_field_aligned_tensor
    public :: pullback_field_aligned_tensor_jvp
    public :: pullback_field_aligned_tensor_vjp
    public :: evaluate_tensor_power_split
    public :: evaluate_tensor_power_split_jvp
    public :: evaluate_tensor_power_split_vjp
    public :: spherical_harmonic
    public :: spherical_harmonic_theta_derivative
    public :: spherical_harmonic_phi_derivative
    public :: spherical_harmonic_product_coefficient
    public :: toroidal_p
    public :: toroidal_q
    public :: toroidal_p_derivative
    public :: toroidal_q_derivative
    public :: toroidal_p_second_derivative
    public :: toroidal_q_second_derivative
    public :: apply_nestor_fourier_response_map
    public :: apply_nestor_fourier_response_map_jvp
    public :: apply_nestor_fourier_response_map_vjp
    public :: evaluate_nestor_fourier_reciprocity_defect
    public :: compute_interface_scalar_jump_average
    public :: compute_interface_vector_traces
    public :: assemble_surface_delta_load
    public :: assemble_surface_delta_load_jvp
    public :: assemble_surface_delta_load_vjp
    public :: assemble_surface_vector_delta_load
    public :: assemble_surface_vector_delta_load_jvp
    public :: assemble_surface_vector_delta_load_vjp
    public :: assemble_interface_jump_penalty
    public :: assemble_interface_jump_penalty_jvp
    public :: assemble_interface_jump_penalty_vjp
    public :: assemble_symmetric_nitsche_interface
    public :: assemble_symmetric_nitsche_interface_jvp
    public :: assemble_symmetric_nitsche_interface_vjp
    public :: assemble_scalar_sipg_interface
    public :: assemble_scalar_sipg_interface_jvp
    public :: assemble_scalar_sipg_interface_vjp
    public :: assemble_vector_jump_penalty
    public :: assemble_vector_jump_penalty_jvp
    public :: assemble_vector_jump_penalty_vjp
    public :: assemble_vector_sipg_interface
    public :: assemble_vector_sipg_interface_jvp
    public :: assemble_vector_sipg_interface_vjp
    public :: assemble_hdg_static_condensation
    public :: assemble_hdg_static_condensation_jvp
    public :: assemble_hdg_static_condensation_vjp
    public :: assemble_hdg_global_skeleton
    public :: assemble_hdg_global_skeleton_jvp
    public :: assemble_hdg_global_skeleton_vjp
    public :: assemble_hdg_global_skeleton_csc
    public :: assemble_hdg_global_skeleton_csc_jvp
    public :: assemble_hdg_global_skeleton_csc_vjp
    public :: assemble_feec_exact_sequence
    public :: assemble_feec_exact_sequence_jvp
    public :: assemble_feec_exact_sequence_vjp
    public :: assemble_enriched_feec_sequence
    public :: assemble_enriched_feec_sequence_jvp
    public :: assemble_enriched_feec_sequence_vjp
    public :: assemble_broken_feec_sequence
    public :: assemble_broken_feec_sequence_jvp
    public :: assemble_broken_feec_sequence_vjp
    public :: assemble_fitted_trace_constraint
    public :: assemble_fitted_trace_constraint_jvp
    public :: assemble_fitted_trace_constraint_vjp
    public :: assemble_feec_commuting_projection
    public :: assemble_feec_commuting_projection_jvp
    public :: assemble_feec_commuting_projection_vjp
    public :: assemble_scalar_numerical_flux
    public :: assemble_scalar_numerical_flux_jvp
    public :: assemble_scalar_numerical_flux_vjp
    public :: NUMERICAL_FLUX_CENTRAL
    public :: NUMERICAL_FLUX_UPWIND
    public :: NUMERICAL_FLUX_LAX_FRIEDRICHS
    public :: assemble_vector_numerical_flux
    public :: assemble_vector_numerical_flux_jvp
    public :: assemble_vector_numerical_flux_vjp
    public :: assemble_vector_entropy_stable_flux
    public :: assemble_vector_entropy_stable_flux_jvp
    public :: assemble_vector_entropy_stable_flux_vjp
    public :: assemble_mortar_trace_coupling
    public :: assemble_mortar_trace_coupling_jvp
    public :: assemble_mortar_trace_coupling_vjp
    public :: assemble_geometry_mortar_trace_coupling
    public :: assemble_geometry_mortar_trace_coupling_jvp
    public :: assemble_geometry_mortar_trace_coupling_vjp
    public :: assemble_geometry_mortar_component_coupling
    public :: assemble_geometry_mortar_component_coupling_jvp
    public :: assemble_geometry_mortar_component_coupling_vjp
    public :: sample_physical_surface_geometry
    public :: sample_physical_surface_geometry_jvp
    public :: sample_physical_surface_geometry_vjp
    public :: evaluate_surface_vector_trace
    public :: evaluate_surface_vector_trace_jvp
    public :: evaluate_surface_vector_trace_vjp
    public :: assemble_fci_boundary_patch_mortar
    public :: assemble_fci_boundary_patch_mortar_jvp
    public :: assemble_fci_boundary_patch_mortar_vjp
    public :: evaluate_cgl_pressure_divergence
    public :: evaluate_cgl_pressure_divergence_jvp
    public :: evaluate_cgl_pressure_divergence_vjp
    public :: interchange_sample_set_t
    public :: initialize_interchange_samples
    public :: validate_interchange_samples
    public :: compare_interchange_samples
    public :: compare_interchange_samples_jvp
    public :: compare_interchange_samples_vjp
    public :: complex_interchange_sample_set_t
    public :: initialize_complex_interchange_samples
    public :: validate_complex_interchange_samples
    public :: compare_complex_interchange_samples
    public :: compare_complex_interchange_samples_jvp
    public :: compare_complex_interchange_samples_vjp
    public :: evaluate_weighted_complex_error
    public :: evaluate_weighted_complex_error_jvp
    public :: evaluate_weighted_complex_error_vjp
    public :: evaluate_weighted_reflection_coefficient
    public :: evaluate_weighted_reflection_coefficient_jvp
    public :: evaluate_weighted_reflection_coefficient_vjp
    public :: assemble_singular_layer_matching
    public :: assemble_singular_layer_matching_jvp
    public :: assemble_singular_layer_matching_vjp
    public :: OBJECTIVE_METADATA_KIND_EQUATION
    public :: OBJECTIVE_METADATA_KIND_OBJECTIVE
    public :: OBJECTIVE_METADATA_KIND_CONSTRAINT
    public :: OBJECTIVE_METADATA_UNSET_ID
    public :: equation_objective_metadata_t
    public :: initialize_equation_objective_metadata
    public :: validate_equation_objective_metadata
    public :: clear_equation_objective_metadata
    public :: evaluate_equation_objective_merit
    public :: evaluate_equation_objective_merit_jvp
    public :: evaluate_equation_objective_merit_vjp
    public :: evaluate_surface_shape_objective
    public :: evaluate_surface_shape_objective_jvp
    public :: evaluate_surface_shape_objective_vjp
    public :: evaluate_surface_integral_constraint
    public :: evaluate_surface_integral_constraint_jvp
    public :: evaluate_surface_integral_constraint_vjp
    public :: near_axis_diagnostic_metadata_t
    public :: evaluate_boozer_like_rotational_transform
    public :: evaluate_boozer_like_rotational_transform_jvp
    public :: evaluate_boozer_like_rotational_transform_vjp
    public :: evaluate_near_axis_diagnostic_metadata
    public :: evaluate_near_axis_diagnostic_metadata_jvp
    public :: evaluate_near_axis_diagnostic_metadata_vjp
    public :: CRITICAL_CONTOUR_LIMITER
    public :: CRITICAL_CONTOUR_NONE
    public :: CRITICAL_CONTOUR_SEPARATRIX
    public :: CRITICAL_POINT_DEGENERATE
    public :: CRITICAL_POINT_O_POINT
    public :: CRITICAL_POINT_REGULAR
    public :: CRITICAL_POINT_X_POINT
    public :: critical_point_metadata_t
    public :: evaluate_critical_point_metadata
    public :: evaluate_critical_point_metadata_jvp
    public :: evaluate_critical_point_metadata_vjp
    public :: validate_critical_point_metadata
    public :: linear_response_cross_metadata_t
    public :: initialize_linear_response_cross_metadata
    public :: validate_linear_response_cross_metadata
    public :: assemble_linear_response_eigen_cross_residual
    public :: assemble_linear_response_eigen_cross_residual_jvp
    public :: assemble_linear_response_eigen_cross_residual_vjp
    public :: equilibrium_interchange_t
    public :: equilibrium_normalization_t
    public :: initialize_equilibrium_interchange
    public :: validate_equilibrium_interchange
    public :: build_equilibrium_interchange_sample_set
    public :: oracle_manifest_schema_magic
    public :: oracle_manifest_schema_version
    public :: oracle_manifest_t
    public :: oracle_normalization_t
    public :: oracle_tolerance_t
    public :: oracle_timing_t
    public :: initialize_oracle_manifest
    public :: validate_oracle_manifest
    public :: read_oracle_manifest
    public :: write_oracle_manifest
    public :: linear_response_interchange_t
    public :: initialize_linear_response_interchange
    public :: validate_linear_response_interchange
    public :: evaluate_linear_response_diagnostics
    public :: assemble_linear_response_residual
    public :: assemble_linear_response_residual_jvp
    public :: assemble_linear_response_residual_vjp
    public :: assemble_linear_response_operator
    public :: assemble_linear_response_operator_jvp
    public :: assemble_linear_response_operator_vjp
    public :: LINEAR_RECIPROCITY_UNDECLARED
    public :: LINEAR_RECIPROCITY_TRANSPOSE
    public :: LINEAR_PASSIVITY_UNDECLARED
    public :: LINEAR_PASSIVITY_HERMITIAN_NONNEGATIVE
    public :: linear_perturbation_metadata_t
    public :: initialize_linear_perturbation_metadata
    public :: validate_linear_perturbation_metadata
    public :: assemble_linear_perturbation_operator
    public :: assemble_linear_perturbation_operator_jvp
    public :: assemble_linear_perturbation_operator_vjp
    public :: RESISTIVE_MHD_AMPERE
    public :: RESISTIVE_MHD_ANISOTROPIC_TRANSPORT
    public :: RESISTIVE_MHD_BLOCK_COUNT
    public :: RESISTIVE_MHD_FARADAY
    public :: RESISTIVE_MHD_FREE_BOUNDARY
    public :: RESISTIVE_MHD_MOMENTUM
    public :: RESISTIVE_MHD_PRESSURE
    public :: RESISTIVE_MHD_TENSOR
    public :: RESISTIVE_MHD_WALL
    public :: nonlinear_resistive_mhd_energy_ledger_t
    public :: assemble_nonlinear_resistive_mhd_residual
    public :: assemble_nonlinear_resistive_mhd_residual_jvp
    public :: assemble_nonlinear_resistive_mhd_residual_vjp
    public :: resistive_mhd_branch_history_t
    public :: resistive_mhd_branch_diagnostics_t
    public :: initialize_resistive_mhd_branch_history
    public :: validate_resistive_mhd_branch_history
    public :: evaluate_resistive_mhd_branch_diagnostics
    public :: compare_resistive_mhd_branch_histories
    public :: evaluate_resistive_mhd_branch_path_metric
    public :: evaluate_resistive_mhd_branch_path_metric_jvp
    public :: evaluate_resistive_mhd_branch_path_metric_vjp
    public :: linear_response_schema_magic
    public :: read_linear_response_interchange
    public :: write_linear_response_interchange
    public :: assemble_generalized_eigen_residual
    public :: assemble_generalized_eigen_residual_jvp
    public :: assemble_generalized_eigen_residual_vjp
    public :: assemble_beltrami_residual
    public :: assemble_beltrami_residual_jvp
    public :: assemble_beltrami_residual_vjp
    public :: beltrami_parity_t
    public :: compare_beltrami_two_region_residual
    public :: validate_beltrami_parity
    public :: validate_beltrami_resonance
    public :: beltrami_shell_parity_t
    public :: compare_beltrami_shell_residual
    public :: validate_beltrami_shell_parity
    public :: assemble_coupled_field_residual
    public :: assemble_coupled_field_residual_jvp
    public :: assemble_coupled_field_residual_vjp
    public :: assemble_block_2x2_residual
    public :: assemble_block_2x2_residual_jvp
    public :: assemble_block_2x2_residual_vjp
    public :: assemble_block_graph_residual
    public :: assemble_block_graph_residual_jvp
    public :: assemble_block_graph_residual_vjp
    public :: assemble_block_graph_csc
    public :: retained_field_split_t, retained_complex_field_split_t
    public :: factor_retained_field_split, factor_retained_complex_field_split
    public :: apply_retained_field_split, apply_retained_complex_field_split
    public :: apply_retained_field_split_jvp, apply_retained_complex_field_split_jvp
    public :: apply_retained_field_split_vjp, apply_retained_complex_field_split_vjp
    public :: free_retained_field_split, free_retained_complex_field_split
    public :: assemble_retained_coupled_schur
    public :: assemble_retained_coupled_schur_jvp
    public :: assemble_retained_coupled_schur_vjp
    public :: assemble_complex_block_graph_residual
    public :: assemble_complex_block_graph_residual_jvp
    public :: assemble_complex_block_graph_residual_vjp
    public :: boundary_operator_contract_t
    public :: initialize_boundary_operator_contract
    public :: initialize_boundary_operator_trace_metadata
    public :: validate_boundary_operator_contract
    public :: boundary_operator_parity_t
    public :: compare_boundary_operator_parity
    public :: compare_boundary_operator_parity_jvp
    public :: compare_boundary_operator_parity_vjp
    public :: evaluate_surface_trace_parity_ledger
    public :: evaluate_surface_trace_parity_ledger_jvp
    public :: evaluate_surface_trace_parity_ledger_vjp
    public :: validate_boundary_operator_parity
    public :: evaluate_weighted_boundary_response_diagnostics
    public :: evaluate_source_trace_map
    public :: evaluate_source_trace_map_jvp
    public :: evaluate_source_trace_map_vjp
    public :: evaluate_weighted_source_trace_reciprocity_defect
    public :: larger_domain_parity_t
    public :: compare_larger_domain_solution
    public :: compare_larger_domain_solution_jvp
    public :: validate_larger_domain_parity
    public :: BOUNDARY_OPERATOR_BACKEND_FEM
    public :: BOUNDARY_OPERATOR_BACKEND_BEM
    public :: BOUNDARY_OPERATOR_BACKEND_DTN
    public :: BOUNDARY_OPERATOR_BACKEND_PML
    public :: BOUNDARY_OPERATOR_BACKEND_NESTOR
    public :: BOUNDARY_OPERATOR_BACKEND_BIEST
    public :: BOUNDARY_OPERATOR_BACKEND_VIRTUAL_CASING
    public :: BOUNDARY_OPERATOR_BACKEND_USER
    public :: BOUNDARY_OPERATOR_TRACE_CHANNEL_SCALAR
    public :: BOUNDARY_OPERATOR_TRACE_CHANNEL_NORMAL
    public :: BOUNDARY_OPERATOR_TRACE_CHANNEL_TANGENTIAL
    public :: BOUNDARY_OPERATOR_TRACE_CHANNEL_MIXED
    public :: assemble_boundary_trace_residual
    public :: assemble_boundary_trace_residual_jvp
    public :: assemble_boundary_trace_residual_vjp
    public :: assemble_free_boundary_port_residual
    public :: assemble_free_boundary_port_residual_jvp
    public :: assemble_free_boundary_port_residual_vjp
    public :: assemble_free_boundary_source_response
    public :: assemble_free_boundary_source_response_jvp
    public :: assemble_free_boundary_source_response_vjp
    public :: apply_toroidal_modal_convolution
    public :: apply_toroidal_modal_convolution_jvp
    public :: apply_toroidal_modal_convolution_vjp
    public :: assemble_complex_coupled_field_residual
    public :: assemble_complex_coupled_field_residual_jvp
    public :: assemble_complex_coupled_field_residual_vjp
    public :: complex_low_rank_matrix_t
    public :: initialize_complex_low_rank_matrix
    public :: validate_complex_low_rank_matrix
    public :: compress_complex_matrix_cross
    public :: materialize_complex_low_rank_matrix
    public :: apply_complex_low_rank_matrix
    public :: apply_complex_low_rank_matrix_jvp
    public :: apply_complex_low_rank_matrix_vjp
    public :: advance_mixed_wave_wall_midpoint
    public :: advance_mixed_wave_wall_midpoint_jvp
    public :: advance_mixed_wave_wall_midpoint_vjp
    public :: evaluate_mixed_wave_wall_energy_balance
    public :: evaluate_batched_vector_enrichment_differential_3d
    public :: evaluate_batched_vector_enrichment_differential_3d_jvp
    public :: evaluate_batched_vector_enrichment_differential_3d_vjp
    public :: assemble_complex_boundary_trace_residual
    public :: assemble_complex_boundary_trace_residual_jvp
    public :: assemble_complex_boundary_trace_residual_vjp
    public :: assemble_mixed_elasticity_residual
    public :: assemble_mixed_elasticity_residual_jvp
    public :: assemble_mixed_elasticity_residual_vjp
    public :: assemble_elasticity_symmetry_constraint
    public :: assemble_elasticity_symmetry_constraint_jvp
    public :: assemble_elasticity_symmetry_constraint_vjp
    public :: assemble_symplectic_map_defect
    public :: assemble_symplectic_map_defect_jvp
    public :: assemble_symplectic_map_defect_vjp
    public :: assemble_glued_feec_sequence
    public :: assemble_glued_feec_sequence_jvp
    public :: assemble_glued_feec_sequence_vjp
    public :: assemble_glued_feec_sequence_csc
    public :: assemble_glued_feec_sequence_csc_jvp
    public :: assemble_glued_feec_sequence_csc_vjp
    public :: assemble_glued_feec_sequence_csc_compositions
    public :: assemble_glued_feec_sequence_csc_compositions_jvp
    public :: assemble_glued_feec_sequence_csc_compositions_vjp
    public :: build_multipatch_signed_dof_map

    ! Plotting interface
    public :: plot

end module fortfem_api
