module fortfem_api
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
    use fortfem_assembly_full_vector_arbitrary_order_2d, only: &
        assemble_triangle_bdm_div_mass_csc, &
        assemble_triangle_bdm_div_mass_element, &
        assemble_triangle_nedelec_second_curl_mass_csc, &
        assemble_triangle_nedelec_second_curl_mass_element
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
        assemble_tetra_nedelec_pml_csc, &
        assemble_tetra_nedelec_pml_csc_jvp, &
        assemble_tetra_nedelec_pml_csc_vjp, &
        assemble_tetra_nedelec_vector_load, &
        assemble_tetra_nedelec_vector_load_samples, &
        assemble_tetra_nedelec_vector_load_samples_jvp, &
        assemble_tetra_nedelec_vector_load_samples_vjp, &
        assemble_tetra_nedelec_weighted_csc
    use fortfem_kinds
    use fortfem_bspline_polar, only: &
        build_bspline_polar_feec_2d_extractions, &
        build_bspline_polar_feec_2d_operators, &
        build_bspline_polar_h1_extraction, evaluate_periodic_bspline_basis
    use fortfem_bspline_multipatch, only: &
        BSPLINE_FACE_X_MAX, BSPLINE_FACE_X_MIN, BSPLINE_FACE_Y_MAX, &
        BSPLINE_FACE_Y_MIN, BSPLINE_FACE_Z_MAX, BSPLINE_FACE_Z_MIN, &
        build_bspline_feec_2d_interface_dofs, &
        build_bspline_feec_2d_two_patch_maps, &
        build_bspline_feec_3d_interface_dofs, &
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
    use fortfem_torus_curved_panel, only: evaluate_torus_curved_panel
    use fortfem_sphere_surface_mesh, only: generate_sphere_surface_mesh
    use fortfem_sphere_curved_panel, only: &
        evaluate_sphere_curved_panel, invert_sphere_curved_panel
    use fortfem_maxwell_sphere_curved_rwg, only: &
        assemble_maxwell_sphere_curved_efie_rwg_3d, &
        assemble_maxwell_sphere_curved_efie_imaginary_rwg_3d, &
        assemble_maxwell_sphere_curved_efie_bc_imaginary_3d, &
        assemble_maxwell_sphere_curved_regularized_cfie_rwg_3d, &
        assemble_maxwell_sphere_curved_regularized_cfie_rhs_rwg_3d, &
        assemble_maxwell_sphere_curved_mfie_exterior_trace_rwg_rbc_3d, &
        assemble_maxwell_sphere_curved_mfie_rwg_rbc_3d, &
        assemble_maxwell_sphere_curved_plane_wave_rhs_rwg_3d, &
        assemble_maxwell_sphere_curved_potential_operators_rwg_3d, &
        assemble_maxwell_sphere_curved_rwg_mass_matrix, &
        assemble_maxwell_sphere_curved_rwg_rbc_pairing, &
        assemble_maxwell_sphere_curved_vector_potential_rwg_3d, &
        evaluate_maxwell_sphere_curved_far_field_rwg_3d, &
        evaluate_maxwell_sphere_curved_localized_rwg_basis, &
        evaluate_maxwell_sphere_curved_magnetic_field_rwg_3d, &
        evaluate_maxwell_sphere_curved_rwg_basis, &
        integrate_maxwell_sphere_curved_adjacent_rwg_pair_3d, &
        integrate_maxwell_sphere_curved_coincident_rwg_pair_3d, &
        solve_maxwell_pec_sphere_curved_efie_rwg_3d, &
        solve_maxwell_pec_sphere_curved_regularized_cfie_rwg_3d
    use fortfem_barycentric_surface_refinement, only: &
        barycentric_refine_surface_mesh, &
        barycentric_refine_torus_surface_mesh
    use fortfem_maxwell_localized_rwg_surface, only: &
        evaluate_maxwell_localized_rwg_basis
    use fortfem_maxwell_bc_surface, only: &
        assemble_maxwell_rwg_rbc_pairing, build_maxwell_bc_transformation
    use fortfem_maxwell_magnetic_rwg_3d, only: &
        evaluate_maxwell_magnetic_field_rwg_3d
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
    use fortfem_maxwell_torus_curved_rwg, only: &
        assemble_maxwell_torus_curved_efie_rwg_3d, &
        assemble_maxwell_torus_curved_efie_imaginary_rwg_3d, &
        assemble_maxwell_torus_curved_efie_bc_imaginary_3d, &
        assemble_maxwell_torus_curved_mfie_offset_trace_rwg_rbc_3d, &
        assemble_maxwell_torus_curved_mfie_rwg_rbc_3d, &
        assemble_maxwell_torus_curved_plane_wave_rhs_rwg_3d, &
        assemble_maxwell_torus_curved_plane_wave_rhs_bc_3d, &
        assemble_maxwell_torus_curved_potential_operators_rwg_3d, &
        assemble_maxwell_torus_curved_regularized_cfie_rwg_3d, &
        assemble_maxwell_torus_curved_rwg_mass_matrix, &
        assemble_maxwell_torus_curved_rwg_rbc_pairing, &
        evaluate_maxwell_torus_curved_far_field_rwg_3d, &
        evaluate_maxwell_torus_curved_localized_rwg_basis, &
        evaluate_maxwell_torus_curved_rwg_basis, &
        integrate_maxwell_torus_curved_adjacent_rwg_pair_3d, &
        integrate_maxwell_torus_curved_coincident_rwg_pair_3d, &
        solve_maxwell_pec_torus_curved_regularized_cfie_rwg_3d, &
        solve_maxwell_pec_torus_curved_regularized_cfie_rwg_multiple_3d
    use fortfem_toroidal_coordinates, only: cartesian_to_toroidal, &
        toroidal_point_to_cartesian, toroidal_vector_to_cartesian
    use fortfem_laplace_representation_3d, only: &
        evaluate_laplace_representation_triangles_3d, &
        evaluate_laplace_representation_torus_curved_3d
    use fortfem_adaptive_surface_bem, only: &
        estimate_helmholtz_p0_two_level_residual_3d, &
        estimate_laplace_p0_two_level_residual_3d, mark_bem_dorfler, &
        refine_surface_mesh_marked
    use fortfem_laplace_torus_curved_bem_3d, only: &
        assemble_laplace_torus_curved_calderon_3d, &
        assemble_laplace_torus_curved_dtn_3d, &
        solve_laplace_bem_dtn_torus_curved_3d
    use fortfem_laplace_torus_curved_fem_bem_coupling_3d, only: &
        assemble_laplace_fem_bem_costabel_torus_curved_3d, &
        solve_laplace_fem_bem_costabel_torus_curved_3d
    use fortfem_laplace_galerkin_3d, only: &
        assemble_laplace_calderon_p1_p0_3d, &
        assemble_laplace_single_layer_p0_3d, &
        assemble_laplace_single_layer_p0_adaptive_3d, &
        solve_laplace_dirichlet_p0_3d
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
    use fortfem_helmholtz_torus_curved_bem_3d, only: &
        assemble_helmholtz_torus_curved_calderon_3d, &
        assemble_helmholtz_torus_curved_dtn_3d, &
        solve_helmholtz_bem_dtn_torus_curved_3d
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
    use fortfem_scalar_helmholtz_pml_slab_1d, only: &
        solve_scalar_helmholtz_pml_slab_1d
    use fortfem_scalar_helmholtz_pml_2d, only: &
        solve_scalar_helmholtz_pml_p1_2d
    use fortfem_scalar_helmholtz_pml_3d, only: &
        solve_scalar_helmholtz_pml_p1_3d
    use fortfem_toroidal_poisson_dtn, only: &
        evaluate_toroidal_harmonic_p, evaluate_toroidal_ampere_field_p, &
        toroidal_poisson_exterior_dtn_p
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
        initialize_triangle_nedelec_first_kind, triangle_nedelec_dof_count, &
        triangle_nedelec_first_kind_t
    use fortfem_triangle_nedelec_second_kind, only: &
        assignment(=), evaluate_triangle_nedelec_second_kind, &
        initialize_triangle_nedelec_second_kind, &
        triangle_nedelec_second_kind_dof_count, &
        triangle_nedelec_second_kind_t
    use fortfem_tetra_nedelec_first_order, only: &
        evaluate_tetra_nedelec_first_order
    use fortfem_tetra_nedelec_arbitrary_order, only: &
        assignment(=), evaluate_tetra_nedelec_first_kind, &
        initialize_tetra_nedelec_first_kind, tetra_nedelec_dof_count, &
        tetra_nedelec_first_kind_t
    use fortfem_tetra_discontinuous_arbitrary_order, only: &
        assignment(=), evaluate_tetra_discontinuous, &
        initialize_tetra_discontinuous, tetra_discontinuous_dof_count, &
        tetra_discontinuous_t
    use fortfem_tetra_discontinuous_projection, only: &
        project_physical_tetra_discontinuous
    use fortfem_tetra_lagrange_arbitrary_order, only: &
        assignment(=), evaluate_tetra_lagrange, initialize_tetra_lagrange, &
        tetra_lagrange_barycentric_indices, tetra_lagrange_dof_count, &
        tetra_lagrange_nodes, tetra_lagrange_t
    use fortfem_tetra_lagrange_global_dof_map, only: &
        build_tetra_lagrange_dof_map
    use fortfem_tetra_feec_operators, only: &
        build_tetra_discrete_curl, build_tetra_discrete_gradient
    use fortfem_tetra_rt_arbitrary_order, only: &
        assignment(=), evaluate_tetra_rt, initialize_tetra_rt, &
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
        interpolate_reference_tetra_nedelec
    use fortfem_tetra_rt_interpolation, only: &
        interpolate_physical_tetra_rt, interpolate_reference_tetra_rt
    use fortfem_triangle_rt_arbitrary_order, only: &
        assignment(=), evaluate_triangle_raviart_thomas, &
        initialize_triangle_raviart_thomas, triangle_rt_basis_t, &
        triangle_rt_dof_count
    use fortfem_triangle_bdm_arbitrary_order, only: &
        assignment(=), evaluate_triangle_bdm, initialize_triangle_bdm, &
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
        evaluate_triangle_nedelec_interpolant, &
        evaluate_triangle_nedelec_second_kind_interpolant, &
        evaluate_triangle_rt_interpolant, interpolate_triangle_bdm, &
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
        sparse_direct_solve_factored_jvp, sparse_direct_solve_factored_vjp
    use fortfem_mixed_poisson_2d, only: &
        solve_mixed_poisson_rt, solve_mixed_poisson_rt0, &
        solve_symbolic_mixed_poisson_rt
    use fortfem_magnetic_box_3d, only: solve_magnetic_box_3d
    ! Public arbitrary-order H(curl) solve, including optional homogeneous PEC.
    use fortfem_tetra_nedelec_solver_3d, only: &
        solve_tetra_nedelec_curl_mass, solve_tetra_nedelec_pml, &
        solve_tetra_nedelec_weighted_curl_mass
    ! Public arbitrary-order H(div) solve, including optional zero normal trace.
    use fortfem_tetra_rt_solver_3d, only: solve_tetra_rt_div_mass
    use fortfem_tetra_mixed_poisson_3d, only: &
        assemble_tetra_dg_source_load_samples, &
        assemble_tetra_dg_source_load_samples_jvp, &
        assemble_tetra_dg_source_load_samples_vjp, &
        solve_symbolic_tetra_mixed_poisson_rt
    use fortfem_tetra_lagrange_solver_3d, only: &
        assignment(=), evaluate_tetra_lagrange_solution, &
        evaluate_tetra_lagrange_solution_prepared, &
        initialize_tetra_lagrange_solution_evaluator, &
        solve_tetra_lagrange_diffusion_reaction, &
        solve_tetra_lagrange_poisson, tetra_lagrange_solution_evaluator_t
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
    public :: generate_sphere_surface_mesh
    public :: evaluate_sphere_curved_panel
    public :: invert_sphere_curved_panel
    public :: evaluate_maxwell_sphere_curved_rwg_basis
    public :: assemble_maxwell_sphere_curved_efie_rwg_3d
    public :: assemble_maxwell_sphere_curved_efie_imaginary_rwg_3d
    public :: assemble_maxwell_sphere_curved_efie_bc_imaginary_3d
    public :: assemble_maxwell_sphere_curved_regularized_cfie_rwg_3d
    public :: assemble_maxwell_sphere_curved_regularized_cfie_rhs_rwg_3d
    public :: assemble_maxwell_sphere_curved_mfie_exterior_trace_rwg_rbc_3d
    public :: assemble_maxwell_sphere_curved_mfie_rwg_rbc_3d
    public :: solve_maxwell_pec_sphere_curved_efie_rwg_3d
    public :: solve_maxwell_pec_sphere_curved_regularized_cfie_rwg_3d
    public :: assemble_maxwell_sphere_curved_rwg_mass_matrix
    public :: assemble_maxwell_sphere_curved_rwg_rbc_pairing
    public :: assemble_maxwell_sphere_curved_plane_wave_rhs_rwg_3d
    public :: assemble_maxwell_sphere_curved_potential_operators_rwg_3d
    public :: assemble_maxwell_sphere_curved_vector_potential_rwg_3d
    public :: evaluate_maxwell_sphere_curved_far_field_rwg_3d
    public :: evaluate_maxwell_sphere_curved_localized_rwg_basis
    public :: evaluate_maxwell_sphere_curved_magnetic_field_rwg_3d
    public :: integrate_maxwell_sphere_curved_coincident_rwg_pair_3d
    public :: integrate_maxwell_sphere_curved_adjacent_rwg_pair_3d
    public :: barycentric_refine_surface_mesh
    public :: barycentric_refine_torus_surface_mesh
    public :: evaluate_maxwell_localized_rwg_basis
    public :: build_maxwell_bc_transformation
    public :: assemble_maxwell_rwg_rbc_pairing
    public :: evaluate_maxwell_magnetic_field_rwg_3d
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
    public :: assemble_maxwell_torus_curved_rwg_rbc_pairing
    public :: assemble_maxwell_torus_curved_efie_rwg_3d
    public :: assemble_maxwell_torus_curved_efie_imaginary_rwg_3d
    public :: assemble_maxwell_torus_curved_efie_bc_imaginary_3d
    public :: assemble_maxwell_torus_curved_mfie_offset_trace_rwg_rbc_3d
    public :: assemble_maxwell_torus_curved_mfie_rwg_rbc_3d
    public :: assemble_maxwell_torus_curved_plane_wave_rhs_rwg_3d
    public :: assemble_maxwell_torus_curved_plane_wave_rhs_bc_3d
    public :: assemble_maxwell_torus_curved_potential_operators_rwg_3d
    public :: assemble_maxwell_torus_curved_regularized_cfie_rwg_3d
    public :: evaluate_maxwell_torus_curved_far_field_rwg_3d
    public :: evaluate_maxwell_torus_curved_localized_rwg_basis
    public :: evaluate_maxwell_torus_curved_rwg_basis
    public :: integrate_maxwell_torus_curved_adjacent_rwg_pair_3d
    public :: integrate_maxwell_torus_curved_coincident_rwg_pair_3d
    public :: solve_maxwell_pec_torus_curved_regularized_cfie_rwg_3d
    public :: solve_maxwell_pec_torus_curved_regularized_cfie_rwg_multiple_3d
    public :: cartesian_to_toroidal
    public :: toroidal_point_to_cartesian
    public :: toroidal_vector_to_cartesian
    public :: evaluate_laplace_representation_triangles_3d
    public :: evaluate_laplace_representation_torus_curved_3d
    public :: assemble_laplace_torus_curved_calderon_3d
    public :: assemble_laplace_torus_curved_dtn_3d
    public :: solve_laplace_bem_dtn_torus_curved_3d
    public :: assemble_laplace_fem_bem_costabel_torus_curved_3d
    public :: solve_laplace_fem_bem_costabel_torus_curved_3d
    public :: assemble_laplace_single_layer_p0_3d
    public :: assemble_laplace_single_layer_p0_adaptive_3d
    public :: assemble_laplace_calderon_p1_p0_3d
    public :: solve_laplace_dirichlet_p0_3d
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
    public :: assemble_helmholtz_torus_curved_calderon_3d
    public :: assemble_helmholtz_torus_curved_dtn_3d
    public :: solve_helmholtz_bem_dtn_torus_curved_3d
    public :: assemble_helmholtz_fem_bem_costabel_torus_curved_3d
    public :: solve_helmholtz_fem_bem_costabel_torus_curved_3d
    public :: assemble_helmholtz_single_layer_p0_3d
    public :: assemble_helmholtz_calderon_p1_p0_3d
    public :: assemble_helmholtz_fem_bem_costabel_3d
    public :: assemble_helmholtz_single_layer_p0_adaptive_3d
    public :: assemble_helmholtz_double_layer_p0_3d
    public :: evaluate_helmholtz_cfie_p0_3d
    public :: solve_helmholtz_cfie_p0_3d
    public :: solve_helmholtz_dirichlet_p0_3d
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
    public :: build_bspline_feec_2d_two_patch_maps
    public :: build_bspline_feec_2d_two_patch_operators_csc
    public :: build_bspline_feec_3d_interface_dofs
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
    public :: solve_tetra_nedelec_pml
    public :: solve_tetra_nedelec_weighted_curl_mass
    public :: solve_tetra_rt_div_mass
    public :: solve_symbolic_tetra_mixed_poisson_rt
    public :: assemble_tetra_dg_source_load_samples
    public :: assemble_tetra_dg_source_load_samples_jvp
    public :: assemble_tetra_dg_source_load_samples_vjp
    public :: evaluate_tetra_lagrange_solution
    public :: evaluate_tetra_lagrange_solution_prepared
    public :: initialize_tetra_lagrange_solution_evaluator
    public :: solve_tetra_lagrange_diffusion_reaction
    public :: solve_tetra_lagrange_poisson
    public :: solve_tetra_lagrange_state
    public :: solve_tetra_lagrange_state_jvp
    public :: solve_tetra_lagrange_state_vjp
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
    public :: circular_helmholtz_dtn_eigenvalue
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
    public :: build_cartesian_pml_element_stretch
    public :: build_cartesian_pml_element_stretch_jvp
    public :: build_cartesian_pml_element_stretch_vjp
    public :: solve_scalar_helmholtz_pml_slab_1d
    public :: solve_scalar_helmholtz_pml_p1_2d
    public :: solve_scalar_helmholtz_pml_p1_3d
    public :: evaluate_toroidal_harmonic_p
    public :: evaluate_toroidal_ampere_field_p
    public :: toroidal_poisson_exterior_dtn_p
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
    public :: initialize_triangle_nedelec_first_kind
    public :: triangle_nedelec_dof_count
    public :: triangle_nedelec_first_kind_t
    public :: evaluate_triangle_nedelec_second_kind
    public :: initialize_triangle_nedelec_second_kind
    public :: triangle_nedelec_second_kind_dof_count
    public :: triangle_nedelec_second_kind_t
    public :: evaluate_tetra_nedelec_first_order
    public :: evaluate_tetra_nedelec_first_kind
    public :: initialize_tetra_nedelec_first_kind
    public :: tetra_nedelec_dof_count
    public :: tetra_nedelec_first_kind_t
    public :: evaluate_tetra_discontinuous
    public :: initialize_tetra_discontinuous
    public :: tetra_discontinuous_dof_count
    public :: tetra_discontinuous_t
    public :: project_physical_tetra_discontinuous
    public :: evaluate_tetra_lagrange
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
    public :: evaluate_tetra_rt
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
    public :: evaluate_triangle_bdm
    public :: initialize_triangle_bdm
    public :: triangle_bdm_basis_t
    public :: triangle_bdm_dof_count
    public :: evaluate_triangle_raviart_thomas
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
    public :: evaluate_triangle_rt_interpolant
    public :: interpolate_triangle_rt
    public :: evaluate_triangle_bdm_interpolant
    public :: evaluate_triangle_nedelec_second_kind_interpolant
    public :: interpolate_triangle_bdm
    public :: interpolate_triangle_nedelec_second_kind
    public :: assemble_helmholtz_adjoint_double_layer_constant
    public :: assemble_helmholtz_double_layer_constant
    public :: assemble_helmholtz_double_layer_mixed_linear
    public :: assemble_helmholtz_hypersingular_linear
    public :: assemble_helmholtz_single_layer_constant
    public :: assemble_triangle_bdm_div_mass_csc
    public :: assemble_triangle_bdm_div_mass_element
    public :: assemble_triangle_nedelec_second_curl_mass_csc
    public :: assemble_triangle_nedelec_second_curl_mass_element
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
    public :: sparse_direct_solve_factored
    public :: sparse_direct_solve_factored_jvp
    public :: sparse_direct_solve_factored_vjp, sparse_direct_free

    ! Advanced solver types and functions
    public :: solver_options_t, solver_stats_t
    public :: solver_options
    public :: cg_solve, pcg_solve, bicgstab_solve, gmres_solve
    public :: jacobi_preconditioner, ilu_preconditioner

    ! Plotting interface
    public :: plot

end module fortfem_api
