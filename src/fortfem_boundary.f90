module fortfem_boundary
    !! Direct boundary-facing facade.
    !!
    !! Mesh boundary geometry remains available here for compatibility.  The
    !! open-boundary primitives are re-exported from their canonical modules
    !! so clients can depend on this small boundary surface without importing
    !! the umbrella ``fortfem_api`` (or any parity/comparison implementation).
    use fortfem_boundary_operator_contract, only: &
        BOUNDARY_OPERATOR_BACKEND_BEM, &
        BOUNDARY_OPERATOR_BACKEND_DTN, &
        BOUNDARY_OPERATOR_BACKEND_FEM, &
        BOUNDARY_OPERATOR_BACKEND_NESTOR, &
        BOUNDARY_OPERATOR_BACKEND_BIEST, &
        BOUNDARY_OPERATOR_BACKEND_PML, &
        BOUNDARY_OPERATOR_BACKEND_USER, &
        BOUNDARY_OPERATOR_BACKEND_VIRTUAL_CASING, &
        BOUNDARY_OPERATOR_TRACE_CHANNEL_MIXED, &
        BOUNDARY_OPERATOR_TRACE_CHANNEL_NORMAL, &
        BOUNDARY_OPERATOR_TRACE_CHANNEL_SCALAR, &
        BOUNDARY_OPERATOR_TRACE_CHANNEL_TANGENTIAL, &
        boundary_operator_contract_t, &
        initialize_boundary_operator_contract, &
        initialize_boundary_operator_trace_metadata, &
        validate_boundary_operator_contract
    use fortfem_circular_dtn_2d, only: &
        apply_circular_helmholtz_dtn, circular_helmholtz_dtn_eigenvalue
    use fortfem_circular_dtn_2d_ad, only: &
        apply_circular_helmholtz_dtn_jvp, apply_circular_helmholtz_dtn_vjp, &
        circular_helmholtz_dtn_eigenvalue_jvp, &
        circular_helmholtz_dtn_eigenvalue_vjp
    use fortfem_curved_acoustic_displacement_ntd_2d, only: &
        apply_curved_acoustic_displacement_ntd_2d, &
        assemble_curved_acoustic_displacement_ntd_form_2d
    use fortfem_kinds, only: dp
    use fortfem_laplace_boundary_operators_2d, only: &
        assemble_laplace_double_layer_mixed_linear, &
        assemble_laplace_hypersingular_linear, &
        assemble_laplace_single_layer_constant
    use fortfem_helmholtz_boundary_operators_2d, only: &
        assemble_helmholtz_adjoint_double_layer_constant, &
        assemble_helmholtz_double_layer_constant, &
        assemble_helmholtz_double_layer_mixed_linear, &
        assemble_helmholtz_hypersingular_linear, &
        assemble_helmholtz_single_layer_constant
    use fortfem_helmholtz_exterior_2d, only: &
        evaluate_helmholtz_combined_potential_adaptive_constant, &
        evaluate_helmholtz_combined_potential_constant, &
        solve_helmholtz_cfie_constant
    use fortfem_free_boundary_port_residual, only: &
        assemble_free_boundary_port_residual, &
        assemble_free_boundary_port_residual_jvp, &
        assemble_free_boundary_port_residual_vjp
    use fortfem_complex_boundary_trace_residual, only: &
        assemble_complex_boundary_trace_residual, &
        assemble_complex_boundary_trace_residual_jvp, &
        assemble_complex_boundary_trace_residual_vjp
    use fortfem_adaptive_surface_bem, only: &
        estimate_helmholtz_p0_two_level_residual_3d, &
        estimate_laplace_p0_two_level_residual_3d, &
        mark_bem_dorfler, refine_surface_mesh_marked
    use fortfem_elasticity_planar_acoustic_dtn_2d, only: &
        solve_elasticity_planar_acoustic_dtn_p1
    use fortfem_elasticity_curved_acoustic_ntd_2d, only: &
        solve_elasticity_curved_acoustic_ntd_p1
    use fortfem_laplace_symmetric_coupling_2d, only: &
        solve_laplace_symmetric_coupling_p1_p0
    use fortfem_planar_helmholtz_dtn, only: &
        apply_planar_helmholtz_dtn, &
        apply_planar_helmholtz_dtn_jvp, &
        apply_planar_helmholtz_dtn_vjp, &
        assemble_planar_helmholtz_dtn_form, &
        assemble_planar_helmholtz_dtn_form_jvp, &
        assemble_planar_helmholtz_dtn_form_vjp
    use fortfem_scalar_helmholtz_planar_dtn_2d, only: &
        solve_scalar_helmholtz_planar_dtn_p1, &
        solve_scalar_helmholtz_planar_dtn_p1_jvp, &
        solve_scalar_helmholtz_planar_dtn_p1_vjp
    use fortfem_spherical_helmholtz_dtn, only: &
        apply_spherical_helmholtz_dtn, &
        apply_spherical_helmholtz_dtn_jvp, &
        apply_spherical_helmholtz_dtn_vjp, &
        spherical_helmholtz_dtn_eigenvalue, &
        spherical_helmholtz_dtn_eigenvalue_jvp, &
        spherical_helmholtz_dtn_eigenvalue_vjp
    use fortfem_planar_maxwell_dtn, only: &
        apply_planar_maxwell_dtn, assemble_planar_maxwell_dtn_form, &
        assemble_planar_maxwell_dtn_form_jvp, &
        assemble_planar_maxwell_dtn_form_vjp
    use fortfem_planar_nedelec_maxwell_dtn, only: &
        assemble_planar_nedelec_maxwell_dtn_form, &
        build_planar_nedelec_trace_sampling, &
        pullback_planar_maxwell_dtn_form, &
        pullback_planar_maxwell_dtn_form_jvp, &
        pullback_planar_maxwell_dtn_form_vjp
    use fortfem_scalar_helmholtz_pml_slab_1d, only: &
        assemble_scalar_helmholtz_pml_slab_1d_matrix, &
        solve_scalar_helmholtz_pml_slab_1d
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
    use fortfem_maxwell_fem_bem_coupling_3d, only: &
        assemble_maxwell_fem_bem_boundary_matrix_3d, &
        assemble_maxwell_fem_bem_system_3d, &
        assemble_maxwell_fem_bem_torus_curved_system_3d, &
        assemble_maxwell_rwg_nedelec_coupling_3d, &
        solve_maxwell_fem_bem_system_3d, &
        solve_maxwell_fem_bem_torus_curved_system_3d
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
    use fortfem_maxwell_fem_bem_state_ad, only: &
        solve_maxwell_fem_bem_linear_state, &
        solve_maxwell_fem_bem_linear_state_jvp, &
        solve_maxwell_fem_bem_linear_state_vjp
    use fortfem_wall_response_condensation, only: &
        condense_wall_response_blocks, &
        condense_wall_response_blocks_jvp, &
        condense_wall_response_blocks_vjp
    use fortfem_maxwell_torus_curved_rwg, only: &
        evaluate_maxwell_torus_curved_far_field_rwg_3d, &
        evaluate_maxwell_torus_curved_magnetic_field_rwg_3d, &
        solve_maxwell_pec_torus_curved_regularized_cfie_rwg_multiple_3d
    use fortfem_helmholtz_galerkin_3d, only: &
        assemble_helmholtz_single_layer_p0_adaptive_3d, &
        assemble_helmholtz_single_layer_p0_3d, &
        evaluate_helmholtz_cfie_p0_3d, solve_helmholtz_cfie_p0_3d, &
        solve_helmholtz_dirichlet_p0_3d
    use fortfem_helmholtz_panel_pair_3d, only: &
        assemble_helmholtz_single_layer_p0_3d_jvp, &
        assemble_helmholtz_single_layer_p0_3d_vjp
    use fortfem_helmholtz_hierarchical_3d, only: &
        apply_helmholtz_cfie_p0_hierarchical_3d, &
        apply_helmholtz_single_layer_p0_hierarchical_3d, &
        solve_helmholtz_cfie_p0_hierarchical_3d, &
        solve_helmholtz_dirichlet_p0_hierarchical_3d
    use fortfem_helmholtz_representation_3d, only: &
        evaluate_helmholtz_representation_triangles_3d, &
        evaluate_helmholtz_representation_torus_curved_3d
    use fortfem_laplace_galerkin_3d, only: &
        assemble_laplace_single_layer_p0_3d, &
        assemble_laplace_single_layer_p0_adaptive_3d, &
        solve_laplace_dirichlet_p0_3d
    use fortfem_laplace_hierarchical_3d, only: &
        apply_laplace_single_layer_p0_hierarchical_3d
    use fortfem_laplace_representation_3d, only: &
        evaluate_laplace_representation_triangles_3d, &
        evaluate_laplace_representation_torus_curved_3d
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
    use fortfem_torus_curved_panel, only: evaluate_torus_curved_panel
    use fortfem_helmholtz_representation_3d, only: &
        evaluate_helmholtz_representation_torus_curved_3d
    use fortfem_helmholtz_representation_ad_3d, only: &
        evaluate_helmholtz_representation_torus_curved_3d_geometry_jvp, &
        evaluate_helmholtz_representation_torus_curved_3d_geometry_vjp
    use fortfem_laplace_torus_curved_bem_3d, only: &
        assemble_laplace_torus_curved_calderon_3d, &
        assemble_laplace_torus_curved_dtn_3d, &
        solve_laplace_bem_dtn_torus_curved_3d
    use fortfem_laplace_torus_curved_bem_ad_3d, only: &
        assemble_laplace_torus_curved_dtn_3d_geometry_jvp, &
        assemble_laplace_torus_curved_dtn_3d_geometry_vjp
    use fortfem_laplace_torus_curved_fem_bem_coupling_3d, only: &
        solve_laplace_fem_bem_costabel_torus_curved_3d
    implicit none

    private
    public :: boundary_t
    public :: BOUNDARY_OPERATOR_BACKEND_BEM
    public :: BOUNDARY_OPERATOR_BACKEND_DTN
    public :: BOUNDARY_OPERATOR_BACKEND_FEM
    public :: BOUNDARY_OPERATOR_BACKEND_NESTOR
    public :: BOUNDARY_OPERATOR_BACKEND_BIEST
    public :: BOUNDARY_OPERATOR_BACKEND_PML
    public :: BOUNDARY_OPERATOR_BACKEND_USER
    public :: BOUNDARY_OPERATOR_BACKEND_VIRTUAL_CASING
    public :: BOUNDARY_OPERATOR_TRACE_CHANNEL_MIXED
    public :: BOUNDARY_OPERATOR_TRACE_CHANNEL_NORMAL
    public :: BOUNDARY_OPERATOR_TRACE_CHANNEL_SCALAR
    public :: BOUNDARY_OPERATOR_TRACE_CHANNEL_TANGENTIAL
    public :: boundary_operator_contract_t
    public :: initialize_boundary_operator_contract
    public :: initialize_boundary_operator_trace_metadata
    public :: validate_boundary_operator_contract
    public :: apply_planar_helmholtz_dtn
    public :: apply_planar_helmholtz_dtn_jvp
    public :: apply_planar_helmholtz_dtn_vjp
    public :: assemble_planar_helmholtz_dtn_form
    public :: assemble_planar_helmholtz_dtn_form_jvp
    public :: assemble_planar_helmholtz_dtn_form_vjp
    public :: solve_scalar_helmholtz_planar_dtn_p1
    public :: solve_scalar_helmholtz_planar_dtn_p1_jvp
    public :: solve_scalar_helmholtz_planar_dtn_p1_vjp
    public :: apply_spherical_helmholtz_dtn
    public :: apply_spherical_helmholtz_dtn_jvp
    public :: apply_spherical_helmholtz_dtn_vjp
    public :: spherical_helmholtz_dtn_eigenvalue
    public :: spherical_helmholtz_dtn_eigenvalue_jvp
    public :: spherical_helmholtz_dtn_eigenvalue_vjp
    public :: apply_circular_helmholtz_dtn
    public :: apply_circular_helmholtz_dtn_jvp
    public :: apply_circular_helmholtz_dtn_vjp
    public :: circular_helmholtz_dtn_eigenvalue
    public :: circular_helmholtz_dtn_eigenvalue_jvp
    public :: circular_helmholtz_dtn_eigenvalue_vjp
    public :: apply_curved_acoustic_displacement_ntd_2d
    public :: assemble_curved_acoustic_displacement_ntd_form_2d
    public :: estimate_helmholtz_p0_two_level_residual_3d
    public :: estimate_laplace_p0_two_level_residual_3d
    public :: mark_bem_dorfler
    public :: refine_surface_mesh_marked
    public :: solve_elasticity_planar_acoustic_dtn_p1
    public :: solve_elasticity_curved_acoustic_ntd_p1
    public :: apply_planar_maxwell_dtn
    public :: assemble_planar_maxwell_dtn_form
    public :: assemble_planar_maxwell_dtn_form_jvp
    public :: assemble_planar_maxwell_dtn_form_vjp
    public :: assemble_planar_nedelec_maxwell_dtn_form
    public :: build_planar_nedelec_trace_sampling
    public :: pullback_planar_maxwell_dtn_form
    public :: pullback_planar_maxwell_dtn_form_jvp
    public :: pullback_planar_maxwell_dtn_form_vjp
    public :: assemble_scalar_helmholtz_pml_slab_1d_matrix
    public :: solve_scalar_helmholtz_pml_slab_1d
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
    public :: assemble_laplace_single_layer_constant
    public :: assemble_laplace_double_layer_mixed_linear
    public :: assemble_laplace_hypersingular_linear
    public :: assemble_helmholtz_double_layer_constant
    public :: assemble_helmholtz_adjoint_double_layer_constant
    public :: assemble_helmholtz_double_layer_mixed_linear
    public :: assemble_helmholtz_hypersingular_linear
    public :: assemble_helmholtz_single_layer_constant
    public :: evaluate_helmholtz_combined_potential_constant
    public :: evaluate_helmholtz_combined_potential_adaptive_constant
    public :: solve_helmholtz_cfie_constant
    public :: solve_laplace_symmetric_coupling_p1_p0
    public :: assemble_free_boundary_port_residual
    public :: assemble_free_boundary_port_residual_jvp
    public :: assemble_free_boundary_port_residual_vjp
    public :: assemble_complex_boundary_trace_residual
    public :: assemble_complex_boundary_trace_residual_jvp
    public :: assemble_complex_boundary_trace_residual_vjp
    public :: assemble_maxwell_fem_bem_boundary_matrix_3d
    public :: assemble_maxwell_fem_bem_system_3d
    public :: assemble_maxwell_fem_bem_torus_curved_system_3d
    public :: assemble_maxwell_rwg_nedelec_coupling_3d
    public :: assemble_maxwell_torus_curved_dtn_rwg_3d
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
    public :: condense_wall_response_blocks
    public :: condense_wall_response_blocks_jvp
    public :: condense_wall_response_blocks_vjp
    public :: assemble_helmholtz_single_layer_p0_3d
    public :: assemble_helmholtz_single_layer_p0_adaptive_3d
    public :: assemble_helmholtz_single_layer_p0_3d_jvp
    public :: assemble_helmholtz_single_layer_p0_3d_vjp
    public :: assemble_laplace_single_layer_p0_3d
    public :: assemble_laplace_single_layer_p0_adaptive_3d
    public :: apply_helmholtz_cfie_p0_hierarchical_3d
    public :: apply_helmholtz_single_layer_p0_hierarchical_3d
    public :: apply_laplace_single_layer_p0_hierarchical_3d
    public :: evaluate_helmholtz_cfie_p0_3d
    public :: evaluate_helmholtz_representation_triangles_3d
    public :: evaluate_helmholtz_representation_torus_curved_3d
    public :: evaluate_laplace_representation_triangles_3d
    public :: evaluate_laplace_representation_torus_curved_3d
    public :: evaluate_maxwell_torus_curved_far_field_rwg_3d
    public :: evaluate_maxwell_torus_curved_magnetic_field_rwg_3d
    public :: solve_helmholtz_cfie_p0_3d
    public :: solve_helmholtz_cfie_p0_hierarchical_3d
    public :: solve_helmholtz_dirichlet_p0_3d
    public :: solve_helmholtz_dirichlet_p0_hierarchical_3d
    public :: solve_laplace_dirichlet_p0_3d
    public :: solve_maxwell_fem_bem_torus_curved_system_3d
    public :: solve_maxwell_fem_bem_system_3d
    public :: solve_maxwell_fem_bem_linear_state
    public :: solve_maxwell_fem_bem_linear_state_jvp
    public :: solve_maxwell_fem_bem_linear_state_vjp
    public :: solve_maxwell_pec_torus_curved_regularized_cfie_rwg_multiple_3d
    public :: evaluate_torus_curved_panel
    public :: evaluate_helmholtz_representation_torus_curved_3d
    public :: evaluate_helmholtz_representation_torus_curved_3d_geometry_jvp
    public :: evaluate_helmholtz_representation_torus_curved_3d_geometry_vjp
    public :: assemble_helmholtz_torus_curved_calderon_3d
    public :: assemble_helmholtz_torus_curved_dtn_3d
    public :: assemble_helmholtz_torus_curved_dtn_3d_geometry_jvp
    public :: assemble_helmholtz_torus_curved_dtn_3d_geometry_vjp
    public :: assemble_helmholtz_fem_bem_costabel_torus_curved_3d
    public :: solve_helmholtz_bem_dtn_torus_curved_3d
    public :: solve_helmholtz_fem_bem_costabel_torus_curved_3d
    public :: solve_laplace_bem_dtn_torus_curved_3d
    public :: assemble_laplace_torus_curved_calderon_3d
    public :: assemble_laplace_torus_curved_dtn_3d
    public :: assemble_laplace_torus_curved_dtn_3d_geometry_jvp
    public :: assemble_laplace_torus_curved_dtn_3d_geometry_vjp
    public :: solve_laplace_fem_bem_costabel_torus_curved_3d

    ! Boundary type for defining domains
    type :: boundary_t
        integer :: n_points = 0
        real(dp), allocatable :: points(:,:) ! (2, n_points)
        integer, allocatable :: labels(:) ! (n_points-1) segment labels
        logical :: is_closed = .false.
    end type boundary_t

end module fortfem_boundary
