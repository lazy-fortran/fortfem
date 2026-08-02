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
    use fortfem_tetra_duffy_quadrature, only: tetra_duffy_quadrature
    use fortfem_tetra_nedelec_global_dof_map, only: build_tetra_nedelec_dof_map
    use fortfem_tetra_piola_maps, only: map_tetra_nedelec_covariant
    use fortfem_assembly_bspline_polar_2d, only: &
        assemble_bspline_polar_h1_operator_csc, &
        assemble_bspline_polar_hcurl_operator_csc, &
        assemble_bspline_polar_l2_mass_csc, restrict_bspline_polar_operator_csc
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
        assemble_tensor_diffusion_matrix_3d, &
        assemble_tensor_diffusion_matrix_3d_jvp, &
        assemble_tensor_diffusion_matrix_3d_vjp
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
    public :: tetra_duffy_quadrature
    public :: assemble_bspline_polar_h1_operator_csc
    public :: assemble_bspline_polar_hcurl_operator_csc
    public :: assemble_bspline_polar_l2_mass_csc
    public :: restrict_bspline_polar_operator_csc
    public :: assemble_tensor_diffusion_matrix_3d
    public :: assemble_tensor_diffusion_matrix_3d_jvp
    public :: assemble_tensor_diffusion_matrix_3d_vjp
    public :: validate_tree_cotree_gauge
    public :: beltrami_parity_t
    public :: compare_beltrami_two_region_residual
    public :: validate_beltrami_parity
    public :: validate_beltrami_resonance
    public :: beltrami_shell_parity_t
    public :: compare_beltrami_shell_residual
    public :: validate_beltrami_shell_parity

end module fortfem_feec
