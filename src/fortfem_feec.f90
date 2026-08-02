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
    public :: diagnose_tree_cotree_iga_invariance
    public :: evaluate_tetra_nedelec_first_kind
    public :: evaluate_tetra_nedelec_first_kind_jvp
    public :: evaluate_tetra_nedelec_first_kind_vjp
    public :: initialize_tetra_nedelec_first_kind
    public :: interpolate_reference_tetra_nedelec
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
    public :: validate_tree_cotree_gauge

end module fortfem_feec
