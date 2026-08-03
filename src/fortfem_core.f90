module fortfem_core
    !! Small canonical facade for foundational topology and geometry.
    !!
    !! This module deliberately re-exports existing implementations instead of
    !! wrapping or duplicating them.  The broader `fortfem_api` remains
    !! available for compatibility; new clients can depend on this low-risk
    !! core surface without importing application-specific operators.
    use fortfem_cell_complex, only: &
        cell_complex_betti_numbers, cell_complex_cocycle_basis, &
        cell_complex_cohomology_cocycle_basis, cell_complex_cycle_basis, &
        cell_complex_homology_cycle_basis, cell_complex_harmonic_one_forms, &
        cell_complex_euler_characteristic, cell_complex_t, &
        initialize_cell_complex, quotient_cell_complex, validate_cell_complex
    use fortfem_physical_surface_geometry, only: &
        sample_physical_surface_geometry, &
        sample_physical_surface_geometry_jvp, &
        sample_physical_surface_geometry_vjp
    use fortfem_boundary_region_graph, only: &
        boundary_region_graph_t, initialize_boundary_region_graph, &
        validate_boundary_region_graph, boundary_region_graph_incidence, &
        boundary_region_graph_components, boundary_region_graph_cycle_basis, &
        boundary_region_graph_interface_samples, &
        boundary_region_graph_interface_metadata
    use fortfem_region_interface_graph, only: &
        region_interface_graph_t, initialize_region_interface_graph, &
        validate_region_interface_graph, region_interface_graph_incidence, &
        region_interface_graph_components, region_interface_graph_cycle_basis
    use fortfem_physical_trace_ownership, only: &
        physical_trace_ownership_t, initialize_physical_trace_ownership, &
        validate_physical_trace_ownership, physical_trace_ownership_maps, &
        physical_trace_ownership_dimension, physical_trace_ownership_point_count, &
        physical_trace_ownership_rank, compare_physical_trace_coordinates
    use fortfem_physical_trace_reconciliation, only: &
        physical_trace_reconciliation_t, &
        initialize_physical_trace_reconciliation, &
        validate_physical_trace_reconciliation, &
        physical_trace_reconciliation_maps, &
        reconcile_physical_trace_values, &
        reconcile_physical_trace_values_jvp, &
        reconcile_physical_trace_values_vjp
    use fortfem_physical_trace_owner_selection, only: &
        physical_trace_owner_selection_t, initialize_physical_trace_owner_selection, &
        validate_physical_trace_owner_selection, physical_trace_owner_selection_maps, &
        gather_physical_trace_values, gather_physical_trace_values_jvp, &
        gather_physical_trace_values_vjp
    use fortfem_mpi_trace_exchange, only: &
        mpi_trace_exchange_schedule_t, initialize_mpi_trace_exchange_schedule, &
        validate_mpi_trace_exchange_schedule, mpi_trace_exchange_schedule_maps, &
        pack_mpi_trace_exchange, pack_mpi_trace_exchange_jvp, pack_mpi_trace_exchange_vjp, &
        unpack_mpi_trace_exchange, unpack_mpi_trace_exchange_jvp, unpack_mpi_trace_exchange_vjp
    use fortfem_patch_graph_trace_contraction, only: &
        assemble_patch_graph_trace_contraction, &
        assemble_patch_graph_trace_contraction_jvp, &
        assemble_patch_graph_trace_contraction_vjp
    use fortfem_toroidal_coordinates, only: &
        cartesian_to_toroidal, cartesian_to_toroidal_jvp, &
        cartesian_to_toroidal_vjp, toroidal_point_to_cartesian, &
        toroidal_point_to_cartesian_jvp, toroidal_point_to_cartesian_vjp, &
        toroidal_vector_to_cartesian, toroidal_vector_to_cartesian_jvp, &
        toroidal_vector_to_cartesian_vjp
    use fortfem_tetra_affine_map, only: invert_tetra_affine_map, &
        invert_tetra_affine_map_jvp, invert_tetra_affine_map_vjp
    use fortfem_triangle_affine_map, only: &
        invert_triangle_affine_map, invert_triangle_affine_map_jvp, &
        invert_triangle_affine_map_vjp
    use fortfem_axis_regular_radial_basis, only: &
        evaluate_axis_regular_radial_basis, &
        evaluate_axis_regular_radial_basis_jvp, &
        evaluate_axis_regular_radial_basis_vjp
    use fortfem_nested_surface_geometry, only: &
        evaluate_nested_surface_geometry, &
        evaluate_nested_surface_geometry_jvp, &
        evaluate_nested_surface_geometry_vjp, &
        evaluate_nested_surface_geometry_coordinate_jvp, &
        evaluate_nested_surface_geometry_coordinate_vjp
    use fortfem_sphere_surface_mesh, only: generate_sphere_surface_mesh
    use fortfem_sphere_curved_panel, only: &
        evaluate_sphere_curved_panel, evaluate_sphere_curved_panel_jvp, &
        evaluate_sphere_curved_panel_vjp, invert_sphere_curved_panel
    use fortfem_solid_torus_tetra_mesh, only: generate_solid_torus_tetra_mesh
    use fortfem_structured_tetra_box_mesh, only: &
        generate_structured_tetra_box_mesh
    use fortfem_torus_surface_mesh, only: generate_torus_surface_mesh
    use fortfem_torus_curved_panel, only: &
        evaluate_torus_curved_panel, evaluate_torus_curved_panel_jvp, &
        evaluate_torus_curved_panel_vjp
    use fortfem_surface_triangle_geometry_3d, only: &
        evaluate_surface_triangle_geometry_3d, &
        evaluate_surface_triangle_geometry_3d_jvp, &
        evaluate_surface_triangle_geometry_3d_vjp
    use fortfem_surface_triangle_measures_3d, only: &
        assemble_surface_triangle_measures_3d, &
        assemble_surface_triangle_measures_3d_jvp, &
        assemble_surface_triangle_measures_3d_vjp
    use fortfem_barycentric_surface_refinement, only: &
        barycentric_refine_torus_surface_mesh, &
        barycentric_refine_torus_surface_mesh_jvp, &
        barycentric_refine_torus_surface_mesh_vjp
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
    use fortfem_geometry_mortar_trace_coupling, only: &
        assemble_geometry_mortar_trace_coupling, &
        assemble_geometry_mortar_trace_coupling_jvp, &
        assemble_geometry_mortar_trace_coupling_vjp
    use fortfem_api_types, only: mesh_t, function_space_t, &
        vector_function_space_t, dirichlet_bc_t
    use fortfem_api_mesh, only: circle_boundary, rectangle_boundary, &
        line_segment, arc_segment, l_shape_boundary, unit_square_mesh, &
        rectangle_mesh, unit_disk_mesh, mesh_from_boundary, mesh_from_domain, &
        structured_quad_mesh, refine_uniform, refine_adaptive, &
        compute_gradient_indicators, find_triangle_edges
    implicit none
    private

    public :: cartesian_to_toroidal
    public :: cartesian_to_toroidal_jvp
    public :: cartesian_to_toroidal_vjp
    public :: cell_complex_betti_numbers
    public :: cell_complex_cocycle_basis
    public :: cell_complex_cohomology_cocycle_basis
    public :: cell_complex_cycle_basis
    public :: cell_complex_euler_characteristic
    public :: cell_complex_homology_cycle_basis
    public :: cell_complex_harmonic_one_forms
    public :: cell_complex_t
    public :: generate_solid_torus_tetra_mesh
    public :: generate_sphere_surface_mesh
    public :: evaluate_sphere_curved_panel
    public :: evaluate_sphere_curved_panel_jvp
    public :: evaluate_sphere_curved_panel_vjp
    public :: invert_sphere_curved_panel
    public :: generate_structured_tetra_box_mesh
    public :: generate_torus_surface_mesh
    public :: barycentric_refine_torus_surface_mesh
    public :: barycentric_refine_torus_surface_mesh_jvp
    public :: barycentric_refine_torus_surface_mesh_vjp
    public :: evaluate_torus_curved_panel
    public :: evaluate_torus_curved_panel_jvp
    public :: evaluate_torus_curved_panel_vjp
    public :: evaluate_surface_triangle_geometry_3d
    public :: evaluate_surface_triangle_geometry_3d_jvp
    public :: evaluate_surface_triangle_geometry_3d_vjp
    public :: assemble_surface_triangle_measures_3d
    public :: assemble_surface_triangle_measures_3d_jvp
    public :: assemble_surface_triangle_measures_3d_vjp
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
    public :: mesh_t
    public :: function_space_t
    public :: vector_function_space_t
    public :: dirichlet_bc_t
    public :: unit_square_mesh
    public :: rectangle_mesh
    public :: unit_disk_mesh
    public :: structured_quad_mesh
    public :: circle_boundary
    public :: rectangle_boundary
    public :: line_segment
    public :: arc_segment
    public :: l_shape_boundary
    public :: mesh_from_boundary
    public :: mesh_from_domain
    public :: refine_uniform
    public :: refine_adaptive
    public :: compute_gradient_indicators
    public :: find_triangle_edges
    public :: sample_physical_surface_geometry
    public :: sample_physical_surface_geometry_jvp
    public :: sample_physical_surface_geometry_vjp
    public :: evaluate_axis_regular_radial_basis
    public :: evaluate_axis_regular_radial_basis_jvp
    public :: evaluate_axis_regular_radial_basis_vjp
    public :: evaluate_nested_surface_geometry
    public :: evaluate_nested_surface_geometry_jvp
    public :: evaluate_nested_surface_geometry_vjp
    public :: evaluate_nested_surface_geometry_coordinate_jvp
    public :: evaluate_nested_surface_geometry_coordinate_vjp
    public :: initialize_cell_complex
    public :: invert_tetra_affine_map
    public :: invert_tetra_affine_map_jvp
    public :: invert_tetra_affine_map_vjp
    public :: invert_triangle_affine_map
    public :: invert_triangle_affine_map_jvp
    public :: invert_triangle_affine_map_vjp
    public :: toroidal_point_to_cartesian
    public :: toroidal_point_to_cartesian_jvp
    public :: toroidal_point_to_cartesian_vjp
    public :: toroidal_vector_to_cartesian
    public :: toroidal_vector_to_cartesian_jvp
    public :: toroidal_vector_to_cartesian_vjp
    public :: validate_cell_complex
    public :: quotient_cell_complex
    public :: boundary_region_graph_t
    public :: initialize_boundary_region_graph
    public :: validate_boundary_region_graph
    public :: boundary_region_graph_incidence
    public :: boundary_region_graph_components
    public :: boundary_region_graph_cycle_basis
    public :: boundary_region_graph_interface_samples
    public :: boundary_region_graph_interface_metadata
    public :: region_interface_graph_t
    public :: initialize_region_interface_graph
    public :: validate_region_interface_graph
    public :: region_interface_graph_incidence
    public :: region_interface_graph_components
    public :: region_interface_graph_cycle_basis
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
    public :: mpi_trace_exchange_schedule_t
    public :: initialize_mpi_trace_exchange_schedule
    public :: validate_mpi_trace_exchange_schedule
    public :: mpi_trace_exchange_schedule_maps
    public :: pack_mpi_trace_exchange
    public :: pack_mpi_trace_exchange_jvp
    public :: pack_mpi_trace_exchange_vjp
    public :: unpack_mpi_trace_exchange
    public :: unpack_mpi_trace_exchange_jvp
    public :: unpack_mpi_trace_exchange_vjp
    public :: assemble_patch_graph_trace_contraction
    public :: assemble_patch_graph_trace_contraction_jvp
    public :: assemble_patch_graph_trace_contraction_vjp
    public :: gather_physical_trace_values
    public :: gather_physical_trace_values_jvp
    public :: gather_physical_trace_values_vjp
    public :: assemble_geometry_mortar_trace_coupling
    public :: assemble_geometry_mortar_trace_coupling_jvp
    public :: assemble_geometry_mortar_trace_coupling_vjp

end module fortfem_core
