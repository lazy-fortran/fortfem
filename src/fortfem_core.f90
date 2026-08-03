module fortfem_core
    !! Small canonical facade for foundational topology and geometry.
    !!
    !! This module deliberately re-exports existing implementations instead of
    !! wrapping or duplicating them.  The broader `fortfem_api` remains
    !! available for compatibility; new clients can depend on this low-risk
    !! core surface without importing application-specific operators.
    use fortfem_cell_complex, only: &
        cell_complex_betti_numbers, cell_complex_euler_characteristic, &
        cell_complex_t, initialize_cell_complex, validate_cell_complex
    use fortfem_boundary_region_graph, only: &
        boundary_region_graph_t, initialize_boundary_region_graph, &
        validate_boundary_region_graph, boundary_region_graph_incidence, &
        boundary_region_graph_components, boundary_region_graph_cycle_basis, &
        boundary_region_graph_interface_samples, &
        boundary_region_graph_interface_metadata
    use fortfem_physical_trace_ownership, only: &
        physical_trace_ownership_t, initialize_physical_trace_ownership, &
        validate_physical_trace_ownership, physical_trace_ownership_maps, &
        physical_trace_ownership_dimension, physical_trace_ownership_point_count, &
        physical_trace_ownership_rank, compare_physical_trace_coordinates
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
    use fortfem_nested_surface_geometry, only: evaluate_nested_surface_geometry
    use fortfem_sphere_surface_mesh, only: generate_sphere_surface_mesh
    use fortfem_solid_torus_tetra_mesh, only: generate_solid_torus_tetra_mesh
    use fortfem_structured_tetra_box_mesh, only: &
        generate_structured_tetra_box_mesh
    use fortfem_torus_surface_mesh, only: generate_torus_surface_mesh
    use fortfem_torus_curved_panel, only: evaluate_torus_curved_panel
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
    use fortfem_api_types, only: mesh_t, function_space_t, &
        vector_function_space_t, dirichlet_bc_t
    use fortfem_api_mesh, only: circle_boundary, rectangle_boundary, &
        line_segment, arc_segment, l_shape_boundary, unit_square_mesh, &
        rectangle_mesh, unit_disk_mesh, mesh_from_boundary, structured_quad_mesh
    implicit none
    private

    public :: cartesian_to_toroidal
    public :: cartesian_to_toroidal_jvp
    public :: cartesian_to_toroidal_vjp
    public :: cell_complex_betti_numbers
    public :: cell_complex_euler_characteristic
    public :: cell_complex_t
    public :: generate_solid_torus_tetra_mesh
    public :: generate_sphere_surface_mesh
    public :: generate_structured_tetra_box_mesh
    public :: generate_torus_surface_mesh
    public :: evaluate_torus_curved_panel
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
    public :: evaluate_axis_regular_radial_basis
    public :: evaluate_axis_regular_radial_basis_jvp
    public :: evaluate_axis_regular_radial_basis_vjp
    public :: evaluate_nested_surface_geometry
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
    public :: boundary_region_graph_t
    public :: initialize_boundary_region_graph
    public :: validate_boundary_region_graph
    public :: boundary_region_graph_incidence
    public :: boundary_region_graph_components
    public :: boundary_region_graph_cycle_basis
    public :: boundary_region_graph_interface_samples
    public :: boundary_region_graph_interface_metadata
    public :: physical_trace_ownership_t
    public :: initialize_physical_trace_ownership
    public :: validate_physical_trace_ownership
    public :: physical_trace_ownership_maps
    public :: physical_trace_ownership_dimension
    public :: physical_trace_ownership_point_count
    public :: physical_trace_ownership_rank
    public :: compare_physical_trace_coordinates

end module fortfem_core
