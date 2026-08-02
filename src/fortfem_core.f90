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
    use fortfem_toroidal_coordinates, only: &
        cartesian_to_toroidal, cartesian_to_toroidal_jvp, &
        cartesian_to_toroidal_vjp, toroidal_point_to_cartesian, &
        toroidal_point_to_cartesian_jvp, toroidal_point_to_cartesian_vjp
    use fortfem_tetra_affine_map, only: invert_tetra_affine_map
    use fortfem_sphere_surface_mesh, only: generate_sphere_surface_mesh
    use fortfem_solid_torus_tetra_mesh, only: generate_solid_torus_tetra_mesh
    use fortfem_torus_surface_mesh, only: generate_torus_surface_mesh
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
    public :: generate_torus_surface_mesh
    public :: initialize_cell_complex
    public :: invert_tetra_affine_map
    public :: toroidal_point_to_cartesian
    public :: toroidal_point_to_cartesian_jvp
    public :: toroidal_point_to_cartesian_vjp
    public :: validate_cell_complex

end module fortfem_core
