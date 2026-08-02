---
title: Example Source Code
---

# Example Source Code and Plots

Every executable example shipped with FortFEM is listed here.
Pages contain the complete source, usage command, and generated plots.

- [simple_poisson](simple_poisson.html) - This example demonstrates solving the classic Poisson equation using FortFEM's FEniCS-inspired API.
- [minimal_mesh_example](minimal_mesh_example.html) - A simple demonstration of FortFEM's mesh generation capabilities.
- [plot_mesh](plot_mesh.html) - This example demonstrates FortFEM's mesh visualization capabilities, showing how to create and plot finite element meshes at different refinement levels.
- [plotting](plotting.html) - This comprehensive example showcases FortFEM's complete plotting and visualization capabilities using multiple colormap options and plot types.
- [mixed_poisson](mixed_poisson.html) - This example is the first vector finite-element problem after the scalar
- [anisotropic_tensor_diffusion](anisotropic_tensor_diffusion.html) - This example solves a manufactured two-dimensional diffusion problem with a
- [xfem_interface_solution](xfem_interface_solution.html) - This physical-first gallery fixture constructs a manufactured scalar and
- [regularized_surface_current_layer](regularized_surface_current_layer.html) - This physical-first slab example resolves a tangential surface current
- [sheet_current_parity](sheet_current_parity.html) - This physical-first slab gallery fixture compares the same tangential sheet
- [curl_curl](curl_curl.html) - This example solves a manufactured three-dimensional curl-curl problem with
- [cgl_pressure_tensor](cgl_pressure_tensor.html) - This manufactured profile evaluates the gyrotropic CGL pressure tensor
- [field_aligned_flux](field_aligned_flux.html) - This manufactured profile evaluates the generated pointwise constitutive
- [fci_parallel_diffusion](fci_parallel_diffusion.html) - This small fixture exercises FortFEM's field-coordinate-independent support
- [fci_sol_gallery](fci_sol_gallery.html) - This neutral FCI/SOL fixture traces three helical field lines on a toroidal
- [fci_quadrilateral_geometry](fci_quadrilateral_geometry.html) - This small geometry fixture computes positive, fixed-topology areas for three
- [fci_polygon_geometry](fci_polygon_geometry.html) - This fixture computes fixed-topology straight, quadratic, cubic, quartic, quintic, and sextic
- [mixed_acoustic_wave](mixed_acoustic_wave.html) - This fixture is a small physical wave problem for the common mixed
- [mixed_wave_3d_structure](mixed_wave_3d_structure.html) - This neutral gallery fixture advances three manufactured Cartesian oscillator
- [mixed_elasticity_wave](mixed_elasticity_wave.html) - This example is a manufactured mixed elastic-wave modal problem. It uses the
- [mixed_wave_wall](mixed_wave_wall.html) - This small structure-preserving example couples a first-order mixed wave port
- [nonlinear_resistive_mhd_gallery](nonlinear_resistive_mhd_gallery.html) - This bounded manufactured fixture is the physical-first gallery seed for
- [tetra_h1_poisson](tetra_h1_poisson.html) - This example solves
- [tetra_mixed_poisson](tetra_mixed_poisson.html) - This example extends the symbolic Darcy form from triangles to a unit cube
- [tetra_nedelec_p_convergence](tetra_nedelec_p_convergence.html) - This example interpolates the manufactured curl-free vector field.  The first
- [circular_dtn_modes](circular_dtn_modes.html) - Boundary data on the unit circle:
- [helmholtz_open_boundary_comparison](helmholtz_open_boundary_comparison.html) - This example compares three truncations of the same one-dimensional outgoing
- [acoustic_fem_dtn](acoustic_fem_dtn.html) - This example couples complex, time-harmonic P1 elasticity to an outgoing
- [maxwell_open_boundary_comparison](maxwell_open_boundary_comparison.html) - This example evaluates the biperiodic planar Maxwell capacity operator on
- [fortfem_mesh_benchmark](fortfem_mesh_benchmark.html) - This example benchmarks FortFEM's mesh generation performance against FreeFEM.
- [solver_benchmark](solver_benchmark.html) - This example solves a manufactured Poisson system, shows the physical
- [biro_tree_cotree_benchmark](biro_tree_cotree_benchmark.html) - This small manufactured case reproduces the structural part of the direct
- [team13_neutral_benchmark](team13_neutral_benchmark.html) - This is a small, solution-first foundation fixture shaped like the TEAM 13
- [interoperability_benchmarks](interoperability_benchmarks.html) - This example reads the neutral CSV records produced by the isolated FEniCSx,
- [linear_perturbation_response](linear_perturbation_response.html) - This closure-neutral fixture composes the seven public linear-perturbation
- [eulerian_island_gallery](eulerian_island_gallery.html) - This small gallery samples the analytic slab-island flux
- [singular_layer_matching_gallery](singular_layer_matching_gallery.html) - This bounded MHD-08/MHD-09 foundation gallery samples analytical complex inner
- [direct_force_campaign_gallery](direct_force_campaign_gallery.html) - This small, closure-neutral campaign exercises FortFEM's public direct
- [laplace_bem_circle_spectrum](laplace_bem_circle_spectrum.html) - Executable FortFEM laplace_bem_circle_spectrum.f90 example.
- [laplace_exterior_bem_circle](laplace_exterior_bem_circle.html) - This example solves the unbounded exterior Dirichlet problem on the unit
- [helmholtz_bem_circle_spectrum](helmholtz_bem_circle_spectrum.html) - Executable FortFEM helmholtz_bem_circle_spectrum.f90 example.
- [helmholtz_cfie_circle](helmholtz_cfie_circle.html) - Executable FortFEM helmholtz_cfie_circle.f90 example.
- [laplace_symmetric_transmission](laplace_symmetric_transmission.html) - Manufactured transmission pair:
- [curved_acoustic_ntd](curved_acoustic_ntd.html) - This example maps the normal displacement of an arbitrary closed polygonal
- [bem_sphere_3d](bem_sphere_3d.html) - This example solves the unit-sphere Laplace Dirichlet problem with the P0
- [adaptive_bem_prolate](adaptive_bem_prolate.html) - This example solves the exterior electrostatic Dirichlet problem with constant
- [adaptive_helmholtz_bem_sphere](adaptive_helmholtz_bem_sphere.html) - This example applies the solve-estimate-mark-refine loop to the complex P0
- [magnetic_curvilinear_metrics](magnetic_curvilinear_metrics.html) - This example evaluates the metric-to-coefficient transform used by the
- [toroidal_analytical](toroidal_analytical.html) - This example evaluates the real \(n=2,m=1\) toroidal harmonic
- [nested_surface_solution](nested_surface_solution.html) - This physical-first gallery fixture evaluates FortFEM's neutral nested-surface
- [iga_shape_sensitivity](iga_shape_sensitivity.html) - This example differentiates a complete scalar isogeometric state problem,
- [maxwell_mesh_adjoint](maxwell_mesh_adjoint.html) - This example benchmarks analytical forward and reverse products for
- [iga_polar_feec](iga_polar_feec.html) - This example constructs the Type-1 polar spline de Rham sequence
- [iga_jorek_flux](iga_jorek_flux.html) - This advanced example evolves the ideal poloidal magnetic-flux subflow
- [maxwell_torus_curved_scattering](maxwell_torus_curved_scattering.html) - This example solves time-harmonic perfect-electric-conductor scattering on an
- [maxwell_torus_fem_bem_solution](maxwell_torus_fem_bem_solution.html) - This fixture solves the neutral curl--curl FEM--BEM system on a small exact

[← Back to examples](../index.html)
