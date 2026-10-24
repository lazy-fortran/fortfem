# Reuse the C API mesh/topology/solver modules and the tensor-line product jets.
# Keep one implementation of each shared source in the native target graph.
add_library(fortfem_meridian STATIC
    src/assembly/assembly_meridian_feec_n.f90
    src/assembly/assembly_meridian_fourier_apply.f90
    src/assembly/assembly_meridian_fourier_curl.f90
    src/assembly/assembly_meridian_fourier_equilibration.f90
    src/assembly/assembly_meridian_fourier_scalar.f90
    src/bem/axisymmetric_gs_boundary_operators_2d.f90
    src/bem/axisymmetric_gs_kernel.f90
    src/bem/meridian_magnetic_coupling.f90
    src/elements/polynomial_candidate_jets.f90
    src/elements/triangle_feec_operators.f90
    src/elements/triangle_lagrange_arbitrary_order.f90
    src/elements/triangle_nedelec_arbitrary_order.f90
    src/elements/triangle_piola_maps.f90
    src/elements/triangle_rt_arbitrary_order.f90
    src/estimators/equilibrated_flux_2d.f90
    src/estimators/equilibrated_flux_topology.f90
    src/estimators/equilibrated_hcurl_meridian_n.f90
    src/estimators/verified_bernstein_majorant.f90
    src/generated/fortfem_equilibrated_flux_products.f90
    src/generated/fortfem_mapped_oscillation_products.f90
    src/generated/fortfem_meridian_fourier_equilibration_products.f90
    src/generated/fortfem_meridian_hcurl_equilibration_products.f90
    src/generated/fortfem_triangle_level_geometry.f90
    src/generated/fortfem_triangle_piola_products.f90
    src/generated/fortfem_verified_bernstein_products.f90
    src/geometry/meridian_curved_mesh.f90
    src/quadrature/triangle_duffy_quadrature.f90
    src/quadrature/triangle_level_quadrature.f90)
set_target_properties(fortfem_meridian PROPERTIES
    Fortran_MODULE_DIRECTORY "${PROJECT_BINARY_DIR}/meridian_modules")
target_include_directories(fortfem_meridian PUBLIC
    "${PROJECT_BINARY_DIR}/meridian_modules"
    "${PROJECT_BINARY_DIR}/modules")
target_link_libraries(fortfem_meridian PUBLIC
    fortfem_capi fortfem_tensor_lines fortnum fortsparse OpenMP::OpenMP_Fortran)

add_library(fortfem_meridian_test_support STATIC
    test/equilibrated_test_mesh.f90
    test/meridian_test_fields.f90
    test/meridian_test_meshes.f90)
set_target_properties(fortfem_meridian_test_support PROPERTIES
    Fortran_MODULE_DIRECTORY "${PROJECT_BINARY_DIR}/meridian_test_modules")
target_include_directories(fortfem_meridian_test_support PUBLIC
    "${PROJECT_BINARY_DIR}/meridian_test_modules")
target_link_libraries(fortfem_meridian_test_support PUBLIC
    fortfem_meridian fortfem_test_check)

foreach(test_name IN ITEMS
        test_equilibrated_hcurl_meridian_n
        test_equilibrated_hcurl_meridian_sheet
        test_meridian_cylinder_chart
        test_meridian_edge_remainder
        test_meridian_feec_cohomology
        test_meridian_feec_exact_sequence
        test_meridian_feec_gs_n0
        test_meridian_feec_manufactured
        test_meridian_feec_threads
        test_meridian_fourier_apply
        test_meridian_fourier_curl
        test_meridian_fourier_equilibration
        test_meridian_fourier_scalar
        test_meridian_hermitian_forms
        test_meridian_magnetic_coupling
        test_triangle_level_conditioning
        test_triangle_level_context
        test_triangle_level_quadrature
        test_verified_bernstein_majorant
        test_verified_mapped_majorant)
    add_executable(${test_name} test/${test_name}.f90)
    set_target_properties(${test_name} PROPERTIES
        Fortran_MODULE_DIRECTORY
            "${PROJECT_BINARY_DIR}/meridian_test_modules/${test_name}")
    target_link_libraries(${test_name} PRIVATE fortfem_meridian_test_support)
    add_test(NAME ${test_name} COMMAND ${test_name})
    set_tests_properties(${test_name} PROPERTIES
        LABELS "meridian" PROCESSORS 8
        ENVIRONMENT "OMP_NUM_THREADS=8;OPENBLAS_NUM_THREADS=1")
endforeach()
