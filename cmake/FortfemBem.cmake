# Native closure for caller-owned parametric P2/P1 Laplace BEM surfaces.
add_library(fortfem_bem STATIC
    src/utils/kinds.f90
    src/generated/fortfem_reference_basis_products.f90
    src/generated/fortfem_parametric_bem_products.f90
    src/generated/fortfem_bem_high_order_products.f90
    src/quadrature/triangle_duffy_quadrature.f90
    src/estimators/equilibrated_flux_topology.f90
    src/bem/parametric_laplace_bem_3d.f90
    src/bem/bem_circulant.f90
    src/bem/parametric_bem_high_order.f90)
add_library(fortfem::bem ALIAS fortfem_bem)
set_target_properties(fortfem_bem PROPERTIES
    Fortran_MODULE_DIRECTORY "${PROJECT_BINARY_DIR}/bem_modules")
target_include_directories(fortfem_bem PUBLIC
    "$<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/bem_modules>")
target_link_libraries(fortfem_bem PUBLIC fortnum OpenMP::OpenMP_Fortran)

if(FORTFEM_BUILD_TESTING)
    add_library(fortfem_bem_test_geometry STATIC
        src/generated/fortfem_sphere_curved_panel_products.f90
        src/generated/fortfem_torus_curved_panel_products.f90
        src/mesh/sphere_curved_panel.f90
        src/mesh/torus_curved_panel.f90
        src/mesh/sphere_surface_mesh.f90
        src/mesh/torus_surface_mesh.f90
        test/parametric_bem_test_geometry.f90)
    set_target_properties(fortfem_bem_test_geometry PROPERTIES
        Fortran_MODULE_DIRECTORY "${PROJECT_BINARY_DIR}/bem_test_modules")
    target_include_directories(fortfem_bem_test_geometry PUBLIC
        "${PROJECT_BINARY_DIR}/bem_test_modules")
    target_link_libraries(fortfem_bem_test_geometry PUBLIC fortfem_bem)

    foreach(test_name IN ITEMS parametric_bem_edge_policy
            parametric_bem_high_order parametric_bem_estimator
            parametric_bem_pair_accuracy bem_shared_edge_affine
            bem_shared_edge_charts parametric_bem_generated)
        add_executable(test_${test_name} test/test_${test_name}.f90)
        target_link_libraries(test_${test_name} PRIVATE
            fortfem_bem_test_geometry fortfem_test_check)
        add_test(NAME test_${test_name} COMMAND test_${test_name})
        set_tests_properties(test_${test_name} PROPERTIES
            LABELS "bem;parametric"
            ENVIRONMENT "OMP_NUM_THREADS=4;OPENBLAS_NUM_THREADS=1")
    endforeach()

    add_executable(test_parametric_bem_circulant
        test/parametric_bem_circulant_test_cases.f90
        test/test_parametric_bem_circulant.f90)
    target_link_libraries(test_parametric_bem_circulant PRIVATE
        fortfem_bem_test_geometry fortfem_test_check)
    add_test(NAME test_parametric_bem_circulant COMMAND test_parametric_bem_circulant)
    set_tests_properties(test_parametric_bem_circulant PROPERTIES
        LABELS "bem;parametric"
        ENVIRONMENT "OMP_NUM_THREADS=4;OPENBLAS_NUM_THREADS=1")
endif()
