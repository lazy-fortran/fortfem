cmake_minimum_required(VERSION 3.22)
project(fortfem_triangle_level_oracles LANGUAGES Fortran)
# Standalone source consumer: the public FortFEM CMake target exposes the C API,
# while Fortran clients such as KIN compile their selected source graph.
get_filename_component(FORTFEM_DIR "${CMAKE_CURRENT_LIST_DIR}/../../.." ABSOLUTE)
set(FORTNUM_DIR "" CACHE PATH "Existing pinned FortNum source checkout")
if(NOT EXISTS "${FORTNUM_DIR}/src/roots/fortnum_roots.f90")
    message(FATAL_ERROR "Set FORTNUM_DIR to an existing pinned FortNum checkout")
endif()
set(CMAKE_Fortran_STANDARD 2018)
set(CMAKE_Fortran_STANDARD_REQUIRED ON)
enable_testing()
add_library(triangle_level_sources STATIC
    "${FORTNUM_DIR}/src/fortnum_kinds.f90"
    "${FORTNUM_DIR}/src/fortnum_status.f90"
    "${FORTNUM_DIR}/src/roots/fortnum_roots.f90"
    "${FORTNUM_DIR}/src/quadrature/fortnum_quadrature.f90"
    "${FORTFEM_DIR}/src/utils/kinds.f90"
    "${FORTFEM_DIR}/src/generated/fortfem_reference_basis_products.f90"
    "${FORTFEM_DIR}/src/generated/fortfem_triangle_level_geometry.f90"
    "${FORTFEM_DIR}/src/quadrature/triangle_level_quadrature.f90"
    "${FORTFEM_DIR}/test/check.f90")
set_target_properties(triangle_level_sources PROPERTIES
    Fortran_MODULE_DIRECTORY "${PROJECT_BINARY_DIR}/modules")
target_include_directories(triangle_level_sources PUBLIC "${PROJECT_BINARY_DIR}/modules")
if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    target_compile_options(triangle_level_sources PRIVATE
        "$<$<CONFIG:Debug>:-fcheck=all;-ffpe-trap=invalid,zero,overflow>")
endif()
foreach(oracle triangle_level_quadrature triangle_level_conditioning)
    add_executable(test_${oracle} "${FORTFEM_DIR}/test/test_${oracle}.f90")
    target_link_libraries(test_${oracle} PRIVATE triangle_level_sources)
    add_test(NAME ${oracle} COMMAND test_${oracle})
endforeach()
