cmake_minimum_required(VERSION 3.22)

project(fortfem VERSION 0.1.0 LANGUAGES C Fortran)

include(FetchContent)
include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

option(FORTFEM_FETCH_DEPENDENCIES
    "Fetch pinned dependencies when targets are not already available" ON)
option(FORTFEM_ENABLE_ENZYME
    "Enable optional Flang/Enzyme element-kernel tournaments" OFF)
option(FORTFEM_ENZYME_REQUIRED
    "Fail configuration when the optional Enzyme toolchain is unavailable" OFF)
option(FORTFEM_TRIANGLE_COMPAT_ONLY
    "Build only the dependency-free triangle_compat mesher library" OFF)

set(CMAKE_POSITION_INDEPENDENT_CODE ON)
set(CMAKE_Fortran_STANDARD 2018)
set(CMAKE_Fortran_STANDARD_REQUIRED ON)

file(GLOB FORTFEM_TRIANGLE_COMPAT_SOURCES CONFIGURE_DEPENDS
    src/triangle_compat/*.f90)
add_library(fortfem_triangle_compat ${FORTFEM_TRIANGLE_COMPAT_SOURCES})
add_library(fortfem::triangle_compat ALIAS fortfem_triangle_compat)
set_target_properties(fortfem_triangle_compat PROPERTIES
    EXPORT_NAME triangle_compat
    Fortran_MODULE_DIRECTORY "${PROJECT_BINARY_DIR}/triangle_compat_modules")
target_include_directories(fortfem_triangle_compat PUBLIC
    $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/triangle_compat_modules>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/fortfem>)

if(FORTFEM_TRIANGLE_COMPAT_ONLY)
    return()
endif()

if(NOT TARGET fortsparse)
    if(NOT FORTFEM_FETCH_DEPENDENCIES)
        message(FATAL_ERROR
            "fortsparse target is required when dependency fetching is disabled")
    endif()
    set(FORTSPARSE_BUILD_TESTING OFF CACHE BOOL "" FORCE)
    set(FORTSPARSE_INSTALL_HELPER OFF CACHE BOOL "" FORCE)
    set(FORTSPARSE_USE_SYSTEM_DEPS ON CACHE BOOL "" FORCE)
    FetchContent_Declare(fortsparse
        GIT_REPOSITORY https://github.com/lazy-fortran/fortsparse.git
        GIT_TAG 3c19c978b914edde02ff58704feec0830008c162
        GIT_SHALLOW FALSE)
    FetchContent_MakeAvailable(fortsparse)
endif()

set(FORTFEM_CAPI_SOURCES
    src/assembly/assembly_mixed_2d.f90
    src/assembly/assembly_nedelec_2d.f90
    src/basis/basis_edge_2d.f90
    src/basis/basis_rt_2d.f90
    src/bindings/fortfem_capi_mesh.f90
    src/bindings/fortfem_capi_solver.f90
    src/fortfem_boundary.f90
    src/interpolation/nedelec_field_2d.f90
    src/interpolation/rt_field_2d.f90
    src/mesh/mesh_2d.f90
    src/mesh/mesh_extension_2d.f90
    src/mesh/triangle_io.f90
    src/quadrature/gauss_quadrature_2d.f90
    src/solvers/fortfem_sparse_direct.f90
    src/triangulation/bowyer_watson.f90
    src/triangulation/constrained_delaunay.f90
    src/triangulation/delaunay_types.f90
    src/triangulation/geometric_predicates.f90
    src/triangulation/mesh_refinement.f90
    src/triangulation/point_location.f90
    src/triangulation/robust_predicates.f90
    src/triangulation/triangulation_fortran.f90
    src/triangulation/triangulator.f90
    src/utils/kinds.f90)

add_library(fortfem_capi SHARED ${FORTFEM_CAPI_SOURCES})
add_library(fortfem::capi ALIAS fortfem_capi)

set_target_properties(fortfem_capi PROPERTIES
    EXPORT_NAME capi
    OUTPUT_NAME fortfem
    Fortran_MODULE_DIRECTORY "${PROJECT_BINARY_DIR}/modules")
target_include_directories(fortfem_capi
    PUBLIC
        $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
    PRIVATE "${PROJECT_BINARY_DIR}/modules")
target_link_libraries(fortfem_capi PRIVATE fortsparse)

install(TARGETS fortfem_capi fortfem_triangle_compat
    EXPORT fortfemTargets
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
install(FILES include/fortfem.h DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

configure_package_config_file(
    cmake/fortfemConfig.cmake.in
    "${PROJECT_BINARY_DIR}/fortfemConfig.cmake"
    INSTALL_DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/fortfem")
write_basic_package_version_file(
    "${PROJECT_BINARY_DIR}/fortfemConfigVersion.cmake"
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY SameMajorVersion)
install(FILES
    "${PROJECT_BINARY_DIR}/fortfemConfig.cmake"
    "${PROJECT_BINARY_DIR}/fortfemConfigVersion.cmake"
    DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/fortfem")
install(EXPORT fortfemTargets
    FILE fortfemTargets.cmake
    NAMESPACE fortfem::
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/fortfem)

if(FORTFEM_ENABLE_ENZYME)
    # Reuse FortNum's tested Flang -> LLVM IR -> Enzyme -> Flang pipeline.
    # Fetching the pinned source does not build FortNum; only its generic CMake
    # fixture support is included.
    FetchContent_Declare(fortnum_enzyme_support
        GIT_REPOSITORY https://github.com/lazy-fortran/fortnum.git
        GIT_TAG 609d7b0a095e006128c08502c1a30242e2889344
        GIT_SHALLOW FALSE
        SOURCE_SUBDIR cmake/enzyme-support)
    FetchContent_MakeAvailable(fortnum_enzyme_support)
    set(FORTNUM_ENABLE_ENZYME ON)
    set(FORTNUM_ENZYME_REQUIRED ${FORTFEM_ENZYME_REQUIRED})
    include("${fortnum_enzyme_support_SOURCE_DIR}/cmake/FortnumEnzyme.cmake")
    include(
        "${fortnum_enzyme_support_SOURCE_DIR}/cmake/FortnumAddEnzymeTest.cmake")
    enable_testing()
    fortnum_add_enzyme_test(
        NAME fortfem_enzyme_h1_geometry_products
        SOURCES
            src/generated/fortfem_bspline_h1_geometry_products.f90
            cmake/enzyme/h1_geometry_products.f90
        MODE forward
        OPT_LEVEL O2)
endif()
