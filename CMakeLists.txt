cmake_minimum_required(VERSION 3.22)

project(fortfem VERSION 0.1.0 LANGUAGES C Fortran)

include(FetchContent)
include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

option(FORTFEM_FETCH_DEPENDENCIES
    "Fetch pinned dependencies when targets are not already available" ON)

set(CMAKE_POSITION_INDEPENDENT_CODE ON)
set(CMAKE_Fortran_STANDARD 2018)
set(CMAKE_Fortran_STANDARD_REQUIRED ON)

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
        GIT_TAG 23f3f6c62a1632dccd41ddde57957eb322a34f7d
        GIT_SHALLOW FALSE)
    FetchContent_MakeAvailable(fortsparse)
endif()

set(FORTFEM_CAPI_SOURCES
    src/basis/basis_rt_2d.f90
    src/bindings/fortfem_capi_mesh.f90
    src/bindings/fortfem_capi_solver.f90
    src/fortfem_boundary.f90
    src/interpolation/rt_field_2d.f90
    src/mesh/mesh_2d.f90
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

install(TARGETS fortfem_capi
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
