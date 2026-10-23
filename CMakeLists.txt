cmake_minimum_required(VERSION 3.22)

project(fortfem VERSION 0.1.0 LANGUAGES C Fortran)

# gfortran < 13 mis-optimizes an allocatable left-hand side of `x = matmul(...)`
# at -O1 and above: the inlined matmul does not reallocate x and writes past
# the heap block (minimal reproducer: gfortran 12.2 -O2, valgrind). Calling the
# library matmul avoids it; newer compilers are not affected.
if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU"
   AND CMAKE_Fortran_COMPILER_VERSION VERSION_LESS 13)
    add_compile_options($<$<COMPILE_LANGUAGE:Fortran>:-finline-matmul-limit=0>)
endif()

include(FetchContent)
include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

option(FORTFEM_FETCH_DEPENDENCIES
    "Fetch pinned dependencies when targets are not already available" ON)
option(FORTFEM_BUILD_TESTING "Build focused C API and sparse adapter tests" OFF)
option(FORTFEM_NATIVE_SPARSE_ONLY
    "Build with project-owned sparse solvers and no external sparse libraries" ON)
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

option(FORTFEM_ELEMENT_KERNELS_ONLY
    "Build the dependency-free batched element and CSR interface only" OFF)
option(FORTFEM_ENABLE_CUDA "Build batched CUDA element kernels" OFF)
set(FORTFEM_CUDA_ARCHITECTURES "" CACHE STRING "Explicit CUDA architecture list")
find_package(OpenMP REQUIRED COMPONENTS Fortran)
add_library(fortfem_element_kernels STATIC
    src/utils/kinds.f90
    src/generated/fortfem_element_precision_products.f90
    src/assembly/element_batch.f90
    src/assembly/element_to_csr.f90)
add_library(fortfem::element_kernels ALIAS fortfem_element_kernels)
set_target_properties(fortfem_element_kernels PROPERTIES
    Fortran_MODULE_DIRECTORY "${PROJECT_BINARY_DIR}/element_modules")
target_include_directories(fortfem_element_kernels PUBLIC
    $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/element_modules>
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>)
target_link_libraries(fortfem_element_kernels PUBLIC OpenMP::OpenMP_Fortran)
if(FORTFEM_BUILD_TESTING)
    enable_testing()
    add_executable(test_element_precision_batch
        test/check.f90 test/test_element_precision_batch.f90)
    target_link_libraries(test_element_precision_batch PRIVATE
        fortfem_element_kernels)
    add_test(NAME fortfem_element_precision_batch
        COMMAND test_element_precision_batch)
    set_tests_properties(fortfem_element_precision_batch PROPERTIES
        LABELS "element;fast" ENVIRONMENT "OMP_NUM_THREADS=4")
endif()
# The laboratory line contract is also available in the dependency-free build.
add_library(fortfem_tensor_lines STATIC
    src/generated/fortfem_lagrange_product_jets.f90
    src/discretisation/tensor_line_tables.f90
    src/discretisation/gll_line.f90
    src/discretisation/tensor_mesh.f90
    src/discretisation/tensor_contractions.f90)
set_target_properties(fortfem_tensor_lines PROPERTIES
    Fortran_MODULE_DIRECTORY "${PROJECT_BINARY_DIR}/tensor_modules")
target_include_directories(fortfem_tensor_lines PUBLIC
    $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/tensor_modules>)
target_link_libraries(fortfem_tensor_lines PUBLIC OpenMP::OpenMP_Fortran)
if(FORTFEM_BUILD_TESTING)
    add_executable(test_lab_gll_line test/check.f90 test/test_lab_gll_line.f90)
    target_link_libraries(test_lab_gll_line PRIVATE fortfem_tensor_lines)
    add_test(NAME test_lab_gll_line COMMAND test_lab_gll_line)
    add_executable(test_lab_tensor_mesh test/test_lab_tensor_mesh.f90)
    target_link_libraries(test_lab_tensor_mesh PRIVATE fortfem_tensor_lines)
    add_test(NAME test_lab_tensor_mesh COMMAND test_lab_tensor_mesh)
    set_tests_properties(test_lab_gll_line PROPERTIES
        LABELS "lab;fast" ENVIRONMENT "OMP_NUM_THREADS=4")
endif()
add_library(fortfem_tensor_operator STATIC
    src/generated/fortfem_tensor_metric_products.f90
    src/discretisation/line_mesh.f90
    src/discretisation/tensor_operator.f90)
set_target_properties(fortfem_tensor_operator PROPERTIES
    Fortran_MODULE_DIRECTORY "${PROJECT_BINARY_DIR}/tensor_modules")
target_include_directories(fortfem_tensor_operator PUBLIC
    $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/tensor_modules>)
target_link_libraries(fortfem_tensor_operator PUBLIC
    fortfem_tensor_lines fortfem_element_kernels)
if(FORTFEM_BUILD_TESTING)
    add_executable(test_lab_tensor_apply test/test_lab_tensor_apply.f90)
    target_link_libraries(test_lab_tensor_apply PRIVATE fortfem_tensor_operator)
    add_test(NAME test_lab_tensor_apply COMMAND test_lab_tensor_apply)
    set_tests_properties(test_lab_tensor_apply PROPERTIES
        LABELS "lab;fast" ENVIRONMENT "OMP_NUM_THREADS=4")
endif()
option(FORTFEM_BUILD_BENCHMARKS "Build element scaling benchmark" OFF)
if(FORTFEM_BUILD_BENCHMARKS)
    add_executable(fortfem_element_precision_scaling
        benchmark/fortfem_element_precision_scaling.f90)
    target_link_libraries(fortfem_element_precision_scaling PRIVATE
        fortfem_element_kernels)
endif()
if(FORTFEM_ENABLE_CUDA)
    if(NOT FORTFEM_CUDA_ARCHITECTURES)
        message(FATAL_ERROR "FORTFEM_CUDA_ARCHITECTURES must be explicit")
    endif()
    set(CMAKE_CUDA_ARCHITECTURES ${FORTFEM_CUDA_ARCHITECTURES})
    enable_language(CUDA)
    add_library(fortfem_element_cuda STATIC src/cuda/element_batch.cu)
    target_include_directories(fortfem_element_cuda PUBLIC include
        PRIVATE src/generated)
    set_target_properties(fortfem_element_cuda PROPERTIES
        CUDA_STANDARD 17 CUDA_STANDARD_REQUIRED ON)
    if(FORTFEM_BUILD_TESTING)
        add_executable(test_element_precision_cuda
            test/test_element_precision_cuda.cu)
        target_link_libraries(test_element_precision_cuda PRIVATE
            fortfem_element_cuda)
        add_test(NAME fortfem_element_precision_cuda
            COMMAND test_element_precision_cuda)
        set_tests_properties(fortfem_element_precision_cuda PROPERTIES
            LABELS "element;cuda" SKIP_RETURN_CODE 77)
    endif()
endif()
if(FORTFEM_ELEMENT_KERNELS_ONLY)
    return()
endif()

if(NOT TARGET fortsparse)
    if(NOT FORTFEM_FETCH_DEPENDENCIES)
        message(FATAL_ERROR
            "fortsparse target is required when dependency fetching is disabled")
    endif()
    if(FORTFEM_NATIVE_SPARSE_ONLY)
        set(FORTSPARSE_ENABLE_SUPERLU OFF CACHE BOOL "" FORCE)
        set(FORTSPARSE_ENABLE_UMFPACK_IPC OFF CACHE BOOL "" FORCE)
        set(FORTSPARSE_ENABLE_MUMPS OFF CACHE BOOL "" FORCE)
        set(FORTSPARSE_ENABLE_NATIVE_CUDA OFF CACHE BOOL "" FORCE)
    endif()
    set(FORTSPARSE_BUILD_TESTING OFF CACHE BOOL "" FORCE)
    set(FORTSPARSE_INSTALL_HELPER OFF CACHE BOOL "" FORCE)
    set(FORTSPARSE_USE_SYSTEM_DEPS ON CACHE BOOL "" FORCE)
    FetchContent_Declare(fortsparse
        GIT_REPOSITORY https://github.com/lazy-fortran/fortsparse.git
        GIT_TAG 7f16e38aaf5ebdfa56b0b8b4f5f57c0c31de212c
        GIT_SHALLOW FALSE)
    FetchContent_MakeAvailable(fortsparse)
endif()

set(FORTFEM_CAPI_SOURCES
    src/generated/fortfem_reference_basis_products.f90
    src/assembly/assembly_mixed_2d.f90
    src/assembly/assembly_nedelec_2d.f90
    src/basis/basis_edge_2d.f90
    src/basis/basis_rt_2d.f90
    src/bindings/fortfem_capi_mesh.f90
    src/bindings/fortfem_capi_solver.f90
    src/boundary/boundary_types.f90
    src/interpolation/nedelec_field_2d.f90
    src/interpolation/rt_field_2d.f90
    src/mesh/mesh_2d.f90
    src/mesh/mesh_extension_2d.f90
    src/mesh/triangle_io.f90
    src/quadrature/gauss_quadrature_2d.f90
    src/solvers/fortfem_sparse_direct.f90
    src/topology/tree_cotree_gauge.f90
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

add_library(fortfem_tensor_solver STATIC src/discretisation/tensor_solver.f90)
set_target_properties(fortfem_tensor_solver PROPERTIES
    Fortran_MODULE_DIRECTORY "${PROJECT_BINARY_DIR}/tensor_modules")
target_include_directories(fortfem_tensor_solver PUBLIC
    $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/tensor_modules>
    PRIVATE "${PROJECT_BINARY_DIR}/modules")
target_link_libraries(fortfem_tensor_solver PUBLIC
    fortfem_tensor_operator fortfem_capi fortsparse)
if(FORTFEM_BUILD_TESTING)
    add_executable(test_lab_tensor_solver test/test_lab_tensor_solver.f90)
    target_link_libraries(test_lab_tensor_solver PRIVATE fortfem_tensor_solver)
    add_test(NAME test_lab_tensor_solver COMMAND test_lab_tensor_solver)
    set_tests_properties(test_lab_tensor_solver PROPERTIES
        LABELS "lab" ENVIRONMENT "OMP_NUM_THREADS=4")
endif()

if(FORTFEM_BUILD_TESTING)
    add_executable(test_lab_tensor_mms test/test_lab_tensor_mms.f90)
    target_link_libraries(test_lab_tensor_mms PRIVATE fortfem_tensor_solver)
    add_test(NAME test_lab_tensor_mms COMMAND test_lab_tensor_mms)
    set_tests_properties(test_lab_tensor_mms PROPERTIES
        LABELS "lab;convergence" ENVIRONMENT "OMP_NUM_THREADS=4")
endif()
if(FORTFEM_BUILD_BENCHMARKS)
    add_executable(fortfem_lab_mf1 benchmark/lab/mf1/tensor_benchmark.f90)
    target_link_libraries(fortfem_lab_mf1 PRIVATE fortfem_tensor_solver)
endif()

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
        GIT_TAG 76e27115347c59e26cf80d73f578327675348f50
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

if(FORTFEM_BUILD_TESTING)
    enable_testing()
    add_library(fortfem_test_check STATIC test/check.f90)
    set_target_properties(fortfem_test_check PROPERTIES
        Fortran_MODULE_DIRECTORY "${PROJECT_BINARY_DIR}/test_modules")
    target_include_directories(fortfem_test_check PUBLIC
        "${PROJECT_BINARY_DIR}/test_modules")
    foreach(test_name IN ITEMS sparse_direct_backend sparse_constrained_derivatives
            sparse_direct_derivatives sparse_direct_complex_derivatives
            capi_solver)
        add_executable(test_${test_name} test/test_${test_name}.f90)
        target_include_directories(test_${test_name} PRIVATE
            "${PROJECT_BINARY_DIR}/modules")
        target_link_libraries(test_${test_name} PRIVATE
            fortfem_capi fortfem_test_check fortsparse)
        add_test(NAME fortfem_${test_name} COMMAND test_${test_name})
        set_tests_properties(fortfem_${test_name} PROPERTIES
            LABELS "sparse;fast" ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1")
    endforeach()
endif()
