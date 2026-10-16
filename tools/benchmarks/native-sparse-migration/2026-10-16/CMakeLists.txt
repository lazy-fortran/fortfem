cmake_minimum_required(VERSION 3.22)
project(fortfem_native_comparison LANGUAGES C Fortran)
set(FORTFEM_BUILD_TESTING OFF CACHE BOOL "" FORCE)
add_subdirectory("${FORTFEM_SOURCE}" fortfem)
enable_testing()
add_library(comparison_check STATIC "${DRIVER_SOURCE}/test/check.f90")
set_target_properties(comparison_check PROPERTIES
    Fortran_MODULE_DIRECTORY "${CMAKE_BINARY_DIR}/comparison_modules")
target_include_directories(comparison_check PUBLIC
    "${CMAKE_BINARY_DIR}/comparison_modules")
foreach(case IN ITEMS sparse_direct_backend sparse_direct_derivatives
        sparse_direct_complex_derivatives sparse_constrained_derivatives capi_solver)
    add_executable(test_${case} "${DRIVER_SOURCE}/test/test_${case}.f90")
    target_include_directories(test_${case} PRIVATE
        "${CMAKE_BINARY_DIR}/fortfem/modules")
    target_link_libraries(test_${case} PRIVATE fortfem_capi comparison_check fortsparse)
    add_test(NAME ${case} COMMAND test_${case})
    set_tests_properties(${case} PROPERTIES
        ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1")
endforeach()
