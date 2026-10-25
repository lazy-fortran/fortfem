# Neutral callback coarse cache; the caller supplies P and the frozen action.
find_package(LAPACK REQUIRED)
add_library(fortfem_callback_galerkin STATIC src/solvers/callback_galerkin.f90)
add_library(fortfem::callback_galerkin ALIAS fortfem_callback_galerkin)
set_target_properties(fortfem_callback_galerkin PROPERTIES
    Fortran_MODULE_DIRECTORY "${PROJECT_BINARY_DIR}/galerkin_modules")
target_include_directories(fortfem_callback_galerkin PUBLIC
    $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/galerkin_modules>)
target_link_libraries(fortfem_callback_galerkin PUBLIC
    LAPACK::LAPACK OpenMP::OpenMP_Fortran)
if(FORTFEM_BUILD_TESTING)
    add_executable(test_callback_galerkin test/test_callback_galerkin.f90)
    target_link_libraries(test_callback_galerkin PRIVATE
        fortfem_callback_galerkin fortfem_test_check)
    add_test(NAME fortfem_callback_galerkin COMMAND test_callback_galerkin)
    set_tests_properties(fortfem_callback_galerkin PROPERTIES
        LABELS "solver;fast"
        ENVIRONMENT "OMP_NUM_THREADS=4;OPENBLAS_NUM_THREADS=1")
endif()
