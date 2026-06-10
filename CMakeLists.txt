cmake_minimum_required(VERSION 3.20)
project(fortfem VERSION 0.1.0 LANGUAGES Fortran)

# CMake is an additional build option; fpm (fpm.toml) remains the primary
# build tool and both build the same source tree.

option(FORTFEM_EXAMPLES "Build example programs" OFF)
option(FORTFEM_TRIANGLE_COMPAT_ONLY
    "Build only the dependency-free triangle_compat mesher library" OFF)
include(CTest)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

set(FORTFEM_MODULE_DIR ${CMAKE_CURRENT_BINARY_DIR}/include)

# The Triangle-compatible mesher has no external dependencies, so embedding
# projects can fetch fortfem and build just this target.
file(GLOB FORTFEM_TRIANGLE_COMPAT_SOURCES CONFIGURE_DEPENDS
    src/triangle_compat/*.f90)
add_library(fortfem_triangle_compat ${FORTFEM_TRIANGLE_COMPAT_SOURCES})
set_target_properties(fortfem_triangle_compat PROPERTIES
    Fortran_MODULE_DIRECTORY ${FORTFEM_MODULE_DIR}
)
target_include_directories(fortfem_triangle_compat PUBLIC
    $<BUILD_INTERFACE:${FORTFEM_MODULE_DIR}>
    $<INSTALL_INTERFACE:include/fortfem>
)

if(FORTFEM_TRIANGLE_COMPAT_ONLY)
    return()
endif()

find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)

# SuiteSparse libraries declared via link = [...] in fpm.toml. Homebrew keeps
# them in a keg-only prefix on macOS.
set(FORTFEM_SUITESPARSE_HINTS
    /opt/homebrew/opt/suite-sparse/lib
    /usr/local/opt/suite-sparse/lib
)
set(FORTFEM_SUITESPARSE_LIBRARIES)
foreach(_name umfpack amd cholmod suitesparseconfig)
    string(TOUPPER ${_name} _upper)
    find_library(${_upper}_LIBRARY ${_name} HINTS ${FORTFEM_SUITESPARSE_HINTS})
    if(NOT ${_upper}_LIBRARY)
        message(FATAL_ERROR "SuiteSparse library '${_name}' not found")
    endif()
    list(APPEND FORTFEM_SUITESPARSE_LIBRARIES ${${_upper}_LIBRARY})
endforeach()

# fortplot is a dependency of the core library (fortfem_api re-exports the
# plot interface), so it is always fetched. fpm.toml pins no revision, hence
# GIT_TAG main. For offline builds point FETCHCONTENT_SOURCE_DIR_FORTPLOT at
# a local clone.
include(FetchContent)
set(FORTPLOT_BUILD_FULL ON)
FetchContent_Declare(fortplot
    GIT_REPOSITORY https://github.com/lazy-fortran/fortplot.git
    GIT_TAG main
)
FetchContent_MakeAvailable(fortplot)

file(GLOB_RECURSE FORTFEM_SOURCES CONFIGURE_DEPENDS src/*.f90)
list(REMOVE_ITEM FORTFEM_SOURCES ${FORTFEM_TRIANGLE_COMPAT_SOURCES})
add_library(fortfem ${FORTFEM_SOURCES})
set_target_properties(fortfem PROPERTIES
    Fortran_MODULE_DIRECTORY ${FORTFEM_MODULE_DIR}
)
target_include_directories(fortfem PUBLIC
    $<BUILD_INTERFACE:${FORTFEM_MODULE_DIR}>
    $<INSTALL_INTERFACE:include/fortfem>
)
target_link_libraries(fortfem PUBLIC
    fortfem_triangle_compat
    fortplot::fortplot
    ${FORTFEM_SUITESPARSE_LIBRARIES}
    LAPACK::LAPACK
    BLAS::BLAS
)

add_executable(fortfem_app app/main.f90)
set_target_properties(fortfem_app PROPERTIES OUTPUT_NAME fortfem)
target_link_libraries(fortfem_app PRIVATE fortfem)

add_executable(triangle_compat_cli app/triangle_compat_cli.f90)
target_link_libraries(triangle_compat_cli PRIVATE fortfem)

include(GNUInstallDirs)
install(TARGETS fortfem fortfem_triangle_compat
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
)
install(DIRECTORY ${FORTFEM_MODULE_DIR}/
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/fortfem
    FILES_MATCHING PATTERN "*.mod"
)

if(BUILD_TESTING)
    # test/check.f90 and test/tc_compat_check.f90 are support modules
    # shared by the test programs.
    add_library(fortfem_check STATIC test/check.f90 test/tc_compat_check.f90)
    set_target_properties(fortfem_check PROPERTIES
        Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/test_modules
    )
    target_include_directories(fortfem_check INTERFACE
        ${CMAKE_BINARY_DIR}/test_modules
    )
    target_link_libraries(fortfem_check PRIVATE fortfem)

    # Tests write plot output to the relative path build/, as under fpm.
    set(FORTFEM_TEST_WORKDIR ${CMAKE_BINARY_DIR}/test_run)
    file(MAKE_DIRECTORY ${FORTFEM_TEST_WORKDIR}/build)

    file(GLOB FORTFEM_TEST_SOURCES CONFIGURE_DEPENDS test/*.f90)
    list(REMOVE_ITEM FORTFEM_TEST_SOURCES
        ${CMAKE_CURRENT_SOURCE_DIR}/test/check.f90
        ${CMAKE_CURRENT_SOURCE_DIR}/test/tc_compat_check.f90
        # run_all_tests shells out to fpm; it is not a ctest case.
        ${CMAKE_CURRENT_SOURCE_DIR}/test/run_all_tests.f90
    )
    foreach(_src ${FORTFEM_TEST_SOURCES})
        get_filename_component(_name ${_src} NAME_WE)
        add_executable(${_name} ${_src})
        target_link_libraries(${_name} PRIVATE fortfem fortfem_check)
        add_test(NAME ${_name}
            COMMAND ${_name}
            WORKING_DIRECTORY ${FORTFEM_TEST_WORKDIR}
        )
    endforeach()

    # The Triangle-compat tests read test/data/ via paths relative to the
    # project root, as under fpm.
    set_tests_properties(
        test_triangle_compat_golden
        test_triangle_compat_oracle
        PROPERTIES WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    )

    # These compare against the external Triangle binary and download it from
    # netlib when /tmp/triangle_bin is missing; without network or a working
    # build they warn and pass on the FortFEM-only assertions.
    set_tests_properties(
        test_triangle_comparison
        test_triangle_hole_comparison
        test_triangle_intersection_comparison
        test_triangle_quality_comparison
        test_lshape_hole_refinement_comparison
        PROPERTIES LABELS network
    )
endif()

if(FORTFEM_EXAMPLES)
    file(GLOB FORTFEM_EXAMPLE_DIRS RELATIVE ${CMAKE_CURRENT_SOURCE_DIR}/example
        ${CMAKE_CURRENT_SOURCE_DIR}/example/*)
    foreach(_dir ${FORTFEM_EXAMPLE_DIRS})
        file(GLOB _example_sources
            ${CMAKE_CURRENT_SOURCE_DIR}/example/${_dir}/*.f90)
        if(_example_sources)
            add_executable(example_${_dir} ${_example_sources})
            target_link_libraries(example_${_dir} PRIVATE fortfem)
        endif()
    endforeach()
endif()
