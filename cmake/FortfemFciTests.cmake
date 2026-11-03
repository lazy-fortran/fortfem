add_library(fortfem_fci_transport_test_support STATIC
    "${PROJECT_SOURCE_DIR}/test/fci_transport_reference.f90")
set_target_properties(fortfem_fci_transport_test_support PROPERTIES
    Fortran_MODULE_DIRECTORY "${PROJECT_BINARY_DIR}/fci_transport_test_modules")
target_include_directories(fortfem_fci_transport_test_support PUBLIC
    "${PROJECT_BINARY_DIR}/fci_transport_test_modules")
target_link_libraries(fortfem_fci_transport_test_support PUBLIC fortfem::api)

file(GLOB fci_test_sources CONFIGURE_DEPENDS
    "${PROJECT_SOURCE_DIR}/test/test_fci_*.f90")
foreach(test_source IN LISTS fci_test_sources)
    get_filename_component(test_name "${test_source}" NAME_WE)
    add_executable(${test_name} "${test_source}")
    if(test_name STREQUAL "test_fci_transport_reference")
        target_link_libraries(${test_name} PRIVATE
            fortfem_fci_transport_test_support fortfem_test_check)
    else()
        target_link_libraries(${test_name} PRIVATE fortfem::api fortfem_test_check)
    endif()
    if(test_name STREQUAL "test_fci_transport_pollution")
        target_link_libraries(${test_name} PRIVATE fortfem_tensor_lines)
    endif()
    if(test_name STREQUAL "test_fci_transport_pollution")
        add_test(NAME ${test_name} COMMAND ${test_name} diagnostics
            "${PROJECT_BINARY_DIR}/pollution_checks.csv")
    else()
        add_test(NAME ${test_name} COMMAND ${test_name})
    endif()
    set_tests_properties(${test_name} PROPERTIES
        LABELS "fci" PROCESSORS 1
        ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1")
    if(test_name MATCHES "^test_fci_transport_(reference|pollution)$")
        set_property(TEST ${test_name} APPEND PROPERTY LABELS fast)
    endif()
endforeach()

add_test(NAME test_fci_transport_pollution_figure_report
    COMMAND test_fci_transport_pollution paper
        "${PROJECT_BINARY_DIR}/pollution_published.csv")
set_tests_properties(test_fci_transport_pollution_figure_report PROPERTIES
    LABELS "fci;benchmark" PROCESSORS 1
    ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1")

add_test(NAME test_fci_transport_pollution_isotropic_gate
    COMMAND test_fci_transport_pollution isotropic
        "${PROJECT_BINARY_DIR}/pollution_isotropic.csv")
set_tests_properties(test_fci_transport_pollution_isotropic_gate PROPERTIES
    LABELS "fci;acceptance;fast" PROCESSORS 1
    ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1")

add_test(NAME test_fci_transport_pollution_parallel_gate
    COMMAND test_fci_transport_pollution parallel
        "${PROJECT_BINARY_DIR}/pollution_parallel.csv")
set_tests_properties(test_fci_transport_pollution_parallel_gate PROPERTIES
    LABELS "fci;acceptance;fast" PROCESSORS 1
    ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1")

foreach(degree RANGE 1 4)
    add_test(NAME test_fci_transport_pollution_scaling_p${degree}
        COMMAND test_fci_transport_pollution scaling ${degree}
            "${PROJECT_BINARY_DIR}/pollution_scaling_p${degree}.csv")
    set_tests_properties(test_fci_transport_pollution_scaling_p${degree} PROPERTIES
        LABELS "fci;benchmark" PROCESSORS 1
        ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1")
endforeach()

add_test(NAME test_fci_transport_reference_threads
    COMMAND test_fci_transport_reference)
set_tests_properties(test_fci_transport_reference_threads PROPERTIES
    LABELS "fci;fast" PROCESSORS 8
    ENVIRONMENT "OMP_NUM_THREADS=8;OPENBLAS_NUM_THREADS=1")
foreach(method IN ITEMS reference direct)
    foreach(threads IN ITEMS 1 8)
        set(name test_fci_transport_${method}_cost_${threads})
        add_test(NAME ${name} COMMAND test_fci_transport_reference ${method})
        set_tests_properties(${name} PROPERTIES
            LABELS "fci;benchmark" PROCESSORS ${threads}
            ENVIRONMENT "OMP_NUM_THREADS=${threads};OPENBLAS_NUM_THREADS=1")
    endforeach()
endforeach()

add_test(NAME test_fci_transport_pollution_sparse_gate
    COMMAND test_fci_transport_pollution sparse)
set_tests_properties(test_fci_transport_pollution_sparse_gate PROPERTIES
    LABELS "fci;acceptance;fast" PROCESSORS 1
    ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1")

# Keep one small measured row per candidate in the fast selection.
foreach(method IN ITEMS h1 symmetric fci)
    set(degrees 1)
    if(method STREQUAL "h1")
        set(degrees 1 2 3 4)
    endif()
    foreach(degree IN LISTS degrees)
        set(name test_fci_transport_pollution_cost_row_${method}_p${degree})
        add_test(NAME ${name} COMMAND ${CMAKE_COMMAND}
            "-Dexecutable=$<TARGET_FILE:test_fci_transport_pollution>"
            "-Dmethod=${method}" "-Ddegree=${degree}"
            "-Doutput=${PROJECT_BINARY_DIR}/pollution_cost_row_${method}_p${degree}.csv"
            -P "${PROJECT_SOURCE_DIR}/cmake/CheckFciPollutionCostRow.cmake")
        set_tests_properties(${name} PROPERTIES
            LABELS "fci;acceptance;fast" PROCESSORS 1
            RESOURCE_LOCK anisotropic_pollution_timing
            ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1")
    endforeach()
endforeach()

foreach(method IN ITEMS h1 symmetric fci)
    set(degrees 1)
    if(method STREQUAL "h1")
        set(degrees 1 2 3 4)
    endif()
    foreach(degree IN LISTS degrees)
        foreach(power IN ITEMS 6 8 10)
            set(name test_fci_transport_pollution_cost_${method}_p${degree}_r${power})
            add_test(NAME ${name} COMMAND test_fci_transport_pollution cost
                sovinec ${method} ${degree} ${power}
                "${PROJECT_BINARY_DIR}/pollution_cost_${method}_p${degree}_r${power}.csv")
            set_tests_properties(${name} PROPERTIES
                LABELS "fci;benchmark" PROCESSORS 1
                RESOURCE_LOCK anisotropic_pollution_timing
                ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1")
        endforeach()
    endforeach()
endforeach()

foreach(method IN ITEMS symmetric fci)
    foreach(power IN ITEMS 6 8 10)
        set(name test_fci_transport_pollution_circle_${method}_r${power})
        add_test(NAME ${name} COMMAND test_fci_transport_pollution cost
            circle ${method} 1 ${power}
            "${PROJECT_BINARY_DIR}/pollution_circle_${method}_r${power}.csv")
        set_tests_properties(${name} PROPERTIES
            LABELS "fci;benchmark" PROCESSORS 1
            RESOURCE_LOCK anisotropic_pollution_timing
            ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1")
    endforeach()
endforeach()

add_test(NAME test_fci_transport_pollution_candidate_parallel
    COMMAND test_fci_transport_pollution parallel_candidates
        "${PROJECT_BINARY_DIR}/pollution_candidate_parallel.csv")
add_test(NAME test_fci_transport_pollution_manufactured_source
    COMMAND test_fci_transport_pollution manufactured_source)
set_tests_properties(test_fci_transport_pollution_candidate_parallel
    test_fci_transport_pollution_manufactured_source PROPERTIES
    LABELS "fci;acceptance;fast" PROCESSORS 1
    ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1")

foreach(problem IN ITEMS sovinec circle manufactured)
    foreach(method IN ITEMS h1 fci symmetric)
        if(NOT problem STREQUAL "manufactured" AND method STREQUAL "h1")
            continue()
        endif()
        set(degrees 1)
        if(method STREQUAL "h1")
            set(degrees 1 2 3 4)
        endif()
        foreach(degree IN LISTS degrees)
            foreach(power IN ITEMS 6 8 10)
                set(name test_fci_transport_pollution_refine_${problem}_${method}_p${degree}_r${power})
                add_test(NAME ${name} COMMAND test_fci_transport_pollution refine
                    ${problem} ${method} ${degree} ${power}
                    "${PROJECT_BINARY_DIR}/pollution_refine_${problem}_${method}_p${degree}_r${power}.csv")
                set_tests_properties(${name} PROPERTIES
                    LABELS "fci;benchmark" PROCESSORS 1
                    RESOURCE_LOCK anisotropic_pollution_timing
                    ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1")
            endforeach()
        endforeach()
    endforeach()
endforeach()

foreach(power IN ITEMS 6 8 10)
    set(name test_fci_transport_pollution_cap_manufactured_h1_p1_r${power})
    add_test(NAME ${name} COMMAND test_fci_transport_pollution refine_cap
        manufactured h1 1 ${power}
        "${PROJECT_BINARY_DIR}/pollution_cap_manufactured_h1_p1_r${power}.csv")
    set_tests_properties(${name} PROPERTIES LABELS "fci;benchmark" PROCESSORS 1
        RESOURCE_LOCK anisotropic_pollution_timing
        ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1")
endforeach()
set(name test_fci_transport_pollution_cap_manufactured_h1_p2_r10)
add_test(NAME ${name} COMMAND test_fci_transport_pollution refine_cap
    manufactured h1 2 10
    "${PROJECT_BINARY_DIR}/pollution_cap_manufactured_h1_p2_r10.csv")
set_tests_properties(${name} PROPERTIES LABELS "fci;benchmark" PROCESSORS 1
    RESOURCE_LOCK anisotropic_pollution_timing
    ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1")

add_test(NAME test_fci_transport_pollution_circle_interpolation
    COMMAND test_fci_transport_pollution circle_interpolation)
set_tests_properties(test_fci_transport_pollution_circle_interpolation PROPERTIES
    LABELS "fci;acceptance;fast" PROCESSORS 1
    ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1")
