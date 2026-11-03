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
    if(test_name MATCHES "^test_fci_transport_(reference|nested)$")
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
    LABELS "fci;acceptance" PROCESSORS 1
    ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1")

add_test(NAME test_fci_transport_pollution_parallel_gate
    COMMAND test_fci_transport_pollution parallel
        "${PROJECT_BINARY_DIR}/pollution_parallel.csv")
set_tests_properties(test_fci_transport_pollution_parallel_gate PROPERTIES
    LABELS "fci;acceptance" PROCESSORS 1
    ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1")

foreach(degree RANGE 1 4)
    add_test(NAME test_fci_transport_pollution_scaling_p${degree}
        COMMAND test_fci_transport_pollution scaling ${degree}
            "${PROJECT_BINARY_DIR}/pollution_scaling_p${degree}.csv")
    set_tests_properties(test_fci_transport_pollution_scaling_p${degree} PROPERTIES
        LABELS "fci;benchmark" PROCESSORS 1
        ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1")
endforeach()

add_test(NAME test_fci_transport_nested_accuracy_gate
    COMMAND test_fci_transport_nested gate)
set_tests_properties(test_fci_transport_nested_accuracy_gate PROPERTIES
    LABELS "fci;acceptance" PROCESSORS 1
    ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1")

add_test(NAME test_fci_transport_reference_threads
    COMMAND test_fci_transport_reference)
set_tests_properties(test_fci_transport_reference_threads PROPERTIES
    LABELS "fci" PROCESSORS 8
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
    LABELS "fci;acceptance" PROCESSORS 1
    ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1")

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
