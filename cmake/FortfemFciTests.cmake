file(GLOB fci_test_sources CONFIGURE_DEPENDS
    "${PROJECT_SOURCE_DIR}/test/test_fci_*.f90")
foreach(test_source IN LISTS fci_test_sources)
    get_filename_component(test_name "${test_source}" NAME_WE)
    add_executable(${test_name} "${test_source}")
    target_link_libraries(${test_name} PRIVATE fortfem::api fortfem_test_check)
    add_test(NAME ${test_name} COMMAND ${test_name})
    set_tests_properties(${test_name} PROPERTIES
        LABELS "fci" PROCESSORS 1
        ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1")
endforeach()

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
