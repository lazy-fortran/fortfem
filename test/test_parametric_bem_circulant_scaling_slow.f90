program test_parametric_bem_circulant_scaling_slow
    use parametric_bem_circulant_test_cases, only: run_circulant_case
    implicit none
    call run_circulant_case(12, .true.)
end program
