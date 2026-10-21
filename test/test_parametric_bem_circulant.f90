program test_parametric_bem_circulant
    use parametric_bem_circulant_test_cases, only: run_circulant_case
    implicit none
    call run_circulant_case(4, .false.)
end program
