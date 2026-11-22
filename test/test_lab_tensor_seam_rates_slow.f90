program test_lab_tensor_seam_rates_slow
    use lab_tensor_seam_cases, only: run_tensor_seam_slow_part
    implicit none
    call run_tensor_seam_slow_part(3)
end program
