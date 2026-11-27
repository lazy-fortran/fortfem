program test_lab_tensor_seam_reproduction_p4_slow
    use lab_tensor_seam_cases, only: run_tensor_seam_slow_part
    implicit none
    call run_tensor_seam_slow_part(4)
end program
