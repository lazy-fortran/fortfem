program test_lab_tensor_seam_slow
    use lab_tensor_seam_cases, only: run_tensor_seam
    implicit none
    call run_tensor_seam(.true.)
end program
