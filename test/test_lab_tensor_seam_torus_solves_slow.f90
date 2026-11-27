program test_lab_tensor_seam_torus_solves_slow
    use lab_tensor_seam_cases, only: run_tensor_seam_slow_part
    implicit none
    call run_tensor_seam_slow_part(8)
end program
