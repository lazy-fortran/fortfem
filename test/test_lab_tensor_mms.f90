program test_lab_tensor_mms
    use lab_tensor_mms_ladder, only: run_tensor_mms
    implicit none
    ! Standard gate: full 2D ladders and 3D levels 0..1, p = 2..4. The full 3D
    ! ladder (about 6 CPU s more) is test_lab_tensor_mms_slow.
    call run_tensor_mms(1, 4)
end program
