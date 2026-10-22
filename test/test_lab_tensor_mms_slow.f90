program test_lab_tensor_mms_slow
    use lab_tensor_mms_ladder, only: run_tensor_mms
    implicit none
    call run_tensor_mms(2, 6)
end program
