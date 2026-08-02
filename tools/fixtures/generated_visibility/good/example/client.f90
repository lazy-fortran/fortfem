program generated_visibility_client
    use fortfem_api, only: evaluate_kernel
    implicit none
    real :: value
    call evaluate_kernel(value)
end program generated_visibility_client
