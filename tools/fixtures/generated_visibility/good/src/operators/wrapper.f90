module fortfem_kernel_wrapper
    use fortfem_generated_kernel, only: generated_kernel
    implicit none
    private
    public :: evaluate_kernel
contains
    subroutine evaluate_kernel(value)
        real, intent(out) :: value
        call generated_kernel(value)
    end subroutine evaluate_kernel
end module fortfem_kernel_wrapper
