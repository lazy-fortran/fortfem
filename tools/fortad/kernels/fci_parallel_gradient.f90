! Extracted from src/generated/fortfem_fci_parallel_products.f90 by
! tools/fortad/extract.py. It is the primal fortsym generated
! its products from, with the offload directives removed:
! they are the consumer's choice, not part of what the
! function means.

pure subroutine generated_fci_parallel_gradient(forward_value, upper_field, backward_value, &
        lower_field, line_length, gradient)
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    real(dp), intent(in) :: forward_value, upper_field, backward_value, lower_field, line_length
    real(dp), intent(out) :: gradient

    gradient = (-backward_value*lower_field + forward_value*upper_field)/line_length

end subroutine generated_fci_parallel_gradient
