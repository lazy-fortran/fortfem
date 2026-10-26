! Extracted from src/generated/fortfem_bem_high_order_products.f90 by
! tools/fortad/extract.py. It is the primal fortsym generated
! its products from, with the offload directives removed:
! they are the consumer's choice, not part of what the
! function means.

pure subroutine bem_trace_product(kernel, jx, jy, test, trial, value)
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    real(dp), intent(in) :: kernel, jx, jy, test, trial
    real(dp), intent(out) :: value

    value = jx*jy*kernel*test*trial

end subroutine bem_trace_product
