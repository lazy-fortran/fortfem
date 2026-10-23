! Extracted from src/generated/fortfem_harmonic_field_products.f90 by
! tools/fortad/extract.py. It is the primal fortsym generated
! its products from, with the offload directives removed:
! they are the consumer's choice, not part of what the
! function means.

pure subroutine field_inner_product(weight, a1, a2, a3, b1, b2, b3, value)
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    real(dp), intent(in) :: weight, a1, a2, a3, b1, b2, b3
    real(dp), intent(out) :: value

    value = weight*(a1*b1 + a2*b2 + a3*b3)

end subroutine field_inner_product
