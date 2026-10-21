! Extracted from src/generated/fortfem_bem_high_order_products.f90 by
! tools/fortad/extract.py. It is the primal fortsym generated
! its products from, with the offload directives removed:
! they are the consumer's choice, not part of what the
! function means.

pure subroutine bem_hyper_product(green, tx, ty, tz, sx, sy, sz, value)
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    real(dp), intent(in) :: green, tx, ty, tz, sx, sy, sz
    real(dp), intent(out) :: value

    value = green*(sx*tx + sy*ty + sz*tz)

end subroutine bem_hyper_product
