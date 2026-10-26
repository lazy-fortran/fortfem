! Extracted from src/generated/fortfem_exterior_harmonic_products.f90 by
! tools/fortad/extract.py. It is the primal fortsym generated
! its products from, with the offload directives removed:
! they are the consumer's choice, not part of what the
! function means.

pure subroutine exterior_gram_density(r1x, r1y, r1z, r2x, r2y, r2z, nx, ny, nz, value)
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    real(dp), intent(in) :: r1x, r1y, r1z, r2x, r2y, r2z, nx, ny, nz
    real(dp), intent(out) :: value

    value = -(nx*(r1y*r2z - r1z*r2y) + ny*(-r1x*r2z + r1z*r2x) + nz*(r1x*r2y - r1y*r2x))

end subroutine exterior_gram_density
