! Extracted from src/generated/fortfem_exterior_harmonic_products.f90 by
! tools/fortad/extract.py. It is the primal fortsym generated
! its products from, with the offload directives removed:
! they are the consumer's choice, not part of what the
! function means.

pure subroutine current_segment(r1x, r1y, r1z, r2x, r2y, r2z, ax, ay, az, bx, by, bz)
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    real(dp), intent(in) :: r1x, r1y, r1z, r2x, r2y, r2z
    real(dp), intent(out) :: ax, ay, az, bx, by, bz
    real(dp) :: t1, t2, t3, t4, t5, t6, t7, t8, t9

    t1 = sqrt(r1x**2 + r1y**2 + r1z**2)
    t2 = sqrt(r2x**2 + r2y**2 + r2z**2)
    t3 = r1x - r2x
    t4 = r1y - r2y
    t5 = r1z - r2z
    t6 = sqrt(t3**2 + t4**2 + t5**2)
    t7 = log((t1 + t2 + t6)/(t1 + t2 - t6))
    t8 = t1*t2*(r1x*r2x + r1y*r2y + r1z*r2z + t1*t2)
    t9 = t1 + t2
    ax = t7*t3/t6
    ay = t7*t4/t6
    az = t7*t5/t6
    bx = t9*(r1y*r2z - r1z*r2y)/t8
    by = t9*(-r1x*r2z + r1z*r2x)/t8
    bz = t9*(r1x*r2y - r1y*r2x)/t8

end subroutine current_segment
