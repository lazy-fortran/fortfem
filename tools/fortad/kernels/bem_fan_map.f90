! Extracted from src/generated/fortfem_bem_high_order_products.f90 by
! tools/fortad/extract.py. It is the primal fortsym generated
! its products from, with the offload directives removed:
! they are the consumer's choice, not part of what the
! function means.

pure subroutine bem_fan_map(apex1, apex2, e1u, e1v, e2u, e2v, radius, angle, determinant, uv1, uv2, &
        weight)
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    real(dp), intent(in) :: apex1, apex2, e1u, e1v, e2u, e2v, radius, angle, determinant
    real(dp), intent(out) :: uv1, uv2, weight
    real(dp) :: t1, t2

    t1 = radius**2
    t2 = -angle + 1
    uv1 = apex1 + t1*(angle*e2u + e1u*t2)
    uv2 = apex2 + t1*(angle*e2v + e1v*t2)
    weight = determinant*radius**3*2

end subroutine bem_fan_map
