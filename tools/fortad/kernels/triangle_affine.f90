! Extracted from src/generated/fortfem_affine_map_products.f90 by
! tools/fortad/extract.py. It is the primal fortsym generated
! its products from, with the offload directives removed:
! they are the consumer's choice, not part of what the
! function means.

pure subroutine generated_triangle_affine(jacobian, relative, reference)
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    real(dp), intent(in) :: jacobian(2,2), relative(2)
    real(dp), intent(out) :: reference(2)
    real(dp) :: t1

    t1 = jacobian(1,1)*jacobian(2,2) - jacobian(1,2)*jacobian(2,1)
    reference(1) = -jacobian(1,2)*relative(2)/t1 + jacobian(2,2)*relative(1)/t1
    reference(2) = jacobian(1,1)*relative(2)/t1 - jacobian(2,1)*relative(1)/t1

end subroutine generated_triangle_affine
