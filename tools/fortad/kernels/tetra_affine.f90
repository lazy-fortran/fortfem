! Extracted from src/generated/fortfem_affine_map_products.f90 by
! tools/fortad/extract.py. It is the primal fortsym generated
! its products from, with the offload directives removed:
! they are the consumer's choice, not part of what the
! function means.

pure subroutine generated_tetra_affine(jacobian, relative, reference)
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    real(dp), intent(in) :: jacobian(3,3), relative(3)
    real(dp), intent(out) :: reference(3)

    reference(1) = relative(1)*(jacobian(2,2)*jacobian(3,3) - jacobian(2,3)*jacobian(3,2))/ &
        (jacobian(1,1)*jacobian(2,2)*jacobian(3,3) - jacobian(1,1)*jacobian(2,3)*jacobian(3,2) - jacobian(1, &
        2)*jacobian(2,1)*jacobian(3,3) + jacobian(1,2)*jacobian(2,3)*jacobian(3,1) + jacobian(1,3)* &
        jacobian(2,1)*jacobian(3,2) - jacobian(1,3)*jacobian(2,2)*jacobian(3,1)) - relative(2)*(jacobian(1, &
        2)*jacobian(3,3) - jacobian(1,3)*jacobian(3,2))/(jacobian(1,1)*jacobian(2,2)*jacobian(3,3) - &
        jacobian(1,1)*jacobian(2,3)*jacobian(3,2) - jacobian(1,2)*jacobian(2,1)*jacobian(3,3) + jacobian(1, &
        2)*jacobian(2,3)*jacobian(3,1) + jacobian(1,3)*jacobian(2,1)*jacobian(3,2) - jacobian(1,3)* &
        jacobian(2,2)*jacobian(3,1)) + relative(3)*(jacobian(1,2)*jacobian(2,3) - jacobian(1,3)*jacobian(2, &
        2))/(jacobian(1,1)*jacobian(2,2)*jacobian(3,3) - jacobian(1,1)*jacobian(2,3)*jacobian(3,2) - &
        jacobian(1,2)*jacobian(2,1)*jacobian(3,3) + jacobian(1,2)*jacobian(2,3)*jacobian(3,1) + jacobian(1, &
        3)*jacobian(2,1)*jacobian(3,2) - jacobian(1,3)*jacobian(2,2)*jacobian(3,1))
    reference(2) = -relative(1)*(jacobian(2,1)*jacobian(3,3) - jacobian(2,3)*jacobian(3,1))/ &
        (jacobian(1,1)*jacobian(2,2)*jacobian(3,3) - jacobian(1,1)*jacobian(2,3)*jacobian(3,2) - jacobian(1, &
        2)*jacobian(2,1)*jacobian(3,3) + jacobian(1,2)*jacobian(2,3)*jacobian(3,1) + jacobian(1,3)* &
        jacobian(2,1)*jacobian(3,2) - jacobian(1,3)*jacobian(2,2)*jacobian(3,1)) + relative(2)*(jacobian(1, &
        1)*jacobian(3,3) - jacobian(1,3)*jacobian(3,1))/(jacobian(1,1)*jacobian(2,2)*jacobian(3,3) - &
        jacobian(1,1)*jacobian(2,3)*jacobian(3,2) - jacobian(1,2)*jacobian(2,1)*jacobian(3,3) + jacobian(1, &
        2)*jacobian(2,3)*jacobian(3,1) + jacobian(1,3)*jacobian(2,1)*jacobian(3,2) - jacobian(1,3)* &
        jacobian(2,2)*jacobian(3,1)) - relative(3)*(jacobian(1,1)*jacobian(2,3) - jacobian(1,3)*jacobian(2, &
        1))/(jacobian(1,1)*jacobian(2,2)*jacobian(3,3) - jacobian(1,1)*jacobian(2,3)*jacobian(3,2) - &
        jacobian(1,2)*jacobian(2,1)*jacobian(3,3) + jacobian(1,2)*jacobian(2,3)*jacobian(3,1) + jacobian(1, &
        3)*jacobian(2,1)*jacobian(3,2) - jacobian(1,3)*jacobian(2,2)*jacobian(3,1))
    reference(3) = relative(1)*(jacobian(2,1)*jacobian(3,2) - jacobian(2,2)*jacobian(3,1))/ &
        (jacobian(1,1)*jacobian(2,2)*jacobian(3,3) - jacobian(1,1)*jacobian(2,3)*jacobian(3,2) - jacobian(1, &
        2)*jacobian(2,1)*jacobian(3,3) + jacobian(1,2)*jacobian(2,3)*jacobian(3,1) + jacobian(1,3)* &
        jacobian(2,1)*jacobian(3,2) - jacobian(1,3)*jacobian(2,2)*jacobian(3,1)) - relative(2)*(jacobian(1, &
        1)*jacobian(3,2) - jacobian(1,2)*jacobian(3,1))/(jacobian(1,1)*jacobian(2,2)*jacobian(3,3) - &
        jacobian(1,1)*jacobian(2,3)*jacobian(3,2) - jacobian(1,2)*jacobian(2,1)*jacobian(3,3) + jacobian(1, &
        2)*jacobian(2,3)*jacobian(3,1) + jacobian(1,3)*jacobian(2,1)*jacobian(3,2) - jacobian(1,3)* &
        jacobian(2,2)*jacobian(3,1)) + relative(3)*(jacobian(1,1)*jacobian(2,2) - jacobian(1,2)*jacobian(2, &
        1))/(jacobian(1,1)*jacobian(2,2)*jacobian(3,3) - jacobian(1,1)*jacobian(2,3)*jacobian(3,2) - &
        jacobian(1,2)*jacobian(2,1)*jacobian(3,3) + jacobian(1,2)*jacobian(2,3)*jacobian(3,1) + jacobian(1, &
        3)*jacobian(2,1)*jacobian(3,2) - jacobian(1,3)*jacobian(2,2)*jacobian(3,1))

end subroutine generated_tetra_affine
