! Extracted from src/generated/fortfem_surface_shape_objective_products.f90 by
! tools/fortad/extract.py. It is the primal fortsym generated
! its products from, with the offload directives removed:
! they are the consumer's choice, not part of what the
! function means.

pure subroutine generated_surface_shape_objective_contribution(candidate, target, weight, &
        contribution)
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    real(dp), intent(in) :: candidate, target, weight
    real(dp), intent(out) :: contribution

    contribution = weight*(candidate - target)**2*5.0000000000000000E-001_dp

end subroutine generated_surface_shape_objective_contribution
