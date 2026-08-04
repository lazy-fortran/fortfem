! Extracted from src/generated/fortfem_regularized_surface_current_products.f90 by
! tools/fortad/extract.py. It is the primal fortsym generated
! its products from, with the offload directives removed:
! they are the consumer's choice, not part of what the
! function means.

pure subroutine generated_regularized_surface_current(signed_distance, sheet_current, thickness, &
        inverse_sqrt_pi, volume_current)
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    real(dp), intent(in) :: signed_distance, sheet_current, thickness, inverse_sqrt_pi
    real(dp), intent(out) :: volume_current

    volume_current = inverse_sqrt_pi*sheet_current*exp(-signed_distance**2/thickness**2)/thickness

end subroutine generated_regularized_surface_current
