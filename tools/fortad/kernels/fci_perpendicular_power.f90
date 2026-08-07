! Extracted from src/generated/fortfem_fci_power_flux_products.f90 by
! tools/fortad/extract.py. It is the primal fortsym generated
! its products from, with the offload directives removed:
! they are the consumer's choice, not part of what the
! function means.

pure subroutine generated_fci_perpendicular_power(field, action, canonical_volume, &
        perpendicular_power)
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    real(dp), intent(in) :: field, action, canonical_volume
    real(dp), intent(out) :: perpendicular_power

    perpendicular_power = action*canonical_volume*field

end subroutine generated_fci_perpendicular_power
