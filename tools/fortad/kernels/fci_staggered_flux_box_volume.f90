! Extracted from src/generated/fortfem_fci_support_volume_products.f90 by
! tools/fortad/extract.py. It is the primal fortsym generated
! its products from, with the offload directives removed:
! they are the consumer's choice, not part of what the
! function means.

pure subroutine generated_fci_staggered_flux_box_volume(forward_flux_expansion, &
        backward_flux_expansion, base_cell_area, toroidal_field, staggered_volume)
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    real(dp), intent(in) :: forward_flux_expansion, backward_flux_expansion, base_cell_area, &
        toroidal_field
    real(dp), intent(out) :: staggered_volume

    staggered_volume = base_cell_area*toroidal_field*(backward_flux_expansion + &
        forward_flux_expansion)

end subroutine generated_fci_staggered_flux_box_volume
