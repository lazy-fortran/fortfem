! Extracted from src/generated/fortfem_fci_power_flux_products.f90 by
! tools/fortad/extract.py. It is the primal fortsym generated
! its products from, with the offload directives removed:
! they are the consumer's choice, not part of what the
! function means.

pure subroutine generated_fci_parallel_flux_power(gradient, coefficient, staggered_volume, &
        parallel_flux, parallel_power)
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    real(dp), intent(in) :: gradient, coefficient, staggered_volume
    real(dp), intent(out) :: parallel_flux, parallel_power

    parallel_flux = -coefficient*gradient
    parallel_power = -coefficient*staggered_volume*gradient**2

end subroutine generated_fci_parallel_flux_power
