! Extracted from src/generated/fortfem_fci_parallel_products.f90 by
! tools/fortad/extract.py. It is the primal fortsym generated
! its products from, with the offload directives removed:
! they are the consumer's choice, not part of what the
! function means.

pure subroutine generated_fci_parallel_diffusion(forward_value, upper_field, backward_value, &
        lower_field, line_length, parallel_coefficient, canonical_volume, staggered_volume, &
        lower_contribution, upper_contribution)
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    real(dp), intent(in) :: forward_value, upper_field, backward_value, lower_field, line_length, &
        parallel_coefficient, canonical_volume, staggered_volume
    real(dp), intent(out) :: lower_contribution, upper_contribution
    real(dp) :: t1

    t1 = -backward_value*lower_field + forward_value*upper_field
    lower_contribution = backward_value*parallel_coefficient*staggered_volume*t1/canonical_volume/ &
        line_length**2
    upper_contribution = -forward_value*parallel_coefficient*staggered_volume*t1/canonical_volume/ &
        line_length**2

end subroutine generated_fci_parallel_diffusion
