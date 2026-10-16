! Extracted from src/generated/fortfem_field_aligned_flux_products.f90 by
! tools/fortad/extract.py. It is the primal fortsym generated
! its products from, with the offload directives removed:
! they are the consumer's choice, not part of what the
! function means.

pure subroutine generated_field_aligned_flux(parallel_coefficient, perpendicular_coefficient, &
        direction_1, direction_2, direction_3, gradient_1, gradient_2, gradient_3, flux_1, flux_2, flux_3)
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    real(dp), intent(in) :: parallel_coefficient, perpendicular_coefficient, direction_1, &
        direction_2, direction_3, gradient_1, gradient_2, gradient_3
    real(dp), intent(out) :: flux_1, flux_2, flux_3
    real(dp) :: t1, t2

    t1 = parallel_coefficient - perpendicular_coefficient
    t2 = direction_1*gradient_1 + direction_2*gradient_2 + direction_3*gradient_3
    flux_1 = direction_1*t1*t2 + gradient_1*perpendicular_coefficient
    flux_2 = direction_2*t1*t2 + gradient_2*perpendicular_coefficient
    flux_3 = direction_3*t1*t2 + gradient_3*perpendicular_coefficient

end subroutine generated_field_aligned_flux
