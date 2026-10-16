! Extracted from src/generated/fortfem_field_aligned_hall_products.f90 by
! tools/fortad/extract.py. It is the primal fortsym generated
! its products from, with the offload directives removed:
! they are the consumer's choice, not part of what the
! function means.

pure subroutine generated_field_aligned_hall(hall_coefficient, direction_1, direction_2, &
        direction_3, hall_direction_1, hall_direction_2, hall_direction_3)
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    real(dp), intent(in) :: hall_coefficient, direction_1, direction_2, direction_3
    real(dp), intent(out) :: hall_direction_1, hall_direction_2, hall_direction_3

    hall_direction_1 = direction_1*hall_coefficient
    hall_direction_2 = direction_2*hall_coefficient
    hall_direction_3 = direction_3*hall_coefficient

end subroutine generated_field_aligned_hall
