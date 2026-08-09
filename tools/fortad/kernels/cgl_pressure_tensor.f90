! Extracted from src/generated/fortfem_cgl_pressure_tensor_products.f90 by
! tools/fortad/extract.py. It is the primal fortsym generated
! its products from, with the offload directives removed:
! they are the consumer's choice, not part of what the
! function means.

pure subroutine generated_cgl_pressure_tensor(p_parallel, p_perpendicular, direction_1, direction_2, &
        direction_3, pressure_11, pressure_22, pressure_33, pressure_12, pressure_13, pressure_23)
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    real(dp), intent(in) :: p_parallel, p_perpendicular, direction_1, direction_2, direction_3
    real(dp), intent(out) :: pressure_11, pressure_22, pressure_33, pressure_12, pressure_13, &
        pressure_23
    real(dp) :: t1

    t1 = p_parallel - p_perpendicular
    pressure_11 = p_perpendicular + direction_1**2*t1
    pressure_22 = p_perpendicular + direction_2**2*t1
    pressure_33 = p_perpendicular + direction_3**2*t1
    pressure_12 = direction_1*direction_2*t1
    pressure_13 = direction_1*direction_3*t1
    pressure_23 = direction_2*direction_3*t1

end subroutine generated_cgl_pressure_tensor
