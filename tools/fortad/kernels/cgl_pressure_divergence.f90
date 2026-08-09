! Extracted from src/generated/fortfem_cgl_pressure_divergence_products.f90 by
! tools/fortad/extract.py. It is the primal fortsym generated
! its products from, with the offload directives removed:
! they are the consumer's choice, not part of what the
! function means.

pure subroutine generated_cgl_pressure_divergence(p_parallel, p_perpendicular, direction_1, &
        direction_2, direction_3, parallel_gradient_1, parallel_gradient_2, parallel_gradient_3, &
        perpendicular_gradient_1, perpendicular_gradient_2, perpendicular_gradient_3, direction_gradient_11, &
        direction_gradient_12, direction_gradient_13, direction_gradient_21, direction_gradient_22, &
        direction_gradient_23, direction_gradient_31, direction_gradient_32, direction_gradient_33, force_1, &
        force_2, force_3)
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    real(dp), intent(in) :: p_parallel, p_perpendicular, direction_1, direction_2, direction_3, &
        parallel_gradient_1, parallel_gradient_2, parallel_gradient_3, perpendicular_gradient_1, &
        perpendicular_gradient_2, perpendicular_gradient_3, direction_gradient_11, direction_gradient_12, &
        direction_gradient_13, direction_gradient_21, direction_gradient_22, direction_gradient_23, &
        direction_gradient_31, direction_gradient_32, direction_gradient_33
    real(dp), intent(out) :: force_1, force_2, force_3
    real(dp) :: t1, t2, t3, t4

    t1 = parallel_gradient_2 - perpendicular_gradient_2
    t2 = parallel_gradient_3 - perpendicular_gradient_3
    t3 = p_parallel - p_perpendicular
    t4 = parallel_gradient_1 - perpendicular_gradient_1
    force_1 = perpendicular_gradient_1 + direction_1*direction_2*t1 + direction_1*direction_3*t2 + &
        direction_1*direction_gradient_11*t3*2 + direction_1**2*t4 + t3*(direction_1*direction_gradient_22 + &
        direction_2*direction_gradient_12) + t3*(direction_1*direction_gradient_33 + direction_3* &
        direction_gradient_13)
    force_2 = perpendicular_gradient_2 + direction_1*direction_2*t4 + direction_2*direction_3*t2 + &
        direction_2*direction_gradient_22*t3*2 + direction_2**2*t1 + t3*(direction_1*direction_gradient_21 + &
        direction_2*direction_gradient_11) + t3*(direction_2*direction_gradient_33 + direction_3* &
        direction_gradient_23)
    force_3 = perpendicular_gradient_3 + direction_1*direction_3*t4 + direction_2*direction_3*t1 + &
        direction_3*direction_gradient_33*t3*2 + direction_3**2*t2 + t3*(direction_1*direction_gradient_31 + &
        direction_3*direction_gradient_11) + t3*(direction_2*direction_gradient_32 + direction_3* &
        direction_gradient_22)

end subroutine generated_cgl_pressure_divergence
