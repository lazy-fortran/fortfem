! Extracted from src/generated/fortfem_laplace_bem_kernel_products.f90 by
! tools/fortad/extract.py. It is the primal fortsym generated
! its products from, with the offload directives removed:
! they are the consumer's choice, not part of what the
! function means.

pure subroutine generated_laplace_single_layer_integrand(first_point_1, first_point_2, &
        first_point_3, second_point_1, second_point_2, second_point_3, first_jacobian, second_jacobian, &
        kernel_scale, value)
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    real(dp), intent(in) :: first_point_1, first_point_2, first_point_3, second_point_1, &
        second_point_2, second_point_3, first_jacobian, second_jacobian, kernel_scale
    real(dp), intent(out) :: value
    real(dp) :: t1, t2, t3

    t1 = first_point_1 - second_point_1
    t2 = first_point_2 - second_point_2
    t3 = first_point_3 - second_point_3
    value = first_jacobian*kernel_scale*second_jacobian/sqrt(t1*t1 + t2*t2 + t3*t3)

end subroutine generated_laplace_single_layer_integrand
