! Extracted from src/generated/fortfem_toroidal_poisson_products.f90 by
! tools/fortad/extract.py. It is the primal fortsym generated
! its products from, with the offload directives removed:
! they are the consumer's choice, not part of what the
! function means.

pure subroutine generated_toroidal_poisson_products(degree, order, scale, eta, theta, phi, radial, &
        radial_derivative, harmonic, field_1, field_2, field_3, dtn_value, normal_derivative)
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    real(dp), intent(in) :: degree, order, scale, eta, theta, phi, radial, radial_derivative
    real(dp), intent(out) :: harmonic, field_1, field_2, field_3, dtn_value, normal_derivative
    real(dp) :: t1, t2, t3, t4, t5, t6, t7, t8

    t1 = degree*theta
    t2 = cos(t1)
    t3 = order*phi
    t4 = cos(t3)
    t5 = cosh(eta) - cos(theta)
    t6 = sqrt(t5)
    t7 = sinh(eta)
    t8 = radial_derivative*t7/radial + t7/(t5*2)
    harmonic = radial*t2*t4*t6
    field_1 = -t2*t4*t5*(radial*t7/(t6*2) + radial_derivative*t7*t6)/scale
    field_2 = -t4*t5*(-degree*radial*sin(t1)*t6 + radial*t2*sin(theta)/(t6*2))/scale
    field_3 = order*radial*t2*sin(t3)*t6*t5/(scale*t7)
    dtn_value = -t5*t8/scale
    normal_derivative = -radial*t2*t4*t6*t5*t8/scale

end subroutine generated_toroidal_poisson_products
