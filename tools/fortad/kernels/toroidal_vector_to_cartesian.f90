! Extracted from src/generated/fortfem_toroidal_coordinates.f90 by
! tools/fortad/extract.py. It is the primal fortsym generated
! its products from, with the offload directives removed:
! they are the consumer's choice, not part of what the
! function means.

pure subroutine generated_toroidal_vector_to_cartesian(eta, theta, phi, component_1, component_2, &
        component_3, cartesian_1, cartesian_2, cartesian_3)
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    real(dp), intent(in) :: eta, theta, phi, component_1, component_2, component_3
    real(dp), intent(out) :: cartesian_1, cartesian_2, cartesian_3
    real(dp) :: t1, t2, t3, t4, t5, t6, t7

    t1 = cos(phi)
    t2 = cosh(eta)
    t3 = cos(theta)
    t4 = -t3*t2 + 1.0000000000000000E+000_dp
    t5 = sin(theta)
    t6 = sinh(eta)
    t7 = sin(phi)
    cartesian_1 = component_1*t1*t4/(t2 - t3) - component_2*t1*t5*t6/(t2 - t3) - component_3*t7
    cartesian_2 = component_1*t7*t4/(t2 - t3) - component_2*t7*t5*t6/(t2 - t3) + component_3*t1
    cartesian_3 = -component_1*t5*t6/(t2 - t3) - component_2*t4/(t2 - t3)

end subroutine generated_toroidal_vector_to_cartesian
