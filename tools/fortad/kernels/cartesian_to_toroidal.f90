! Extracted from src/generated/fortfem_toroidal_coordinates.f90 by
! tools/fortad/extract.py. It is the primal fortsym generated
! its products from, with the offload directives removed:
! they are the consumer's choice, not part of what the
! function means.

pure subroutine generated_cartesian_to_toroidal(point_1, point_2, point_3, scale, eta, theta, phi)
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    real(dp), intent(in) :: point_1, point_2, point_3, scale
    real(dp), intent(out) :: eta, theta, phi
    real(dp) :: t1, t2

    t1 = sqrt(point_1*point_1 + point_2*point_2)
    t2 = point_3*point_3
    eta = log(((scale + t1)**2 + t2)/((t1 - scale)**2 + t2))*5.0000000000000000E-001_dp
    theta = atan2(point_3*scale*2.0000000000000000E+000_dp, t2 - scale*scale + t1*t1)
    phi = atan2(point_2, point_1)

end subroutine generated_cartesian_to_toroidal
