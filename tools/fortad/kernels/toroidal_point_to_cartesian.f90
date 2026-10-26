! Extracted from src/generated/fortfem_toroidal_coordinates.f90 by
! tools/fortad/extract.py. It is the primal fortsym generated
! its products from, with the offload directives removed:
! they are the consumer's choice, not part of what the
! function means.

pure subroutine generated_toroidal_point_to_cartesian(scale, eta, theta, phi, point)
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    real(dp), intent(in) :: scale, eta, theta, phi
    real(dp), intent(out) :: point(3)
    real(dp) :: t1, t2, t3

    t1 = sinh(eta)
    t2 = cos(theta)
    t3 = cosh(eta) - t2
    point(1) = scale*cos(phi)*t1/t3
    point(2) = scale*sin(phi)*t1/t3
    point(3) = scale*sin(theta)/t3

end subroutine generated_toroidal_point_to_cartesian
