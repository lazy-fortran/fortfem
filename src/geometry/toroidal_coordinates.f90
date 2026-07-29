module fortfem_toroidal_coordinates
    use fortfem_kinds, only: dp
    implicit none
    private

    public :: cartesian_to_toroidal
    public :: toroidal_point_to_cartesian
    public :: toroidal_vector_to_cartesian

contains

    pure subroutine cartesian_to_toroidal(point, scale, eta, theta, phi)
        real(dp), intent(in) :: point(3), scale
        real(dp), intent(out) :: eta, theta, phi
        real(dp) :: cylindrical_radius

        cylindrical_radius = norm2(point(1:2))
        eta = 0.5_dp*log( &
            ((cylindrical_radius + scale)**2 + point(3)**2)/ &
            ((cylindrical_radius - scale)**2 + point(3)**2))
        theta = atan2( &
            2.0_dp*scale*point(3), &
            cylindrical_radius**2 + point(3)**2 - scale**2)
        phi = atan2(point(2), point(1))
    end subroutine cartesian_to_toroidal

    pure subroutine toroidal_point_to_cartesian( &
            scale, eta, theta, phi, point)
        real(dp), intent(in) :: scale, eta, theta, phi
        real(dp), intent(out) :: point(3)
        real(dp) :: denominator

        denominator = cosh(eta) - cos(theta)
        point = [ &
            scale*sinh(eta)*cos(phi)/denominator, &
            scale*sinh(eta)*sin(phi)/denominator, &
            scale*sin(theta)/denominator]
    end subroutine toroidal_point_to_cartesian

    pure subroutine toroidal_vector_to_cartesian( &
            eta, theta, phi, components, cartesian)
        real(dp), intent(in) :: eta, theta, phi, components(3)
        real(dp), intent(out) :: cartesian(3)
        real(dp) :: azimuthal(3), denominator, eta_radial, eta_vertical
        real(dp) :: radial(3), theta_radial, theta_vertical

        denominator = cosh(eta) - cos(theta)
        eta_radial = (1.0_dp - cosh(eta)*cos(theta))/denominator
        eta_vertical = -sinh(eta)*sin(theta)/denominator
        theta_radial = eta_vertical
        theta_vertical = -eta_radial
        radial = [cos(phi), sin(phi), 0.0_dp]
        azimuthal = [-sin(phi), cos(phi), 0.0_dp]
        cartesian = components(1)*( &
            eta_radial*radial + [0.0_dp, 0.0_dp, eta_vertical]) + &
            components(2)*( &
            theta_radial*radial + [0.0_dp, 0.0_dp, theta_vertical]) + &
            components(3)*azimuthal
    end subroutine toroidal_vector_to_cartesian

end module fortfem_toroidal_coordinates
