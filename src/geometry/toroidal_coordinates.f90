module fortfem_toroidal_coordinates
    use fortfem_kinds, only: dp
    use fortfem_generated_toroidal_coordinates, only: &
        generated_toroidal_point_to_cartesian
    use fortfem_generated_toroidal_coordinates_jvp, only: &
        generated_toroidal_point_to_cartesian_jvp
    use fortfem_generated_toroidal_coordinates_vjp, only: &
        generated_toroidal_point_to_cartesian_vjp
    use fortfem_generated_cartesian_to_toroidal, only: &
        generated_cartesian_to_toroidal
    use fortfem_generated_cartesian_to_toroidal_jvp, only: &
        generated_cartesian_to_toroidal_jvp
    use fortfem_generated_cartesian_to_toroidal_vjp, only: &
        generated_cartesian_to_toroidal_vjp
    use fortfem_generated_toroidal_vector_to_cartesian_jvp, only: &
        generated_toroidal_vector_to_cartesian_jvp
    use fortfem_generated_toroidal_vector_to_cartesian_vjp, only: &
        generated_toroidal_vector_to_cartesian_vjp
    implicit none
    private

    public :: cartesian_to_toroidal
    public :: cartesian_to_toroidal_jvp
    public :: cartesian_to_toroidal_vjp
    public :: toroidal_point_to_cartesian
    public :: toroidal_point_to_cartesian_jvp
    public :: toroidal_point_to_cartesian_vjp
    public :: toroidal_vector_to_cartesian
    public :: toroidal_vector_to_cartesian_jvp
    public :: toroidal_vector_to_cartesian_vjp

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
        call generated_toroidal_point_to_cartesian( &
            scale, eta, theta, phi, point)
    end subroutine toroidal_point_to_cartesian

    pure subroutine toroidal_point_to_cartesian_jvp( &
            scale, eta, theta, phi, scale_dot, eta_dot, theta_dot, phi_dot, &
            point_dot, status)
        real(dp), intent(in) :: scale, eta, theta, phi
        real(dp), intent(in) :: scale_dot, eta_dot, theta_dot, phi_dot
        real(dp), intent(out) :: point_dot(3)
        integer, intent(out) :: status
        real(dp) :: denominator

        point_dot = 0.0_dp
        status = 1
        if (scale <= 0.0_dp .or. eta <= 0.0_dp) return
        denominator = cosh(eta) - cos(theta)
        if (denominator <= tiny(1.0_dp)) return
        call generated_toroidal_point_to_cartesian_jvp( &
            scale, eta, theta, phi, scale_dot, eta_dot, theta_dot, phi_dot, &
            point_dot)
        status = 0
    end subroutine toroidal_point_to_cartesian_jvp

    pure subroutine toroidal_point_to_cartesian_vjp( &
            scale, eta, theta, phi, point_bar, scale_bar, eta_bar, theta_bar, &
            phi_bar, status)
        real(dp), intent(in) :: scale, eta, theta, phi, point_bar(3)
        real(dp), intent(out) :: scale_bar, eta_bar, theta_bar, phi_bar
        integer, intent(out) :: status
        real(dp) :: denominator

        scale_bar = 0.0_dp
        eta_bar = 0.0_dp
        theta_bar = 0.0_dp
        phi_bar = 0.0_dp
        status = 1
        if (scale <= 0.0_dp .or. eta <= 0.0_dp) return
        denominator = cosh(eta) - cos(theta)
        if (denominator <= tiny(1.0_dp)) return
        call generated_toroidal_point_to_cartesian_vjp( &
            scale, eta, theta, phi, point_bar(1), point_bar(2), point_bar(3), &
            scale_bar, eta_bar, theta_bar, phi_bar)
        status = 0
    end subroutine toroidal_point_to_cartesian_vjp

    pure subroutine cartesian_to_toroidal_jvp( &
            point, scale, point_dot, scale_dot, eta_dot, theta_dot, phi_dot, &
            status)
        real(dp), intent(in) :: point(3), scale, point_dot(3), scale_dot
        real(dp), intent(out) :: eta_dot, theta_dot, phi_dot
        integer, intent(out) :: status
        real(dp) :: cylindrical_radius

        eta_dot = 0.0_dp
        theta_dot = 0.0_dp
        phi_dot = 0.0_dp
        status = 1
        if (scale <= 0.0_dp) return
        cylindrical_radius = norm2(point(1:2))
        if (cylindrical_radius <= tiny(1.0_dp)) return
        call generated_cartesian_to_toroidal_jvp( &
            point(1), point(2), point(3), scale, point_dot(1), point_dot(2), &
            point_dot(3), scale_dot, eta_dot, theta_dot, phi_dot)
        status = 0
    end subroutine cartesian_to_toroidal_jvp

    pure subroutine cartesian_to_toroidal_vjp( &
            point, scale, eta_bar, theta_bar, phi_bar, point_bar, scale_bar, &
            status)
        real(dp), intent(in) :: point(3), scale, eta_bar, theta_bar, phi_bar
        real(dp), intent(out) :: point_bar(3), scale_bar
        integer, intent(out) :: status
        real(dp) :: cylindrical_radius

        point_bar = 0.0_dp
        scale_bar = 0.0_dp
        status = 1
        if (scale <= 0.0_dp) return
        cylindrical_radius = norm2(point(1:2))
        if (cylindrical_radius <= tiny(1.0_dp)) return
        call generated_cartesian_to_toroidal_vjp( &
            point(1), point(2), point(3), scale, eta_bar, theta_bar, phi_bar, &
            point_bar(1), point_bar(2), point_bar(3), scale_bar)
        status = 0
    end subroutine cartesian_to_toroidal_vjp

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

    pure subroutine toroidal_vector_to_cartesian_jvp( &
            eta, theta, phi, components, eta_dot, theta_dot, phi_dot, &
            components_dot, cartesian_dot, status)
        real(dp), intent(in) :: eta, theta, phi, components(3)
        real(dp), intent(in) :: eta_dot, theta_dot, phi_dot, components_dot(3)
        real(dp), intent(out) :: cartesian_dot(3)
        integer, intent(out) :: status
        real(dp) :: denominator

        cartesian_dot = 0.0_dp
        status = 1
        if (eta <= 0.0_dp) return
        denominator = cosh(eta) - cos(theta)
        if (denominator <= tiny(1.0_dp)) return
        call generated_toroidal_vector_to_cartesian_jvp( &
            eta, theta, phi, components(1), components(2), components(3), &
            eta_dot, theta_dot, phi_dot, components_dot(1), components_dot(2), &
            components_dot(3), cartesian_dot(1), cartesian_dot(2), &
            cartesian_dot(3))
        status = 0
    end subroutine toroidal_vector_to_cartesian_jvp

    pure subroutine toroidal_vector_to_cartesian_vjp( &
            eta, theta, phi, components, cartesian_bar, eta_bar, theta_bar, &
            phi_bar, components_bar, status)
        real(dp), intent(in) :: eta, theta, phi, components(3), cartesian_bar(3)
        real(dp), intent(out) :: eta_bar, theta_bar, phi_bar, components_bar(3)
        integer, intent(out) :: status
        real(dp) :: denominator

        eta_bar = 0.0_dp
        theta_bar = 0.0_dp
        phi_bar = 0.0_dp
        components_bar = 0.0_dp
        status = 1
        if (eta <= 0.0_dp) return
        denominator = cosh(eta) - cos(theta)
        if (denominator <= tiny(1.0_dp)) return
        call generated_toroidal_vector_to_cartesian_vjp( &
            eta, theta, phi, components(1), components(2), components(3), &
            cartesian_bar(1), cartesian_bar(2), cartesian_bar(3), eta_bar, &
            theta_bar, phi_bar, components_bar(1), components_bar(2), &
            components_bar(3))
        status = 0
    end subroutine toroidal_vector_to_cartesian_vjp

end module fortfem_toroidal_coordinates
