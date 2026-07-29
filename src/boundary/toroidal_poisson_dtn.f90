module fortfem_toroidal_poisson_dtn
    ! Scalar toroidal harmonics and the exterior Laplace DtN map.
    !
    ! The computational interior of eta=eta0 is eta>eta0. Its outward normal
    ! points toward decreasing eta, giving d/dn=-(cosh(eta)-cos(theta))/a
    ! times d/deta. The separated radial function is Hobson
    ! P_(n-1/2)^m(cosh eta), consistent with DLMF 14.19.
    use fortfem_kinds, only: dp
    use fortnum_special_toroidal, only: toroidal_p, toroidal_p_derivative
    implicit none
    private

    public :: evaluate_toroidal_harmonic_p
    public :: evaluate_toroidal_ampere_field_p
    public :: toroidal_poisson_exterior_dtn_p

contains

    pure subroutine evaluate_toroidal_harmonic_p( &
            degree_index, order, eta, theta, phi, value, status)
        integer, intent(in) :: degree_index, order
        real(dp), intent(in) :: eta, theta, phi
        real(dp), intent(out) :: value
        integer, intent(out) :: status
        real(dp) :: denominator, radial

        if (degree_index < 0 .or. order < 0 .or. eta <= 0.0_dp) then
            value = 0.0_dp
            status = 1
            return
        end if
        denominator = cosh(eta) - cos(theta)
        radial = toroidal_p(degree_index, order, cosh(eta))
        value = sqrt(denominator)*radial* &
            cos(real(degree_index, dp)*theta)* &
            cos(real(order, dp)*phi)
        status = 0
    end subroutine evaluate_toroidal_harmonic_p

    pure subroutine evaluate_toroidal_ampere_field_p( &
            degree_index, order, scale, eta, theta, phi, field, status)
        integer, intent(in) :: degree_index, order
        real(dp), intent(in) :: scale, eta, theta, phi
        real(dp), intent(out) :: field(3)
        integer, intent(out) :: status
        real(dp) :: denominator, radial, radial_derivative
        real(dp) :: angular_theta, angular_phi
        real(dp) :: derivative_eta, derivative_theta, derivative_phi

        if (degree_index < 0 .or. order < 0 .or. &
            scale <= 0.0_dp .or. eta <= 0.0_dp) then
            field = 0.0_dp
            status = 1
            return
        end if

        denominator = cosh(eta) - cos(theta)
        radial = toroidal_p(degree_index, order, cosh(eta))
        radial_derivative = &
            toroidal_p_derivative(degree_index, order, cosh(eta))
        angular_theta = real(degree_index, dp)*theta
        angular_phi = real(order, dp)*phi

        derivative_eta = ( &
            sinh(eta)*radial/(2.0_dp*sqrt(denominator)) + &
            sqrt(denominator)*radial_derivative*sinh(eta))* &
            cos(angular_theta)*cos(angular_phi)
        derivative_theta = ( &
            sin(theta)*radial*cos(angular_theta)/ &
            (2.0_dp*sqrt(denominator)) - &
            sqrt(denominator)*radial*real(degree_index, dp)* &
            sin(angular_theta))*cos(angular_phi)
        derivative_phi = -sqrt(denominator)*radial* &
            real(order, dp)*cos(angular_theta)*sin(angular_phi)

        field(1) = -denominator*derivative_eta/scale
        field(2) = -denominator*derivative_theta/scale
        field(3) = -denominator*derivative_phi/(scale*sinh(eta))
        status = 0
    end subroutine evaluate_toroidal_ampere_field_p

    pure subroutine toroidal_poisson_exterior_dtn_p( &
            degree_index, order, scale, eta, theta, phi, &
            value, normal_derivative, dtn_value, status)
        integer, intent(in) :: degree_index, order
        real(dp), intent(in) :: scale, eta, theta, phi
        real(dp), intent(out) :: value, normal_derivative, dtn_value
        integer, intent(out) :: status
        real(dp) :: denominator, radial, radial_derivative
        real(dp) :: eta_derivative

        if (scale <= 0.0_dp) then
            value = 0.0_dp
            normal_derivative = 0.0_dp
            dtn_value = 0.0_dp
            status = 1
            return
        end if
        call evaluate_toroidal_harmonic_p( &
            degree_index, order, eta, theta, phi, value, status)
        if (status /= 0) then
            normal_derivative = 0.0_dp
            dtn_value = 0.0_dp
            return
        end if

        denominator = cosh(eta) - cos(theta)
        radial = toroidal_p(degree_index, order, cosh(eta))
        radial_derivative = &
            toroidal_p_derivative(degree_index, order, cosh(eta))
        if (radial == 0.0_dp) then
            normal_derivative = 0.0_dp
            dtn_value = 0.0_dp
            status = 2
            return
        end if
        eta_derivative = sinh(eta)/(2.0_dp*denominator) + &
            sinh(eta)*radial_derivative/radial
        dtn_value = -denominator*eta_derivative/scale
        normal_derivative = dtn_value*value
        status = 0
    end subroutine toroidal_poisson_exterior_dtn_p

end module fortfem_toroidal_poisson_dtn
