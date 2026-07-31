module fortfem_toroidal_poisson_dtn_ad
    !! Analytical products for toroidal P harmonics and the exterior DtN map.
    use fortfem_generated_toroidal_poisson_products, only: &
        generated_toroidal_poisson_products
    use fortfem_generated_toroidal_poisson_products_jvp, only: &
        generated_toroidal_poisson_products_jvp
    use fortfem_generated_toroidal_poisson_products_vjp, only: &
        generated_toroidal_poisson_products_vjp
    use fortfem_kinds, only: dp
    use fortfem_toroidal_poisson_dtn, only: &
        evaluate_toroidal_ampere_field_p_primal => &
        evaluate_toroidal_ampere_field_p, &
        evaluate_toroidal_harmonic_p_primal => &
        evaluate_toroidal_harmonic_p, &
        toroidal_poisson_exterior_dtn_p_primal => &
        toroidal_poisson_exterior_dtn_p
    use fortnum_special_toroidal, only: toroidal_p, toroidal_p_derivative
    implicit none

    private

    public :: evaluate_toroidal_harmonic_p_jvp
    public :: evaluate_toroidal_harmonic_p_vjp
    public :: evaluate_toroidal_ampere_field_p_jvp
    public :: evaluate_toroidal_ampere_field_p_vjp
    public :: toroidal_poisson_exterior_dtn_p_jvp
    public :: toroidal_poisson_exterior_dtn_p_vjp

contains

    pure subroutine evaluate_toroidal_harmonic_p_jvp( &
            degree_index, order, eta, theta, phi, eta_dot, theta_dot, &
            phi_dot, value_dot, status)
        integer, intent(in) :: degree_index, order
        real(dp), intent(in) :: eta, theta, phi, eta_dot, theta_dot, phi_dot
        real(dp), intent(out) :: value_dot
        integer, intent(out) :: status

        real(dp) :: radial, radial_derivative, radial_dot
        real(dp) :: radial_derivative_dot
        real(dp) :: unused(5)

        value_dot = 0.0_dp
        status = 1
        call evaluate_toroidal_harmonic_p_primal( &
            degree_index, order, eta, theta, phi, unused(1), status)
        if (status /= 0) return
        radial = toroidal_p(degree_index, order, cosh(eta))
        radial_derivative = toroidal_p_derivative( &
            degree_index, order, cosh(eta))
        radial_dot = radial_derivative*sinh(eta)*eta_dot
        radial_derivative_dot = toroidal_radial_second_derivative( &
            degree_index, order, eta)*sinh(eta)*eta_dot
        call generated_toroidal_poisson_products_jvp( &
            real(degree_index, dp), real(order, dp), 1.0_dp, eta, theta, phi, &
            radial, radial_derivative, 0.0_dp, 0.0_dp, 0.0_dp, eta_dot, &
            theta_dot, phi_dot, radial_dot, radial_derivative_dot, value_dot, &
            unused(1), unused(2), unused(3), unused(4), unused(5))
        status = 0
    end subroutine evaluate_toroidal_harmonic_p_jvp

    pure subroutine evaluate_toroidal_harmonic_p_vjp( &
            degree_index, order, eta, theta, phi, value_bar, eta_bar, &
            theta_bar, phi_bar, status)
        integer, intent(in) :: degree_index, order
        real(dp), intent(in) :: eta, theta, phi, value_bar
        real(dp), intent(out) :: eta_bar, theta_bar, phi_bar
        integer, intent(out) :: status

        real(dp) :: radial, radial_derivative, radial_second_derivative
        real(dp) :: degree_bar, order_bar, scale_bar, radial_bar
        real(dp) :: radial_derivative_bar
        real(dp) :: unused(5)

        eta_bar = 0.0_dp
        theta_bar = 0.0_dp
        phi_bar = 0.0_dp
        status = 1
        call evaluate_toroidal_harmonic_p_primal( &
            degree_index, order, eta, theta, phi, unused(1), status)
        if (status /= 0) return
        radial = toroidal_p(degree_index, order, cosh(eta))
        radial_derivative = toroidal_p_derivative( &
            degree_index, order, cosh(eta))
        radial_second_derivative = toroidal_radial_second_derivative( &
            degree_index, order, eta)
        call generated_toroidal_poisson_products_vjp( &
            real(degree_index, dp), real(order, dp), 1.0_dp, eta, theta, phi, &
            radial, radial_derivative, value_bar, 0.0_dp, 0.0_dp, 0.0_dp, &
            0.0_dp, 0.0_dp, degree_bar, order_bar, scale_bar, eta_bar, &
            theta_bar, phi_bar, radial_bar, radial_derivative_bar)
        eta_bar = eta_bar + (radial_bar*radial_derivative + &
            radial_derivative_bar*radial_second_derivative)*sinh(eta)
        status = 0
    end subroutine evaluate_toroidal_harmonic_p_vjp

    pure subroutine evaluate_toroidal_ampere_field_p_jvp( &
            degree_index, order, scale, eta, theta, phi, scale_dot, eta_dot, &
            theta_dot, phi_dot, field_dot, status)
        integer, intent(in) :: degree_index, order
        real(dp), intent(in) :: scale, eta, theta, phi
        real(dp), intent(in) :: scale_dot, eta_dot, theta_dot, phi_dot
        real(dp), intent(out) :: field_dot(3)
        integer, intent(out) :: status

        real(dp) :: radial, radial_derivative, radial_dot
        real(dp) :: radial_derivative_dot
        real(dp) :: unused(3)

        field_dot = 0.0_dp
        status = 1
        call evaluate_toroidal_ampere_field_p_primal( &
            degree_index, order, scale, eta, theta, phi, unused, status)
        if (status /= 0) return
        radial = toroidal_p(degree_index, order, cosh(eta))
        radial_derivative = toroidal_p_derivative( &
            degree_index, order, cosh(eta))
        radial_dot = radial_derivative*sinh(eta)*eta_dot
        radial_derivative_dot = toroidal_radial_second_derivative( &
            degree_index, order, eta)*sinh(eta)*eta_dot
        call generated_toroidal_poisson_products_jvp( &
            real(degree_index, dp), real(order, dp), scale, eta, theta, phi, &
            radial, radial_derivative, 0.0_dp, 0.0_dp, scale_dot, eta_dot, &
            theta_dot, phi_dot, radial_dot, radial_derivative_dot, unused(1), &
            field_dot(1), field_dot(2), field_dot(3), unused(2), unused(3))
        status = 0
    end subroutine evaluate_toroidal_ampere_field_p_jvp

    pure subroutine evaluate_toroidal_ampere_field_p_vjp( &
            degree_index, order, scale, eta, theta, phi, field_bar, scale_bar, &
            eta_bar, theta_bar, phi_bar, status)
        integer, intent(in) :: degree_index, order
        real(dp), intent(in) :: scale, eta, theta, phi, field_bar(3)
        real(dp), intent(out) :: scale_bar, eta_bar, theta_bar, phi_bar
        integer, intent(out) :: status

        real(dp) :: radial, radial_derivative, radial_second_derivative
        real(dp) :: degree_bar, order_bar, radial_bar, radial_derivative_bar
        real(dp) :: unused(3)

        scale_bar = 0.0_dp
        eta_bar = 0.0_dp
        theta_bar = 0.0_dp
        phi_bar = 0.0_dp
        status = 1
        call evaluate_toroidal_ampere_field_p_primal( &
            degree_index, order, scale, eta, theta, phi, unused, status)
        if (status /= 0) return
        radial = toroidal_p(degree_index, order, cosh(eta))
        radial_derivative = toroidal_p_derivative( &
            degree_index, order, cosh(eta))
        radial_second_derivative = toroidal_radial_second_derivative( &
            degree_index, order, eta)
        call generated_toroidal_poisson_products_vjp( &
            real(degree_index, dp), real(order, dp), scale, eta, theta, phi, &
            radial, radial_derivative, 0.0_dp, field_bar(1), field_bar(2), &
            field_bar(3), 0.0_dp, 0.0_dp, degree_bar, order_bar, scale_bar, &
            eta_bar, theta_bar, phi_bar, radial_bar, radial_derivative_bar)
        eta_bar = eta_bar + (radial_bar*radial_derivative + &
            radial_derivative_bar*radial_second_derivative)*sinh(eta)
        status = 0
    end subroutine evaluate_toroidal_ampere_field_p_vjp

    pure subroutine toroidal_poisson_exterior_dtn_p_jvp( &
            degree_index, order, scale, eta, theta, phi, scale_dot, eta_dot, &
            theta_dot, phi_dot, value_dot, normal_derivative_dot, &
            dtn_value_dot, status)
        integer, intent(in) :: degree_index, order
        real(dp), intent(in) :: scale, eta, theta, phi
        real(dp), intent(in) :: scale_dot, eta_dot, theta_dot, phi_dot
        real(dp), intent(out) :: value_dot, normal_derivative_dot, dtn_value_dot
        integer, intent(out) :: status

        real(dp) :: radial, radial_derivative, radial_dot
        real(dp) :: radial_derivative_dot
        real(dp) :: unused(3)

        value_dot = 0.0_dp
        normal_derivative_dot = 0.0_dp
        dtn_value_dot = 0.0_dp
        status = 1
        call toroidal_poisson_exterior_dtn_p_primal( &
            degree_index, order, scale, eta, theta, phi, unused(1), &
            unused(2), unused(3), status)
        if (status /= 0) return
        radial = toroidal_p(degree_index, order, cosh(eta))
        radial_derivative = toroidal_p_derivative( &
            degree_index, order, cosh(eta))
        radial_dot = radial_derivative*sinh(eta)*eta_dot
        radial_derivative_dot = toroidal_radial_second_derivative( &
            degree_index, order, eta)*sinh(eta)*eta_dot
        call generated_toroidal_poisson_products_jvp( &
            real(degree_index, dp), real(order, dp), scale, eta, theta, phi, &
            radial, radial_derivative, 0.0_dp, 0.0_dp, scale_dot, eta_dot, &
            theta_dot, phi_dot, radial_dot, radial_derivative_dot, value_dot, &
            unused(1), unused(2), unused(3), dtn_value_dot, &
            normal_derivative_dot)
        status = 0
    end subroutine toroidal_poisson_exterior_dtn_p_jvp

    pure subroutine toroidal_poisson_exterior_dtn_p_vjp( &
            degree_index, order, scale, eta, theta, phi, value_bar, &
            normal_derivative_bar, dtn_value_bar, scale_bar, eta_bar, &
            theta_bar, phi_bar, status)
        integer, intent(in) :: degree_index, order
        real(dp), intent(in) :: scale, eta, theta, phi
        real(dp), intent(in) :: value_bar, normal_derivative_bar, dtn_value_bar
        real(dp), intent(out) :: scale_bar, eta_bar, theta_bar, phi_bar
        integer, intent(out) :: status

        real(dp) :: radial, radial_derivative, radial_second_derivative
        real(dp) :: degree_bar, order_bar, radial_bar, radial_derivative_bar

        scale_bar = 0.0_dp
        eta_bar = 0.0_dp
        theta_bar = 0.0_dp
        phi_bar = 0.0_dp
        status = 1
        call toroidal_poisson_exterior_dtn_p_primal( &
            degree_index, order, scale, eta, theta, phi, radial_bar, &
            radial_derivative_bar, radial_second_derivative, status)
        if (status /= 0) return
        radial = toroidal_p(degree_index, order, cosh(eta))
        radial_derivative = toroidal_p_derivative( &
            degree_index, order, cosh(eta))
        radial_second_derivative = toroidal_radial_second_derivative( &
            degree_index, order, eta)
        call generated_toroidal_poisson_products_vjp( &
            real(degree_index, dp), real(order, dp), scale, eta, theta, phi, &
            radial, radial_derivative, value_bar, 0.0_dp, 0.0_dp, 0.0_dp, &
            dtn_value_bar, normal_derivative_bar, degree_bar, order_bar, &
            scale_bar, eta_bar, theta_bar, phi_bar, radial_bar, &
            radial_derivative_bar)
        eta_bar = eta_bar + (radial_bar*radial_derivative + &
            radial_derivative_bar*radial_second_derivative)*sinh(eta)
        status = 0
    end subroutine toroidal_poisson_exterior_dtn_p_vjp

    pure function toroidal_radial_second_derivative( &
            degree_index, order, eta) result(value)
        integer, intent(in) :: degree_index, order
        real(dp), intent(in) :: eta
        real(dp) :: value
        real(dp) :: argument, denominator, radial, radial_derivative
        real(dp) :: degree, order_real

        argument = cosh(eta)
        denominator = argument*argument - 1.0_dp
        radial = toroidal_p(degree_index, order, argument)
        radial_derivative = toroidal_p_derivative( &
            degree_index, order, argument)
        degree = real(degree_index, dp) - 0.5_dp
        order_real = real(order, dp)
        value = (degree*(degree + 1.0_dp) + order_real**2/denominator)* &
            radial/denominator - 2.0_dp*argument*radial_derivative/denominator
    end function toroidal_radial_second_derivative

end module fortfem_toroidal_poisson_dtn_ad
