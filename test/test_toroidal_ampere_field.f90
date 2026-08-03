program test_toroidal_ampere_field
    use check, only: check_condition, check_summary
    use fortfem_fourier, only: evaluate_toroidal_ampere_field_p
    use fortfem_kinds, only: dp
    implicit none

    real(dp), parameter :: scale = 1.7_dp
    real(dp), parameter :: eta = &
        1.3169578969248167086250463473079684_dp
    real(dp), parameter :: theta = 0.6_dp
    real(dp), parameter :: phi = 0.4_dp
    real(dp), parameter :: reference(3) = [ &
        -2.5743996635144626758125703168794369_dp, &
        5.4760717158563791944410421313196825_dp, &
        0.27257717557660818930919485899476534_dp]
    real(dp) :: field(3)
    real(dp) :: curl_field(3)
    integer :: status
    logical :: passed

    call evaluate_toroidal_ampere_field_p( &
        2, 1, scale, eta, theta, phi, field, status)
    passed = status == 0 .and. &
        maxval(abs(field - reference)) < 1.2e-11_dp
    call check_condition(passed, &
        "Toroidal Ampere field matches a 50-digit mpmath gradient oracle")
    call numerical_curl(curl_field)
    passed = passed .and. maxval(abs(curl_field)) < 2.0e-8_dp
    call check_condition(maxval(abs(curl_field)) < 2.0e-8_dp, &
        "Toroidal Ampere field satisfies curl H = 0")

    call evaluate_toroidal_ampere_field_p( &
        2, 1, -1.0_dp, eta, theta, phi, field, status)
    call check_condition(status /= 0, &
        "Toroidal Ampere field rejects a negative scale")

    call check_summary("Toroidal Ampere analytical field")
    if (.not. passed) error stop 1

contains

    subroutine numerical_curl(curl_value)
        real(dp), intent(out) :: curl_value(3)
        real(dp), parameter :: step = 2.0e-5_dp
        real(dp) :: denominator, scale_factor(3)

        denominator = cosh(eta) - cos(theta)
        scale_factor = [scale/denominator, scale/denominator, &
            scale*sinh(eta)/denominator]
        curl_value(1) = ( &
            (weighted(3, eta, theta + step, phi) - &
            weighted(3, eta, theta - step, phi)) - &
            (weighted(2, eta, theta, phi + step) - &
            weighted(2, eta, theta, phi - step)))/ &
            (2.0_dp*step*scale_factor(2)*scale_factor(3))
        curl_value(2) = ( &
            (weighted(1, eta, theta, phi + step) - &
            weighted(1, eta, theta, phi - step)) - &
            (weighted(3, eta + step, theta, phi) - &
            weighted(3, eta - step, theta, phi)))/ &
            (2.0_dp*step*scale_factor(3)*scale_factor(1))
        curl_value(3) = ( &
            (weighted(2, eta + step, theta, phi) - &
            weighted(2, eta - step, theta, phi)) - &
            (weighted(1, eta, theta + step, phi) - &
            weighted(1, eta, theta - step, phi)))/ &
            (2.0_dp*step*scale_factor(1)*scale_factor(2))
    end subroutine numerical_curl

    function weighted(component, eta_value, theta_value, phi_value) &
            result(value)
        integer, intent(in) :: component
        real(dp), intent(in) :: eta_value, theta_value, phi_value
        real(dp) :: value
        real(dp) :: local_field(3), denominator, scale_factor(3)
        integer :: local_status

        call evaluate_toroidal_ampere_field_p( &
            2, 1, scale, eta_value, theta_value, phi_value, &
            local_field, local_status)
        if (local_status /= 0) error stop "Unexpected toroidal field status"
        denominator = cosh(eta_value) - cos(theta_value)
        scale_factor = [scale/denominator, scale/denominator, &
            scale*sinh(eta_value)/denominator]
        value = scale_factor(component)*local_field(component)
    end function weighted

end program test_toroidal_ampere_field
