program test_magnetic_curvilinear_coefficients_2d
    use check, only: check_condition, check_summary
    use fortfem_api, only: &
        scalar_reluctivity_curvilinear_fourier_coefficients
    use fortfem_kinds, only: dp
    implicit none

    real(dp), parameter :: pi = acos(-1.0_dp)
    real(dp) :: curl_weight, expected(2, 2), mass_tensor(2, 2)
    real(dp) :: metric(3, 3), reluctivity
    real(dp) :: alpha, beta, gamma, radius, theta, period
    integer :: status
    logical :: all_passed

    all_passed = .true.
    reluctivity = 2.75_dp

    period = 3.5_dp
    metric = 0.0_dp
    metric(1, 1) = 1.0_dp
    metric(2, 2) = 1.0_dp
    metric(3, 3) = (period / (2.0_dp * pi))**2
    call scalar_reluctivity_curvilinear_fourier_coefficients( &
        metric, reluctivity, curl_weight, mass_tensor, status)
    expected = 0.0_dp
    expected(1, 1) = reluctivity * 2.0_dp * pi / period
    expected(2, 2) = expected(1, 1)
    call record_condition(status == 0, &
        "Cartesian periodic metric is accepted")
    call record_condition(abs( &
        curl_weight - reluctivity * period / (2.0_dp * pi)) < 2.0e-14_dp, &
        "Periodic Cartesian metric yields its exact curl coefficient")
    call record_condition(maxval(abs(mass_tensor - expected)) < 2.0e-14_dp, &
        "Periodic Cartesian metric yields its exact transverse tensor")

    radius = 1.7_dp
    metric = 0.0_dp
    metric(1, 1) = 1.0_dp
    metric(2, 2) = 1.0_dp
    metric(3, 3) = radius**2
    call scalar_reluctivity_curvilinear_fourier_coefficients( &
        metric, reluctivity, curl_weight, mass_tensor, status)
    expected = 0.0_dp
    expected(1, 1) = reluctivity / radius
    expected(2, 2) = expected(1, 1)
    call record_condition(abs(curl_weight - reluctivity * radius) < &
        2.0e-14_dp, "Cylindrical metric yields R times reluctivity")
    call record_condition(maxval(abs(mass_tensor - expected)) < 2.0e-14_dp, &
        "Cylindrical metric yields reluctivity divided by R")

    radius = 2.3_dp
    theta = 0.8_dp
    metric = 0.0_dp
    metric(1, 1) = 1.0_dp
    metric(2, 2) = radius**2
    metric(3, 3) = radius**2 * sin(theta)**2
    call scalar_reluctivity_curvilinear_fourier_coefficients( &
        metric, reluctivity, curl_weight, mass_tensor, status)
    expected = 0.0_dp
    expected(1, 1) = reluctivity / sin(theta)
    expected(2, 2) = reluctivity / (radius**2 * sin(theta))
    call record_condition(abs(curl_weight - reluctivity * sin(theta)) < &
        2.0e-14_dp, "Spherical metric yields the paper curl coefficient")
    call record_condition(maxval(abs(mass_tensor - expected)) < 2.0e-14_dp, &
        "Spherical metric yields the paper transverse tensor")

    alpha = 0.4_dp
    beta = 1.3_dp
    gamma = 0.7_dp
    metric = 0.0_dp
    metric(1, 1) = 1.0_dp
    metric(1, 2) = alpha
    metric(2, 1) = alpha
    metric(2, 2) = alpha**2 + beta**2
    metric(3, 3) = gamma**2
    call scalar_reluctivity_curvilinear_fourier_coefficients( &
        metric, reluctivity, curl_weight, mass_tensor, status)
    expected(1, 1) = &
        reluctivity * (alpha**2 + beta**2) / (beta * gamma)
    expected(1, 2) = -reluctivity * alpha / (beta * gamma)
    expected(2, 1) = expected(1, 2)
    expected(2, 2) = reluctivity / (beta * gamma)
    call record_condition(abs(curl_weight - reluctivity * gamma / beta) < &
        2.0e-14_dp, "Sheared chart yields its exact curl density")
    call record_condition(maxval(abs(mass_tensor - expected)) < 2.0e-14_dp, &
        "Sheared chart yields its independently inverted metric tensor")

    metric(2, 2) = alpha**2
    call scalar_reluctivity_curvilinear_fourier_coefficients( &
        metric, reluctivity, curl_weight, mass_tensor, status)
    call record_condition(status /= 0, &
        "Singular transverse metrics are rejected")

    metric = 0.0_dp
    metric(1, 1) = -1.0_dp
    metric(2, 2) = -1.0_dp
    metric(3, 3) = 1.0_dp
    call scalar_reluctivity_curvilinear_fourier_coefficients( &
        metric, reluctivity, curl_weight, mass_tensor, status)
    call record_condition(status /= 0, &
        "Negative-definite transverse metrics are rejected")

    call check_summary("Magnetic curvilinear coefficients")
    if (.not. all_passed) error stop 1

contains

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_magnetic_curvilinear_coefficients_2d
