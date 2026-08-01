program test_fortnum_spherical_dependency
    ! Integration oracle for the FortNum revision consumed by FortFEM.
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use fortfem_api, only: spherical_harmonic, &
        spherical_harmonic_product_coefficient, &
        spherical_harmonic_theta_derivative, spherical_harmonic_phi_derivative
    implicit none

    real(dp), parameter :: pi = 3.1415926535897932384626433832795_dp
    real(dp), parameter :: expected = 1.0_dp/sqrt(4.0_dp*pi)
    complex(dp) :: value, theta_derivative, phi_derivative
    real(dp) :: product_coefficient

    value = spherical_harmonic(0, 0, 0.7_dp, 1.2_dp)
    if (abs(value - cmplx(expected, 0.0_dp, dp)) > 3.0e-13_dp) then
        write (error_unit, "(a,2es24.16)") &
            "FortNum spherical-harmonic pin failed: ", real(value, dp), &
            aimag(value)
        error stop 1
    end if
    theta_derivative = spherical_harmonic_theta_derivative(0, 0, 0.7_dp, 1.2_dp)
    phi_derivative = spherical_harmonic_phi_derivative(0, 0, 0.7_dp, 1.2_dp)
    if (abs(theta_derivative) > 3.0e-13_dp .or. &
        abs(phi_derivative) > 3.0e-13_dp) then
        write (error_unit, "(a,4es24.16)") &
            "FortNum spherical-harmonic derivative pin failed: ", &
            real(theta_derivative, dp), aimag(theta_derivative), &
            real(phi_derivative, dp), aimag(phi_derivative)
        error stop 1
    end if
    product_coefficient = spherical_harmonic_product_coefficient( &
        1, 0, 1, 0, 0, 0)
    if (abs(product_coefficient - expected) > 3.0e-13_dp) then
        write (error_unit, "(a,es24.16)") &
            "FortNum spherical product pin failed: ", product_coefficient
        error stop 1
    end if
    write (*, "(a)") "PASS"
end program test_fortnum_spherical_dependency
