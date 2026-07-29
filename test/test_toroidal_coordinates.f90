program test_toroidal_coordinates
    use check, only: check_condition, check_summary
    use fortfem_api, only: cartesian_to_toroidal, &
        toroidal_point_to_cartesian, toroidal_vector_to_cartesian
    use fortfem_kinds, only: dp
    implicit none

    real(dp), parameter :: scale = 1.7_dp
    real(dp), parameter :: eta_expected = 0.8_dp
    real(dp), parameter :: theta_expected = -1.1_dp
    real(dp), parameter :: phi_expected = 0.6_dp
    real(dp) :: cartesian(3), eta, phi, theta
    real(dp) :: cartesian_field(3), components(3), theta_field(3)
    logical :: all_passed

    all_passed = .true.
    call toroidal_point_to_cartesian( &
        scale, eta_expected, theta_expected, phi_expected, cartesian)
    call cartesian_to_toroidal(cartesian, scale, eta, theta, phi)
    call record_condition(abs(eta - eta_expected) < 1.0e-13_dp, &
        "Cartesian inverse recovers toroidal eta")
    call record_condition(abs(theta - theta_expected) < 1.0e-13_dp, &
        "Cartesian inverse recovers toroidal theta")
    call record_condition(abs(phi - phi_expected) < 1.0e-13_dp, &
        "Cartesian inverse recovers toroidal phi")

    components = [1.0_dp, 0.0_dp, 0.0_dp]
    call toroidal_vector_to_cartesian( &
        eta_expected, theta_expected, phi_expected, components, &
        cartesian_field)
    call record_condition(abs(norm2(cartesian_field) - 1.0_dp) < &
        1.0e-13_dp, "Toroidal eta unit vector remains normalized")
    components = [0.0_dp, 1.0_dp, 0.0_dp]
    call toroidal_vector_to_cartesian( &
        eta_expected, theta_expected, phi_expected, components, theta_field)
    call record_condition(abs(norm2(theta_field) - 1.0_dp) < 1.0e-13_dp, &
        "Toroidal theta unit vector remains normalized")
    call record_condition(abs(dot_product(cartesian_field, theta_field)) < &
        1.0e-13_dp, "Toroidal eta and theta unit vectors are orthogonal")

    call check_summary("Toroidal coordinates")
    if (.not. all_passed) error stop 1

contains

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(*), intent(in) :: description

        call check_condition(condition, description)
        all_passed = all_passed .and. condition
    end subroutine record_condition
end program test_toroidal_coordinates
