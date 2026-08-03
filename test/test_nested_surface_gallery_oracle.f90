program test_nested_surface_gallery_oracle
    use check, only: check_condition, check_summary
    use fortfem_core, only: &
        evaluate_axis_regular_radial_basis, evaluate_nested_surface_geometry
    use fortfem_fourier_mode_registry, only: &
        fourier_mode_registry_t, initialize_fourier_mode_registry
    use fortfem_kinds, only: dp
    use fortsparse, only: fortsparse_status_t
    implicit none

    real(dp), parameter :: pi = acos(-1.0_dp)
    real(dp), parameter :: major_radius = 2.4_dp
    real(dp), parameter :: minor_radius = 0.7_dp
    real(dp), parameter :: tolerance = 3.0e-12_dp
    integer, parameter :: sample_count = 6
    type(fourier_mode_registry_t) :: registry
    type(fortsparse_status_t) :: status
    complex(dp) :: coefficients(3, 2), radial_coefficients(2)
    complex(dp) :: radial_values(sample_count)
    real(dp) :: rho(sample_count), theta(sample_count), zeta(sample_count)
    real(dp) :: mapped(3, sample_count), physical(3, sample_count)
    real(dp) :: mapped_jacobian(3, 3, sample_count)
    real(dp) :: physical_jacobian(3, 3, sample_count)
    real(dp) :: metric(3, 3, sample_count), volume(sample_count)
    real(dp) :: expected_mapped(3, sample_count)
    real(dp) :: expected_physical(3, sample_count)
    real(dp) :: expected_radial(sample_count), radius
    integer :: sample
    logical :: all_passed

    all_passed = .true.
    call initialize_fourier_mode_registry( &
        registry, [0, 1], [0, 0], 1, 0.0_dp, 0.0_dp, .false., &
        radial_powers=[0, 1], status=status)
    call record_condition(status%code == 0, &
        "gallery Fourier/radial registry initializes")

    coefficients = cmplx(0.0_dp, 0.0_dp, dp)
    coefficients(1, 1) = cmplx(major_radius, 0.0_dp, dp)
    coefficients(1, 2) = cmplx(minor_radius, 0.0_dp, dp)
    coefficients(2, 2) = cmplx(0.0_dp, -minor_radius, dp)
    rho = [1.0_dp, 0.72_dp, 0.45_dp, 1.0_dp, 1.0_dp, 1.0_dp]
    theta = [0.0_dp, 0.5_dp*pi, pi, 1.5_dp*pi, 0.0_dp, 2.0_dp*pi]
    zeta = [0.0_dp, pi/3.0_dp, 0.5_dp*pi, pi, 0.0_dp, 2.0_dp*pi]

    call evaluate_nested_surface_geometry( &
        registry, coefficients, rho, theta, zeta, mapped, physical, &
        mapped_jacobian, physical_jacobian, metric, volume, status)
    call record_condition(status%code == 0, &
        "gallery circular-torus embedding evaluates")

    do sample = 1, sample_count
        radius = major_radius + minor_radius*rho(sample)*cos(theta(sample))
        expected_mapped(:, sample) = [ &
            radius, minor_radius*rho(sample)*sin(theta(sample)), 0.0_dp]
        expected_physical(:, sample) = [ &
            radius*cos(zeta(sample)), radius*sin(zeta(sample)), &
            minor_radius*rho(sample)*sin(theta(sample))]
    end do
    call record_condition( &
        maxval(abs(mapped - expected_mapped)) < tolerance, &
        "gallery mapped coordinates match the analytical circular torus")
    call record_condition( &
        maxval(abs(physical - expected_physical)) < tolerance, &
        "gallery Cartesian coordinates match an independent torus oracle")
    call record_condition( &
        maxval(abs(physical(:, 5) - physical(:, 6))) < tolerance, &
        "gallery surface closes at the doubly-periodic seam")
    call record_condition( &
        maxval(abs(volume + minor_radius**2*rho*expected_mapped(1, :))) < &
        tolerance, &
        "gallery volume Jacobian matches the analytical coordinate orientation")

    radial_coefficients = [ &
        cmplx(0.25_dp, 0.0_dp, dp), cmplx(0.75_dp, 0.0_dp, dp)]
    call evaluate_axis_regular_radial_basis( &
        0, [0, 2], radial_coefficients, rho, radial_values, status)
    expected_radial = 0.25_dp + 0.75_dp*rho**2
    call record_condition(status%code == 0, &
        "gallery axis-regular radial profile evaluates")
    call record_condition( &
        maxval(abs(real(radial_values, dp) - expected_radial)) < tolerance, &
        "gallery radial profile matches its independent polynomial oracle")

    call check_summary("nested-surface gallery oracle")
    if (.not. all_passed) error stop 1

contains

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: description

        call check_condition(condition, description)
        all_passed = all_passed .and. condition
    end subroutine record_condition

end program test_nested_surface_gallery_oracle
