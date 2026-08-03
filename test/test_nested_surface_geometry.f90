program test_nested_surface_geometry
    use check, only: check_condition, check_summary
    use fortfem_core, only: &
        evaluate_nested_surface_geometry, &
        evaluate_nested_surface_geometry_jvp, &
        evaluate_nested_surface_geometry_vjp, &
        evaluate_nested_surface_geometry_coordinate_jvp, &
        evaluate_nested_surface_geometry_coordinate_vjp
    use fortfem_fourier, only: &
        fourier_mode_registry_t, initialize_fourier_mode_registry
    use fortfem_kinds, only: dp
    use fortsparse, only: fortsparse_status_t
    implicit none

    real(dp), parameter :: tolerance = 2.0e-11_dp
    real(dp), parameter :: derivative_tolerance = 2.0e-8_dp
    real(dp), parameter :: step = 1.0e-6_dp
    integer, parameter :: sample_count = 2
    integer, parameter :: mode_count = 3
    type(fourier_mode_registry_t) :: registry
    type(fortsparse_status_t) :: status
    complex(dp) :: coefficients(3, mode_count), coefficients_dot(3, mode_count)
    complex(dp) :: coefficients_bar(3, mode_count)
    real(dp) :: rho(sample_count), theta(sample_count), zeta(sample_count)
    real(dp) :: mapped(3, sample_count), physical(3, sample_count)
    real(dp) :: mapped_jacobian(3, 3, sample_count)
    real(dp) :: physical_jacobian(3, 3, sample_count)
    real(dp) :: metric(3, 3, sample_count), volume(sample_count)
    real(dp) :: mapped_dot(3, sample_count), physical_dot(3, sample_count)
    real(dp) :: mapped_jacobian_dot(3, 3, sample_count)
    real(dp) :: physical_jacobian_dot(3, 3, sample_count)
    real(dp) :: metric_dot(3, 3, sample_count), volume_dot(sample_count)
    real(dp) :: rho_dot(sample_count), theta_dot(sample_count), zeta_dot(sample_count)
    real(dp) :: coordinate_mapped_dot(3, sample_count)
    real(dp) :: coordinate_physical_dot(3, sample_count)
    real(dp) :: coordinate_rho_bar(sample_count)
    real(dp) :: coordinate_theta_bar(sample_count)
    real(dp) :: coordinate_zeta_bar(sample_count)
    real(dp) :: mapped_plus(3, sample_count), physical_plus(3, sample_count)
    real(dp) :: mapped_minus(3, sample_count), physical_minus(3, sample_count)
    real(dp) :: mapped_jacobian_plus(3, 3, sample_count)
    real(dp) :: mapped_jacobian_minus(3, 3, sample_count)
    real(dp) :: physical_jacobian_plus(3, 3, sample_count)
    real(dp) :: physical_jacobian_minus(3, 3, sample_count)
    real(dp) :: metric_plus(3, 3, sample_count), metric_minus(3, 3, sample_count)
    real(dp) :: volume_plus(sample_count), volume_minus(sample_count)
    real(dp) :: mapped_bar(3, sample_count), physical_bar(3, sample_count)
    real(dp) :: mapped_jacobian_bar(3, 3, sample_count)
    real(dp) :: physical_jacobian_bar(3, 3, sample_count)
    real(dp) :: metric_bar(3, 3, sample_count), volume_bar(sample_count)
    real(dp) :: expected_mapped(3), basis(3), forward_pairing, reverse_pairing
    real(dp) :: bad_rho(1), bad_mapped(3, 1), bad_physical(3, 1)
    real(dp) :: bad_mapped_jacobian(3, 3, 1), bad_physical_jacobian(3, 3, 1)
    real(dp) :: bad_metric(3, 3, 1), bad_volume(1)
    real(dp) :: bad_coordinate_dot(1), bad_coordinate_mapped_dot(3, 1)
    real(dp) :: bad_coordinate_physical_dot(3, 1)
    real(dp) :: bad_coordinate_bar(3, 1)
    integer :: status_code, mode
    logical :: all_passed

    all_passed = .true.
    call initialize_fourier_mode_registry( &
        registry, [0, 1, 0], [0, 0, 1], 2, 0.0_dp, 0.0_dp, .false., &
        radial_powers=[0, 1, 0], status=status)
    call record_condition(status%code == 0, &
        "nested-surface test registry initializes")

    coefficients = cmplx(0.0_dp, 0.0_dp, dp)
    coefficients(:, 1) = [ &
        cmplx(2.5_dp, 0.0_dp, dp), cmplx(0.3_dp, 0.0_dp, dp), &
        cmplx(0.0_dp, 0.0_dp, dp)]
    coefficients(:, 2) = [ &
        cmplx(0.4_dp, -0.1_dp, dp), cmplx(0.2_dp, 0.05_dp, dp), &
        cmplx(0.1_dp, 0.0_dp, dp)]
    coefficients(:, 3) = [ &
        cmplx(0.05_dp, 0.02_dp, dp), cmplx(0.0_dp, 0.0_dp, dp), &
        cmplx(0.03_dp, -0.01_dp, dp)]
    rho = [0.0_dp, 0.4_dp]
    theta = [0.7_dp, 0.4_dp]
    zeta = [0.3_dp, 0.8_dp]

    call evaluate_nested_surface_geometry( &
        registry, coefficients, rho, theta, zeta, mapped, physical, &
        mapped_jacobian, physical_jacobian, metric, volume, status)
    call record_condition(status%code == 0, &
        "nested-surface evaluator accepts finite modal data")
    call evaluate_mode_oracle( &
        coefficients, rho(2), theta(2), zeta(2), basis, expected_mapped)
    call record_condition(maxval(abs(mapped(:, 2) - expected_mapped)) < tolerance, &
        "inverse R/Z/lambda coordinates match independent Fourier oracle")
    call record_condition(maxval(abs(metric - transpose_metric(metric))) < tolerance, &
        "physical metric is symmetric")
    call record_condition(all(volume == volume), &
        "volume Jacobian diagnostic is finite")

    call evaluate_nested_surface_geometry( &
        registry, coefficients, rho, theta, zeta + [0.0_dp, acos(-1.0_dp)], &
        mapped_plus, physical_plus, mapped_jacobian_plus, physical_jacobian_plus, &
        metric_plus, volume_plus, status)
    call record_condition(status%code == 0, &
        "field-period shifted samples are accepted")
    call record_condition(maxval(abs(mapped_plus(:, 2) - mapped(:, 2))) < tolerance, &
        "mapped coordinates respect the two-field-period zeta period")
    call record_condition( &
        maxval(abs(physical_plus(1:2, 2) + physical(1:2, 2))) < tolerance .and. &
        abs(physical_plus(3, 2) - physical(3, 2)) < tolerance, &
        "Cartesian coordinates rotate by the two-field-period zeta period")
    call record_condition(maxval(abs(mapped(:, 1) - mapped_at_theta( &
        registry, coefficients, rho(1), theta(1) + 1.1_dp, zeta(1)))) < tolerance, &
        "axis radial powers remove theta dependence at rho=0")

    coefficients_dot = reshape([ &
        cmplx(0.01_dp, -0.02_dp, dp), cmplx(-0.03_dp, 0.01_dp, dp), &
        cmplx(0.02_dp, 0.04_dp, dp), cmplx(-0.01_dp, 0.03_dp, dp), &
        cmplx(0.05_dp, -0.02_dp, dp), cmplx(-0.04_dp, -0.01_dp, dp), &
        cmplx(0.03_dp, 0.02_dp, dp), cmplx(0.02_dp, -0.05_dp, dp), &
        cmplx(-0.02_dp, 0.03_dp, dp)], [3, mode_count])
    call evaluate_nested_surface_geometry_jvp( &
        registry, coefficients, rho, theta, zeta, coefficients_dot, &
        mapped_dot, physical_dot, mapped_jacobian_dot, physical_jacobian_dot, &
        metric_dot, volume_dot, status)
    call record_condition(status%code == 0, &
        "nested-surface coefficient JVP succeeds")
    call evaluate_nested_surface_geometry( &
        registry, coefficients + step*coefficients_dot, rho, theta, zeta, &
        mapped_plus, physical_plus, mapped_jacobian_plus, physical_jacobian_plus, &
        metric_plus, volume_plus, status)
    call evaluate_nested_surface_geometry( &
        registry, coefficients - step*coefficients_dot, rho, theta, zeta, &
        mapped_minus, physical_minus, mapped_jacobian_minus, physical_jacobian_minus, &
        metric_minus, volume_minus, status)
    call record_condition(maxval(abs(mapped_dot - &
        (mapped_plus - mapped_minus)/(2.0_dp*step))) < derivative_tolerance, &
        "inverse-coordinate JVP matches central differences")
    call record_condition(maxval(abs(physical_dot - &
        (physical_plus - physical_minus)/(2.0_dp*step))) < derivative_tolerance, &
        "Cartesian-coordinate JVP matches central differences")
    call record_condition(maxval(abs(metric_dot - &
        (metric_plus - metric_minus)/(2.0_dp*step))) < derivative_tolerance, &
        "metric JVP matches central differences")
    call record_condition(maxval(abs(volume_dot - &
        (volume_plus - volume_minus)/(2.0_dp*step))) < derivative_tolerance, &
        "volume JVP matches central differences")

    mapped_bar = reshape([(0.01_dp*real(mode, dp), mode=1, 6)], [3, sample_count])
    physical_bar = reshape([(0.02_dp*real(mode, dp), mode=1, 6)], [3, sample_count])
    mapped_jacobian_bar = reshape([(0.003_dp*real(mode, dp), mode=1, 18)], &
        [3, 3, sample_count])
    physical_jacobian_bar = reshape([(0.004_dp*real(mode, dp), mode=1, 18)], &
        [3, 3, sample_count])
    metric_bar = reshape([(0.005_dp*real(mode, dp), mode=1, 18)], &
        [3, 3, sample_count])
    volume_bar = [0.07_dp, -0.08_dp]
    call evaluate_nested_surface_geometry_vjp( &
        registry, coefficients, rho, theta, zeta, mapped_bar, physical_bar, &
        mapped_jacobian_bar, physical_jacobian_bar, metric_bar, volume_bar, &
        coefficients_bar, status)
    call record_condition(status%code == 0, &
        "nested-surface coefficient VJP succeeds")
    forward_pairing = sum(mapped_bar*mapped_dot) + sum(physical_bar*physical_dot) + &
        sum(mapped_jacobian_bar*mapped_jacobian_dot) + &
        sum(physical_jacobian_bar*physical_jacobian_dot) + &
        sum(metric_bar*metric_dot) + sum(volume_bar*volume_dot)
    reverse_pairing = real(sum(conjg(coefficients_bar)*coefficients_dot), dp)
    call record_condition( &
        abs(forward_pairing - reverse_pairing) < derivative_tolerance, &
        "nested-surface coefficient VJP satisfies the dot-product oracle")

    rho_dot = [0.0_dp, 0.08_dp]
    theta_dot = [0.02_dp, -0.03_dp]
    zeta_dot = [-0.04_dp, 0.05_dp]
    call evaluate_nested_surface_geometry_coordinate_jvp( &
        registry, coefficients, rho, theta, zeta, rho_dot, theta_dot, zeta_dot, &
        coordinate_mapped_dot, coordinate_physical_dot, status)
    call record_condition(status%code == 0, &
        "nested-surface coordinate JVP succeeds")
    call evaluate_nested_surface_geometry( &
        registry, coefficients, rho + step*rho_dot, theta + step*theta_dot, &
        zeta + step*zeta_dot, mapped_plus, physical_plus, mapped_jacobian_plus, &
        physical_jacobian_plus, metric_plus, volume_plus, status)
    call evaluate_nested_surface_geometry( &
        registry, coefficients, rho - step*rho_dot, theta - step*theta_dot, &
        zeta - step*zeta_dot, mapped_minus, physical_minus, mapped_jacobian_minus, &
        physical_jacobian_minus, metric_minus, volume_minus, status)
    call record_condition(maxval(abs(coordinate_mapped_dot - &
        (mapped_plus - mapped_minus)/(2.0_dp*step))) < derivative_tolerance, &
        "mapped-coordinate JVP matches central differences")
    call record_condition(maxval(abs(coordinate_physical_dot - &
        (physical_plus - physical_minus)/(2.0_dp*step))) < derivative_tolerance, &
        "Cartesian-coordinate JVP matches central differences")

    mapped_bar = reshape([(0.013_dp*real(mode, dp), mode=1, 6)], &
        [3, sample_count])
    physical_bar = reshape([(0.017_dp*real(mode, dp), mode=1, 6)], &
        [3, sample_count])
    call evaluate_nested_surface_geometry_coordinate_vjp( &
        registry, coefficients, rho, theta, zeta, mapped_bar, physical_bar, &
        coordinate_rho_bar, coordinate_theta_bar, coordinate_zeta_bar, status)
    call record_condition(status%code == 0, &
        "nested-surface coordinate VJP succeeds")
    forward_pairing = sum(mapped_bar*coordinate_mapped_dot) + &
        sum(physical_bar*coordinate_physical_dot)
    reverse_pairing = sum(coordinate_rho_bar*rho_dot) + &
        sum(coordinate_theta_bar*theta_dot) + sum(coordinate_zeta_bar*zeta_dot)
    call record_condition( &
        abs(forward_pairing - reverse_pairing) < derivative_tolerance, &
        "nested-surface coordinate VJP satisfies the dot-product oracle")

    bad_coordinate_dot = [0.1_dp]
    call evaluate_nested_surface_geometry_coordinate_jvp( &
        registry, coefficients, rho, theta, zeta, bad_coordinate_dot, &
        theta_dot, zeta_dot, bad_coordinate_mapped_dot, &
        bad_coordinate_physical_dot, status)
    call record_condition(status%code /= 0, &
        "nested-surface coordinate JVP rejects an axis-length mismatch")
    call evaluate_nested_surface_geometry_coordinate_vjp( &
        registry, coefficients, rho, theta, zeta, bad_coordinate_bar, &
        physical_bar, coordinate_rho_bar, coordinate_theta_bar, &
        coordinate_zeta_bar, status)
    call record_condition(status%code /= 0, &
        "nested-surface coordinate VJP rejects a sample-shape mismatch")

    bad_rho = [-0.1_dp]
    call evaluate_nested_surface_geometry( &
        registry, coefficients, bad_rho, theta(:1), zeta(:1), bad_mapped, &
        bad_physical, bad_mapped_jacobian, bad_physical_jacobian, bad_metric, &
        bad_volume, status)
    status_code = status%code
    call record_condition(status_code /= 0, &
        "nested-surface evaluator rejects negative radial coordinates")

    call check_summary("nested-surface geometry")
    if (.not. all_passed) error stop 1

contains

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: description

        call check_condition(condition, description)
        all_passed = all_passed .and. condition
    end subroutine record_condition

    subroutine evaluate_mode_oracle(coefficients, radius, poloidal, toroidal, &
            basis, fields)
        complex(dp), intent(in) :: coefficients(:, :)
        real(dp), intent(in) :: radius, poloidal, toroidal
        real(dp), intent(out) :: basis(3), fields(3)
        complex(dp) :: mode_basis(3)

        mode_basis(1) = cmplx(1.0_dp, 0.0_dp, dp)
        mode_basis(2) = radius*exp(cmplx(0.0_dp, poloidal, dp))
        mode_basis(3) = exp(cmplx(0.0_dp, 2.0_dp*toroidal, dp))
        basis = real(mode_basis, dp)
        fields = [ &
            real(sum(coefficients(1, :)*mode_basis), dp), &
            real(sum(coefficients(2, :)*mode_basis), dp), &
            real(sum(coefficients(3, :)*mode_basis), dp)]
    end subroutine evaluate_mode_oracle

    function mapped_at_theta(registry, coefficients, radius, poloidal, toroidal) &
            result(value)
        type(fourier_mode_registry_t), intent(in) :: registry
        complex(dp), intent(in) :: coefficients(:, :)
        real(dp), intent(in) :: radius, poloidal, toroidal
        real(dp) :: value(3), mapped_samples(3, 1), dummy_physical(3, 1)
        real(dp) :: jacobian(3, 3, 1), dummy_jacobian(3, 3, 1)
        real(dp) :: dummy_metric(3, 3, 1), dummy_volume(1)

        call evaluate_nested_surface_geometry( &
            registry, coefficients, [radius], [poloidal], [toroidal], mapped_samples, &
            dummy_physical, jacobian, dummy_jacobian, dummy_metric, &
            dummy_volume, status)
        value = mapped_samples(:, 1)
    end function mapped_at_theta

    pure function transpose_metric(values) result(transposed)
        real(dp), intent(in) :: values(:, :, :)
        real(dp) :: transposed(size(values, 2), size(values, 1), size(values, 3))
        integer :: sample

        do sample = 1, size(values, 3)
            transposed(:, :, sample) = transpose(values(:, :, sample))
        end do
    end function transpose_metric

end program test_nested_surface_geometry
