program test_fourier_mode_energy
    use check, only: check_condition, check_summary
    use fortfem_fourier, only: &
        assemble_fourier_mode_energy, assemble_fourier_mode_energy_jvp, &
        assemble_fourier_mode_energy_vjp, &
        fourier_coefficients_conjugate_symmetric, &
        fourier_mode_registry_t, initialize_fourier_mode_registry
    use fortfem_kinds, only: dp
    use fortsparse, only: fortsparse_status_t
    implicit none

    integer, parameter :: component_count = 2, point_count = 2, mode_count = 3
    integer, parameter :: poloidal(mode_count) = [-1, 0, 1]
    integer, parameter :: toroidal(mode_count) = [1, 0, -1]
    integer, parameter :: radial_power(mode_count) = [1, 0, 1]
    real(dp), parameter :: normalization(mode_count) = [1.1_dp, 0.8_dp, 1.1_dp]
    real(dp), parameter :: point_weights(point_count) = [0.7_dp, 1.2_dp]
    real(dp), parameter :: mode_weights(mode_count) = [1.3_dp, 0.9_dp, 1.3_dp]
    real(dp), parameter :: point_weights_dot(point_count) = [0.2_dp, -0.15_dp]
    real(dp), parameter :: mode_weights_dot(mode_count) = [-0.1_dp, 0.08_dp, 0.12_dp]
    real(dp), parameter :: eps = 1.0e-7_dp
    type(fourier_mode_registry_t) :: registry
    type(fortsparse_status_t) :: status
    complex(dp) :: coefficients(component_count, point_count, mode_count)
    complex(dp) :: coefficients_dot(component_count, point_count, mode_count)
    complex(dp) :: coefficients_bar(component_count, point_count, mode_count)
    real(dp) :: mode_energy(mode_count), mode_energy_dot(mode_count)
    real(dp) :: mode_energy_plus(mode_count), mode_energy_minus(mode_count)
    real(dp) :: mode_energy_bar(mode_count)
    real(dp) :: total_energy, total_energy_dot, total_energy_plus, total_energy_minus
    real(dp) :: total_energy_bar, point_weights_bar(point_count)
    real(dp) :: mode_weights_bar(mode_count), oracle(mode_count), oracle_dot(mode_count)
    real(dp) :: residual, lhs, rhs
    logical :: symmetric
    integer :: component, point, mode

    call initialize_fourier_mode_registry( &
        registry, poloidal, toroidal, 2, 0.15_dp, -0.25_dp, .true., &
        radial_power, normalization, status)
    if (status%code /= 0) error stop 1
    do mode = 1, mode_count
        do point = 1, point_count
            do component = 1, component_count
                coefficients(component, point, mode) = cmplx( &
                    0.15_dp*real(component + 2*point + mode, dp), &
                    -0.04_dp*real(component + point + mode, dp), dp)
                coefficients_dot(component, point, mode) = cmplx( &
                    -0.03_dp*real(component + point, dp), &
                    0.02_dp*real(mode + 2*point, dp), dp)
            end do
        end do
    end do
    coefficients(:, :, 2) = cmplx(real(coefficients(:, :, 2), dp), 0.0_dp, dp)
    coefficients(:, :, 3) = conjg(coefficients(:, :, 1))

    symmetric = fourier_coefficients_conjugate_symmetric( &
        registry, coefficients, 1.0e-14_dp, residual, status)
    call check_condition(status%code == 0 .and. symmetric .and. residual < 1.0e-14_dp, &
        "Fourier modal coefficients satisfy real-packed conjugacy")
    coefficients(1, 1, 3) = coefficients(1, 1, 3) + cmplx(0.1_dp, 0.0_dp, dp)
    symmetric = fourier_coefficients_conjugate_symmetric( &
        registry, coefficients, 1.0e-14_dp, residual, status)
    call check_condition(status%code == 0 .and. .not. symmetric .and. residual > 0.05_dp, &
        "Fourier modal symmetry diagnostic reports a broken pair")
    coefficients(:, :, 3) = conjg(coefficients(:, :, 1))

    do mode = 1, mode_count
        oracle(mode) = 0.0_dp
        oracle_dot(mode) = 0.0_dp
        do point = 1, point_count
            do component = 1, component_count
                oracle(mode) = oracle(mode) + 0.5_dp*mode_weights(mode)* &
                    point_weights(point)*abs(coefficients(component, point, mode))**2
                oracle_dot(mode) = oracle_dot(mode) + 0.5_dp*point_weights(point)* &
                    mode_weights_dot(mode)*abs(coefficients(component, point, mode))**2 + &
                    0.5_dp*mode_weights(mode)*point_weights_dot(point)* &
                    abs(coefficients(component, point, mode))**2 + &
                    mode_weights(mode)*point_weights(point)*real(conjg( &
                    coefficients(component, point, mode))*coefficients_dot( &
                    component, point, mode), dp)
            end do
        end do
    end do
    call assemble_fourier_mode_energy( &
        registry, coefficients, point_weights, mode_weights, mode_energy, &
        total_energy, status)
    call check_condition(status%code == 0 .and. &
        maxval(abs(mode_energy - oracle)) < 1.0e-14_dp .and. &
        abs(total_energy - sum(oracle)) < 1.0e-14_dp, &
        "weighted Fourier modal energy matches the independent oracle")
    call assemble_fourier_mode_energy_jvp( &
        registry, coefficients, coefficients_dot, point_weights, point_weights_dot, &
        mode_weights, mode_weights_dot, mode_energy_dot, total_energy_dot, status)
    call assemble_fourier_mode_energy( &
        registry, coefficients + eps*coefficients_dot, &
        point_weights + eps*point_weights_dot, mode_weights + eps*mode_weights_dot, &
        mode_energy_plus, total_energy_plus, status)
    call assemble_fourier_mode_energy( &
        registry, coefficients - eps*coefficients_dot, &
        point_weights - eps*point_weights_dot, mode_weights - eps*mode_weights_dot, &
        mode_energy_minus, total_energy_minus, status)
    call check_condition(status%code == 0 .and. &
        maxval(abs(mode_energy_dot - oracle_dot)) < 1.0e-14_dp .and. &
        maxval(abs((mode_energy_plus - mode_energy_minus)/(2.0_dp*eps) - &
        mode_energy_dot)) < 1.0e-8_dp .and. &
        abs((total_energy_plus - total_energy_minus)/(2.0_dp*eps) - &
        total_energy_dot) < 1.0e-8_dp, &
        "Fourier modal energy JVP matches the product rule and differences")

    mode_energy_bar = [0.4_dp, -0.2_dp, 0.7_dp]
    total_energy_bar = -0.35_dp
    call assemble_fourier_mode_energy_vjp( &
        registry, coefficients, point_weights, mode_weights, mode_energy_bar, &
        total_energy_bar, coefficients_bar, point_weights_bar, mode_weights_bar, &
        status)
    lhs = dot_product(mode_energy_bar, mode_energy_dot) + &
        total_energy_bar*total_energy_dot
    rhs = real(sum(conjg(coefficients_bar)*coefficients_dot), dp) + &
        dot_product(point_weights_bar, point_weights_dot) + &
        dot_product(mode_weights_bar, mode_weights_dot)
    call check_condition(status%code == 0 .and. abs(lhs - rhs) < 1.0e-13_dp, &
        "Fourier modal energy VJP satisfies the real adjoint identity")
    call check_summary("Fourier modal energy and symmetry")
end program test_fourier_mode_energy
