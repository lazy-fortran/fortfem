program test_fourier_mode_linear_operator
    use, intrinsic :: ieee_arithmetic, only: ieee_quiet_nan, ieee_value
    use check, only: check_condition, check_summary
    use fortfem_fourier, only: &
        apply_fourier_mode_linear_operator, &
        apply_fourier_mode_linear_operator_jvp, &
        apply_fourier_mode_linear_operator_vjp
    use fortfem_kinds, only: dp
    use fortsparse, only: fortsparse_status_t
    implicit none

    integer, parameter :: output_count = 3, input_count = 2, mode_count = 4
    real(dp), parameter :: eps = 1.0e-7_dp
    complex(dp), allocatable :: operator(:, :, :), operator_dot(:, :, :)
    complex(dp), allocatable :: coefficients(:, :), coefficients_dot(:, :)
    complex(dp), allocatable :: output(:, :), output_dot(:, :)
    complex(dp), allocatable :: output_plus(:, :), output_minus(:, :)
    complex(dp), allocatable :: output_bar(:, :)
    complex(dp), allocatable :: operator_bar(:, :, :), coefficients_bar(:, :)
    complex(dp), allocatable :: oracle(:, :), oracle_dot(:, :)
    complex(dp), allocatable :: bad_output(:, :), nonfinite_operator(:, :, :)
    type(fortsparse_status_t) :: status
    integer :: input_component, output_component, mode
    real(dp) :: lhs, rhs

    allocate(operator(output_count, input_count, mode_count), &
        operator_dot(output_count, input_count, mode_count), &
        coefficients(input_count, mode_count), &
        coefficients_dot(input_count, mode_count), &
        output(output_count, mode_count), output_dot(output_count, mode_count), &
        output_plus(output_count, mode_count), output_minus(output_count, mode_count), &
        output_bar(output_count, mode_count), &
        operator_bar(output_count, input_count, mode_count), &
        coefficients_bar(input_count, mode_count), &
        oracle(output_count, mode_count), oracle_dot(output_count, mode_count), &
        bad_output(1, mode_count), &
        nonfinite_operator(output_count, input_count, mode_count))

    do mode = 1, mode_count
        do input_component = 1, input_count
            coefficients(input_component, mode) = cmplx( &
                0.17_dp*real(input_component + mode, dp), &
                -0.09_dp*real(2*input_component + mode, dp), dp)
            coefficients_dot(input_component, mode) = cmplx( &
                -0.03_dp*real(mode, dp), &
                0.04_dp*real(input_component, dp), dp)
            do output_component = 1, output_count
                operator(output_component, input_component, mode) = cmplx( &
                    0.08_dp*real(output_component + input_component + mode, dp), &
                    -0.02_dp*real(output_component + 2*mode, dp), dp)
                operator_dot(output_component, input_component, mode) = cmplx( &
                    -0.01_dp*real(output_component + mode, dp), &
                    0.025_dp*real(input_component, dp), dp)
            end do
        end do
    end do

    oracle = cmplx(0.0_dp, 0.0_dp, dp)
    oracle_dot = cmplx(0.0_dp, 0.0_dp, dp)
    do mode = 1, mode_count
        do output_component = 1, output_count
            do input_component = 1, input_count
                oracle(output_component, mode) = &
                    oracle(output_component, mode) + &
                    operator(output_component, input_component, mode)* &
                    coefficients(input_component, mode)
                oracle_dot(output_component, mode) = &
                    oracle_dot(output_component, mode) + &
                    operator_dot(output_component, input_component, mode)* &
                    coefficients(input_component, mode) + &
                    operator(output_component, input_component, mode)* &
                    coefficients_dot(input_component, mode)
            end do
        end do
    end do

    call apply_fourier_mode_linear_operator( &
        operator, coefficients, output, status)
    call check_condition(status%code == 0 .and. &
        maxval(abs(output - oracle)) < 1.0e-13_dp, &
        "mode-diagonal Fourier operator matches independent oracle")

    call apply_fourier_mode_linear_operator_jvp( &
        operator, coefficients, operator_dot, coefficients_dot, output_dot, status)
    call apply_fourier_mode_linear_operator( &
        operator + eps*operator_dot, coefficients + eps*coefficients_dot, &
        output_plus, status)
    call apply_fourier_mode_linear_operator( &
        operator - eps*operator_dot, coefficients - eps*coefficients_dot, &
        output_minus, status)
    call check_condition(status%code == 0 .and. &
        maxval(abs(output_dot - oracle_dot)) < 1.0e-13_dp .and. &
        maxval(abs((output_plus - output_minus)/(2.0_dp*eps) - output_dot)) &
        < 1.0e-8_dp, "mode-diagonal Fourier operator JVP matches product rule")

    do mode = 1, mode_count
        do output_component = 1, output_count
            output_bar(output_component, mode) = cmplx( &
                0.05_dp*real(output_component + mode, dp), &
                -0.03_dp*real(2*output_component + mode, dp), dp)
        end do
    end do
    call apply_fourier_mode_linear_operator_vjp( &
        operator, coefficients, output_bar, operator_bar, coefficients_bar, status)
    lhs = real(sum(conjg(output_bar)*oracle_dot), dp)
    rhs = real(sum(conjg(operator_bar)*operator_dot) + &
        sum(conjg(coefficients_bar)*coefficients_dot), dp)
    call check_condition(status%code == 0 .and. abs(lhs - rhs) < 1.0e-12_dp, &
        "mode-diagonal Fourier operator VJP satisfies the complex dot identity")

    call apply_fourier_mode_linear_operator( &
        operator, coefficients, bad_output, status)
    call check_condition(status%code /= 0, &
        "mode-diagonal Fourier operator rejects an incompatible output shape")

    nonfinite_operator = operator
    nonfinite_operator(1, 1, 1) = cmplx( &
        ieee_value(0.0_dp, ieee_quiet_nan), 0.0_dp, dp)
    call apply_fourier_mode_linear_operator( &
        nonfinite_operator, coefficients, output, status)
    call check_condition(status%code /= 0, &
        "mode-diagonal Fourier operator rejects nonfinite input")

    call check_summary("mode-diagonal Fourier linear operator")
end program test_fourier_mode_linear_operator
