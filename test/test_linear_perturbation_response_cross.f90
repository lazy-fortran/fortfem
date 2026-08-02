program test_linear_perturbation_response_cross
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use check, only: check_condition, check_summary
    use fortfem_interop, only: &
        assemble_linear_perturbation_operator, &
        assemble_linear_perturbation_operator_jvp, &
        assemble_linear_response_residual, &
        assemble_linear_response_residual_jvp, &
        evaluate_linear_response_diagnostics
    implicit none

    integer, parameter :: n = 3
    complex(dp), parameter :: imaginary_unit = cmplx(0.0_dp, 1.0_dp, dp)
    complex(dp) :: lorentz(n, n), pressure(n, n), inertia(n, n)
    complex(dp) :: vacuum(n, n), wall(n, n), resistive(n, n), singular(n, n)
    complex(dp) :: operator(n, n), expected(n, n), operator_dot(n, n)
    complex(dp) :: expected_dot(n, n), state(n), state_dot(n), source(n)
    complex(dp) :: source_dot(n), residual(n), residual_dot(n)
    complex(dp) :: response(n, n), frequency, frequency_dot
    real(dp) :: reciprocity_error, passivity_margin
    integer :: column, row, status
    logical :: all_passed

    all_passed = .true.
    call build_blocks( &
        lorentz, pressure, inertia, vacuum, wall, resistive, singular)
    frequency = cmplx(0.72_dp, 0.0_dp, dp)
    frequency_dot = cmplx(0.13_dp, 0.0_dp, dp)
    state = [ &
        cmplx(0.8_dp, -0.2_dp, dp), cmplx(-0.3_dp, 0.5_dp, dp), &
        cmplx(0.25_dp, 0.15_dp, dp)]
    state_dot = [ &
        cmplx(-0.1_dp, 0.04_dp, dp), cmplx(0.03_dp, -0.08_dp, dp), &
        cmplx(0.06_dp, 0.02_dp, dp)]

    call assemble_linear_perturbation_operator( &
        lorentz, pressure, inertia, vacuum, wall, resistive, singular, &
        frequency, operator, status)
    call record_condition(status == 0, "seven-block operator composes")
    do column = 1, n
        do row = 1, n
            expected(row, column) = lorentz(row, column) + &
                pressure(row, column) + vacuum(row, column) + &
                wall(row, column) + singular(row, column) - &
                frequency*frequency*inertia(row, column) + &
                imaginary_unit*frequency*resistive(row, column)
        end do
    end do
    call record_condition(maxval(abs(operator - expected)) < 2.0e-14_dp, &
        "cross fixture matches an independent elementwise operator oracle")

    source = cmplx(0.0_dp, 0.0_dp, dp)
    do column = 1, n
        do row = 1, n
            source(row) = source(row) + expected(row, column)*state(column)
        end do
    end do
    call assemble_linear_response_residual( &
        operator, state, source, residual, status)
    call record_condition(status == 0 .and. maxval(abs(residual)) < 2.0e-14_dp, &
        "manufactured complex state satisfies the forced residual")

    call assemble_linear_perturbation_operator_jvp( &
        lorentz, pressure, inertia, vacuum, wall, resistive, singular, &
        frequency, 0.0_dp*lorentz, 0.0_dp*pressure, 0.0_dp*inertia, &
        0.0_dp*vacuum, 0.0_dp*wall, 0.0_dp*resistive, &
        0.0_dp*singular, frequency_dot, operator_dot, status)
    expected_dot = -2.0_dp*frequency*frequency_dot*inertia + &
        imaginary_unit*frequency_dot*resistive
    call record_condition(maxval(abs(operator_dot - expected_dot)) < 2.0e-14_dp, &
        "frequency JVP crosses the seven-block composition")
    source_dot = matmul(expected_dot, state) + matmul(expected, state_dot)
    call assemble_linear_response_residual_jvp( &
        operator, state, source, operator_dot, state_dot, source_dot, &
        residual_dot, status)
    call record_condition( &
        status == 0 .and. maxval(abs(residual_dot)) < 3.0e-14_dp, &
        "composed response JVP preserves the manufactured zero residual")

    response = cmplx(0.0_dp, 0.0_dp, dp)
    response(1, 1) = cmplx(1.2_dp, 0.10_dp, dp)
    response(2, 2) = cmplx(1.0_dp, 0.08_dp, dp)
    response(3, 3) = cmplx(0.9_dp, 0.06_dp, dp)
    response(1, 2) = cmplx(0.08_dp, 0.01_dp, dp)
    response(2, 1) = response(1, 2)
    response(2, 3) = cmplx(0.05_dp, 0.0_dp, dp)
    response(3, 2) = response(2, 3)
    call evaluate_linear_response_diagnostics( &
        response, reciprocity_error, passivity_margin, status)
    call record_condition(status == 0 .and. reciprocity_error == 0.0_dp, &
        "transpose-reciprocal response is identified")
    call record_condition(passivity_margin > 0.8_dp, &
        "passive response has a positive certified margin")

    call check_summary("linear perturbation response cross fixture")
    if (.not. all_passed) error stop 1

contains

    subroutine build_blocks( &
            lorentz, pressure, inertia, vacuum, wall, resistive, singular)
        complex(dp), intent(out) :: lorentz(:, :), pressure(:, :)
        complex(dp), intent(out) :: inertia(:, :), vacuum(:, :), wall(:, :)
        complex(dp), intent(out) :: resistive(:, :), singular(:, :)

        lorentz = reshape([ &
            3.2_dp, -0.25_dp, 0.08_dp, &
            -0.25_dp, 2.8_dp, -0.18_dp, &
            0.08_dp, -0.18_dp, 2.5_dp], shape(lorentz))
        pressure = reshape([ &
            0.55_dp, 0.04_dp, 0.0_dp, &
            0.04_dp, 0.45_dp, 0.03_dp, &
            0.0_dp, 0.03_dp, 0.35_dp], shape(pressure))
        inertia = cmplx(0.0_dp, 0.0_dp, dp)
        inertia(1, 1) = 1.0_dp
        inertia(2, 2) = 0.85_dp
        inertia(3, 3) = 0.7_dp
        vacuum = cmplx(0.0_dp, 0.0_dp, dp)
        vacuum(1, 1) = 0.22_dp
        vacuum(2, 2) = 0.18_dp
        vacuum(3, 3) = 0.15_dp
        wall = cmplx(0.0_dp, 0.0_dp, dp)
        wall(1, 1) = 0.09_dp
        wall(2, 2) = 0.07_dp
        wall(3, 3) = 0.06_dp
        resistive = cmplx(0.0_dp, 0.0_dp, dp)
        resistive(1, 1) = 0.32_dp
        resistive(2, 2) = 0.26_dp
        resistive(3, 3) = 0.21_dp
        singular = cmplx(0.0_dp, 0.0_dp, dp)
        singular(1, 1) = 0.11_dp
        singular(2, 2) = 0.08_dp
        singular(3, 3) = 0.05_dp
    end subroutine build_blocks

    subroutine record_condition(condition, message)
        logical, intent(in) :: condition
        character(*), intent(in) :: message

        call check_condition(condition, message)
        all_passed = all_passed .and. condition
    end subroutine record_condition

end program test_linear_perturbation_response_cross
