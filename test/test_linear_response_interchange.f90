program test_linear_response_interchange
    use, intrinsic :: iso_fortran_env, only: real64
    use check, only: check_condition, check_summary
    use fortfem_interop, only: &
        evaluate_linear_response_diagnostics, &
        assemble_linear_response_operator, &
        assemble_linear_response_operator_jvp, &
        assemble_linear_response_operator_vjp, &
        assemble_linear_response_residual, &
        assemble_linear_response_residual_jvp, &
        assemble_linear_response_residual_vjp, &
        initialize_linear_response_interchange, &
        linear_response_interchange_t, validate_linear_response_interchange
    implicit none

    integer, parameter :: dp = real64
    integer, parameter :: n = 3, mode_count = 3, response_count = 2
    complex(dp), parameter :: frequency = cmplx(0.7_dp, -0.2_dp, dp)
    complex(dp) :: equilibrium(n, n), inertia(n, n), resistive(n, n)
    complex(dp) :: vacuum(n, n), wall(n, n), response(response_count, response_count)
    complex(dp) :: operator(n, n), operator_dot(n, n), operator_bar(n, n)
    complex(dp) :: state(n), state_dot(n), state_bar(n)
    complex(dp) :: source(n), source_dot(n), source_bar(n)
    complex(dp) :: residual(n), residual_dot(n), residual_bar(n)
    complex(dp) :: residual_plus(n), operator_bar_residual(n, n)
    complex(dp) :: equilibrium_dot(n, n), inertia_dot(n, n)
    complex(dp) :: resistive_dot(n, n), vacuum_dot(n, n), wall_dot(n, n)
    complex(dp) :: operator_plus(n, n), operator_minus(n, n)
    complex(dp), allocatable :: equilibrium_bar(:, :), inertia_bar(:, :)
    complex(dp), allocatable :: resistive_bar(:, :), vacuum_bar(:, :), wall_bar(:, :)
    complex(dp) :: frequency_dot, frequency_bar
    integer :: mode_numbers(2, mode_count), status
    character(len=32) :: mode_labels(mode_count), response_labels(response_count)
    type(linear_response_interchange_t) :: data, copy, invalid
    real(dp) :: epsilon, finite_difference_error, lhs, rhs
    real(dp) :: reciprocity_error, passivity_bound
    logical :: all_passed

    all_passed = .true.
    call build_blocks(equilibrium, inertia, resistive, vacuum, wall, response)
    mode_numbers = reshape([0, 1, 1, -1, 2, 0], shape(mode_numbers))
    mode_labels = [character(len=32) :: "m0_n1", "m1_nminus1", "m2_n0"]
    response_labels = [character(len=32) :: "vacuum_trace", "wall_trace"]

    call initialize_linear_response_interchange( &
        data, frequency, mode_numbers, mode_labels, equilibrium, inertia, &
        resistive, vacuum, wall, response, response_labels, 2.5_dp, status)
    call record_condition(status == 0, "linear response data initializes")
    call record_condition(validate_linear_response_interchange(data, status), &
        "linear response data validates")
    call record_condition(data%state_count == n .and. &
        data%mode_count == mode_count .and. &
        data%response_count == response_count, "response dimensions are retained")
    call evaluate_linear_response_diagnostics( &
        data%response_matrix, reciprocity_error, passivity_bound, status)
    call record_condition(status == 0 .and. reciprocity_error > 1.0e-2_dp, &
        "response reciprocity defect is reported")
    call record_condition(passivity_bound > 0.4_dp, &
        "response Hermitian passivity lower bound is positive")

    copy = data
    call record_condition(validate_linear_response_interchange(copy, status), &
        "linear response assignment preserves a valid deep copy")
    call record_condition(all(copy%equilibrium_block == data%equilibrium_block), &
        "linear response assignment copies blocks")

    call assemble_linear_response_operator( &
        data%equilibrium_block, data%inertia_block, data%resistive_block, &
        data%vacuum_block, data%wall_block, data%frequency, operator, status)
    call record_condition(status == 0, "linear response operator assembles")
    call record_condition(maxval(abs(operator - expected_operator( &
        equilibrium, inertia, resistive, vacuum, wall, frequency))) < 1.0e-14_dp, &
        "linear response block signs and frequency factors are correct")

    state = [ &
        cmplx(0.3_dp, -0.4_dp, dp), cmplx(-0.2_dp, 0.1_dp, dp), &
        cmplx(0.5_dp, 0.2_dp, dp)]
    source = [ &
        cmplx(0.7_dp, 0.3_dp, dp), cmplx(-0.1_dp, 0.6_dp, dp), &
        cmplx(0.2_dp, -0.5_dp, dp)]
    call assemble_linear_response_residual( &
        operator, state, source, residual, status)
    call record_condition(status == 0 .and. maxval(abs(residual - &
        (matmul(operator, state) - source))) < 1.0e-14_dp, &
        "linear response residual matches the matrix action")

    equilibrium_dot = cmplx(0.0_dp, 0.0_dp, dp)
    inertia_dot = cmplx(0.0_dp, 0.0_dp, dp)
    resistive_dot = cmplx(0.0_dp, 0.0_dp, dp)
    vacuum_dot = cmplx(0.0_dp, 0.0_dp, dp)
    wall_dot = cmplx(0.0_dp, 0.0_dp, dp)
    equilibrium_dot(1, 2) = cmplx(0.11_dp, -0.07_dp, dp)
    inertia_dot(2, 2) = cmplx(-0.05_dp, 0.03_dp, dp)
    resistive_dot(3, 1) = cmplx(0.04_dp, 0.09_dp, dp)
    vacuum_dot(1, 3) = cmplx(-0.02_dp, 0.08_dp, dp)
    wall_dot(2, 1) = cmplx(0.06_dp, -0.01_dp, dp)
    frequency_dot = cmplx(0.03_dp, -0.04_dp, dp)
    call assemble_linear_response_operator_jvp( &
        equilibrium, inertia, resistive, vacuum, wall, frequency, &
        equilibrium_dot, inertia_dot, resistive_dot, vacuum_dot, wall_dot, &
        frequency_dot, operator_dot, status)
    call record_condition(status == 0, "linear response operator JVP assembles")
    epsilon = 1.0e-7_dp
    call assemble_linear_response_operator( &
        equilibrium + epsilon*equilibrium_dot, inertia + epsilon*inertia_dot, &
        resistive + epsilon*resistive_dot, vacuum + epsilon*vacuum_dot, &
        wall + epsilon*wall_dot, frequency + epsilon*frequency_dot, &
        operator_plus, status)
    call assemble_linear_response_operator( &
        equilibrium - epsilon*equilibrium_dot, inertia - epsilon*inertia_dot, &
        resistive - epsilon*resistive_dot, vacuum - epsilon*vacuum_dot, &
        wall - epsilon*wall_dot, frequency - epsilon*frequency_dot, &
        operator_minus, status)
    finite_difference_error = maxval(abs(operator_dot - &
        (operator_plus - operator_minus)/(2.0_dp*epsilon)))
    call record_condition(finite_difference_error < 2.0e-8_dp, &
        "linear response operator JVP matches central differences")

    state_dot = [ &
        cmplx(-0.02_dp, 0.07_dp, dp), cmplx(0.06_dp, -0.01_dp, dp), &
        cmplx(0.03_dp, 0.04_dp, dp)]
    source_dot = [ &
        cmplx(0.01_dp, -0.03_dp, dp), cmplx(-0.05_dp, 0.02_dp, dp), &
        cmplx(0.04_dp, 0.01_dp, dp)]
    call assemble_linear_response_residual_jvp( &
        operator, state, source, operator_dot, state_dot, source_dot, &
        residual_dot, status)
    call record_condition(status == 0, "linear response residual JVP assembles")
    call assemble_linear_response_residual( &
        operator + epsilon*operator_dot, state + epsilon*state_dot, &
        source + epsilon*source_dot, residual_plus, status)
    finite_difference_error = maxval(abs(residual_dot - &
        (residual_plus - residual)/epsilon))
    call record_condition(finite_difference_error < 2.0e-8_dp, &
        "linear response residual JVP matches a forward difference")

    residual_bar = [ &
        cmplx(0.2_dp, -0.4_dp, dp), cmplx(-0.3_dp, 0.1_dp, dp), &
        cmplx(0.5_dp, 0.2_dp, dp)]
    call assemble_linear_response_residual_vjp( &
        operator, state, source, residual_bar, operator_bar_residual, &
        state_bar, source_bar, status)
    call record_condition(status == 0, "linear response residual VJP assembles")
    lhs = real(sum(conjg(residual_bar)*residual_dot), dp)
    rhs = real(sum(conjg(operator_bar_residual)*operator_dot) + &
        sum(conjg(state_bar)*state_dot) + sum(conjg(source_bar)*source_dot), dp)
    call record_condition(abs(lhs - rhs) < 2.0e-13_dp, &
        "linear response residual VJP satisfies the real complex adjoint identity")

    operator_bar = reshape([ &
        cmplx(0.2_dp, -0.1_dp, dp), cmplx(-0.3_dp, 0.4_dp, dp), &
        cmplx(0.5_dp, 0.2_dp, dp), cmplx(-0.1_dp, 0.3_dp, dp), &
        cmplx(0.7_dp, -0.2_dp, dp), cmplx(-0.4_dp, 0.6_dp, dp), &
        cmplx(0.11_dp, 0.31_dp, dp), cmplx(-0.21_dp, 0.17_dp, dp), &
        cmplx(0.09_dp, -0.13_dp, dp)], shape(operator_bar))
    call assemble_linear_response_operator_vjp( &
        inertia, resistive, frequency, operator_bar, equilibrium_bar, &
        inertia_bar, resistive_bar, vacuum_bar, wall_bar, frequency_bar, status)
    call record_condition(status == 0, "linear response operator VJP assembles")
    lhs = real(sum(conjg(operator_bar)*operator_dot), dp)
    rhs = real(sum(conjg(equilibrium_bar)*equilibrium_dot) + &
        sum(conjg(inertia_bar)*inertia_dot) + &
        sum(conjg(resistive_bar)*resistive_dot) + &
        sum(conjg(vacuum_bar)*vacuum_dot) + &
        sum(conjg(wall_bar)*wall_dot) + conjg(frequency_bar)*frequency_dot, dp)
    call record_condition(abs(lhs - rhs) < 2.0e-13_dp, &
        "linear response operator VJP satisfies the real complex adjoint identity")

    invalid = data
    invalid%mode_numbers(:, 2) = invalid%mode_numbers(:, 1)
    call record_condition(.not. validate_linear_response_interchange(invalid, status), &
        "duplicate Fourier modes are rejected")
    call check_summary("linear response interchange")
    if (.not. all_passed) error stop 1

contains

    subroutine build_blocks(equilibrium, inertia, resistive, vacuum, wall, response)
        complex(dp), intent(out) :: equilibrium(:, :), inertia(:, :)
        complex(dp), intent(out) :: resistive(:, :), vacuum(:, :), wall(:, :)
        complex(dp), intent(out) :: response(:, :)

        equilibrium = cmplx(0.0_dp, 0.0_dp, dp)
        inertia = cmplx(0.0_dp, 0.0_dp, dp)
        resistive = cmplx(0.0_dp, 0.0_dp, dp)
        vacuum = cmplx(0.0_dp, 0.0_dp, dp)
        wall = cmplx(0.0_dp, 0.0_dp, dp)
        equilibrium = reshape([ &
            cmplx(2.0_dp, 0.1_dp, dp), cmplx(0.2_dp, -0.3_dp, dp), cmplx(0.1_dp, 0.4_dp, dp), &
            cmplx(-0.1_dp, 0.2_dp, dp), cmplx(1.5_dp, -0.2_dp, dp), cmplx(0.3_dp, 0.1_dp, dp), &
            cmplx(0.4_dp, -0.1_dp, dp), cmplx(0.2_dp, 0.3_dp, dp), cmplx(1.2_dp, 0.0_dp, dp)], shape(equilibrium))
        inertia = diag([cmplx(1.0_dp, 0.0_dp, dp), cmplx(1.2_dp, 0.1_dp, dp), cmplx(0.8_dp, -0.2_dp, dp)])
        resistive = diag([cmplx(0.4_dp, 0.0_dp, dp), cmplx(0.2_dp, 0.1_dp, dp), cmplx(0.3_dp, -0.1_dp, dp)])
        vacuum = cmplx(0.0_dp, 0.0_dp, dp)
        vacuum(1, 1) = cmplx(0.3_dp, -0.2_dp, dp)
        vacuum(2, 3) = cmplx(-0.1_dp, 0.3_dp, dp)
        wall = cmplx(0.0_dp, 0.0_dp, dp)
        wall(2, 2) = cmplx(0.15_dp, 0.05_dp, dp)
        wall(3, 1) = cmplx(0.07_dp, -0.04_dp, dp)
        response = reshape([ &
            cmplx(1.0_dp, 0.0_dp, dp), cmplx(0.1_dp, -0.2_dp, dp), &
            cmplx(-0.3_dp, 0.4_dp, dp), cmplx(0.8_dp, 0.1_dp, dp)], shape(response))
    end subroutine build_blocks

    pure function diag(values) result(matrix)
        complex(dp), intent(in) :: values(:)
        complex(dp) :: matrix(size(values), size(values))
        integer :: index

        matrix = cmplx(0.0_dp, 0.0_dp, dp)
        do index = 1, size(values)
            matrix(index, index) = values(index)
        end do
    end function diag

    pure function expected_operator(equilibrium, inertia, resistive, vacuum, wall, frequency) result(matrix)
        complex(dp), intent(in) :: equilibrium(:, :), inertia(:, :), resistive(:, :)
        complex(dp), intent(in) :: vacuum(:, :), wall(:, :), frequency
        complex(dp) :: matrix(size(equilibrium, 1), size(equilibrium, 2))

        matrix = equilibrium - frequency**2*inertia + &
            cmplx(0.0_dp, 1.0_dp, dp)*frequency*resistive + vacuum + wall
    end function expected_operator

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_linear_response_interchange
