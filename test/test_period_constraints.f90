program test_period_constraints
    use, intrinsic :: iso_fortran_env, only: real64
    use check, only: check_condition, check_summary
    use fortfem_period_constraints, only: &
        assemble_period_constraints, assemble_period_constraints_jvp, &
        assemble_period_constraints_vjp
    implicit none

    integer, parameter :: dp = real64
    integer, parameter :: edge_count = 4, period_count = 2
    real(dp) :: cycles(edge_count, period_count)
    real(dp) :: cycles_dot(edge_count, period_count)
    real(dp) :: cycles_bar(edge_count, period_count)
    complex(dp) :: edge_values(edge_count), edge_values_dot(edge_count)
    complex(dp) :: target(period_count), target_dot(period_count)
    complex(dp) :: residual(period_count), residual_dot(period_count)
    complex(dp) :: residual_bar(period_count), edge_values_bar(edge_count)
    complex(dp) :: target_bar(period_count), residual_plus(period_count)
    integer :: status
    real(dp) :: epsilon, finite_difference_error, lhs, rhs
    logical :: all_passed

    all_passed = .true.
    cycles = reshape([ &
        1.0_dp, -1.0_dp, 0.0_dp, 0.0_dp, &
        0.0_dp, 0.0_dp, 1.0_dp, -1.0_dp], shape(cycles))
    edge_values = [ &
        cmplx(0.7_dp, -0.2_dp, dp), cmplx(-0.1_dp, 0.4_dp, dp), &
        cmplx(0.3_dp, 0.6_dp, dp), cmplx(-0.5_dp, 0.1_dp, dp)]
    target = [cmplx(0.9_dp, -0.6_dp, dp), cmplx(0.8_dp, 0.5_dp, dp)]

    call assemble_period_constraints(cycles, edge_values, target, residual, status)
    call record_condition(status == 0, "cycle-period residual assembles")
    call record_condition(maxval(abs(residual - &
        (matmul(transpose(cycles), edge_values) - target))) < 1.0e-14_dp, &
        "cycle-period residual matches the matrix oracle")

    cycles_dot = reshape([ &
        0.03_dp, -0.02_dp, 0.04_dp, 0.01_dp, &
        -0.01_dp, 0.05_dp, 0.02_dp, -0.03_dp], shape(cycles_dot))
    edge_values_dot = [ &
        cmplx(0.04_dp, 0.02_dp, dp), cmplx(-0.01_dp, 0.03_dp, dp), &
        cmplx(0.02_dp, -0.05_dp, dp), cmplx(0.06_dp, 0.01_dp, dp)]
    target_dot = [cmplx(-0.03_dp, 0.02_dp, dp), cmplx(0.01_dp, -0.04_dp, dp)]
    call assemble_period_constraints_jvp( &
        cycles, edge_values, target, cycles_dot, edge_values_dot, target_dot, &
        residual_dot, status)
    call record_condition(status == 0, "cycle-period JVP assembles")
    epsilon = 1.0e-7_dp
    call assemble_period_constraints( &
        cycles + epsilon*cycles_dot, edge_values + epsilon*edge_values_dot, &
        target + epsilon*target_dot, residual_plus, status)
    finite_difference_error = maxval(abs(residual_dot - &
        (residual_plus - residual)/epsilon))
    call record_condition(finite_difference_error < 2.0e-8_dp, &
        "cycle-period JVP matches a forward difference")

    residual_bar = [ &
        cmplx(0.2_dp, -0.4_dp, dp), cmplx(-0.3_dp, 0.1_dp, dp)]
    call assemble_period_constraints_vjp( &
        cycles, edge_values, target, residual_bar, cycles_bar, edge_values_bar, &
        target_bar, status)
    call record_condition(status == 0, "cycle-period VJP assembles")
    lhs = real(sum(conjg(residual_bar)*residual_dot), dp)
    rhs = sum(cycles_bar*cycles_dot) + &
        real(sum(conjg(edge_values_bar)*edge_values_dot) + &
        sum(conjg(target_bar)*target_dot), dp)
    call record_condition(abs(lhs - rhs) < 2.0e-13_dp, &
        "cycle-period VJP satisfies the real complex adjoint identity")

    call assemble_period_constraints(cycles(:edge_count - 1, :), edge_values, &
        target, residual, status)
    call record_condition(status /= 0, "incompatible cycle dimensions are rejected")

    call check_summary("cycle-period constraints")
    if (.not. all_passed) error stop 1

contains

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_period_constraints
