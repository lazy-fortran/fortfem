program test_fitted_trace_constraint
    use check, only: check_condition, check_summary
    use fortfem_api, only: assemble_fitted_trace_constraint, &
        assemble_fitted_trace_constraint_jvp, assemble_fitted_trace_constraint_vjp
    use fortfem_kinds, only: dp
    use fortsparse, only: fortsparse_status_t
    implicit none

    integer, parameter :: quadrature_count = 3, multiplier_count = 2
    integer, parameter :: plus_count = 2, minus_count = 1
    real(dp), parameter :: epsilon_fd = 1.0e-7_dp
    real(dp), parameter :: multiplier_trace(quadrature_count, multiplier_count) = &
        reshape([0.5_dp, -0.2_dp, 0.7_dp, 0.3_dp, -0.4_dp, 0.9_dp], &
        [quadrature_count, multiplier_count])
    real(dp), parameter :: plus_trace(quadrature_count, plus_count) = reshape([ &
        1.0_dp, 0.2_dp, -0.3_dp, 0.4_dp, -0.5_dp, 0.6_dp], &
        [quadrature_count, plus_count])
    real(dp), parameter :: minus_trace(quadrature_count, minus_count) = &
        reshape([0.8_dp, -0.1_dp, 0.5_dp], [quadrature_count, minus_count])
    real(dp), parameter :: surface_weights(quadrature_count) = [0.7_dp, 1.1_dp, 0.9_dp]
    real(dp), parameter :: multiplier_trace_dot(quadrature_count, multiplier_count) = &
        reshape([0.1_dp, -0.3_dp, 0.5_dp, -0.7_dp, 0.9_dp, -1.1_dp], &
        [quadrature_count, multiplier_count])
    real(dp), parameter :: plus_trace_dot(quadrature_count, plus_count) = reshape([ &
        -0.2_dp, 0.4_dp, -0.6_dp, 0.8_dp, -1.0_dp, 1.2_dp], &
        [quadrature_count, plus_count])
    real(dp), parameter :: minus_trace_dot(quadrature_count, minus_count) = &
        reshape([-0.25_dp, 0.35_dp, -0.45_dp], [quadrature_count, minus_count])
    real(dp), parameter :: surface_weights_dot(quadrature_count) = [0.2_dp, -0.1_dp, 0.3_dp]
    real(dp) :: matrix(multiplier_count, plus_count + minus_count)
    real(dp) :: matrix_dot(size(matrix, 1), size(matrix, 2))
    real(dp) :: expected(size(matrix, 1), size(matrix, 2))
    real(dp) :: matrix_plus(size(matrix, 1), size(matrix, 2))
    real(dp) :: matrix_minus(size(matrix, 1), size(matrix, 2))
    real(dp) :: multiplier_trace_bar(size(multiplier_trace, 1), &
        size(multiplier_trace, 2))
    real(dp) :: plus_trace_bar(size(plus_trace, 1), size(plus_trace, 2))
    real(dp) :: minus_trace_bar(size(minus_trace, 1), size(minus_trace, 2))
    real(dp) :: surface_weights_bar(size(surface_weights))
    real(dp) :: matrix_bar(size(matrix, 1), size(matrix, 2))
    real(dp) :: lhs, rhs
    real(dp) :: weights_plus(size(surface_weights)), weights_minus(size(surface_weights))
    real(dp) :: multiplier_plus(size(multiplier_trace, 1), size(multiplier_trace, 2))
    real(dp) :: multiplier_minus(size(multiplier_trace, 1), size(multiplier_trace, 2))
    real(dp) :: plus_plus(size(plus_trace, 1), size(plus_trace, 2))
    real(dp) :: plus_minus(size(plus_trace, 1), size(plus_trace, 2))
    real(dp) :: minus_plus(size(minus_trace, 1), size(minus_trace, 2))
    real(dp) :: minus_minus(size(minus_trace, 1), size(minus_trace, 2))
    real(dp) :: invalid_weights(size(surface_weights))
    type(fortsparse_status_t) :: status
    integer :: quadrature, multiplier, column
    logical :: all_passed

    all_passed = .true.
    expected = 0.0_dp
    do quadrature = 1, quadrature_count
        do multiplier = 1, multiplier_count
            do column = 1, plus_count
                expected(multiplier, column) = expected(multiplier, column) + &
                    surface_weights(quadrature)*multiplier_trace( &
                    quadrature, multiplier)*plus_trace(quadrature, column)
            end do
            expected(multiplier, plus_count + 1) = &
                expected(multiplier, plus_count + 1) - surface_weights(quadrature)* &
                multiplier_trace(quadrature, multiplier)*minus_trace(quadrature, 1)
        end do
    end do
    call assemble_fitted_trace_constraint( &
        multiplier_trace, plus_trace, minus_trace, surface_weights, matrix, status)
    call record_condition(status%code == 0 .and. maxval(abs(matrix - expected)) < &
        1.0e-14_dp, "fitted trace constraint matches the signed cross-mass oracle")

    call assemble_fitted_trace_constraint_jvp( &
        multiplier_trace, plus_trace, minus_trace, surface_weights, &
        multiplier_trace_dot, plus_trace_dot, minus_trace_dot, surface_weights_dot, &
        matrix_dot, status)
    call record_condition(status%code == 0, &
        "fitted trace constraint accepts all fixed-topology increments")
    multiplier_plus = multiplier_trace + epsilon_fd*multiplier_trace_dot
    multiplier_minus = multiplier_trace - epsilon_fd*multiplier_trace_dot
    plus_plus = plus_trace + epsilon_fd*plus_trace_dot
    plus_minus = plus_trace - epsilon_fd*plus_trace_dot
    minus_plus = minus_trace + epsilon_fd*minus_trace_dot
    minus_minus = minus_trace - epsilon_fd*minus_trace_dot
    weights_plus = surface_weights + epsilon_fd*surface_weights_dot
    weights_minus = surface_weights - epsilon_fd*surface_weights_dot
    call assemble_fitted_trace_constraint( &
        multiplier_plus, plus_plus, minus_plus, weights_plus, matrix_plus, status)
    call assemble_fitted_trace_constraint( &
        multiplier_minus, plus_minus, minus_minus, weights_minus, matrix_minus, status)
    call record_condition(maxval(abs(matrix_dot - &
        (matrix_plus - matrix_minus)/(2.0_dp*epsilon_fd))) < 2.0e-8_dp, &
        "fitted trace constraint JVP matches independent finite differences")

    matrix_bar = reshape([(0.13_dp*real(column, dp), &
        column = 1, size(matrix))], shape(matrix))
    call assemble_fitted_trace_constraint_vjp( &
        multiplier_trace, plus_trace, minus_trace, surface_weights, matrix_bar, &
        multiplier_trace_bar, plus_trace_bar, minus_trace_bar, &
        surface_weights_bar, status)
    lhs = sum(matrix_bar*matrix_dot)
    rhs = sum(multiplier_trace_bar*multiplier_trace_dot) + &
        sum(plus_trace_bar*plus_trace_dot) + sum(minus_trace_bar*minus_trace_dot) + &
        sum(surface_weights_bar*surface_weights_dot)
    call record_condition(status%code == 0 .and. abs(lhs - rhs) < 1.0e-12_dp, &
        "fitted trace constraint VJP satisfies the real dot-product identity")

    invalid_weights = surface_weights
    invalid_weights(2) = 0.0_dp
    call assemble_fitted_trace_constraint( &
        multiplier_trace, plus_trace, minus_trace, invalid_weights, matrix, status)
    call record_condition(status%code /= 0, &
        "fitted trace constraint rejects non-positive surface weights")
    call check_summary("fitted trace constraint")
    if (.not. all_passed) error stop 1

contains

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_fitted_trace_constraint
