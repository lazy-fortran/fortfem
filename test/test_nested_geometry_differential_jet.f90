program test_nested_geometry_differential_jet
    use check, only: check_condition, check_summary
    use fortfem_kinds, only: dp
    use fortfem_nested_geometry_differential_jet, only: &
        evaluate_nested_geometry_differential_jet, &
        evaluate_nested_geometry_differential_jet_jvp, &
        evaluate_nested_geometry_differential_jet_vjp, &
        evaluate_nested_geometry_polynomial_jet, &
        validate_nested_geometry_axis_flags
    use fortsparse, only: FORTSPARSE_OK, fortsparse_status_t
    implicit none

    integer, parameter :: n = 2
    real(dp), parameter :: eps = 2.0e-7_dp
    real(dp), parameter :: tol = 3.0e-7_dp
    real(dp) :: constant(3), linear(3, 3), quadratic(3, 3, 3), cubic(3, 3, 3, 3)
    real(dp) :: coordinates(3, n), coordinates_dot(3, n)
    real(dp) :: value(3, n), jacobian(3, 3, n), hessian(3, 3, 3, n)
    real(dp) :: third(3, 3, 3, 3, n)
    real(dp) :: metric(3, 3, n), metric_gradient(3, 3, 3, n)
    real(dp) :: metric_hessian(3, 3, 3, 3, n)
    real(dp) :: determinant(n), determinant_gradient(3, n), determinant_hessian(3, 3, n)
    real(dp) :: inverse(3, 3, n), inverse_gradient(3, 3, 3, n)
    real(dp) :: value_dot(3, n), jacobian_dot(3, 3, n), hessian_dot(3, 3, 3, n)
    real(dp) :: third_dot(3, 3, 3, 3, n)
    real(dp) :: metric_dot(3, 3, n), metric_gradient_dot(3, 3, 3, n)
    real(dp) :: metric_hessian_dot(3, 3, 3, 3, n)
    real(dp) :: determinant_dot(n), determinant_gradient_dot(3, n)
    real(dp) :: determinant_hessian_dot(3, 3, n)
    real(dp) :: inverse_dot(3, 3, n), inverse_gradient_dot(3, 3, 3, n)
    real(dp) :: value_plus(3, n), jacobian_plus(3, 3, n), hessian_plus(3, 3, 3, n)
    real(dp) :: third_plus(3, 3, 3, 3, n), metric_plus(3, 3, n)
    real(dp) :: metric_gradient_plus(3, 3, 3, n), metric_hessian_plus(3, 3, 3, 3, n)
    real(dp) :: determinant_plus(n), determinant_gradient_plus(3, n)
    real(dp) :: determinant_hessian_plus(3, 3, n), inverse_plus(3, 3, n)
    real(dp) :: inverse_gradient_plus(3, 3, 3, n)
    real(dp) :: value_minus(3, n), jacobian_minus(3, 3, n), hessian_minus(3, 3, 3, n)
    real(dp) :: third_minus(3, 3, 3, 3, n), metric_minus(3, 3, n)
    real(dp) :: metric_gradient_minus(3, 3, 3, n), metric_hessian_minus(3, 3, 3, 3, n)
    real(dp) :: determinant_minus(n), determinant_gradient_minus(3, n)
    real(dp) :: determinant_hessian_minus(3, 3, n), inverse_minus(3, 3, n)
    real(dp) :: inverse_gradient_minus(3, 3, 3, n)
    real(dp) :: metric_bar(3, 3, n), metric_gradient_bar(3, 3, 3, n)
    real(dp) :: metric_hessian_bar(3, 3, 3, 3, n), determinant_bar(n)
    real(dp) :: determinant_gradient_bar(3, n), determinant_hessian_bar(3, 3, n)
    real(dp) :: inverse_bar(3, 3, n), inverse_gradient_bar(3, 3, 3, n)
    real(dp) :: value_bar(3, n), jacobian_bar(3, 3, n), hessian_bar(3, 3, 3, n)
    real(dp) :: third_bar(3, 3, 3, 3, n), lhs, rhs
    real(dp) :: singular_jacobian(3, 3, n)
    logical :: axis_limit(n), bad_axis(1), all_passed
    type(fortsparse_status_t) :: status
    integer :: i

    all_passed = .true.
    constant = [0.3_dp, -0.2_dp, 0.7_dp]
    linear = reshape([1.2_dp, 0.1_dp, -0.1_dp, 0.2_dp, 1.1_dp, 0.05_dp, &
                      0.0_dp, -0.2_dp, 0.9_dp], [3, 3])
    quadratic = 0.0_dp
    quadratic(1, 1, 1) = 0.2_dp
    quadratic(2, 2, 3) = -0.15_dp
    quadratic(3, 3, 2) = 0.1_dp
    cubic = 0.0_dp
    cubic(1, 1, 2, 3) = 0.04_dp
    cubic(2, 2, 2, 1) = -0.03_dp
    cubic(3, 3, 1, 2) = 0.02_dp
    coordinates = reshape([0.2_dp, -0.1_dp, 0.3_dp, 0.4_dp, 0.1_dp, -0.2_dp], [3, n])
    coordinates_dot = reshape([0.03_dp, -0.02_dp, 0.01_dp, -0.04_dp, 0.05_dp, 0.02_dp], [3, n])
    axis_limit = [.false., .false.]
    call evaluate_nested_geometry_polynomial_jet( &
        constant, linear, quadratic, cubic, coordinates, axis_limit, value, jacobian, &
        hessian, third, metric, metric_gradient, metric_hessian, determinant, &
        determinant_gradient, determinant_hessian, inverse, inverse_gradient, status)
    call record(status%code == FORTSPARSE_OK, "polynomial jet evaluates")
    call record(maxval(abs(metric - transpose_metric(metric))) < 2.0e-13_dp, &
                "metric is symmetric")
call record(abs(determinant(1) - determinant3_oracle(jacobian(:, :, 1))) < 2.0e-13_dp, &
                "determinant matches independent 3-by-3 oracle")

    value_dot = 0.0_dp
    jacobian_dot = 0.0_dp
    hessian_dot = reshape([(0.001_dp*real(i, dp), i=1, 54)], [3, 3, 3, n])
    third_dot = reshape([(0.0003_dp*real(i, dp), i=1, 162)], [3, 3, 3, 3, n])
    jacobian_dot = reshape([(0.002_dp*real(i, dp), i=1, 18)], [3, 3, n])
    call evaluate_nested_geometry_differential_jet_jvp( &
        value, jacobian, hessian, third, axis_limit, value_dot, jacobian_dot, &
        hessian_dot, third_dot, metric_dot, metric_gradient_dot, metric_hessian_dot, &
      determinant_dot, determinant_gradient_dot, determinant_hessian_dot, inverse_dot, &
        inverse_gradient_dot, status)
    call record(status%code == FORTSPARSE_OK, "jet JVP evaluates")
    call evaluate_nested_geometry_differential_jet( &
        value + eps*value_dot, jacobian + eps*jacobian_dot, hessian + eps*hessian_dot, &
        third + eps*third_dot, axis_limit, metric_plus, metric_gradient_plus, &
        metric_hessian_plus, determinant_plus, determinant_gradient_plus, &
        determinant_hessian_plus, inverse_plus, inverse_gradient_plus, status)
    call evaluate_nested_geometry_differential_jet( &
        value - eps*value_dot, jacobian - eps*jacobian_dot, hessian - eps*hessian_dot, &
        third - eps*third_dot, axis_limit, metric_minus, metric_gradient_minus, &
        metric_hessian_minus, determinant_minus, determinant_gradient_minus, &
        determinant_hessian_minus, inverse_minus, inverse_gradient_minus, status)
call record(maxval(abs(metric_dot - (metric_plus - metric_minus)/(2.0_dp*eps))) < tol, &
                "metric JVP matches finite differences")
    call record(maxval(abs(determinant_hessian_dot - &
          (determinant_hessian_plus - determinant_hessian_minus)/(2.0_dp*eps))) < tol, &
                "determinant Hessian JVP matches finite differences")

    metric_bar = reshape([(0.01_dp*real(i, dp), i=1, 18)], [3, 3, n])
    metric_gradient_bar = reshape([(0.002_dp*real(i, dp), i=1, 54)], [3, 3, 3, n])
    metric_hessian_bar = reshape([(0.0003_dp*real(i, dp), i=1, 162)], [3, 3, 3, 3, n])
    determinant_bar = [0.3_dp, -0.2_dp]
    determinant_gradient_bar = reshape([(0.004_dp*real(i, dp), i=1, 6)], [3, n])
    determinant_hessian_bar = reshape([(0.0007_dp*real(i, dp), i=1, 18)], [3, 3, n])
    inverse_bar = reshape([(0.005_dp*real(i, dp), i=1, 18)], [3, 3, n])
    inverse_gradient_bar = reshape([(0.0006_dp*real(i, dp), i=1, 54)], [3, 3, 3, n])
    call evaluate_nested_geometry_differential_jet_vjp( &
        value, jacobian, hessian, third, axis_limit, metric_bar, metric_gradient_bar, &
        metric_hessian_bar, determinant_bar, determinant_gradient_bar, &
        determinant_hessian_bar, inverse_bar, inverse_gradient_bar, value_bar, &
        jacobian_bar, hessian_bar, third_bar, status)
    lhs = sum(metric_bar*metric_dot) + sum(metric_gradient_bar*metric_gradient_dot) + &
   sum(metric_hessian_bar*metric_hessian_dot) + sum(determinant_bar*determinant_dot) + &
          sum(determinant_gradient_bar*determinant_gradient_dot) + &
 sum(determinant_hessian_bar*determinant_hessian_dot) + sum(inverse_bar*inverse_dot) + &
          sum(inverse_gradient_bar*inverse_gradient_dot)
    rhs = sum(jacobian_bar*jacobian_dot) + sum(hessian_bar*hessian_dot) + sum(third_bar*third_dot)
    call record(abs(lhs-rhs) < 3.0e-7_dp, "jet VJP satisfies independent dot-product oracle")

    bad_axis = [.false.]
    call validate_nested_geometry_axis_flags(bad_axis, n, status)
  call record(status%code /= FORTSPARSE_OK, "axis-limit flags reject a sample mismatch")
    singular_jacobian = jacobian
    singular_jacobian(:, :, 1) = 0.0_dp
    call evaluate_nested_geometry_differential_jet( &
        value, singular_jacobian, hessian, third, axis_limit, metric, metric_gradient, &
      metric_hessian, determinant, determinant_gradient, determinant_hessian, inverse, &
        inverse_gradient, status)
    call record(status%code /= FORTSPARSE_OK, "singular non-axis jet is rejected")
    axis_limit(1) = .true.
    call evaluate_nested_geometry_differential_jet( &
        value, singular_jacobian, hessian, third, axis_limit, metric, metric_gradient, &
      metric_hessian, determinant, determinant_gradient, determinant_hessian, inverse, &
        inverse_gradient, status)
    call record(status%code == FORTSPARSE_OK .and. maxval(abs(inverse(:, :, 1))) == 0.0_dp, &
                "caller-marked axis limit avoids hidden inverse division")
    call check_summary("nested geometry differential jet")
    if (.not. all_passed) error stop 1

contains

    subroutine record(condition, description)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: description
        call check_condition(condition, description)
        all_passed = all_passed .and. condition
    end subroutine record

    pure function transpose_metric(input) result(output)
        real(dp), intent(in) :: input(:, :, :)
        real(dp) :: output(size(input, 2), size(input, 1), size(input, 3))
        integer :: sample
        do sample = 1, size(input, 3)
            output(:, :, sample) = transpose(input(:, :, sample))
        end do
    end function transpose_metric

    pure real(dp) function determinant3_oracle(matrix) result(value)
        real(dp), intent(in) :: matrix(3, 3)
        value = matrix(1, 1)*(matrix(2, 2)*matrix(3, 3) - matrix(2, 3)*matrix(3, 2)) - &
                matrix(1, 2)*(matrix(2, 1)*matrix(3, 3) - matrix(2, 3)*matrix(3, 1)) + &
                matrix(1, 3)*(matrix(2, 1)*matrix(3, 2) - matrix(2, 2)*matrix(3, 1))
    end function determinant3_oracle

end program test_nested_geometry_differential_jet
