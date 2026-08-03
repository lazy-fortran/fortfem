program test_complex_geometry_mortar_component_coupling
    use, intrinsic :: ieee_arithmetic, only: ieee_quiet_nan, ieee_value
    use check, only: check_condition, check_summary
    use fortfem_complex_geometry_mortar_component_coupling, only: &
        assemble_complex_geometry_mortar_component_coupling, &
        assemble_complex_geometry_mortar_component_coupling_jvp, &
        assemble_complex_geometry_mortar_component_coupling_vjp
    use fortfem_kinds, only: dp
    use fortsparse, only: fortsparse_status_t
    implicit none

    integer, parameter :: component_count = 2, quadrature_count = 2
    integer, parameter :: test_count = 2, trial_count = 2
    real(dp), parameter :: step = 1.0e-7_dp
    complex(dp) :: test_trace(component_count, quadrature_count, test_count)
    complex(dp) :: trial_trace(component_count, quadrature_count, trial_count)
    complex(dp) :: metric(component_count, component_count, quadrature_count)
    complex(dp) :: test_dot(component_count, quadrature_count, test_count)
    complex(dp) :: trial_dot(component_count, quadrature_count, trial_count)
    complex(dp) :: metric_dot(component_count, component_count, quadrature_count)
    real(dp) :: weights(quadrature_count), jacobian(quadrature_count)
    real(dp) :: weights_dot(quadrature_count), jacobian_dot(quadrature_count)
    complex(dp) :: matrix(test_count, trial_count), expected(test_count, trial_count)
    real(dp) :: physical_weights(quadrature_count)
    complex(dp) :: matrix_dot(test_count, trial_count)
    real(dp) :: physical_weights_dot(quadrature_count)
    complex(dp) :: matrix_plus(test_count, trial_count)
    complex(dp) :: matrix_minus(test_count, trial_count)
    real(dp) :: physical_plus(quadrature_count), physical_minus(quadrature_count)
    complex(dp) :: test_plus(component_count, quadrature_count, test_count)
    complex(dp) :: test_minus(component_count, quadrature_count, test_count)
    complex(dp) :: trial_plus(component_count, quadrature_count, trial_count)
    complex(dp) :: trial_minus(component_count, quadrature_count, trial_count)
    complex(dp) :: metric_plus(component_count, component_count, quadrature_count)
    complex(dp) :: metric_minus(component_count, component_count, quadrature_count)
    real(dp) :: weights_plus(quadrature_count), weights_minus(quadrature_count)
    real(dp) :: jacobian_plus(quadrature_count), jacobian_minus(quadrature_count)
    complex(dp) :: matrix_bar(test_count, trial_count)
    real(dp) :: physical_bar(quadrature_count)
    complex(dp) :: test_bar(component_count, quadrature_count, test_count)
    complex(dp) :: trial_bar(component_count, quadrature_count, trial_count)
    complex(dp) :: metric_bar(component_count, component_count, quadrature_count)
    real(dp) :: weights_bar(quadrature_count), jacobian_bar(quadrature_count)
    complex(dp) :: bad_metric(component_count, component_count, quadrature_count)
    complex(dp) :: bad_test_dot(component_count, quadrature_count, test_count)
    complex(dp) :: short_matrix(test_count - 1, trial_count)
    real(dp) :: lhs, rhs
    logical :: all_passed
    type(fortsparse_status_t) :: status

    all_passed = .true.
    test_trace = cmplx(reshape([ &
        0.8_dp, -0.2_dp, 0.3_dp, 0.1_dp, 0.4_dp, -0.7_dp, -0.5_dp, &
        0.6_dp], shape(test_trace)), reshape([ &
        -0.1_dp, 0.4_dp, 0.2_dp, -0.3_dp, 0.5_dp, 0.1_dp, -0.2_dp, &
        0.3_dp], shape(test_trace)), dp)
    trial_trace = cmplx(reshape([ &
        0.5_dp, 0.4_dp, -0.2_dp, 0.6_dp, -0.3_dp, 0.1_dp, 0.7_dp, &
        -0.8_dp], shape(trial_trace)), reshape([ &
        0.3_dp, -0.2_dp, 0.5_dp, 0.1_dp, -0.4_dp, 0.6_dp, -0.1_dp, &
        0.2_dp], shape(trial_trace)), dp)
    metric = cmplx(reshape([ &
        1.0_dp, 0.2_dp, -0.1_dp, 1.3_dp, 0.7_dp, -0.4_dp, 0.5_dp, &
        1.1_dp], shape(metric)), reshape([ &
        0.1_dp, -0.3_dp, 0.2_dp, 0.4_dp, -0.2_dp, 0.5_dp, 0.3_dp, &
        -0.1_dp], shape(metric)), dp)
    test_dot = cmplx(reshape([ &
        0.03_dp, -0.02_dp, 0.05_dp, -0.04_dp, 0.01_dp, 0.02_dp, &
        -0.01_dp, 0.06_dp], shape(test_dot)), reshape([ &
        -0.02_dp, 0.01_dp, 0.04_dp, 0.03_dp, -0.05_dp, 0.02_dp, &
        0.06_dp, -0.01_dp], shape(test_dot)), dp)
    trial_dot = cmplx(reshape([ &
        -0.01_dp, 0.06_dp, 0.02_dp, 0.04_dp, -0.03_dp, 0.05_dp, &
        0.02_dp, -0.04_dp], shape(trial_dot)), reshape([ &
        0.03_dp, -0.01_dp, 0.02_dp, -0.05_dp, 0.04_dp, 0.01_dp, &
        -0.02_dp, 0.06_dp], shape(trial_dot)), dp)
    metric_dot = cmplx(reshape([ &
        0.02_dp, -0.03_dp, 0.04_dp, 0.01_dp, -0.05_dp, 0.06_dp, &
        0.01_dp, -0.02_dp], shape(metric_dot)), reshape([ &
        -0.01_dp, 0.05_dp, -0.04_dp, 0.02_dp, 0.03_dp, -0.02_dp, &
        0.06_dp, 0.01_dp], shape(metric_dot)), dp)
    weights = [1.5_dp, 0.8_dp]
    jacobian = [1.2_dp, 2.0_dp]
    weights_dot = [0.2_dp, -0.1_dp]
    jacobian_dot = [0.04_dp, 0.05_dp]

    call independent_cross_mass( &
        test_trace, trial_trace, metric, weights, jacobian, expected)
    call assemble_complex_geometry_mortar_component_coupling( &
        test_trace, trial_trace, weights, jacobian, metric, matrix, &
        physical_weights, status)
    call record(status%code == 0 .and. &
        maxval(abs(matrix - expected)) < 2.0e-13_dp, &
        "complex component mortar matches an independent nested-loop oracle")
    call record(maxval(abs(physical_weights - weights*jacobian)) < 1.0e-14_dp, &
        "complex component mortar returns real physical weights")

    call assemble_complex_geometry_mortar_component_coupling_jvp( &
        test_trace, trial_trace, weights, jacobian, metric, test_dot, &
        trial_dot, weights_dot, jacobian_dot, metric_dot, matrix_dot, &
        physical_weights_dot, status)
    test_plus = test_trace + step*test_dot
    test_minus = test_trace - step*test_dot
    trial_plus = trial_trace + step*trial_dot
    trial_minus = trial_trace - step*trial_dot
    metric_plus = metric + step*metric_dot
    metric_minus = metric - step*metric_dot
    weights_plus = weights + step*weights_dot
    weights_minus = weights - step*weights_dot
    jacobian_plus = jacobian + step*jacobian_dot
    jacobian_minus = jacobian - step*jacobian_dot
    call assemble_complex_geometry_mortar_component_coupling( &
        test_plus, trial_plus, weights_plus, jacobian_plus, metric_plus, &
        matrix_plus, physical_plus, status)
    call assemble_complex_geometry_mortar_component_coupling( &
        test_minus, trial_minus, weights_minus, jacobian_minus, metric_minus, &
        matrix_minus, physical_minus, status)
    call record(status%code == 0 .and. maxval(abs(matrix_dot - &
        (matrix_plus - matrix_minus)/(2.0_dp*step))) < 2.0e-8_dp, &
        "complex component mortar JVP matches central reassembly")
    call record(maxval(abs(physical_weights_dot - &
        (physical_plus - physical_minus)/(2.0_dp*step))) < 2.0e-8_dp, &
        "complex physical-weight JVP matches central reassembly")

    matrix_bar = cmplx(reshape([0.4_dp, -0.3_dp, 0.2_dp, 0.5_dp], &
        shape(matrix_bar)), reshape([-0.2_dp, 0.6_dp, 0.1_dp, -0.4_dp], &
        shape(matrix_bar)), dp)
    physical_bar = [0.7_dp, -0.2_dp]
    call assemble_complex_geometry_mortar_component_coupling_vjp( &
        test_trace, trial_trace, weights, jacobian, metric, matrix_bar, &
        physical_bar, test_bar, trial_bar, weights_bar, jacobian_bar, &
        metric_bar, status)
    lhs = real(sum(conjg(matrix_bar)*matrix_dot), dp) + &
        sum(physical_bar*physical_weights_dot)
    rhs = real(sum(conjg(test_bar)*test_dot) + &
        sum(conjg(trial_bar)*trial_dot) + &
        sum(conjg(metric_bar)*metric_dot), dp) + &
        sum(weights_bar*weights_dot) + sum(jacobian_bar*jacobian_dot)
    call record(status%code == 0 .and. abs(lhs - rhs) < 3.0e-13_dp, &
        "complex component mortar VJP satisfies the real-part adjoint identity")

    call assemble_complex_geometry_mortar_component_coupling( &
        test_trace, trial_trace, weights, jacobian, metric, short_matrix, &
        physical_weights, status)
    call record(status%code /= 0, &
        "complex component mortar rejects an incompatible output shape")
    bad_metric = metric
    bad_metric(1, 1, 1) = cmplx( &
        real(metric(1, 1, 1), dp), ieee_value(0.0_dp, ieee_quiet_nan), dp)
    call assemble_complex_geometry_mortar_component_coupling( &
        test_trace, trial_trace, weights, jacobian, bad_metric, matrix, &
        physical_weights, status)
    call record(status%code /= 0, &
        "complex component mortar rejects a non-finite imaginary part")
    bad_test_dot = test_dot
    bad_test_dot(1, 1, 1) = cmplx( &
        ieee_value(0.0_dp, ieee_quiet_nan), 0.0_dp, dp)
    call assemble_complex_geometry_mortar_component_coupling_jvp( &
        test_trace, trial_trace, weights, jacobian, metric, bad_test_dot, &
        trial_dot, weights_dot, jacobian_dot, metric_dot, matrix_dot, &
        physical_weights_dot, status)
    call record(status%code /= 0, &
        "complex component mortar JVP rejects a non-finite direction")

    call check_summary("complex geometry mortar component coupling")
    if (.not. all_passed) error stop 1

contains

    subroutine independent_cross_mass( &
            test_values, trial_values, component_metric, reference_weights, &
            surface_jacobian, result)
        complex(dp), intent(in) :: test_values(:, :, :), trial_values(:, :, :)
        complex(dp), intent(in) :: component_metric(:, :, :)
        real(dp), intent(in) :: reference_weights(:), surface_jacobian(:)
        complex(dp), intent(out) :: result(:, :)
        integer :: q, i, j, a, b

        result = cmplx(0.0_dp, 0.0_dp, dp)
        do j = 1, size(trial_values, 3)
            do i = 1, size(test_values, 3)
                do q = 1, size(reference_weights)
                    do b = 1, size(trial_values, 1)
                        do a = 1, size(test_values, 1)
                            result(i, j) = result(i, j) + &
                                reference_weights(q)*surface_jacobian(q)* &
                                test_values(a, q, i)*component_metric(a, b, q)* &
                                trial_values(b, q, j)
                        end do
                    end do
                end do
            end do
        end do
    end subroutine independent_cross_mass

    subroutine record(condition, description)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record

end program test_complex_geometry_mortar_component_coupling
