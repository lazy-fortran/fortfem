program test_force_balance_metric_residual
    use check, only: check_condition, check_summary
    use fortfem_feec, only: &
        assemble_metric_force_balance_residual, &
        assemble_metric_force_balance_residual_jvp, &
        assemble_metric_force_balance_residual_vjp
    use fortfem_kinds, only: dp
    use fortsparse, only: fortsparse_status_t
    implicit none

    integer, parameter :: volume_count = 2, boundary_count = 1, sheet_count = 1
    integer, parameter :: test_count = 2, component_count = 2
    real(dp), parameter :: epsilon_fd = 1.0e-7_dp
    real(dp) :: volume_test(volume_count, test_count, component_count)
    real(dp) :: volume_force(volume_count, component_count)
    real(dp) :: quadrature_weights(volume_count), volume_jacobian(volume_count)
    real(dp) :: boundary_test(boundary_count, test_count, component_count)
    real(dp) :: boundary_force(boundary_count, component_count)
    real(dp) :: boundary_weights(boundary_count)
    real(dp) :: sheet_test(sheet_count, test_count, component_count)
    real(dp) :: sheet_force(sheet_count, component_count)
    real(dp) :: sheet_weights(sheet_count)
    real(dp) :: volume_test_dot(volume_count, test_count, component_count)
    real(dp) :: volume_force_dot(volume_count, component_count)
    real(dp) :: quadrature_weights_dot(volume_count), volume_jacobian_dot(volume_count)
    real(dp) :: boundary_test_dot(boundary_count, test_count, component_count)
    real(dp) :: boundary_force_dot(boundary_count, component_count)
    real(dp) :: boundary_weights_dot(boundary_count)
    real(dp) :: sheet_test_dot(sheet_count, test_count, component_count)
    real(dp) :: sheet_force_dot(sheet_count, component_count)
    real(dp) :: sheet_weights_dot(sheet_count)
    real(dp) :: residual(test_count, component_count)
    real(dp) :: residual_dot(test_count, component_count)
    real(dp) :: residual_plus(test_count, component_count)
    real(dp) :: residual_minus(test_count, component_count)
    real(dp) :: residual_bar(test_count, component_count)
    real(dp) :: volume_test_bar(volume_count, test_count, component_count)
    real(dp) :: volume_force_bar(volume_count, component_count)
    real(dp) :: quadrature_weights_bar(volume_count), volume_jacobian_bar(volume_count)
    real(dp) :: boundary_test_bar(boundary_count, test_count, component_count)
    real(dp) :: boundary_force_bar(boundary_count, component_count)
    real(dp) :: boundary_weights_bar(boundary_count)
    real(dp) :: sheet_test_bar(sheet_count, test_count, component_count)
    real(dp) :: sheet_force_bar(sheet_count, component_count)
    real(dp) :: sheet_weights_bar(sheet_count)
    real(dp) :: expected(test_count, component_count), lhs, rhs
    real(dp) :: q
    type(fortsparse_status_t) :: status
    integer :: sample, test_function, component
    logical :: all_passed

    all_passed = .true.
    volume_test = reshape([ &
        1.0_dp, -0.4_dp, 0.7_dp, 0.2_dp, -0.3_dp, 0.8_dp, &
        0.5_dp, 1.2_dp], shape(volume_test))
    volume_force = reshape([0.8_dp, -1.1_dp, 0.6_dp, 1.4_dp], shape(volume_force))
    quadrature_weights = [0.7_dp, 1.3_dp]
    volume_jacobian = [2.0_dp, 0.5_dp]
    boundary_test = reshape([0.6_dp, -0.2_dp, 1.1_dp, 0.4_dp], shape(boundary_test))
    boundary_force = reshape([0.3_dp, -0.9_dp], shape(boundary_force))
    boundary_weights = [0.8_dp]
    sheet_test = reshape([-0.7_dp, 0.2_dp, 0.5_dp, 1.0_dp], shape(sheet_test))
    sheet_force = reshape([1.2_dp, -0.6_dp], shape(sheet_force))
    sheet_weights = [0.9_dp]
    volume_test_dot = 0.03_dp*reshape( &
        [(real(sample, dp), sample = 1, size(volume_test))], &
        shape(volume_test))
    volume_force_dot = 0.04_dp*reshape( &
        [(real(sample, dp), sample = 1, size(volume_force))], &
        shape(volume_force))
    quadrature_weights_dot = [0.11_dp, -0.08_dp]
    volume_jacobian_dot = [-0.07_dp, 0.12_dp]
    boundary_test_dot = 0.02_dp*reshape( &
        [(real(sample, dp), sample = 1, size(boundary_test))], &
        shape(boundary_test))
    boundary_force_dot = 0.05_dp*reshape( &
        [(real(sample, dp), sample = 1, size(boundary_force))], &
        shape(boundary_force))
    boundary_weights_dot = [-0.06_dp]
    sheet_test_dot = 0.06_dp*reshape( &
        [(real(sample, dp), sample = 1, size(sheet_test))], &
        shape(sheet_test))
    sheet_force_dot = 0.07_dp*reshape( &
        [(real(sample, dp), sample = 1, size(sheet_force))], &
        shape(sheet_force))
    sheet_weights_dot = [0.09_dp]

    call assemble_metric_force_balance_residual( &
        volume_test, volume_force, quadrature_weights, volume_jacobian, &
        boundary_test, boundary_force, boundary_weights, sheet_test, sheet_force, &
        sheet_weights, residual, status)
    expected = 0.0_dp
    do sample = 1, volume_count
        q = quadrature_weights(sample)*volume_jacobian(sample)
        do test_function = 1, test_count
            do component = 1, component_count
                expected(test_function, component) = &
                    expected(test_function, component) + q*volume_test( &
                    sample, test_function, component)*volume_force(sample, component)
            end do
        end do
    end do
    call add_expected(boundary_test, boundary_force, boundary_weights, expected)
    call add_expected(sheet_test, sheet_force, sheet_weights, expected)
    call record_condition( &
        status%code == 0 .and. maxval(abs(residual - expected)) < 1.0e-14_dp, &
        "metric force residual matches independent Jacobian-weighted oracle")

    call assemble_metric_force_balance_residual_jvp( &
        volume_test, volume_force, quadrature_weights, volume_jacobian, boundary_test, &
        boundary_force, boundary_weights, sheet_test, sheet_force, sheet_weights, &
        volume_test_dot, volume_force_dot, quadrature_weights_dot, &
        volume_jacobian_dot, &
        boundary_test_dot, boundary_force_dot, boundary_weights_dot, sheet_test_dot, &
        sheet_force_dot, sheet_weights_dot, residual_dot, status)
    call assemble_metric_force_balance_residual( &
        volume_test + epsilon_fd*volume_test_dot, &
        volume_force + epsilon_fd*volume_force_dot, &
        quadrature_weights + epsilon_fd*quadrature_weights_dot, &
        volume_jacobian + epsilon_fd*volume_jacobian_dot, &
        boundary_test + epsilon_fd*boundary_test_dot, &
        boundary_force + epsilon_fd*boundary_force_dot, &
        boundary_weights + epsilon_fd*boundary_weights_dot, &
        sheet_test + epsilon_fd*sheet_test_dot, &
        sheet_force + epsilon_fd*sheet_force_dot, &
        sheet_weights + epsilon_fd*sheet_weights_dot, &
        residual_plus, status)
    call assemble_metric_force_balance_residual( &
        volume_test - epsilon_fd*volume_test_dot, &
        volume_force - epsilon_fd*volume_force_dot, &
        quadrature_weights - epsilon_fd*quadrature_weights_dot, &
        volume_jacobian - epsilon_fd*volume_jacobian_dot, &
        boundary_test - epsilon_fd*boundary_test_dot, &
        boundary_force - epsilon_fd*boundary_force_dot, &
        boundary_weights - epsilon_fd*boundary_weights_dot, &
        sheet_test - epsilon_fd*sheet_test_dot, &
        sheet_force - epsilon_fd*sheet_force_dot, &
        sheet_weights - epsilon_fd*sheet_weights_dot, &
        residual_minus, status)
    call record_condition(status%code == 0 .and. maxval(abs(residual_dot - &
        (residual_plus - residual_minus)/(2.0_dp*epsilon_fd))) < 2.0e-8_dp, &
        "metric force residual JVP matches independent central difference")

    residual_bar = reshape([0.4_dp, -0.7_dp, 0.2_dp, 0.9_dp], shape(residual_bar))
    call assemble_metric_force_balance_residual_vjp( &
        volume_test, volume_force, quadrature_weights, volume_jacobian, boundary_test, &
        boundary_force, boundary_weights, sheet_test, sheet_force, sheet_weights, &
        residual_bar, &
        volume_test_bar, volume_force_bar, quadrature_weights_bar, &
        volume_jacobian_bar, &
        boundary_test_bar, boundary_force_bar, boundary_weights_bar, sheet_test_bar, &
        sheet_force_bar, sheet_weights_bar, status)
    lhs = sum(residual_bar*residual_dot)
    rhs = sum(volume_test_bar*volume_test_dot) + &
        sum(volume_force_bar*volume_force_dot) + &
        sum(quadrature_weights_bar*quadrature_weights_dot) + &
        sum(volume_jacobian_bar*volume_jacobian_dot) + &
        sum(boundary_test_bar*boundary_test_dot) + &
        sum(boundary_force_bar*boundary_force_dot) + &
        sum(boundary_weights_bar*boundary_weights_dot) + &
        sum(sheet_test_bar*sheet_test_dot) + &
        sum(sheet_force_bar*sheet_force_dot) + sum(sheet_weights_bar*sheet_weights_dot)
    call record_condition(status%code == 0 .and. abs(lhs - rhs) < 1.0e-12_dp, &
        "metric force residual VJP satisfies the real dot-product identity")

    volume_jacobian(1) = 0.0_dp
    call assemble_metric_force_balance_residual( &
        volume_test, volume_force, quadrature_weights, volume_jacobian, boundary_test, &
        boundary_force, boundary_weights, sheet_test, sheet_force, sheet_weights, &
        residual, status)
    call record_condition( &
        status%code /= 0, "metric force residual rejects nonpositive Jacobian")
    call check_summary("metric force-balance residual")
    if (.not. all_passed) error stop 1

contains

    subroutine add_expected(test, force, weights, result)
        real(dp), intent(in) :: test(:, :, :), force(:, :), weights(:)
        real(dp), intent(inout) :: result(:, :)
        integer :: local_sample, local_test, local_component

        do local_sample = 1, size(weights)
            do local_test = 1, test_count
                do local_component = 1, component_count
                    result(local_test, local_component) = &
                        result(local_test, local_component) + weights(local_sample)* &
                        test(local_sample, local_test, local_component)* &
                        force(local_sample, local_component)
                end do
            end do
        end do
    end subroutine add_expected

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_force_balance_metric_residual
