program test_force_balance_residual
    use check, only: check_condition, check_summary
    use fortfem_api, only: assemble_force_balance_residual, &
        assemble_force_balance_residual_jvp, assemble_force_balance_residual_vjp
    use fortfem_kinds, only: dp
    use fortsparse, only: fortsparse_status_t
    implicit none

    integer, parameter :: volume_count = 2, boundary_count = 1, sheet_count = 2
    integer, parameter :: test_count = 2, component_count = 3
    real(dp), parameter :: epsilon = 1.0e-7_dp
    real(dp) :: volume_test(volume_count, test_count, component_count)
    real(dp) :: volume_force(volume_count, component_count)
    real(dp) :: volume_weights(volume_count)
    real(dp) :: boundary_test(boundary_count, test_count, component_count)
    real(dp) :: boundary_force(boundary_count, component_count)
    real(dp) :: boundary_weights(boundary_count)
    real(dp) :: sheet_test(sheet_count, test_count, component_count)
    real(dp) :: sheet_force(sheet_count, component_count)
    real(dp) :: sheet_weights(sheet_count)
    real(dp) :: volume_test_dot(volume_count, test_count, component_count)
    real(dp) :: volume_force_dot(volume_count, component_count)
    real(dp) :: volume_weights_dot(volume_count)
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
    real(dp) :: volume_weights_bar(volume_count)
    real(dp) :: boundary_test_bar(boundary_count, test_count, component_count)
    real(dp) :: boundary_force_bar(boundary_count, component_count)
    real(dp) :: boundary_weights_bar(boundary_count)
    real(dp) :: sheet_test_bar(sheet_count, test_count, component_count)
    real(dp) :: sheet_force_bar(sheet_count, component_count)
    real(dp) :: sheet_weights_bar(sheet_count)
    real(dp) :: residual_reference(test_count, component_count)
    real(dp) :: lhs, rhs
    type(fortsparse_status_t) :: status
    integer :: i
    logical :: all_passed

    all_passed = .true.
    volume_test = reshape([ &
        1.0_dp, -0.5_dp, 0.25_dp, 0.75_dp, 1.25_dp, -1.5_dp, &
        0.5_dp, 0.8_dp, -0.7_dp, 0.4_dp, 1.1_dp, -0.9_dp], &
        shape(volume_test))
    volume_force = reshape([0.8_dp, -1.2_dp, 0.4_dp, 1.3_dp, -0.6_dp, 0.9_dp], &
        shape(volume_force))
    volume_weights = [1.2_dp, 0.7_dp]
    boundary_test = reshape([0.6_dp, -0.2_dp, 1.4_dp, -0.8_dp, 0.3_dp, 0.5_dp], &
        shape(boundary_test))
    boundary_force = reshape([1.1_dp, -0.4_dp, 0.9_dp], shape(boundary_force))
    boundary_weights = [0.9_dp]
    sheet_test = reshape([ &
        0.2_dp, 0.7_dp, -0.3_dp, 1.0_dp, -0.4_dp, 0.6_dp, &
        0.9_dp, -0.5_dp, 0.8_dp, 0.1_dp, 1.2_dp, -0.7_dp], shape(sheet_test))
    sheet_force = reshape([0.3_dp, 1.1_dp, -0.8_dp, -0.6_dp, 0.5_dp, 0.4_dp], &
        shape(sheet_force))
    sheet_weights = [0.4_dp, 1.3_dp]

    volume_test_dot = 0.02_dp*reshape([(real(i, dp), i = 1, size(volume_test))], &
        shape(volume_test))
    volume_force_dot = 0.03_dp*reshape([(real(i, dp), i = 1, size(volume_force))], &
        shape(volume_force))
    volume_weights_dot = [0.11_dp, -0.08_dp]
    boundary_test_dot = 0.04_dp*reshape([(real(i, dp), i = 1, size(boundary_test))], &
        shape(boundary_test))
    boundary_force_dot = 0.05_dp*reshape([(real(i, dp), i = 1, size(boundary_force))], &
        shape(boundary_force))
    boundary_weights_dot = [-0.06_dp]
    sheet_test_dot = 0.06_dp*reshape([(real(i, dp), i = 1, size(sheet_test))], &
        shape(sheet_test))
    sheet_force_dot = 0.07_dp*reshape([(real(i, dp), i = 1, size(sheet_force))], &
        shape(sheet_force))
    sheet_weights_dot = [0.09_dp, -0.12_dp]

    call assemble_force_balance_residual( &
        volume_test, volume_force, volume_weights, boundary_test, boundary_force, &
        boundary_weights, sheet_test, sheet_force, sheet_weights, residual, status)
    call assemble_reference(residual_reference)
    call record_condition(status%code == 0 .and. &
        maxval(abs(residual - residual_reference)) < 1.0e-14_dp, &
        "force-balance composition matches an independent volume/boundary/sheet oracle")

    call assemble_force_balance_residual_jvp( &
        volume_test, volume_force, volume_weights, boundary_test, boundary_force, &
        boundary_weights, sheet_test, sheet_force, sheet_weights, volume_test_dot, &
        volume_force_dot, volume_weights_dot, boundary_test_dot, boundary_force_dot, &
        boundary_weights_dot, sheet_test_dot, sheet_force_dot, sheet_weights_dot, &
        residual_dot, status)
    call assemble_force_balance_residual( &
        volume_test + epsilon*volume_test_dot, volume_force + epsilon*volume_force_dot, &
        volume_weights + epsilon*volume_weights_dot, boundary_test + epsilon*boundary_test_dot, &
        boundary_force + epsilon*boundary_force_dot, &
        boundary_weights + epsilon*boundary_weights_dot, sheet_test + epsilon*sheet_test_dot, &
        sheet_force + epsilon*sheet_force_dot, sheet_weights + epsilon*sheet_weights_dot, &
        residual_plus, status)
    call assemble_force_balance_residual( &
        volume_test - epsilon*volume_test_dot, volume_force - epsilon*volume_force_dot, &
        volume_weights - epsilon*volume_weights_dot, boundary_test - epsilon*boundary_test_dot, &
        boundary_force - epsilon*boundary_force_dot, &
        boundary_weights - epsilon*boundary_weights_dot, sheet_test - epsilon*sheet_test_dot, &
        sheet_force - epsilon*sheet_force_dot, sheet_weights - epsilon*sheet_weights_dot, &
        residual_minus, status)
    call record_condition(maxval(abs(residual_dot - &
        (residual_plus - residual_minus)/(2.0_dp*epsilon))) < 2.0e-8_dp, &
        "force-balance JVP matches an independent central difference")

    residual_bar = reshape([0.4_dp, -0.7_dp, 0.2_dp, 0.9_dp, -0.3_dp, 0.6_dp], &
        shape(residual_bar))
    call assemble_force_balance_residual_vjp( &
        volume_test, volume_force, volume_weights, boundary_test, boundary_force, &
        boundary_weights, sheet_test, sheet_force, sheet_weights, residual_bar, &
        volume_test_bar, volume_force_bar, volume_weights_bar, boundary_test_bar, &
        boundary_force_bar, boundary_weights_bar, sheet_test_bar, sheet_force_bar, &
        sheet_weights_bar, status)
    lhs = sum(residual_bar*residual_dot)
    rhs = sum(volume_test_bar*volume_test_dot) + sum(volume_force_bar*volume_force_dot) + &
        sum(volume_weights_bar*volume_weights_dot) + sum(boundary_test_bar*boundary_test_dot) + &
        sum(boundary_force_bar*boundary_force_dot) + sum(boundary_weights_bar*boundary_weights_dot) + &
        sum(sheet_test_bar*sheet_test_dot) + sum(sheet_force_bar*sheet_force_dot) + &
        sum(sheet_weights_bar*sheet_weights_dot)
    call record_condition(status%code == 0 .and. abs(lhs - rhs) < 1.0e-12_dp, &
        "force-balance VJP satisfies the real dot-product identity")

    volume_weights(1) = -volume_weights(1)
    call assemble_force_balance_residual( &
        volume_test, volume_force, volume_weights, boundary_test, boundary_force, &
        boundary_weights, sheet_test, sheet_force, sheet_weights, residual, status)
    call record_condition(status%code /= 0, "force-balance composition rejects nonpositive measures")
    call check_summary("force-balance residual")
    if (.not. all_passed) error stop 1

contains

    subroutine assemble_reference(reference)
        real(dp), intent(out) :: reference(:, :)
        integer :: sample, test_function, component

        reference = 0.0_dp
        do sample = 1, volume_count
            do test_function = 1, test_count
                do component = 1, component_count
                    reference(test_function, component) = &
                        reference(test_function, component) + volume_weights(sample)* &
                        volume_test(sample, test_function, component)* &
                        volume_force(sample, component)
                end do
            end do
        end do
        call add_reference(boundary_test, boundary_force, boundary_weights, reference)
        call add_reference(sheet_test, sheet_force, sheet_weights, reference)
    end subroutine assemble_reference

    subroutine add_reference(test, force, weights, reference)
        real(dp), intent(in) :: test(:, :, :), force(:, :), weights(:)
        real(dp), intent(inout) :: reference(:, :)
        integer :: sample, test_function, component

        do sample = 1, size(weights)
            do test_function = 1, test_count
                do component = 1, component_count
                    reference(test_function, component) = &
                        reference(test_function, component) + weights(sample)* &
                        test(sample, test_function, component)*force(sample, component)
                end do
            end do
        end do
    end subroutine add_reference

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_force_balance_residual
