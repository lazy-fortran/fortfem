program test_force_balance_validation
    !! Independent contract tests for force-balance array shape handling.
    use check, only: check_condition, check_summary
    use fortfem_api, only: assemble_force_balance_residual, &
        assemble_force_balance_residual_jvp, assemble_force_balance_residual_vjp
    use fortfem_kinds, only: dp
    use fortsparse, only: fortsparse_status_t
    implicit none

    real(dp) :: volume_test(1, 1, 2), volume_force(1, 2), volume_weights(1)
    real(dp) :: boundary_test(1, 1, 2), boundary_force(1, 2), boundary_weights(1)
    real(dp) :: sheet_test(1, 1, 2), sheet_force(1, 2), sheet_weights(1)
    real(dp) :: empty_test(0, 1, 2), empty_force(0, 2), empty_weights(0)
    real(dp) :: bad_force(1, 3), bad_weights(2)
    real(dp) :: volume_test_dot(1, 1, 2), volume_force_dot(1, 2)
    real(dp) :: volume_weights_dot(1)
    real(dp) :: boundary_test_dot(1, 1, 2), boundary_force_dot(1, 2)
    real(dp) :: boundary_weights_dot(1)
    real(dp) :: sheet_test_dot(1, 1, 2), sheet_force_dot(1, 2)
    real(dp) :: sheet_weights_dot(1)
    real(dp) :: bad_boundary_test_dot(2, 1, 2)
    real(dp) :: residual(1, 2), residual_dot(1, 2)
    real(dp) :: residual_bar(1, 2), bad_residual_bar(2, 2)
    real(dp) :: volume_test_bar(1, 1, 2), volume_force_bar(1, 2)
    real(dp) :: volume_weights_bar(1)
    real(dp) :: boundary_test_bar(1, 1, 2), boundary_force_bar(1, 2)
    real(dp) :: boundary_weights_bar(1)
    real(dp) :: sheet_test_bar(1, 1, 2), sheet_force_bar(1, 2)
    real(dp) :: sheet_weights_bar(1)
    type(fortsparse_status_t) :: status
    logical :: all_passed

    all_passed = .true.
    volume_test = reshape([1.0_dp, 2.0_dp], shape(volume_test))
    volume_force = reshape([3.0_dp, 4.0_dp], shape(volume_force))
    volume_weights = [1.0_dp]
    boundary_test = reshape([0.5_dp, -1.0_dp], shape(boundary_test))
    boundary_force = reshape([2.0_dp, -3.0_dp], shape(boundary_force))
    boundary_weights = [0.75_dp]
    sheet_test = reshape([-0.25_dp, 0.8_dp], shape(sheet_test))
    sheet_force = reshape([1.5_dp, 0.25_dp], shape(sheet_force))
    sheet_weights = [0.5_dp]
    residual_bar = reshape([0.3_dp, -0.4_dp], shape(residual_bar))
    volume_test_dot = 0.1_dp
    volume_force_dot = -0.2_dp
    volume_weights_dot = [0.05_dp]
    boundary_test_dot = -0.3_dp
    boundary_force_dot = 0.4_dp
    boundary_weights_dot = [0.1_dp]
    sheet_test_dot = 0.2_dp
    sheet_force_dot = -0.15_dp
    sheet_weights_dot = [-0.08_dp]

    call assemble_force_balance_residual( &
        volume_test, volume_force, volume_weights, empty_test, empty_force, &
        empty_weights, empty_test, empty_force, empty_weights, residual, status)
    call record_condition(status%code == 0 .and. &
        maxval(abs(residual - reshape([3.0_dp, 8.0_dp], shape(residual)))) < &
        1.0e-14_dp, &
        "zero-row boundary and sheet terms are accepted and ignored")

    call assemble_force_balance_residual( &
        volume_test, volume_force, volume_weights, boundary_test, bad_force, &
        boundary_weights, sheet_test, sheet_force, sheet_weights, residual, status)
    call record_condition(status%code /= 0, &
        "force balance rejects a force component-shape mismatch")

    call assemble_force_balance_residual( &
        volume_test, volume_force, volume_weights, boundary_test, boundary_force, &
        bad_weights, sheet_test, sheet_force, sheet_weights, residual, status)
    call record_condition(status%code /= 0, &
        "force balance rejects a quadrature-measure length mismatch")

    call assemble_force_balance_residual_jvp( &
        volume_test, volume_force, volume_weights, boundary_test, boundary_force, &
        boundary_weights, sheet_test, sheet_force, sheet_weights, volume_test_dot, &
        volume_force_dot, volume_weights_dot, bad_boundary_test_dot, &
        boundary_force_dot, boundary_weights_dot, sheet_test_dot, sheet_force_dot, &
        sheet_weights_dot, residual_dot, status)
    call record_condition(status%code /= 0, &
        "force balance JVP rejects an increment shape mismatch")

    bad_residual_bar = 1.0_dp
    call assemble_force_balance_residual_vjp( &
        volume_test, volume_force, volume_weights, boundary_test, boundary_force, &
        boundary_weights, sheet_test, sheet_force, sheet_weights, bad_residual_bar, &
        volume_test_bar, volume_force_bar, volume_weights_bar, boundary_test_bar, &
        boundary_force_bar, boundary_weights_bar, sheet_test_bar, sheet_force_bar, &
        sheet_weights_bar, status)
    call record_condition(status%code /= 0, &
        "force balance VJP rejects a residual-cotangent shape mismatch")

    call assemble_force_balance_residual_vjp( &
        volume_test, volume_force, volume_weights, boundary_test, boundary_force, &
        boundary_weights, sheet_test, sheet_force, sheet_weights, residual_bar, &
        volume_test_bar, volume_force_bar, volume_weights_bar, boundary_test_bar, &
        boundary_force_bar, boundary_weights_bar, sheet_test_bar, sheet_force_bar, &
        sheet_weights_bar, status)
    call record_condition(status%code == 0, &
        "force balance VJP accepts compatible output shapes")

    call check_summary("force-balance validation")
    if (.not. all_passed) error stop 1

contains

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_force_balance_validation
