program test_free_boundary_port_residual
    use check, only: check_condition, check_summary
    use fortfem_api, only: assemble_free_boundary_port_residual, &
        assemble_free_boundary_port_residual_jvp, &
        assemble_free_boundary_port_residual_vjp
    use fortfem_kinds, only: dp
    use fortsparse, only: fortsparse_status_t
    implicit none

    integer, parameter :: quadrature_count = 3, component_count = 2
    real(dp), parameter :: epsilon = 1.0e-7_dp
    real(dp) :: physical_trace(quadrature_count, component_count)
    real(dp) :: external_target(quadrature_count, component_count)
    real(dp) :: sheet_current(quadrature_count, component_count)
    real(dp) :: weights(quadrature_count)
    real(dp) :: physical_trace_dot(quadrature_count, component_count)
    real(dp) :: external_target_dot(quadrature_count, component_count)
    real(dp) :: sheet_current_dot(quadrature_count, component_count)
    real(dp) :: weights_dot(quadrature_count)
    real(dp) :: residual(quadrature_count, component_count)
    real(dp) :: residual_without_sheet(quadrature_count, component_count)
    real(dp) :: residual_dot(quadrature_count, component_count)
    real(dp) :: residual_plus(quadrature_count, component_count)
    real(dp) :: residual_minus(quadrature_count, component_count)
    real(dp) :: residual_bar(quadrature_count, component_count)
    real(dp) :: physical_trace_bar(quadrature_count, component_count)
    real(dp) :: external_target_bar(quadrature_count, component_count)
    real(dp) :: sheet_current_bar(quadrature_count, component_count)
    real(dp) :: weights_bar(quadrature_count)
    real(dp) :: reference(quadrature_count, component_count)
    real(dp) :: reference_without_sheet(quadrature_count, component_count)
    real(dp) :: reference_dot(quadrature_count, component_count)
    real(dp) :: lhs, rhs
    real(dp) :: invalid_residual(1, component_count)
    real(dp) :: zero_trace(0, component_count), zero_target(0, component_count)
    real(dp) :: zero_weights(0), zero_residual(0, component_count)
    real(dp) :: negative_weights(quadrature_count)
    type(fortsparse_status_t) :: status
    integer :: quadrature, component

    physical_trace = reshape([ &
        1.2_dp, -0.7_dp, 0.3_dp, 1.1_dp, -0.4_dp, 0.8_dp], &
        shape(physical_trace))
    external_target = reshape([ &
        0.2_dp, -0.1_dp, 0.5_dp, 0.6_dp, -0.9_dp, 0.4_dp], &
        shape(external_target))
    sheet_current = reshape([ &
        0.05_dp, 0.03_dp, -0.2_dp, 0.1_dp, 0.4_dp, -0.15_dp], &
        shape(sheet_current))
    weights = [1.5_dp, 0.7_dp, 2.1_dp]
    physical_trace_dot = reshape([ &
        -0.2_dp, 0.4_dp, 0.3_dp, -0.5_dp, 0.6_dp, 0.1_dp], &
        shape(physical_trace_dot))
    external_target_dot = reshape([ &
        0.1_dp, -0.3_dp, 0.2_dp, 0.5_dp, -0.4_dp, 0.6_dp], &
        shape(external_target_dot))
    sheet_current_dot = reshape([ &
        -0.05_dp, 0.08_dp, 0.12_dp, -0.07_dp, 0.03_dp, 0.11_dp], &
        shape(sheet_current_dot))
    weights_dot = [0.2_dp, -0.1_dp, 0.3_dp]
    residual_bar = reshape([0.4_dp, -0.6_dp, 0.7_dp, 0.2_dp, -0.3_dp, 0.8_dp], &
        shape(residual_bar))

    reference = 0.0_dp
    reference_without_sheet = 0.0_dp
    reference_dot = 0.0_dp
    do quadrature = 1, quadrature_count
        do component = 1, component_count
            reference(quadrature, component) = weights(quadrature)* &
                (physical_trace(quadrature, component) - &
                external_target(quadrature, component) - &
                sheet_current(quadrature, component))
            reference_without_sheet(quadrature, component) = weights(quadrature)* &
                (physical_trace(quadrature, component) - &
                external_target(quadrature, component))
            reference_dot(quadrature, component) = weights_dot(quadrature)* &
                (physical_trace(quadrature, component) - &
                external_target(quadrature, component) - &
                sheet_current(quadrature, component)) + weights(quadrature)* &
                (physical_trace_dot(quadrature, component) - &
                external_target_dot(quadrature, component) - &
                sheet_current_dot(quadrature, component))
        end do
    end do

    call assemble_free_boundary_port_residual( &
        physical_trace, external_target, weights, residual, status, sheet_current)
    call check_condition(status%code == 0 .and. &
        maxval(abs(residual - reference)) < 1.0e-14_dp, &
        "free-boundary port matches independent physical trace oracle")

    call assemble_free_boundary_port_residual( &
        physical_trace, external_target, weights, residual_without_sheet, status)
    call check_condition(status%code == 0 .and. &
        maxval(abs(residual_without_sheet - reference_without_sheet)) < 1.0e-14_dp, &
        "free-boundary port accepts an omitted sheet-current target")

    call assemble_free_boundary_port_residual_jvp( &
        physical_trace, external_target, weights, physical_trace_dot, &
        external_target_dot, weights_dot, residual_dot, status, sheet_current, &
        sheet_current_dot)
    call assemble_free_boundary_port_residual( &
        physical_trace + epsilon*physical_trace_dot, &
        external_target + epsilon*external_target_dot, weights + epsilon*weights_dot, &
        residual_plus, status, sheet_current + epsilon*sheet_current_dot)
    call assemble_free_boundary_port_residual( &
        physical_trace - epsilon*physical_trace_dot, &
        external_target - epsilon*external_target_dot, weights - epsilon*weights_dot, &
        residual_minus, status, sheet_current - epsilon*sheet_current_dot)
    call check_condition(status%code == 0 .and. &
        maxval(abs(residual_dot - reference_dot)) < 1.0e-14_dp .and. &
        maxval(abs(residual_dot - (residual_plus - residual_minus) / &
        (2.0_dp*epsilon))) < 2.0e-8_dp, &
        "free-boundary port JVP matches oracle and central difference")

    call assemble_free_boundary_port_residual_vjp( &
        physical_trace, external_target, weights, residual_bar, &
        physical_trace_bar, external_target_bar, weights_bar, status, &
        sheet_current, sheet_current_bar)
    lhs = sum(residual_bar*residual_dot)
    rhs = sum(physical_trace_bar*physical_trace_dot) + &
        sum(external_target_bar*external_target_dot) + &
        dot_product(weights_bar, weights_dot) + &
        sum(sheet_current_bar*sheet_current_dot)
    call check_condition(status%code == 0 .and. abs(lhs - rhs) < 1.0e-12_dp, &
        "free-boundary port VJP satisfies the real dot-product identity")

    negative_weights = weights
    negative_weights(2) = -negative_weights(2)
    call assemble_free_boundary_port_residual( &
        physical_trace, external_target, negative_weights, residual, status)
    call check_condition(status%code /= 0, &
        "free-boundary port rejects nonpositive work weights")

    call assemble_free_boundary_port_residual( &
        zero_trace, zero_target, zero_weights, zero_residual, status)
    call check_condition(status%code /= 0, &
        "free-boundary port rejects a zero-row physical trace")

    call assemble_free_boundary_port_residual( &
        physical_trace, external_target, weights, invalid_residual, status)
    call check_condition(status%code /= 0, &
        "free-boundary port rejects an incompatible residual shape")

    call assemble_free_boundary_port_residual_jvp( &
        physical_trace, external_target, weights, physical_trace_dot, &
        external_target_dot, weights_dot, residual_dot, status, sheet_current)
    call check_condition(status%code /= 0, &
        "free-boundary port JVP rejects a missing sheet-current increment")

    call assemble_free_boundary_port_residual_vjp( &
        physical_trace, external_target, weights, residual_bar, &
        physical_trace_bar, external_target_bar, weights_bar, status, sheet_current)
    call check_condition(status%code /= 0, &
        "free-boundary port VJP rejects a missing sheet-current cotangent")

    call check_summary("free-boundary port residual")
end program test_free_boundary_port_residual
