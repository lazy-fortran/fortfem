program test_surface_current_trace_residual
    use check, only: check_condition, check_summary
    use fortfem_boundary, only: &
        assemble_surface_current_trace_residual, &
        assemble_surface_current_trace_residual_jvp, &
        assemble_surface_current_trace_residual_vjp
    use fortfem_kinds, only: dp
    use fortsparse, only: fortsparse_status_t
    implicit none

    real(dp), parameter :: test_basis(2, 2, 3) = reshape([ &
        1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, &
        0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp], [2, 2, 3])
    real(dp), parameter :: trial_basis(2, 3, 3) = reshape([ &
        1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, &
        0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, &
        1.0_dp, -1.0_dp, 0.5_dp, 0.5_dp, 0.0_dp, -0.5_dp], [2, 3, 3])
    real(dp), parameter :: weights(2) = [2.0_dp, 3.0_dp]
    real(dp), parameter :: coefficients(3) = [0.5_dp, -1.0_dp, 2.0_dp]
    real(dp), parameter :: target_current(2, 3) = reshape([ &
        0.2_dp, 0.3_dp, -0.4_dp, 0.5_dp, -0.6_dp, 0.7_dp], [2, 3])
    real(dp), parameter :: test_basis_dot(2, 2, 3) = reshape([ &
        0.1_dp, -0.2_dp, 0.3_dp, -0.4_dp, 0.5_dp, -0.6_dp, &
        0.7_dp, -0.8_dp, 0.9_dp, -1.0_dp, 1.1_dp, -1.2_dp], [2, 2, 3])
    real(dp), parameter :: trial_basis_dot(2, 3, 3) = reshape([ &
        -0.1_dp, 0.2_dp, -0.3_dp, 0.4_dp, -0.5_dp, 0.6_dp, &
        -0.7_dp, 0.8_dp, -0.9_dp, 1.0_dp, -1.1_dp, 1.2_dp, &
        0.2_dp, -0.3_dp, 0.4_dp, -0.5_dp, 0.6_dp, -0.7_dp], [2, 3, 3])
    real(dp), parameter :: weights_dot(2) = [0.2_dp, -0.3_dp]
    real(dp), parameter :: coefficients_dot(3) = [-0.4_dp, 0.5_dp, 0.6_dp]
    real(dp), parameter :: target_current_dot(2, 3) = reshape([ &
        0.3_dp, -0.2_dp, 0.1_dp, -0.4_dp, 0.5_dp, -0.6_dp], [2, 3])
    real(dp), parameter :: residual_bar(2) = [0.8_dp, -0.7_dp]
    real(dp), parameter :: eps = 1.0e-6_dp
    real(dp) :: residual(2), residual_dot(2), residual_plus(2)
    real(dp) :: residual_minus(2), oracle_residual(2), oracle_dot(2)
    real(dp) :: test_basis_bar(2, 2, 3), trial_basis_bar(2, 3, 3)
    real(dp) :: weights_bar(2), coefficients_bar(3), target_current_bar(2, 3)
    real(dp) :: lhs, rhs
    real(dp) :: bad_residual(1)
    type(fortsparse_status_t) :: status
    integer :: quadrature, test_dof, trial_dof
    real(dp) :: current(3), current_dot(3), state(3), state_dot(3)

    oracle_residual = 0.0_dp
    do quadrature = 1, 2
        current = target_current(quadrature, :)
        do trial_dof = 1, 3
            current = current + trial_basis(quadrature, trial_dof, :)* &
                coefficients(trial_dof)
        end do
        do test_dof = 1, 2
            oracle_residual(test_dof) = oracle_residual(test_dof) + &
                weights(quadrature)*dot_product( &
                test_basis(quadrature, test_dof, :), &
                current - 2.0_dp*target_current(quadrature, :))
        end do
    end do
    ! The residual is K_h - target; the expression above starts from target
    ! so subtracting one more target gives the same independent state oracle.
    call assemble_surface_current_trace_residual( &
        test_basis, trial_basis, weights, coefficients, target_current, &
        residual, status)
    call check_condition(status%code == 0, &
        "surface-current trace residual accepts independent trace spaces")
    call check_condition(maxval(abs(residual - oracle_residual)) < 1.0e-14_dp, &
        "surface-current trace residual matches independent state oracle")

    oracle_dot = 0.0_dp
    do quadrature = 1, 2
        current = 0.0_dp
        current_dot = -target_current_dot(quadrature, :)
        do trial_dof = 1, 3
            current = current + trial_basis(quadrature, trial_dof, :)* &
                coefficients(trial_dof)
            current_dot = current_dot + &
                trial_basis_dot(quadrature, trial_dof, :)*coefficients(trial_dof) + &
                trial_basis(quadrature, trial_dof, :)*coefficients_dot(trial_dof)
        end do
        state = current - target_current(quadrature, :)
        state_dot = current_dot
        do test_dof = 1, 2
            oracle_dot(test_dof) = oracle_dot(test_dof) + &
                weights_dot(quadrature)*dot_product( &
                test_basis(quadrature, test_dof, :), state) + &
                weights(quadrature)*dot_product( &
                test_basis_dot(quadrature, test_dof, :), state) + &
                weights(quadrature)*dot_product( &
                test_basis(quadrature, test_dof, :), state_dot)
        end do
    end do
    call assemble_surface_current_trace_residual_jvp( &
        test_basis, trial_basis, weights, coefficients, target_current, &
        test_basis_dot, trial_basis_dot, weights_dot, coefficients_dot, &
        target_current_dot, residual_dot, status)
    call assemble_surface_current_trace_residual( &
        test_basis + eps*test_basis_dot, trial_basis + eps*trial_basis_dot, &
        weights + eps*weights_dot, coefficients + eps*coefficients_dot, &
        target_current + eps*target_current_dot, residual_plus, status)
    call assemble_surface_current_trace_residual( &
        test_basis - eps*test_basis_dot, trial_basis - eps*trial_basis_dot, &
        weights - eps*weights_dot, coefficients - eps*coefficients_dot, &
        target_current - eps*target_current_dot, residual_minus, status)
    call check_condition(status%code == 0 .and. &
        maxval(abs(residual - oracle_residual)) < 1.0e-14_dp .and. &
        maxval(abs(residual_dot - oracle_dot)) < 1.0e-14_dp .and. &
        maxval(abs((residual_plus - residual_minus)/(2.0_dp*eps) - &
        residual_dot)) < 1.0e-8_dp, &
        "surface-current trace JVP matches product rule and finite differences")

    call assemble_surface_current_trace_residual_vjp( &
        test_basis, trial_basis, weights, coefficients, target_current, &
        residual_bar, test_basis_bar, trial_basis_bar, weights_bar, &
        coefficients_bar, target_current_bar, status)
    lhs = dot_product(residual_bar, residual_dot)
    rhs = sum(test_basis_bar*test_basis_dot) + &
        sum(trial_basis_bar*trial_basis_dot) + &
        dot_product(weights_bar, weights_dot) + &
        dot_product(coefficients_bar, coefficients_dot) + &
        sum(target_current_bar*target_current_dot)
    call check_condition(status%code == 0 .and. abs(lhs - rhs) < 1.0e-12_dp, &
        "surface-current trace VJP satisfies the real dot-product identity")

    call assemble_surface_current_trace_residual( &
        test_basis, trial_basis, weights, coefficients, target_current, &
        bad_residual, status)
    call check_condition(status%code /= 0, &
        "surface-current trace residual rejects incompatible output size")

    call check_summary("surface current trace residual")
end program test_surface_current_trace_residual
