program test_total_pressure_balance
    use check, only: check_condition, check_summary
    use fortfem_api, only: &
        assemble_total_pressure_jump, assemble_total_pressure_jump_jvp, &
        assemble_total_pressure_jump_vjp
    use fortfem_kinds, only: dp
    use fortsparse, only: FORTSPARSE_OK, fortsparse_status_t
    implicit none

    real(dp), parameter :: epsilon = 1.0e-7_dp
    real(dp), parameter :: tolerance = 2.0e-8_dp
    real(dp) :: plus_pressure, minus_pressure, plus_field(3), minus_field(3)
    real(dp) :: permeability, target, residual, residual_plus, residual_minus
    real(dp) :: plus_pressure_dot, minus_pressure_dot, plus_field_dot(3)
    real(dp) :: minus_field_dot(3), permeability_dot, target_dot, residual_dot
    real(dp) :: plus_pressure_bar, minus_pressure_bar, plus_field_bar(3)
    real(dp) :: minus_field_bar(3), permeability_bar, target_bar
    real(dp) :: primal_dot, adjoint_dot, expected
    type(fortsparse_status_t) :: status

    plus_pressure = 1.7_dp
    minus_pressure = 0.9_dp
    plus_field = [0.4_dp, -0.7_dp, 1.1_dp]
    minus_field = [-0.2_dp, 0.3_dp, 0.8_dp]
    permeability = 2.5_dp
    target = 0.15_dp
    plus_pressure_dot = -0.2_dp
    minus_pressure_dot = 0.35_dp
    plus_field_dot = [0.1_dp, 0.2_dp, -0.3_dp]
    minus_field_dot = [-0.4_dp, 0.05_dp, 0.25_dp]
    permeability_dot = 0.12_dp
    target_dot = -0.08_dp

    call assemble_total_pressure_jump(plus_pressure, minus_pressure, plus_field, &
        minus_field, permeability, target, residual, status)
    call check_condition(status%code == FORTSPARSE_OK, &
        "total-pressure jump accepts finite traces")
    expected = plus_pressure + 0.5_dp*dot_product(plus_field, plus_field)/ &
        permeability - minus_pressure - 0.5_dp*dot_product(minus_field, minus_field)/ &
        permeability - target
    call check_condition(abs(residual - expected) < tolerance, &
        "total-pressure jump matches independent scalar oracle")

    call assemble_total_pressure_jump_jvp(plus_pressure, minus_pressure, plus_field, &
        minus_field, permeability, target, plus_pressure_dot, minus_pressure_dot, &
        plus_field_dot, minus_field_dot, permeability_dot, target_dot, &
        residual_dot, status)
    call check_condition(status%code == FORTSPARSE_OK, &
        "total-pressure jump JVP accepts finite directions")
    call assemble_total_pressure_jump(plus_pressure + epsilon*plus_pressure_dot, &
        minus_pressure + epsilon*minus_pressure_dot, &
        plus_field + epsilon*plus_field_dot, &
        minus_field + epsilon*minus_field_dot, &
        permeability + epsilon*permeability_dot, &
        target + epsilon*target_dot, residual_plus, status)
    call assemble_total_pressure_jump(plus_pressure - epsilon*plus_pressure_dot, &
        minus_pressure - epsilon*minus_pressure_dot, &
        plus_field - epsilon*plus_field_dot, &
        minus_field - epsilon*minus_field_dot, &
        permeability - epsilon*permeability_dot, &
        target - epsilon*target_dot, residual_minus, status)
    call check_condition(abs(residual_dot - (residual_plus - residual_minus)/ &
        (2.0_dp*epsilon)) < tolerance, &
        "total-pressure jump JVP matches finite difference")

    call assemble_total_pressure_jump_vjp(plus_pressure, minus_pressure, plus_field, &
        minus_field, permeability, target, 1.3_dp, plus_pressure_bar, &
        minus_pressure_bar, &
        plus_field_bar, minus_field_bar, permeability_bar, target_bar, status)
    call check_condition(status%code == FORTSPARSE_OK, &
        "total-pressure jump VJP accepts finite cotangent")
    primal_dot = 1.3_dp*(plus_pressure_dot - minus_pressure_dot + &
        dot_product(plus_field, plus_field_dot)/permeability - &
        dot_product(minus_field, minus_field_dot)/permeability - &
        0.5_dp*(dot_product(plus_field, plus_field) - &
        dot_product(minus_field, minus_field))* &
        permeability_dot/permeability**2 - target_dot)
    adjoint_dot = plus_pressure_bar*plus_pressure_dot + &
        minus_pressure_bar*minus_pressure_dot + &
        dot_product(plus_field_bar, plus_field_dot) + &
        dot_product(minus_field_bar, minus_field_dot) + &
        permeability_bar*permeability_dot + target_bar*target_dot
    call check_condition(abs(primal_dot - adjoint_dot) < tolerance, &
        "total-pressure jump VJP satisfies real adjoint identity")

    call assemble_total_pressure_jump(plus_pressure, minus_pressure, plus_field, &
        minus_field, 0.0_dp, target, residual, status)
    call check_condition(status%code /= FORTSPARSE_OK, &
        "total-pressure jump rejects non-positive permeability")
    call check_summary("total-pressure balance")
end program test_total_pressure_balance
