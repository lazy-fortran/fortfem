program test_heaviside_enrichment
    use check, only: check_condition, check_summary
    use fortfem_api, only: evaluate_shifted_heaviside_enrichment, &
        evaluate_shifted_heaviside_enrichment_jvp, &
        evaluate_shifted_heaviside_enrichment_vjp
    use fortfem_kinds, only: dp
    use fortsparse, only: fortsparse_status_t
    implicit none

    real(dp), parameter :: level_values(4) = [1.0_dp, -2.0_dp, 3.0_dp, -4.0_dp]
    real(dp), parameter :: anchor_values(4) = [-1.0_dp, -2.0_dp, 4.0_dp, 5.0_dp]
    real(dp), parameter :: level_dot(4) = [0.2_dp, -0.3_dp, 0.1_dp, 0.4_dp]
    real(dp), parameter :: anchor_dot(4) = [-0.1_dp, 0.2_dp, 0.3_dp, -0.5_dp]
    real(dp), parameter :: values_bar(4) = [0.4_dp, -0.2_dp, 0.7_dp, 0.1_dp]
    real(dp), parameter :: expected_values(4) = [1.0_dp, 0.0_dp, 0.0_dp, -1.0_dp]
    real(dp) :: values(4), values_dot(4), level_bar(4), anchor_bar(4)
    real(dp) :: lhs, rhs
    type(fortsparse_status_t) :: status

    call evaluate_shifted_heaviside_enrichment( &
        level_values, anchor_values, values, status)
    call check_condition(status%code == 0 .and. &
        maxval(abs(values - expected_values)) < 1.0e-14_dp, &
        "shifted Heaviside enrichment matches the sign oracle")

    call evaluate_shifted_heaviside_enrichment_jvp( &
        level_values, anchor_values, level_dot, anchor_dot, values_dot, status)
    call check_condition(status%code == 0 .and. maxval(abs(values_dot)) < &
        1.0e-14_dp, "fixed-sign Heaviside JVP is zero away from topology events")

    call evaluate_shifted_heaviside_enrichment_vjp( &
        level_values, anchor_values, values_bar, level_bar, anchor_bar, status)
    lhs = sum(values_bar*values_dot)
    rhs = sum(level_bar*level_dot) + sum(anchor_bar*anchor_dot)
    call check_condition(status%code == 0 .and. abs(lhs - rhs) < 1.0e-14_dp .and. &
        maxval(abs(level_bar)) < 1.0e-14_dp .and. &
        maxval(abs(anchor_bar)) < 1.0e-14_dp, &
        "fixed-sign Heaviside VJP satisfies the zero derivative identity")

    call evaluate_shifted_heaviside_enrichment( &
        [0.0_dp, -1.0_dp], [-1.0_dp, 1.0_dp], values(1:2), status)
    call check_condition(status%code /= 0, &
        "Heaviside enrichment rejects a level-set topology event")

    call check_summary("shifted Heaviside enrichment")
end program test_heaviside_enrichment
