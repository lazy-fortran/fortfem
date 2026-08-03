program test_xfem_vector_blending_correction
    use check, only: check_condition, check_summary
    use fortfem_feec, only: evaluate_vector_blending_corrected_enrichment, &
        evaluate_vector_blending_corrected_enrichment_jvp, &
        evaluate_vector_blending_corrected_enrichment_vjp
    use fortfem_kinds, only: dp
    use fortsparse, only: fortsparse_status_t
    implicit none

    real(dp), parameter :: base_values(3, 4) = reshape([ &
        0.2_dp, 0.5_dp, 0.3_dp, 0.4_dp, 0.4_dp, 0.2_dp, &
        0.1_dp, 0.7_dp, 0.2_dp, 0.3_dp, 0.2_dp, 0.5_dp], [3, 4])
    logical, parameter :: enriched_mask(3) = [.false., .true., .false.]
    real(dp), parameter :: enrichment_values(2, 4) = reshape([ &
        2.0_dp, -0.5_dp, -1.0_dp, 0.25_dp, 0.5_dp, 1.5_dp, &
        3.0_dp, -2.0_dp], [2, 4])
    real(dp), parameter :: base_dot(3, 4) = reshape([ &
        0.0_dp, 0.1_dp, -0.2_dp, 0.2_dp, -0.2_dp, 0.1_dp, &
        -0.1_dp, 0.4_dp, 0.0_dp, 0.3_dp, -0.3_dp, 0.2_dp], [3, 4])
    real(dp), parameter :: enrichment_dot(2, 4) = reshape([ &
        0.05_dp, -0.1_dp, 0.2_dp, -0.05_dp, -0.1_dp, 0.3_dp, &
        0.4_dp, -0.2_dp], [2, 4])
    real(dp), parameter :: expected_values(2, 4) = reshape([ &
        1.0_dp, -0.25_dp, -0.4_dp, 0.1_dp, 0.35_dp, 1.05_dp, &
        0.6_dp, -0.4_dp], [2, 4])
    real(dp), parameter :: expected_dot(2, 4) = reshape([ &
        0.225_dp, -0.1_dp, 0.28_dp, -0.07_dp, 0.13_dp, 0.81_dp, &
        -0.82_dp, 0.56_dp], [2, 4])
    real(dp), parameter :: enriched_bar(2, 4) = reshape([ &
        0.4_dp, -0.2_dp, 0.7_dp, 0.1_dp, -0.3_dp, 0.5_dp, &
        -0.6_dp, 0.8_dp], [2, 4])
    real(dp) :: corrected_values(2, 4), corrected_dot(2, 4)
    real(dp) :: base_bar(3, 4), enrichment_bar(2, 4)
    real(dp) :: lhs, rhs
    type(fortsparse_status_t) :: status

    call evaluate_vector_blending_corrected_enrichment( &
        base_values, enriched_mask, enrichment_values, corrected_values, status)
    call check_condition(status%code == 0 .and. &
        maxval(abs(corrected_values - expected_values)) < 1.0e-14_dp, &
        "vector corrected enrichment matches the independent ramp oracle")

    call evaluate_vector_blending_corrected_enrichment_jvp( &
        base_values, enriched_mask, enrichment_values, base_dot, enrichment_dot, &
        corrected_dot, status)
    call check_condition(status%code == 0 .and. &
        maxval(abs(corrected_dot - expected_dot)) < 1.0e-14_dp, &
        "vector corrected enrichment JVP matches the componentwise product rule")

    call evaluate_vector_blending_corrected_enrichment_vjp( &
        base_values, enriched_mask, enrichment_values, enriched_bar, base_bar, &
        enrichment_bar, status)
    lhs = sum(enriched_bar*corrected_dot)
    rhs = sum(base_bar*base_dot) + sum(enrichment_bar*enrichment_dot)
    call check_condition(status%code == 0 .and. abs(lhs - rhs) < 1.0e-14_dp, &
        "vector corrected enrichment VJP satisfies the real dot-product identity")

    call check_summary("XFEM vector blending correction")
end program test_xfem_vector_blending_correction
