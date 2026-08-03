program test_feec_exact_sequence_ad
    use check, only: check_condition, check_summary
    use fortfem_feec, only: assemble_feec_exact_sequence, &
        assemble_feec_exact_sequence_jvp, assemble_feec_exact_sequence_vjp
    use fortfem_kinds, only: dp
    use fortsparse, only: fortsparse_status_t
    implicit none

    real(dp), parameter :: gradient(3, 3) = reshape([ &
        -1.0_dp, 0.0_dp, 1.0_dp, 1.0_dp, -1.0_dp, 0.0_dp, &
        0.0_dp, 1.0_dp, -1.0_dp], [3, 3])
    real(dp), parameter :: curl(2, 3) = reshape([1.0_dp, 2.0_dp, 1.0_dp, &
        2.0_dp, 1.0_dp, 2.0_dp], [2, 3])
    real(dp), parameter :: divergence(1, 2) = reshape([1.0_dp, -0.5_dp], [1, 2])
    real(dp), parameter :: gradient_dot(3, 3) = reshape([ &
        0.1_dp, -0.2_dp, 0.3_dp, -0.4_dp, 0.5_dp, -0.6_dp, &
        0.7_dp, -0.8_dp, 0.9_dp], [3, 3])
    real(dp), parameter :: curl_dot(2, 3) = reshape([ &
        -0.2_dp, 0.3_dp, -0.4_dp, 0.5_dp, -0.6_dp, 0.7_dp], [2, 3])
    real(dp), parameter :: divergence_dot(1, 2) = reshape([0.8_dp, -0.9_dp], [1, 2])
    real(dp), parameter :: curl_gradient_bar(2, 3) = reshape([ &
        0.4_dp, -0.3_dp, 0.2_dp, -0.1_dp, 0.6_dp, -0.5_dp], [2, 3])
    real(dp), parameter :: divergence_curl_bar(1, 3) = reshape([0.7_dp, -0.8_dp, 0.9_dp], [1, 3])
    real(dp) :: curl_gradient(2, 3), divergence_curl(1, 3)
    real(dp) :: curl_gradient_dot(2, 3), divergence_curl_dot(1, 3)
    real(dp) :: gradient_bar(3, 3), curl_bar(2, 3), divergence_bar(1, 2)
    real(dp) :: curl_gradient_expected(2, 3), divergence_curl_expected(1, 3)
    real(dp) :: lhs, rhs
    type(fortsparse_status_t) :: status

    call assemble_feec_exact_sequence( &
        gradient, curl, divergence, curl_gradient, divergence_curl, status)
    call check_condition(status%code == 0 .and. &
        maxval(abs(curl_gradient)) < 1.0e-14_dp .and. &
        maxval(abs(divergence_curl)) < 1.0e-14_dp, &
        "FEEC incidence maps satisfy the independent exact-sequence oracle")

    call assemble_feec_exact_sequence_jvp( &
        gradient, curl, divergence, gradient_dot, curl_dot, divergence_dot, &
        curl_gradient_dot, divergence_curl_dot, status)
    curl_gradient_expected = matmul(curl_dot, gradient) + matmul(curl, gradient_dot)
    divergence_curl_expected = matmul(divergence_dot, curl) + matmul(divergence, curl_dot)
    call check_condition(status%code == 0 .and. &
        maxval(abs(curl_gradient_dot - curl_gradient_expected)) < 1.0e-14_dp .and. &
        maxval(abs(divergence_curl_dot - divergence_curl_expected)) < 1.0e-14_dp, &
        "FEEC exact-sequence JVP differentiates both incidence compositions")

    call assemble_feec_exact_sequence_vjp( &
        gradient, curl, divergence, curl_gradient_bar, divergence_curl_bar, &
        gradient_bar, curl_bar, divergence_bar, status)
    lhs = sum(curl_gradient_bar*curl_gradient_dot) + &
        sum(divergence_curl_bar*divergence_curl_dot)
    rhs = sum(gradient_bar*gradient_dot) + sum(curl_bar*curl_dot) + &
        sum(divergence_bar*divergence_dot)
    call check_condition(status%code == 0 .and. abs(lhs - rhs) < 1.0e-13_dp, &
        "FEEC exact-sequence VJP satisfies the real dot-product identity")

    call assemble_feec_exact_sequence( &
        gradient, curl(:, 1:2), divergence, curl_gradient, divergence_curl, status)
    call check_condition(status%code /= 0, &
        "FEEC exact-sequence diagnostic rejects incompatible incidence dimensions")
    call check_summary("FEEC exact sequence AD")
end program test_feec_exact_sequence_ad
