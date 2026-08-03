program test_enriched_feec_sequence
    use check, only: check_condition, check_summary
    use fortfem_feec, only: assemble_enriched_feec_sequence, &
        assemble_enriched_feec_sequence_jvp, assemble_enriched_feec_sequence_vjp
    use fortfem_kinds, only: dp
    use fortsparse, only: fortsparse_status_t
    implicit none

    integer, parameter :: n0 = 2, n1 = 3, n2 = 2, n3 = 1
    integer, parameter :: e0 = 3, e1 = 2, e2 = 3, e3 = 2
    real(dp), parameter :: gradient(n1, n0) = reshape([ &
        1.0_dp, -2.0_dp, 0.5_dp, 0.25_dp, 1.5_dp, -0.75_dp], [n1, n0])
    real(dp), parameter :: curl(n2, n1) = reshape([ &
        0.2_dp, -0.4_dp, 0.8_dp, 1.1_dp, -0.3_dp, 0.6_dp], [n2, n1])
    real(dp), parameter :: divergence(n3, n2) = reshape([0.7_dp, -1.2_dp], [n3, n2])
    real(dp), parameter :: scalar_map(n0, e0) = reshape([ &
        1.0_dp, -0.2_dp, 0.3_dp, 0.8_dp, -0.7_dp, 1.1_dp], [n0, e0])
    real(dp), parameter :: hcurl_map(n1, e1) = reshape([ &
        0.4_dp, 1.2_dp, -0.8_dp, -0.6_dp, 0.5_dp, 0.9_dp], [n1, e1])
    real(dp), parameter :: hdiv_map(n2, e2) = reshape([ &
        1.1_dp, -0.5_dp, 0.2_dp, 0.7_dp, -0.9_dp, 0.3_dp], [n2, e2])
    real(dp), parameter :: l2_map(n3, e3) = reshape([0.6_dp, -1.4_dp], [n3, e3])
    real(dp), parameter :: gradient_dot(n1, n0) = reshape([ &
        -0.3_dp, 0.5_dp, 0.7_dp, -0.9_dp, 1.1_dp, -1.3_dp], [n1, n0])
    real(dp), parameter :: curl_dot(n2, n1) = reshape([ &
        0.4_dp, -0.6_dp, 0.8_dp, -1.0_dp, 1.2_dp, -1.4_dp], [n2, n1])
    real(dp), parameter :: divergence_dot(n3, n2) = reshape([0.9_dp, -0.7_dp], [n3, n2])
    real(dp), parameter :: scalar_map_dot(n0, e0) = reshape([ &
        0.2_dp, -0.4_dp, 0.6_dp, -0.8_dp, 1.0_dp, -1.2_dp], [n0, e0])
    real(dp), parameter :: hcurl_map_dot(n1, e1) = reshape([ &
        -0.1_dp, 0.3_dp, -0.5_dp, 0.7_dp, -0.9_dp, 1.1_dp], [n1, e1])
    real(dp), parameter :: hdiv_map_dot(n2, e2) = reshape([ &
        0.5_dp, -0.7_dp, 0.9_dp, -1.1_dp, 1.3_dp, -1.5_dp], [n2, e2])
    real(dp), parameter :: l2_map_dot(n3, e3) = reshape([0.8_dp, -1.0_dp], [n3, e3])
    real(dp), parameter :: enriched_gradient_bar(e1, e0) = reshape([ &
        0.2_dp, -0.4_dp, 0.6_dp, -0.8_dp, 1.0_dp, -1.2_dp], [e1, e0])
    real(dp), parameter :: enriched_curl_bar(e2, e1) = reshape([ &
        -0.3_dp, 0.5_dp, -0.7_dp, 0.9_dp, -1.1_dp, 1.3_dp], [e2, e1])
    real(dp), parameter :: enriched_divergence_bar(e3, e2) = reshape([ &
        0.4_dp, -0.6_dp, 0.8_dp, -1.0_dp, 1.2_dp, -1.4_dp], [e3, e2])
    real(dp), parameter :: curl_gradient_bar(e2, e0) = reshape([ &
        0.7_dp, -0.9_dp, 1.1_dp, -1.3_dp, 1.5_dp, -1.7_dp, &
        0.3_dp, -0.5_dp, 0.7_dp], [e2, e0])
    real(dp), parameter :: divergence_curl_bar(e3, e1) = reshape([ &
        -0.8_dp, 1.0_dp, -1.2_dp, 1.4_dp], [e3, e1])
    real(dp) :: enriched_gradient(e1, e0), enriched_curl(e2, e1)
    real(dp) :: enriched_divergence(e3, e2), curl_gradient(e2, e0)
    real(dp) :: divergence_curl(e3, e1)
    real(dp) :: enriched_gradient_dot(e1, e0), enriched_curl_dot(e2, e1)
    real(dp) :: enriched_divergence_dot(e3, e2), curl_gradient_dot(e2, e0)
    real(dp) :: divergence_curl_dot(e3, e1)
    real(dp) :: gradient_bar(n1, n0), curl_bar(n2, n1), divergence_bar(n3, n2)
    real(dp) :: scalar_map_bar(n0, e0), hcurl_map_bar(n1, e1)
    real(dp) :: hdiv_map_bar(n2, e2), l2_map_bar(n3, e3)
    real(dp) :: expected_gradient(e1, e0), expected_curl(e2, e1)
    real(dp) :: expected_divergence(e3, e2), expected_curl_gradient(e2, e0)
    real(dp) :: expected_divergence_curl(e3, e1), lhs, rhs
    type(fortsparse_status_t) :: status

    call assemble_enriched_feec_sequence( &
        gradient, curl, divergence, scalar_map, hcurl_map, hdiv_map, l2_map, &
        enriched_gradient, enriched_curl, enriched_divergence, curl_gradient, &
        divergence_curl, status)
    expected_gradient = matmul(transpose(hcurl_map), matmul(gradient, scalar_map))
    expected_curl = matmul(transpose(hdiv_map), matmul(curl, hcurl_map))
    expected_divergence = matmul(transpose(l2_map), matmul(divergence, hdiv_map))
    expected_curl_gradient = matmul(expected_curl, expected_gradient)
    expected_divergence_curl = matmul(expected_divergence, expected_curl)
    call check_condition(status%code == 0 .and. &
        maxval(abs(enriched_gradient - expected_gradient)) < 1.0e-14_dp .and. &
        maxval(abs(enriched_curl - expected_curl)) < 1.0e-14_dp .and. &
        maxval(abs(enriched_divergence - expected_divergence)) < 1.0e-14_dp .and. &
        maxval(abs(curl_gradient - expected_curl_gradient)) < 1.0e-14_dp .and. &
        maxval(abs(divergence_curl - expected_divergence_curl)) < 1.0e-14_dp, &
        "enriched FEEC value matches an independent dense matrix oracle")

    call assemble_enriched_feec_sequence_jvp( &
        gradient, curl, divergence, scalar_map, hcurl_map, hdiv_map, l2_map, &
        gradient_dot, curl_dot, divergence_dot, scalar_map_dot, hcurl_map_dot, &
        hdiv_map_dot, l2_map_dot, enriched_gradient_dot, enriched_curl_dot, &
        enriched_divergence_dot, curl_gradient_dot, divergence_curl_dot, status)
    expected_gradient = matmul(transpose(hcurl_map_dot), matmul(gradient, scalar_map)) + &
        matmul(transpose(hcurl_map), matmul(gradient_dot, scalar_map)) + &
        matmul(transpose(hcurl_map), matmul(gradient, scalar_map_dot))
    expected_curl = matmul(transpose(hdiv_map_dot), matmul(curl, hcurl_map)) + &
        matmul(transpose(hdiv_map), matmul(curl_dot, hcurl_map)) + &
        matmul(transpose(hdiv_map), matmul(curl, hcurl_map_dot))
    expected_divergence = matmul(transpose(l2_map_dot), matmul(divergence, hdiv_map)) + &
        matmul(transpose(l2_map), matmul(divergence_dot, hdiv_map)) + &
        matmul(transpose(l2_map), matmul(divergence, hdiv_map_dot))
    expected_curl_gradient = matmul(expected_curl, &
        matmul(transpose(hcurl_map), matmul(gradient, scalar_map))) + &
        matmul(matmul(transpose(hdiv_map), matmul(curl, hcurl_map)), expected_gradient)
    expected_divergence_curl = matmul(expected_divergence, &
        matmul(transpose(hdiv_map), matmul(curl, hcurl_map))) + &
        matmul(matmul(transpose(l2_map), matmul(divergence, hdiv_map)), expected_curl)
    call check_condition(status%code == 0 .and. &
        maxval(abs(enriched_gradient_dot - expected_gradient)) < 1.0e-13_dp .and. &
        maxval(abs(enriched_curl_dot - expected_curl)) < 1.0e-13_dp .and. &
        maxval(abs(enriched_divergence_dot - expected_divergence)) < 1.0e-13_dp .and. &
        maxval(abs(curl_gradient_dot - expected_curl_gradient)) < 1.0e-13_dp .and. &
        maxval(abs(divergence_curl_dot - expected_divergence_curl)) < 1.0e-13_dp, &
        "enriched FEEC JVP matches the independent product-rule oracle")

    call assemble_enriched_feec_sequence_vjp( &
        gradient, curl, divergence, scalar_map, hcurl_map, hdiv_map, l2_map, &
        enriched_gradient_bar, enriched_curl_bar, enriched_divergence_bar, &
        curl_gradient_bar, divergence_curl_bar, gradient_bar, curl_bar, &
        divergence_bar, scalar_map_bar, hcurl_map_bar, hdiv_map_bar, l2_map_bar, &
        status)
    lhs = sum(enriched_gradient_bar*enriched_gradient_dot) + &
        sum(enriched_curl_bar*enriched_curl_dot) + &
        sum(enriched_divergence_bar*enriched_divergence_dot) + &
        sum(curl_gradient_bar*curl_gradient_dot) + &
        sum(divergence_curl_bar*divergence_curl_dot)
    rhs = sum(gradient_bar*gradient_dot) + sum(curl_bar*curl_dot) + &
        sum(divergence_bar*divergence_dot) + sum(scalar_map_bar*scalar_map_dot) + &
        sum(hcurl_map_bar*hcurl_map_dot) + sum(hdiv_map_bar*hdiv_map_dot) + &
        sum(l2_map_bar*l2_map_dot)
    call check_condition(status%code == 0 .and. abs(lhs - rhs) < 1.0e-12_dp, &
        "enriched FEEC VJP satisfies the real dot-product identity")

    call assemble_enriched_feec_sequence( &
        gradient, curl, divergence, scalar_map, hcurl_map(:, 1:1), hdiv_map, l2_map, &
        enriched_gradient, enriched_curl, enriched_divergence, curl_gradient, &
        divergence_curl, status)
    call check_condition(status%code /= 0, &
        "enriched FEEC rejects incompatible enriched-map dimensions")
    call check_summary("enriched FEEC sequence")
end program test_enriched_feec_sequence
