program test_feec_facade
    use check, only: check_condition, check_summary
    use fortfem_feec, only: &
        apply_tree_cotree_prolongation, apply_tree_cotree_restriction, &
        assemble_feec_exact_sequence, assemble_feec_exact_sequence_jvp, &
        assemble_feec_exact_sequence_vjp, build_tree_cotree_gauge, &
        tree_cotree_gauge_t
    use fortfem_kinds, only: dp
    use fortsparse, only: fortsparse_status_t
    implicit none

    integer, parameter :: incidence(3, 3) = reshape([ &
        -1, 1, 0, 0, -1, 1, -1, 0, 1], [3, 3])
    real(dp), parameter :: gradient(3, 3) = reshape([ &
        -1.0_dp, 0.0_dp, 1.0_dp, 1.0_dp, -1.0_dp, 0.0_dp, &
        0.0_dp, 1.0_dp, -1.0_dp], [3, 3])
    real(dp), parameter :: curl(2, 3) = reshape([ &
        1.0_dp, 2.0_dp, 1.0_dp, 2.0_dp, 1.0_dp, 2.0_dp], [2, 3])
    real(dp), parameter :: divergence(1, 2) = reshape([ &
        1.0_dp, -0.5_dp], [1, 2])
    real(dp), parameter :: gradient_dot(3, 3) = reshape([ &
        0.1_dp, -0.2_dp, 0.3_dp, -0.4_dp, 0.5_dp, -0.6_dp, &
        0.7_dp, -0.8_dp, 0.9_dp], [3, 3])
    real(dp), parameter :: curl_dot(2, 3) = reshape([ &
        -0.2_dp, 0.3_dp, -0.4_dp, 0.5_dp, -0.6_dp, 0.7_dp], [2, 3])
    real(dp), parameter :: divergence_dot(1, 2) = reshape([ &
        0.8_dp, -0.9_dp], [1, 2])
    real(dp), parameter :: curl_gradient_bar(2, 3) = reshape([ &
        0.4_dp, -0.3_dp, 0.2_dp, -0.1_dp, 0.6_dp, -0.5_dp], [2, 3])
    real(dp), parameter :: divergence_curl_bar(1, 3) = reshape([ &
        0.7_dp, -0.8_dp, 0.9_dp], [1, 3])
    real(dp) :: curl_gradient(2, 3), divergence_curl(1, 3)
    real(dp) :: curl_gradient_dot(2, 3), divergence_curl_dot(1, 3)
    real(dp) :: gradient_bar(3, 3), curl_bar(2, 3), divergence_bar(1, 2)
    real(dp) :: expected_gradient_dot(2, 3), expected_divergence_dot(1, 3)
    real(dp) :: lhs, rhs
    real(dp), allocatable :: reduced(:), prolonged(:)
    type(fortsparse_status_t) :: status
    type(tree_cotree_gauge_t) :: gauge
    integer :: gauge_status

    call assemble_feec_exact_sequence( &
        gradient, curl, divergence, curl_gradient, divergence_curl, status)
    call check_condition(status%code == 0 .and. &
        maxval(abs(curl_gradient)) < 1.0e-14_dp .and. &
        maxval(abs(divergence_curl)) < 1.0e-14_dp, &
        "canonical FEEC facade preserves the de Rham composition invariant")

    call assemble_feec_exact_sequence_jvp( &
        gradient, curl, divergence, gradient_dot, curl_dot, divergence_dot, &
        curl_gradient_dot, divergence_curl_dot, status)
    expected_gradient_dot = matmul(curl_dot, gradient) + &
        matmul(curl, gradient_dot)
    expected_divergence_dot = matmul(divergence_dot, curl) + &
        matmul(divergence, curl_dot)
    call check_condition(status%code == 0 .and. &
        maxval(abs(curl_gradient_dot - expected_gradient_dot)) < 1.0e-14_dp .and. &
        maxval(abs(divergence_curl_dot - expected_divergence_dot)) < 1.0e-14_dp, &
        "canonical FEEC facade exposes the product-rule JVP")

    call assemble_feec_exact_sequence_vjp( &
        gradient, curl, divergence, curl_gradient_bar, divergence_curl_bar, &
        gradient_bar, curl_bar, divergence_bar, status)
    lhs = sum(curl_gradient_bar*curl_gradient_dot) + &
        sum(divergence_curl_bar*divergence_curl_dot)
    rhs = sum(gradient_bar*gradient_dot) + sum(curl_bar*curl_dot) + &
        sum(divergence_bar*divergence_dot)
    call check_condition(status%code == 0 .and. abs(lhs - rhs) < 1.0e-13_dp, &
        "canonical FEEC facade preserves the exact-sequence VJP pairing")

    call build_tree_cotree_gauge(incidence, gauge, gauge_status)
    call apply_tree_cotree_restriction(gauge, [10.0_dp, 20.0_dp, 30.0_dp], &
        reduced, gauge_status)
    call apply_tree_cotree_prolongation(gauge, reduced, prolonged, gauge_status)
    call check_condition(gauge_status == 0 .and. size(reduced) == 1 .and. &
        abs(reduced(1) - 30.0_dp) < 1.0e-14_dp .and. &
        maxval(abs(prolonged - [0.0_dp, 0.0_dp, 30.0_dp])) < 1.0e-14_dp, &
        "canonical FEEC facade exposes a structure-preserving tree gauge")

    call check_summary("FEEC canonical facade")
end program test_feec_facade
