program test_batched_vector_enrichment_differential_3d
    use, intrinsic :: iso_fortran_env, only: real64
    use check, only: check_condition, check_summary
    use fortfem_feec, only: &
        evaluate_batched_vector_enrichment_differential_3d, &
        evaluate_batched_vector_enrichment_differential_3d_jvp, &
        evaluate_batched_vector_enrichment_differential_3d_vjp
    implicit none

    integer, parameter :: dp = real64, basis_count = 2, point_count = 3
    real(dp) :: base(3, basis_count, point_count), base_curl(3, basis_count, point_count)
    real(dp) :: base_div(basis_count, point_count), activation(basis_count, point_count)
    real(dp) :: gradient(3, basis_count, point_count)
    real(dp) :: enriched(3, basis_count, point_count), curl(3, basis_count, point_count)
    real(dp) :: divergence(basis_count, point_count)
    real(dp) :: base_dot(3, basis_count, point_count), curl_dot(3, basis_count, point_count)
    real(dp) :: div_dot(basis_count, point_count), activation_dot(basis_count, point_count)
    real(dp) :: gradient_dot(3, basis_count, point_count)
    real(dp) :: enriched_dot(3, basis_count, point_count), curl_out_dot(3, basis_count, point_count)
    real(dp) :: divergence_dot(basis_count, point_count)
    real(dp) :: enriched_plus(3, basis_count, point_count), curl_plus(3, basis_count, point_count)
    real(dp) :: div_plus(basis_count, point_count), enriched_minus(3, basis_count, point_count)
    real(dp) :: curl_minus(3, basis_count, point_count), div_minus(basis_count, point_count)
    real(dp) :: enriched_bar(3, basis_count, point_count), curl_bar(3, basis_count, point_count)
    real(dp) :: divergence_bar(basis_count, point_count)
    real(dp) :: base_bar(3, basis_count, point_count), base_curl_bar(3, basis_count, point_count)
    real(dp) :: base_div_bar(basis_count, point_count), activation_bar(basis_count, point_count)
    real(dp) :: gradient_bar(3, basis_count, point_count)
    real(dp) :: lhs, rhs, epsilon
    integer :: basis, point, status

    do point = 1, point_count
        do basis = 1, basis_count
            base(:, basis, point) = [0.2_dp*basis, -0.1_dp*point, 0.3_dp + 0.05_dp*basis*point]
            base_curl(:, basis, point) = [0.4_dp - 0.1_dp*point, 0.2_dp*basis, -0.3_dp]
            base_div(basis, point) = 0.15_dp*basis - 0.04_dp*point
            activation(basis, point) = 0.6_dp + 0.05_dp*basis - 0.02_dp*point
            gradient(:, basis, point) = [0.1_dp*basis, -0.2_dp*point, 0.05_dp*(basis + point)]
            base_dot(:, basis, point) = [0.03_dp*point, -0.02_dp*basis, 0.01_dp]
            curl_dot(:, basis, point) = [-0.01_dp, 0.02_dp*point, 0.03_dp*basis]
            div_dot(basis, point) = 0.02_dp*basis
            activation_dot(basis, point) = -0.01_dp*point
            gradient_dot(:, basis, point) = [0.01_dp, 0.02_dp*basis, -0.015_dp*point]
        end do
    end do

    call evaluate_batched_vector_enrichment_differential_3d( &
        base, base_curl, base_div, activation, gradient, enriched, curl, divergence, status)
    call check_condition(status == 0, "batched vector enrichment accepts compatible basis data")
    call check_condition(maxval(abs(divergence)) > 0.0_dp, &
        "batched vector enrichment computes divergence values")
    call evaluate_batched_vector_enrichment_differential_3d_jvp( &
        base, base_curl, base_div, activation, gradient, base_dot, curl_dot, div_dot, &
        activation_dot, gradient_dot, enriched_dot, curl_out_dot, divergence_dot, status)
    epsilon = 1.0e-7_dp
    call evaluate_batched_vector_enrichment_differential_3d( &
        base + epsilon*base_dot, base_curl + epsilon*curl_dot, base_div + epsilon*div_dot, &
        activation + epsilon*activation_dot, gradient + epsilon*gradient_dot, &
        enriched_plus, curl_plus, div_plus, status)
    call evaluate_batched_vector_enrichment_differential_3d( &
        base - epsilon*base_dot, base_curl - epsilon*curl_dot, base_div - epsilon*div_dot, &
        activation - epsilon*activation_dot, gradient - epsilon*gradient_dot, &
        enriched_minus, curl_minus, div_minus, status)
    call check_condition(maxval(abs(enriched_dot - (enriched_plus - enriched_minus)/(2.0_dp*epsilon))) < 2.0e-8_dp .and. &
        maxval(abs(curl_out_dot - (curl_plus - curl_minus)/(2.0_dp*epsilon))) < 2.0e-8_dp .and. &
        maxval(abs(divergence_dot - (div_plus - div_minus)/(2.0_dp*epsilon))) < 2.0e-8_dp, &
        "batched vector enrichment JVP matches central differences")

    enriched_bar = 0.0_dp
    curl_bar = 0.0_dp
    divergence_bar = 0.0_dp
    enriched_bar(1, 1, 1) = 0.3_dp
    enriched_bar(2, 2, 2) = -0.2_dp
    curl_bar(1, 1, 2) = 0.4_dp
    curl_bar(3, 2, 3) = -0.1_dp
    divergence_bar(1, 3) = 0.25_dp
    call evaluate_batched_vector_enrichment_differential_3d_vjp( &
        base, base_curl, base_div, activation, gradient, enriched_bar, curl_bar, &
        divergence_bar, base_bar, base_curl_bar, base_div_bar, activation_bar, gradient_bar, status)
    lhs = sum(enriched_bar*enriched_dot) + sum(curl_bar*curl_out_dot) + sum(divergence_bar*divergence_dot)
    rhs = sum(base_bar*base_dot) + sum(base_curl_bar*curl_dot) + sum(base_div_bar*div_dot) + &
        sum(activation_bar*activation_dot) + sum(gradient_bar*gradient_dot)
    call check_condition(status == 0 .and. abs(lhs - rhs) < 2.0e-12_dp, &
        "batched vector enrichment VJP satisfies the real adjoint identity")
    call check_summary("Batched vector enrichment differential")
end program test_batched_vector_enrichment_differential_3d
