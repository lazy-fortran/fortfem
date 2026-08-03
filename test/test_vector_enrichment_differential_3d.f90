program test_vector_enrichment_differential_3d
    use check, only: check_condition, check_summary
    use fortfem_feec, only: &
        evaluate_vector_enrichment_differential_3d, &
        evaluate_vector_enrichment_differential_3d_jvp, &
        evaluate_vector_enrichment_differential_3d_vjp
    use fortfem_kinds, only: dp
    use fortsparse, only: fortsparse_status_t
    implicit none

    integer, parameter :: point_count = 2
    real(dp), parameter :: eps = 1.0e-7_dp
    real(dp) :: base(3, point_count), base_dot(3, point_count)
    real(dp) :: base_gradient(3, 3, point_count)
    real(dp) :: base_gradient_dot(3, 3, point_count)
    real(dp) :: activation(point_count), activation_dot(point_count)
    real(dp) :: activation_gradient(3, point_count)
    real(dp) :: activation_gradient_dot(3, point_count)
    real(dp) :: enriched(3, point_count), enriched_dot(3, point_count)
    real(dp) :: curl(3, point_count), curl_dot(3, point_count)
    real(dp) :: divergence(point_count), divergence_dot(point_count)
    real(dp) :: enriched_plus(3, point_count), enriched_minus(3, point_count)
    real(dp) :: curl_plus(3, point_count), curl_minus(3, point_count)
    real(dp) :: divergence_plus(point_count), divergence_minus(point_count)
    real(dp) :: enriched_bar(3, point_count), curl_bar(3, point_count)
    real(dp) :: divergence_bar(point_count)
    real(dp) :: base_bar(3, point_count), base_gradient_bar(3, 3, point_count)
    real(dp) :: activation_bar(point_count)
    real(dp) :: activation_gradient_bar(3, point_count)
    real(dp) :: enriched_ref(3, point_count), curl_ref(3, point_count)
    real(dp) :: divergence_ref(point_count), curl_base(3), grad_base(3)
    real(dp) :: lhs, rhs
    real(dp) :: bad_enriched(1, point_count)
    type(fortsparse_status_t) :: status
    integer :: point, component, coordinate

    do point = 1, point_count
        do component = 1, 3
            base(component, point) = 0.2_dp*real(component + point, dp)
            base_dot(component, point) = -0.03_dp*real(component + 2*point, dp)
            activation_gradient(component, point) = &
                0.1_dp*real(component - point, dp)
            activation_gradient_dot(component, point) = &
                -0.02_dp*real(component + point, dp)
            do coordinate = 1, 3
                base_gradient(component, coordinate, point) = &
                    0.04_dp*real(component + coordinate + point, dp)
                base_gradient_dot(component, coordinate, point) = &
                    -0.01_dp*real(component + 2*coordinate + point, dp)
            end do
        end do
        activation(point) = 0.8_dp + 0.1_dp*real(point, dp)
        activation_dot(point) = -0.04_dp*real(point, dp)
    end do

    enriched_ref = 0.0_dp
    curl_ref = 0.0_dp
    divergence_ref = 0.0_dp
    do point = 1, point_count
        enriched_ref(:, point) = activation(point)*base(:, point)
        curl_base = [ &
            base_gradient(3, 2, point) - base_gradient(2, 3, point), &
            base_gradient(1, 3, point) - base_gradient(3, 1, point), &
            base_gradient(2, 1, point) - base_gradient(1, 2, point)]
        curl_ref(:, point) = activation(point)*curl_base + &
            cross3(activation_gradient(:, point), base(:, point))
        grad_base = [base_gradient(1, 1, point), &
            base_gradient(2, 2, point), base_gradient(3, 3, point)]
        divergence_ref(point) = activation(point)*sum(grad_base) + &
            dot_product(activation_gradient(:, point), base(:, point))
    end do
    call evaluate_vector_enrichment_differential_3d( &
        base, base_gradient, activation, activation_gradient, enriched, curl, &
        divergence, status)
    call check_condition(status%code == 0 .and. &
        maxval(abs(enriched - enriched_ref)) < 1.0e-14_dp .and. &
        maxval(abs(curl - curl_ref)) < 1.0e-14_dp .and. &
        maxval(abs(divergence - divergence_ref)) < 1.0e-14_dp, &
        "vector enrichment product rule matches independent curl/divergence oracle")

    call evaluate_vector_enrichment_differential_3d_jvp( &
        base, base_gradient, activation, activation_gradient, base_dot, &
        base_gradient_dot, activation_dot, activation_gradient_dot, &
        enriched_dot, curl_dot, divergence_dot, status)
    call evaluate_vector_enrichment_differential_3d( &
        base + eps*base_dot, base_gradient + eps*base_gradient_dot, &
        activation + eps*activation_dot, &
        activation_gradient + eps*activation_gradient_dot, enriched_plus, &
        curl_plus, divergence_plus, status)
    call evaluate_vector_enrichment_differential_3d( &
        base - eps*base_dot, base_gradient - eps*base_gradient_dot, &
        activation - eps*activation_dot, &
        activation_gradient - eps*activation_gradient_dot, enriched_minus, &
        curl_minus, divergence_minus, status)
    call check_condition(status%code == 0 .and. &
        maxval(abs((enriched_plus - enriched_minus)/(2.0_dp*eps) - enriched_dot)) &
        < 1.0e-8_dp .and. &
        maxval(abs((curl_plus - curl_minus)/(2.0_dp*eps) - curl_dot)) &
        < 1.0e-8_dp .and. &
        maxval(abs((divergence_plus - divergence_minus)/(2.0_dp*eps) - &
        divergence_dot)) &
        < 1.0e-8_dp, "vector enrichment JVP matches central differences")

    do point = 1, point_count
        do component = 1, 3
            enriched_bar(component, point) = &
                0.07_dp*real(component + point, dp)
            curl_bar(component, point) = &
                -0.05_dp*real(component + 2*point, dp)
            activation_gradient_bar(component, point) = 0.0_dp
        end do
        divergence_bar(point) = 0.03_dp*real(point, dp)
    end do
    call evaluate_vector_enrichment_differential_3d_vjp( &
        base, base_gradient, activation, activation_gradient, enriched_bar, &
        curl_bar, divergence_bar, base_bar, base_gradient_bar, activation_bar, &
        activation_gradient_bar, status)
    lhs = sum(enriched_bar*enriched_dot) + sum(curl_bar*curl_dot) + &
        dot_product(divergence_bar, divergence_dot)
    rhs = sum(base_bar*base_dot) + sum(base_gradient_bar*base_gradient_dot) + &
        dot_product(activation_bar, activation_dot) + &
        sum(activation_gradient_bar*activation_gradient_dot)
    call check_condition(status%code == 0 .and. abs(lhs - rhs) < 1.0e-12_dp, &
        "vector enrichment VJP satisfies the real dot-product identity")

    call evaluate_vector_enrichment_differential_3d( &
        base, base_gradient, activation, activation_gradient, bad_enriched, curl, &
        divergence, status)
    call check_condition(status%code /= 0, &
        "vector enrichment rejects an incompatible output shape")
    call check_summary("vector enrichment differential")

contains

    pure function cross3(left, right) result(product)
        real(dp), intent(in) :: left(3), right(3)
        real(dp) :: product(3)

        product = [left(2)*right(3) - left(3)*right(2), &
            left(3)*right(1) - left(1)*right(3), &
            left(1)*right(2) - left(2)*right(1)]
    end function cross3

end program test_vector_enrichment_differential_3d
