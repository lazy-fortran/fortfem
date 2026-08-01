program test_broken_feec_sequence
    use check, only: check_condition, check_summary
    use fortfem_api, only: assemble_broken_feec_sequence, &
        assemble_broken_feec_sequence_jvp, assemble_broken_feec_sequence_vjp
    use fortfem_kinds, only: dp
    use fortsparse, only: fortsparse_status_t
    implicit none

    integer, parameter :: cell_count = 2, scalar_count = 1
    integer, parameter :: hcurl_count = 3, hdiv_count = 2, l2_count = 1
    real(dp), parameter :: epsilon_fd = 1.0e-7_dp
    real(dp), parameter :: local_gradient(hcurl_count, scalar_count, cell_count) = &
        reshape([1.0_dp, 1.0_dp, 1.0_dp, 2.0_dp, -1.0_dp, 1.0_dp], &
        [hcurl_count, scalar_count, cell_count])
    real(dp), parameter :: local_curl(hdiv_count, hcurl_count, cell_count) = &
        reshape([1.0_dp, 2.0_dp, -1.0_dp, -2.0_dp, 0.0_dp, 0.0_dp, &
        1.0_dp, 2.0_dp, 1.0_dp, 2.0_dp, -1.0_dp, -2.0_dp], &
        [hdiv_count, hcurl_count, cell_count])
    real(dp), parameter :: local_divergence(l2_count, hdiv_count, cell_count) = &
        reshape([1.0_dp, -0.5_dp, 1.0_dp, -0.5_dp], &
        [l2_count, hdiv_count, cell_count])
    real(dp), parameter :: local_gradient_dot(hcurl_count, scalar_count, cell_count) = &
        reshape([0.2_dp, -0.3_dp, 0.4_dp, -0.5_dp, 0.6_dp, -0.7_dp], &
        [hcurl_count, scalar_count, cell_count])
    real(dp), parameter :: local_curl_dot(hdiv_count, hcurl_count, cell_count) = &
        reshape([0.1_dp, -0.2_dp, 0.3_dp, -0.4_dp, 0.5_dp, -0.6_dp, &
        -0.7_dp, 0.8_dp, -0.9_dp, 1.0_dp, -1.1_dp, 1.2_dp], &
        [hdiv_count, hcurl_count, cell_count])
    real(dp), parameter :: local_divergence_dot(l2_count, hdiv_count, cell_count) = &
        reshape([0.3_dp, -0.4_dp, 0.5_dp, -0.6_dp], &
        [l2_count, hdiv_count, cell_count])
    real(dp) :: gradient(hcurl_count*cell_count, scalar_count*cell_count)
    real(dp) :: curl(hdiv_count*cell_count, hcurl_count*cell_count)
    real(dp) :: divergence(l2_count*cell_count, hdiv_count*cell_count)
    real(dp) :: curl_gradient(hdiv_count*cell_count, scalar_count*cell_count)
    real(dp) :: divergence_curl(l2_count*cell_count, hcurl_count*cell_count)
    real(dp) :: gradient_dot(size(gradient, 1), size(gradient, 2))
    real(dp) :: curl_dot(size(curl, 1), size(curl, 2))
    real(dp) :: divergence_dot(size(divergence, 1), size(divergence, 2))
    real(dp) :: curl_gradient_dot(size(curl_gradient, 1), size(curl_gradient, 2))
    real(dp) :: divergence_curl_dot(size(divergence_curl, 1), &
        size(divergence_curl, 2))
    real(dp) :: local_gradient_bar(size(local_gradient, 1), &
        size(local_gradient, 2), size(local_gradient, 3))
    real(dp) :: local_curl_bar(size(local_curl, 1), size(local_curl, 2), &
        size(local_curl, 3))
    real(dp) :: local_divergence_bar(size(local_divergence, 1), &
        size(local_divergence, 2), size(local_divergence, 3))
    real(dp) :: gradient_bar(size(gradient, 1), size(gradient, 2))
    real(dp) :: curl_bar(size(curl, 1), size(curl, 2))
    real(dp) :: divergence_bar(size(divergence, 1), size(divergence, 2))
    real(dp) :: curl_gradient_bar(size(curl_gradient, 1), size(curl_gradient, 2))
    real(dp) :: divergence_curl_bar(size(divergence_curl, 1), &
        size(divergence_curl, 2))
    real(dp) :: gradient_plus(size(gradient, 1), size(gradient, 2))
    real(dp) :: curl_plus(size(curl, 1), size(curl, 2))
    real(dp) :: divergence_plus(size(divergence, 1), size(divergence, 2))
    real(dp) :: curl_gradient_plus(size(curl_gradient, 1), size(curl_gradient, 2))
    real(dp) :: divergence_curl_plus(size(divergence_curl, 1), &
        size(divergence_curl, 2))
    real(dp) :: gradient_minus(size(gradient, 1), size(gradient, 2))
    real(dp) :: curl_minus(size(curl, 1), size(curl, 2))
    real(dp) :: divergence_minus(size(divergence, 1), size(divergence, 2))
    real(dp) :: curl_gradient_minus(size(curl_gradient, 1), size(curl_gradient, 2))
    real(dp) :: divergence_curl_minus(size(divergence_curl, 1), &
        size(divergence_curl, 2))
    real(dp) :: local_gradient_plus(size(local_gradient, 1), &
        size(local_gradient, 2), size(local_gradient, 3))
    real(dp) :: local_curl_plus(size(local_curl, 1), size(local_curl, 2), &
        size(local_curl, 3))
    real(dp) :: local_divergence_plus(size(local_divergence, 1), &
        size(local_divergence, 2), size(local_divergence, 3))
    real(dp) :: local_gradient_minus(size(local_gradient, 1), &
        size(local_gradient, 2), size(local_gradient, 3))
    real(dp) :: local_curl_minus(size(local_curl, 1), size(local_curl, 2), &
        size(local_curl, 3))
    real(dp) :: local_divergence_minus(size(local_divergence, 1), &
        size(local_divergence, 2), size(local_divergence, 3))
    real(dp) :: lhs, rhs
    type(fortsparse_status_t) :: status
    integer :: i
    logical :: all_passed

    all_passed = .true.
    call assemble_broken_feec_sequence( &
        local_gradient, local_curl, local_divergence, gradient, curl, divergence, &
        curl_gradient, divergence_curl, status)
    call record_condition(status%code == 0 .and. &
        maxval(abs(curl_gradient)) < 1.0e-14_dp .and. &
        maxval(abs(divergence_curl)) < 1.0e-14_dp, &
        "broken local maps preserve both exact-sequence compositions")
    call record_condition(maxval(abs(gradient(1:3, 2:2))) < 1.0e-14_dp .and. &
        maxval(abs(gradient(4:6, 1:1))) < 1.0e-14_dp, &
        "broken FEEC assembly does not couple independent cells")

    call assemble_broken_feec_sequence_jvp( &
        local_gradient, local_curl, local_divergence, local_gradient_dot, &
        local_curl_dot, local_divergence_dot, gradient_dot, curl_dot, &
        divergence_dot, curl_gradient_dot, divergence_curl_dot, status)
    call record_condition(status%code == 0, &
        "broken FEEC assembly accepts local operator increments")
    local_gradient_plus = local_gradient + epsilon_fd*local_gradient_dot
    local_curl_plus = local_curl + epsilon_fd*local_curl_dot
    local_divergence_plus = local_divergence + epsilon_fd*local_divergence_dot
    local_gradient_minus = local_gradient - epsilon_fd*local_gradient_dot
    local_curl_minus = local_curl - epsilon_fd*local_curl_dot
    local_divergence_minus = local_divergence - epsilon_fd*local_divergence_dot
    call assemble_broken_feec_sequence( &
        local_gradient_plus, local_curl_plus, local_divergence_plus, gradient_plus, &
        curl_plus, divergence_plus, curl_gradient_plus, divergence_curl_plus, status)
    call assemble_broken_feec_sequence( &
        local_gradient_minus, local_curl_minus, local_divergence_minus, gradient_minus, &
        curl_minus, divergence_minus, curl_gradient_minus, divergence_curl_minus, &
        status)
    call record_condition(maxval(abs(gradient_dot - &
        (gradient_plus - gradient_minus)/(2.0_dp*epsilon_fd))) < 2.0e-8_dp .and. &
        maxval(abs(curl_dot - (curl_plus - curl_minus)/(2.0_dp*epsilon_fd))) < &
        2.0e-8_dp .and. maxval(abs(divergence_dot - &
        (divergence_plus - divergence_minus)/(2.0_dp*epsilon_fd))) < 2.0e-8_dp, &
        "broken FEEC map JVP matches independent finite differences")

    gradient_bar = reshape([(0.1_dp*real(i, dp), i = 1, size(gradient))], &
        shape(gradient))
    curl_bar = reshape([(-0.07_dp*real(i, dp), i = 1, size(curl))], shape(curl))
    divergence_bar = reshape([(0.05_dp*real(i, dp), i = 1, size(divergence))], &
        shape(divergence))
    curl_gradient_bar = reshape([(-0.03_dp*real(i, dp), &
        i = 1, size(curl_gradient))], shape(curl_gradient))
    divergence_curl_bar = reshape([(0.02_dp*real(i, dp), &
        i = 1, size(divergence_curl))], shape(divergence_curl))
    call assemble_broken_feec_sequence_vjp( &
        local_gradient, local_curl, local_divergence, gradient_bar, curl_bar, &
        divergence_bar, curl_gradient_bar, divergence_curl_bar, &
        local_gradient_bar, local_curl_bar, local_divergence_bar, status)
    lhs = sum(gradient_bar*gradient_dot) + sum(curl_bar*curl_dot) + &
        sum(divergence_bar*divergence_dot) + &
        sum(curl_gradient_bar*curl_gradient_dot) + &
        sum(divergence_curl_bar*divergence_curl_dot)
    rhs = sum(local_gradient_bar*local_gradient_dot) + &
        sum(local_curl_bar*local_curl_dot) + &
        sum(local_divergence_bar*local_divergence_dot)
    call record_condition(status%code == 0 .and. abs(lhs - rhs) < 1.0e-12_dp, &
        "broken FEEC VJP satisfies the real dot-product identity")

    call assemble_broken_feec_sequence( &
        local_gradient, local_curl(:, :, 1:1), local_divergence, gradient, curl, &
        divergence, curl_gradient, divergence_curl, status)
    call record_condition(status%code /= 0, &
        "broken FEEC assembly rejects inconsistent cell counts")
    call check_summary("broken FEEC sequence")
    if (.not. all_passed) error stop 1

contains

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_broken_feec_sequence
