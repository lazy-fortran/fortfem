program test_glued_feec_sequence
    use check, only: check_condition, check_summary
    use fortfem_api, only: assemble_glued_feec_sequence, &
        assemble_glued_feec_sequence_jvp, assemble_glued_feec_sequence_vjp
    use fortfem_kinds, only: dp
    use fortsparse, only: fortsparse_status_t
    implicit none

    integer, parameter :: cells = 2, n0 = 2, n1 = 2, n2 = 1, n3 = 1
    integer, parameter :: g0 = 3, g1 = 3, g2 = 1, g3 = 2
    real(dp), parameter :: eps = 1.0e-7_dp
    real(dp) :: local_gradient(n1, n0, cells), local_curl(n2, n1, cells)
    real(dp) :: local_divergence(n3, n2, cells)
    real(dp) :: local_gradient_dot(n1, n0, cells), local_curl_dot(n2, n1, cells)
    real(dp) :: local_divergence_dot(n3, n2, cells)
    integer, parameter :: scalar_map(n0, cells) = reshape([1, 2, 2, 3], [n0, cells])
    integer, parameter :: hcurl_map(n1, cells) = reshape([1, 2, 2, -3], [n1, cells])
    integer, parameter :: hdiv_map(n2, cells) = reshape([1, 1], [n2, cells])
    integer, parameter :: l2_map(n3, cells) = reshape([1, 2], [n3, cells])
    real(dp) :: gradient(g1, g0), curl(g2, g1), divergence(g3, g2)
    real(dp) :: curl_gradient(g2, g0), divergence_curl(g3, g1)
    real(dp) :: gradient_dot(g1, g0), curl_dot(g2, g1), divergence_dot(g3, g2)
    real(dp) :: curl_gradient_dot(g2, g0), divergence_curl_dot(g3, g1)
    real(dp) :: gradient_plus(g1, g0), curl_plus(g2, g1), divergence_plus(g3, g2)
    real(dp) :: curl_gradient_plus(g2, g0), divergence_curl_plus(g3, g1)
    real(dp) :: gradient_minus(g1, g0), curl_minus(g2, g1), divergence_minus(g3, g2)
    real(dp) :: curl_gradient_minus(g2, g0), divergence_curl_minus(g3, g1)
    real(dp) :: gradient_bar(g1, g0), curl_bar(g2, g1), divergence_bar(g3, g2)
    real(dp) :: curl_gradient_bar(g2, g0), divergence_curl_bar(g3, g1)
    real(dp) :: local_gradient_bar(n1, n0, cells), local_curl_bar(n2, n1, cells)
    real(dp) :: local_divergence_bar(n3, n2, cells)
    real(dp) :: lhs, rhs
    type(fortsparse_status_t) :: status
    integer :: i

    local_gradient(:, :, 1) = reshape([1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp], [n1, n0])
    local_gradient(:, :, 2) = reshape([-1.0_dp, 0.5_dp, 2.0_dp, -2.5_dp], [n1, n0])
    local_curl(:, :, 1) = reshape([0.5_dp, -1.0_dp], [n2, n1])
    local_curl(:, :, 2) = reshape([1.25_dp, 0.75_dp], [n2, n1])
    local_divergence(:, :, 1) = reshape([1.5_dp], [n3, n2])
    local_divergence(:, :, 2) = reshape([-0.75_dp], [n3, n2])
    local_gradient_dot = reshape([0.2_dp, -0.3_dp, 0.4_dp, -0.5_dp, &
        0.6_dp, -0.7_dp, 0.8_dp, -0.9_dp], [n1, n0, cells])
    local_curl_dot = reshape([0.1_dp, -0.2_dp, 0.3_dp, -0.4_dp], [n2, n1, cells])
    local_divergence_dot = reshape([0.5_dp, -0.6_dp], [n3, n2, cells])

    call assemble_glued_feec_sequence( &
        local_gradient, local_curl, local_divergence, scalar_map, hcurl_map, &
        hdiv_map, l2_map, g0, g1, g2, g3, gradient, curl, divergence, &
        curl_gradient, divergence_curl, status)
    call assemble_reference(local_gradient, local_curl, local_divergence, &
        gradient_plus, curl_plus, divergence_plus, curl_gradient_plus, &
        divergence_curl_plus)
    call check_condition(status%code == 0 .and. &
        maxval(abs(gradient - gradient_plus)) < 1.0e-14_dp .and. &
        maxval(abs(curl - curl_plus)) < 1.0e-14_dp .and. &
        maxval(abs(divergence - divergence_plus)) < 1.0e-14_dp .and. &
        maxval(abs(curl_gradient - curl_gradient_plus)) < 1.0e-14_dp .and. &
        maxval(abs(divergence_curl - divergence_curl_plus)) < 1.0e-14_dp, &
        "signed FEEC gluing matches an independent local-to-global oracle")

    call assemble_glued_feec_sequence_jvp( &
        local_gradient, local_curl, local_divergence, scalar_map, hcurl_map, &
        hdiv_map, l2_map, g0, g1, g2, g3, local_gradient_dot, local_curl_dot, &
        local_divergence_dot, gradient_dot, curl_dot, divergence_dot, &
        curl_gradient_dot, divergence_curl_dot, status)
    call assemble_reference(local_gradient_dot, local_curl_dot, local_divergence_dot, &
        gradient_plus, curl_plus, divergence_plus, curl_gradient_plus, &
        divergence_curl_plus)
    call check_condition(status%code == 0 .and. &
        maxval(abs(gradient_dot - gradient_plus)) < 1.0e-14_dp .and. &
        maxval(abs(curl_dot - curl_plus)) < 1.0e-14_dp .and. &
        maxval(abs(divergence_dot - divergence_plus)) < 1.0e-14_dp .and. &
        maxval(abs(curl_gradient_dot - ( &
        matmul(curl_dot, gradient) + matmul(curl, gradient_dot)))) < 1.0e-14_dp .and. &
        maxval(abs(divergence_curl_dot - ( &
        matmul(divergence_dot, curl) + matmul(divergence, curl_dot)))) < 1.0e-14_dp, &
        "FEEC gluing JVP applies local scatter and composition product rules")

    call assemble_glued_feec_sequence( &
        local_gradient + eps*local_gradient_dot, local_curl + eps*local_curl_dot, &
        local_divergence + eps*local_divergence_dot, scalar_map, hcurl_map, hdiv_map, &
        l2_map, g0, g1, g2, g3, gradient_plus, curl_plus, divergence_plus, &
        curl_gradient_plus, divergence_curl_plus, status)
    call assemble_glued_feec_sequence( &
        local_gradient - eps*local_gradient_dot, local_curl - eps*local_curl_dot, &
        local_divergence - eps*local_divergence_dot, scalar_map, hcurl_map, hdiv_map, &
        l2_map, g0, g1, g2, g3, gradient_minus, curl_minus, divergence_minus, &
        curl_gradient_minus, divergence_curl_minus, status)
    call check_condition(maxval(abs(gradient_dot - &
        (gradient_plus - gradient_minus)/(2.0_dp*eps))) < 2.0e-8_dp .and. &
        maxval(abs(curl_dot - (curl_plus - curl_minus)/(2.0_dp*eps))) < 2.0e-8_dp .and. &
        maxval(abs(divergence_dot - (divergence_plus - divergence_minus)/ &
        (2.0_dp*eps))) < 2.0e-8_dp .and. &
        maxval(abs(curl_gradient_dot - (curl_gradient_plus - curl_gradient_minus)/ &
        (2.0_dp*eps))) < 2.0e-8_dp .and. &
        maxval(abs(divergence_curl_dot - (divergence_curl_plus - &
        divergence_curl_minus)/(2.0_dp*eps))) < 2.0e-8_dp, &
        "FEEC gluing JVP matches an independent central difference")

    gradient_bar = reshape([(0.1_dp*real(i, dp), i = 1, size(gradient))], shape(gradient))
    curl_bar = reshape([0.2_dp, -0.1_dp, 0.4_dp], shape(curl))
    divergence_bar = reshape([-0.3_dp, 0.4_dp], shape(divergence))
    curl_gradient_bar = reshape([0.5_dp, -0.6_dp, 0.7_dp], shape(curl_gradient))
    divergence_curl_bar = reshape([0.8_dp, -0.9_dp, 1.0_dp, -1.1_dp, 1.2_dp, -1.3_dp], &
        shape(divergence_curl))
    call assemble_glued_feec_sequence_vjp( &
        local_gradient, local_curl, local_divergence, scalar_map, hcurl_map, &
        hdiv_map, l2_map, g0, g1, g2, g3, gradient_bar, curl_bar, divergence_bar, &
        curl_gradient_bar, divergence_curl_bar, local_gradient_bar, local_curl_bar, &
        local_divergence_bar, status)
    lhs = sum(gradient_bar*gradient_dot) + sum(curl_bar*curl_dot) + &
        sum(divergence_bar*divergence_dot) + sum(curl_gradient_bar*curl_gradient_dot) + &
        sum(divergence_curl_bar*divergence_curl_dot)
    rhs = sum(local_gradient_bar*local_gradient_dot) + &
        sum(local_curl_bar*local_curl_dot) + sum(local_divergence_bar*local_divergence_dot)
    call check_condition(status%code == 0 .and. abs(lhs - rhs) < 1.0e-12_dp, &
        "FEEC gluing VJP satisfies the real dot-product identity")

    call assemble_glued_feec_sequence( &
        local_gradient, local_curl, local_divergence, reshape([0, 2, 2, 3], [n0, cells]), &
        hcurl_map, hdiv_map, l2_map, g0, g1, g2, g3, gradient, curl, divergence, &
        curl_gradient, divergence_curl, status)
    call check_condition(status%code /= 0, "FEEC gluing rejects zero DOF IDs")
    call check_summary("glued FEEC sequence")

contains

    subroutine assemble_reference(gradient_local, curl_local, divergence_local, &
            gradient_global, curl_global, divergence_global, curl_gradient_global, &
            divergence_curl_global)
        real(dp), intent(in) :: gradient_local(:, :, :), curl_local(:, :, :)
        real(dp), intent(in) :: divergence_local(:, :, :)
        real(dp), intent(out) :: gradient_global(:, :), curl_global(:, :)
        real(dp), intent(out) :: divergence_global(:, :), curl_gradient_global(:, :)
        real(dp), intent(out) :: divergence_curl_global(:, :)
        integer :: cell, row, column, global_row, global_column

        gradient_global = 0.0_dp
        curl_global = 0.0_dp
        divergence_global = 0.0_dp
        do cell = 1, cells
            do row = 1, n1
                global_row = abs(hcurl_map(row, cell))
                do column = 1, n0
                    global_column = abs(scalar_map(column, cell))
                    gradient_global(global_row, global_column) = &
                        gradient_global(global_row, global_column) + &
                        real(sign(1, hcurl_map(row, cell))*sign(1, scalar_map(column, cell)), dp) * &
                        gradient_local(row, column, cell)
                end do
            end do
            do row = 1, n2
                global_row = abs(hdiv_map(row, cell))
                do column = 1, n1
                    global_column = abs(hcurl_map(column, cell))
                    curl_global(global_row, global_column) = &
                        curl_global(global_row, global_column) + &
                        real(sign(1, hdiv_map(row, cell))*sign(1, hcurl_map(column, cell)), dp) * &
                        curl_local(row, column, cell)
                end do
            end do
            do row = 1, n3
                global_row = abs(l2_map(row, cell))
                do column = 1, n2
                    global_column = abs(hdiv_map(column, cell))
                    divergence_global(global_row, global_column) = &
                        divergence_global(global_row, global_column) + &
                        real(sign(1, l2_map(row, cell))*sign(1, hdiv_map(column, cell)), dp) * &
                        divergence_local(row, column, cell)
                end do
            end do
        end do
        curl_gradient_global = matmul(curl_global, gradient_global)
        divergence_curl_global = matmul(divergence_global, curl_global)
    end subroutine assemble_reference

end program test_glued_feec_sequence
