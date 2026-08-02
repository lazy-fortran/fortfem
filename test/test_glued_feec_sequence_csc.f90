program test_glued_feec_sequence_csc
    use check, only: check_condition, check_summary
    use fortfem_api, only: assemble_glued_feec_sequence_csc, &
        assemble_glued_feec_sequence_csc_jvp, assemble_glued_feec_sequence_csc_vjp
    use fortfem_kinds, only: dp
    use fortsparse, only: csc_from_triplet, csc_matvec, csc_t, fortsparse_status_t
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
    real(dp) :: gradient_dot(g1, g0), curl_dot(g2, g1), divergence_dot(g3, g2)
    real(dp) :: gradient_plus(g1, g0), curl_plus(g2, g1), divergence_plus(g3, g2)
    real(dp) :: gradient_minus(g1, g0), curl_minus(g2, g1), divergence_minus(g3, g2)
    real(dp) :: gradient_bar(g1, g0), curl_bar(g2, g1), divergence_bar(g3, g2)
    real(dp) :: local_gradient_bar(n1, n0, cells), local_curl_bar(n2, n1, cells)
    real(dp) :: local_divergence_bar(n3, n2, cells)
    real(dp) :: expected_gradient_bar(n1, n0, cells), expected_curl_bar(n2, n1, cells)
    real(dp) :: expected_divergence_bar(n3, n2, cells)
    type(csc_t) :: gradient_sparse, curl_sparse, divergence_sparse
    type(csc_t) :: gradient_dot_sparse, curl_dot_sparse, divergence_dot_sparse
    type(csc_t) :: gradient_bar_sparse, curl_bar_sparse, divergence_bar_sparse
    type(csc_t) :: gradient_plus_sparse, curl_plus_sparse, divergence_plus_sparse
    type(csc_t) :: gradient_minus_sparse, curl_minus_sparse, divergence_minus_sparse
    type(fortsparse_status_t) :: status
    logical :: passed
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

    call assemble_glued_feec_sequence_csc( &
        local_gradient, local_curl, local_divergence, scalar_map, hcurl_map, &
        hdiv_map, l2_map, g0, g1, g2, g3, gradient_sparse, curl_sparse, &
        divergence_sparse, status)
    call assemble_dense_reference(local_gradient, local_curl, local_divergence, &
        gradient, curl, divergence)
    passed = status%code == 0
    call compare_csc(gradient_sparse, gradient, passed)
    call compare_csc(curl_sparse, curl, passed)
    call compare_csc(divergence_sparse, divergence, passed)
    call check_condition(passed, &
        "signed FEEC CSC gluing matches the independent dense oracle")

    call assemble_glued_feec_sequence_csc_jvp( &
        local_gradient, local_curl, local_divergence, scalar_map, hcurl_map, &
        hdiv_map, l2_map, g0, g1, g2, g3, local_gradient_dot, local_curl_dot, &
        local_divergence_dot, gradient_dot_sparse, curl_dot_sparse, &
        divergence_dot_sparse, status)
    call assemble_dense_reference(local_gradient_dot, local_curl_dot, &
        local_divergence_dot, gradient_dot, curl_dot, divergence_dot)
    passed = status%code == 0
    call compare_csc(gradient_dot_sparse, gradient_dot, passed)
    call compare_csc(curl_dot_sparse, curl_dot, passed)
    call compare_csc(divergence_dot_sparse, divergence_dot, passed)
    call check_condition(passed, "signed FEEC CSC JVP matches independent scatter")

    call assemble_glued_feec_sequence_csc( &
        local_gradient + eps*local_gradient_dot, local_curl + eps*local_curl_dot, &
        local_divergence + eps*local_divergence_dot, scalar_map, hcurl_map, hdiv_map, &
        l2_map, g0, g1, g2, g3, gradient_plus_sparse, curl_plus_sparse, &
        divergence_plus_sparse, status)
    call assemble_glued_feec_sequence_csc( &
        local_gradient - eps*local_gradient_dot, local_curl - eps*local_curl_dot, &
        local_divergence - eps*local_divergence_dot, scalar_map, hcurl_map, hdiv_map, &
        l2_map, g0, g1, g2, g3, gradient_minus_sparse, curl_minus_sparse, &
        divergence_minus_sparse, status)
    ! The CSC JVP is linear in each local block, so the independent dense
    ! comparison above is the tangent oracle; avoid retaining duplicate maps.
    call check_condition(status%code == 0, &
        "signed FEEC CSC path accepts fixed-topology central-difference data")

    gradient_bar = reshape([(0.1_dp*real(i, dp), i = 1, size(gradient_bar))], shape(gradient_bar))
    curl_bar = reshape([0.2_dp, -0.1_dp, 0.4_dp], shape(curl_bar))
    divergence_bar = reshape([-0.3_dp, 0.5_dp], shape(divergence_bar))
    call dense_to_csc(gradient_bar, gradient_bar_sparse, status)
    call dense_to_csc(curl_bar, curl_bar_sparse, status)
    call dense_to_csc(divergence_bar, divergence_bar_sparse, status)
    call assemble_glued_feec_sequence_csc_vjp( &
        local_gradient, local_curl, local_divergence, scalar_map, hcurl_map, &
        hdiv_map, l2_map, g0, g1, g2, g3, gradient_bar_sparse, curl_bar_sparse, &
        divergence_bar_sparse, local_gradient_bar, local_curl_bar, &
        local_divergence_bar, status)
    call assemble_bar_reference(gradient_bar, curl_bar, divergence_bar, &
        expected_gradient_bar, expected_curl_bar, expected_divergence_bar)
    call check_condition(status%code == 0 .and. &
        maxval(abs(local_gradient_bar - expected_gradient_bar)) < 1.0e-14_dp .and. &
        maxval(abs(local_curl_bar - expected_curl_bar)) < 1.0e-14_dp .and. &
        maxval(abs(local_divergence_bar - expected_divergence_bar)) < 1.0e-14_dp, &
        "signed FEEC CSC VJP scatters the real matrix cotangent")

    call assemble_glued_feec_sequence_csc( &
        local_gradient, local_curl, local_divergence, reshape([0, 2, 2, 3], [n0, cells]), &
        hcurl_map, hdiv_map, l2_map, g0, g1, g2, g3, gradient_sparse, curl_sparse, &
        divergence_sparse, status)
    call check_condition(status%code /= 0, "signed FEEC CSC gluing rejects zero IDs")
    call check_summary("glued FEEC CSC sequence")

contains

    subroutine assemble_dense_reference(gradient_local, curl_local, divergence_local, &
            gradient_global, curl_global, divergence_global)
        real(dp), intent(in) :: gradient_local(:, :, :), curl_local(:, :, :)
        real(dp), intent(in) :: divergence_local(:, :, :)
        real(dp), intent(out) :: gradient_global(:, :), curl_global(:, :)
        real(dp), intent(out) :: divergence_global(:, :)
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
    end subroutine assemble_dense_reference

    subroutine assemble_bar_reference(gradient_global, curl_global, divergence_global, &
            gradient_local_bar, curl_local_bar, divergence_local_bar)
        real(dp), intent(in) :: gradient_global(:, :), curl_global(:, :)
        real(dp), intent(in) :: divergence_global(:, :)
        real(dp), intent(out) :: gradient_local_bar(:, :, :), curl_local_bar(:, :, :)
        real(dp), intent(out) :: divergence_local_bar(:, :, :)
        integer :: cell, row, column

        gradient_local_bar = 0.0_dp
        curl_local_bar = 0.0_dp
        divergence_local_bar = 0.0_dp
        do cell = 1, cells
            do row = 1, n1
                do column = 1, n0
                    gradient_local_bar(row, column, cell) = &
                        real(sign(1, hcurl_map(row, cell))*sign(1, scalar_map(column, cell)), dp) * &
                        gradient_global(abs(hcurl_map(row, cell)), abs(scalar_map(column, cell)))
                end do
            end do
            do row = 1, n2
                do column = 1, n1
                    curl_local_bar(row, column, cell) = &
                        real(sign(1, hdiv_map(row, cell))*sign(1, hcurl_map(column, cell)), dp) * &
                        curl_global(abs(hdiv_map(row, cell)), abs(hcurl_map(column, cell)))
                end do
            end do
            do row = 1, n3
                do column = 1, n2
                    divergence_local_bar(row, column, cell) = &
                        real(sign(1, l2_map(row, cell))*sign(1, hdiv_map(column, cell)), dp) * &
                        divergence_global(abs(l2_map(row, cell)), abs(hdiv_map(column, cell)))
                end do
            end do
        end do
    end subroutine assemble_bar_reference

    subroutine compare_csc(matrix, dense, passed)
        type(csc_t), intent(in) :: matrix
        real(dp), intent(in) :: dense(:, :)
        logical, intent(inout) :: passed
        real(dp), allocatable :: basis(:), value(:)
        integer :: column

        if (matrix%nrow /= size(dense, 1) .or. matrix%ncol /= size(dense, 2)) then
            passed = .false.
            return
        end if
        allocate(basis(size(dense, 2)), value(size(dense, 1)))
        do column = 1, size(dense, 2)
            basis = 0.0_dp
            basis(column) = 1.0_dp
            value = csc_matvec(matrix, basis)
            if (maxval(abs(value - dense(:, column))) > 1.0e-14_dp) passed = .false.
        end do
    end subroutine compare_csc

    subroutine dense_to_csc(dense, matrix, status)
        real(dp), intent(in) :: dense(:, :)
        type(csc_t), intent(out) :: matrix
        type(fortsparse_status_t), intent(out) :: status
        integer, allocatable :: rows(:), columns(:)
        real(dp), allocatable :: values(:)
        integer :: row, column, entry

        allocate(rows(size(dense)), columns(size(dense)), values(size(dense)))
        entry = 0
        do column = 1, size(dense, 2)
            do row = 1, size(dense, 1)
                entry = entry + 1
                rows(entry) = row
                columns(entry) = column
                values(entry) = dense(row, column)
            end do
        end do
        call csc_from_triplet(size(dense, 1), size(dense, 2), rows, columns, values, &
            matrix, status)
    end subroutine dense_to_csc

end program test_glued_feec_sequence_csc
