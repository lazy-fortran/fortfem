module fortfem_glued_feec_sequence
    !! Signed local-to-global composition for compatible FEEC maps.
    !!
    !! Each local cell map is accumulated into a caller-owned global numbering
    !! using signed integer IDs. Shared scalar, edge, face, and cell moments
    !! are summed; orientation signs are applied to the relevant rows and
    !! columns. The same primitive therefore covers conforming meshes, broken
    !! DG/HDG numbering, cut-cell spaces, and oriented multipatch IGA maps.
    !! No mesh, metric, trace law, or constitutive model is selected here.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortfem_kinds, only: dp
    use fortsparse, only: fortsparse_status_t, status_set, &
        FORTSPARSE_INVALID_MATRIX, FORTSPARSE_OK
    implicit none
    private

    public :: assemble_glued_feec_sequence
    public :: assemble_glued_feec_sequence_jvp
    public :: assemble_glued_feec_sequence_vjp

contains

    subroutine assemble_glued_feec_sequence( &
            local_gradient, local_curl, local_divergence, scalar_map, hcurl_map, &
            hdiv_map, l2_map, scalar_count, hcurl_count, hdiv_count, l2_count, &
            gradient, curl, divergence, curl_gradient, divergence_curl, status)
        real(dp), intent(in) :: local_gradient(:, :, :), local_curl(:, :, :)
        real(dp), intent(in) :: local_divergence(:, :, :)
        integer, intent(in) :: scalar_map(:, :), hcurl_map(:, :), hdiv_map(:, :)
        integer, intent(in) :: l2_map(:, :)
        integer, intent(in) :: scalar_count, hcurl_count, hdiv_count, l2_count
        real(dp), intent(out) :: gradient(:, :), curl(:, :), divergence(:, :)
        real(dp), intent(out) :: curl_gradient(:, :), divergence_curl(:, :)
        type(fortsparse_status_t), intent(out) :: status

        gradient = 0.0_dp
        curl = 0.0_dp
        divergence = 0.0_dp
        curl_gradient = 0.0_dp
        divergence_curl = 0.0_dp
        call assemble_global_maps( &
            local_gradient, local_curl, local_divergence, scalar_map, hcurl_map, &
            hdiv_map, l2_map, scalar_count, hcurl_count, hdiv_count, l2_count, &
            gradient, curl, divergence, status)
        if (status%code /= FORTSPARSE_OK) return
        curl_gradient = matmul(curl, gradient)
        divergence_curl = matmul(divergence, curl)
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine assemble_glued_feec_sequence

    subroutine assemble_glued_feec_sequence_jvp( &
            local_gradient, local_curl, local_divergence, scalar_map, hcurl_map, &
            hdiv_map, l2_map, scalar_count, hcurl_count, hdiv_count, l2_count, &
            local_gradient_dot, local_curl_dot, local_divergence_dot, gradient_dot, &
            curl_dot, divergence_dot, curl_gradient_dot, divergence_curl_dot, status)
        real(dp), intent(in) :: local_gradient(:, :, :), local_curl(:, :, :)
        real(dp), intent(in) :: local_divergence(:, :, :)
        integer, intent(in) :: scalar_map(:, :), hcurl_map(:, :), hdiv_map(:, :)
        integer, intent(in) :: l2_map(:, :)
        integer, intent(in) :: scalar_count, hcurl_count, hdiv_count, l2_count
        real(dp), intent(in) :: local_gradient_dot(:, :, :), local_curl_dot(:, :, :)
        real(dp), intent(in) :: local_divergence_dot(:, :, :)
        real(dp), intent(out) :: gradient_dot(:, :), curl_dot(:, :), divergence_dot(:, :)
        real(dp), intent(out) :: curl_gradient_dot(:, :), divergence_curl_dot(:, :)
        type(fortsparse_status_t), intent(out) :: status
        real(dp), allocatable :: gradient(:, :), curl(:, :), divergence(:, :)

        gradient_dot = 0.0_dp
        curl_dot = 0.0_dp
        divergence_dot = 0.0_dp
        curl_gradient_dot = 0.0_dp
        divergence_curl_dot = 0.0_dp
        allocate(gradient(size(gradient_dot, 1), size(gradient_dot, 2)))
        allocate(curl(size(curl_dot, 1), size(curl_dot, 2)))
        allocate(divergence(size(divergence_dot, 1), size(divergence_dot, 2)))
        call assemble_global_maps( &
            local_gradient, local_curl, local_divergence, scalar_map, hcurl_map, &
            hdiv_map, l2_map, scalar_count, hcurl_count, hdiv_count, l2_count, &
            gradient, curl, divergence, status)
        if (status%code /= FORTSPARSE_OK) then
            deallocate(gradient, curl, divergence)
            return
        end if
        if (.not. valid_local_direction( &
            local_gradient, local_curl, local_divergence, local_gradient_dot, &
            local_curl_dot, local_divergence_dot)) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "glued FEEC JVP has incompatible local increments")
            deallocate(gradient, curl, divergence)
            return
        end if
        call assemble_global_maps( &
            local_gradient_dot, local_curl_dot, local_divergence_dot, scalar_map, &
            hcurl_map, hdiv_map, l2_map, scalar_count, hcurl_count, hdiv_count, &
            l2_count, gradient_dot, curl_dot, divergence_dot, status)
        if (status%code /= FORTSPARSE_OK) then
            deallocate(gradient, curl, divergence)
            return
        end if
        curl_gradient_dot = matmul(curl_dot, gradient) + matmul(curl, gradient_dot)
        divergence_curl_dot = matmul(divergence_dot, curl) + &
            matmul(divergence, curl_dot)
        deallocate(gradient, curl, divergence)
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine assemble_glued_feec_sequence_jvp

    subroutine assemble_glued_feec_sequence_vjp( &
            local_gradient, local_curl, local_divergence, scalar_map, hcurl_map, &
            hdiv_map, l2_map, scalar_count, hcurl_count, hdiv_count, l2_count, &
            gradient_bar, curl_bar, divergence_bar, curl_gradient_bar, &
            divergence_curl_bar, local_gradient_bar, local_curl_bar, &
            local_divergence_bar, status)
        real(dp), intent(in) :: local_gradient(:, :, :), local_curl(:, :, :)
        real(dp), intent(in) :: local_divergence(:, :, :)
        integer, intent(in) :: scalar_map(:, :), hcurl_map(:, :), hdiv_map(:, :)
        integer, intent(in) :: l2_map(:, :)
        integer, intent(in) :: scalar_count, hcurl_count, hdiv_count, l2_count
        real(dp), intent(in) :: gradient_bar(:, :), curl_bar(:, :), divergence_bar(:, :)
        real(dp), intent(in) :: curl_gradient_bar(:, :), divergence_curl_bar(:, :)
        real(dp), intent(out) :: local_gradient_bar(:, :, :), local_curl_bar(:, :, :)
        real(dp), intent(out) :: local_divergence_bar(:, :, :)
        type(fortsparse_status_t), intent(out) :: status
        real(dp), allocatable :: gradient(:, :), curl(:, :), divergence(:, :)
        real(dp), allocatable :: gradient_total_bar(:, :), curl_total_bar(:, :)
        real(dp), allocatable :: divergence_total_bar(:, :)

        local_gradient_bar = 0.0_dp
        local_curl_bar = 0.0_dp
        local_divergence_bar = 0.0_dp
        allocate(gradient(size(gradient_bar, 1), size(gradient_bar, 2)))
        allocate(curl(size(curl_bar, 1), size(curl_bar, 2)))
        allocate(divergence(size(divergence_bar, 1), size(divergence_bar, 2)))
        call assemble_global_maps( &
            local_gradient, local_curl, local_divergence, scalar_map, hcurl_map, &
            hdiv_map, l2_map, scalar_count, hcurl_count, hdiv_count, l2_count, &
            gradient, curl, divergence, status)
        if (status%code /= FORTSPARSE_OK) then
            deallocate(gradient, curl, divergence)
            return
        end if
        if (.not. valid_output_bars( &
            gradient, curl, divergence, gradient_bar, curl_bar, divergence_bar, &
            curl_gradient_bar, divergence_curl_bar, local_gradient_bar, &
            local_curl_bar, local_divergence_bar)) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "glued FEEC VJP has incompatible cotangents")
            deallocate(gradient, curl, divergence)
            return
        end if
        if (any(shape(local_gradient_bar) /= shape(local_gradient)) .or. &
            any(shape(local_curl_bar) /= shape(local_curl)) .or. &
            any(shape(local_divergence_bar) /= shape(local_divergence))) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "glued FEEC VJP has incompatible local cotangents")
            deallocate(gradient, curl, divergence)
            return
        end if
        allocate(gradient_total_bar(size(gradient, 1), size(gradient, 2)))
        allocate(curl_total_bar(size(curl, 1), size(curl, 2)))
        allocate(divergence_total_bar(size(divergence, 1), size(divergence, 2)))
        gradient_total_bar = gradient_bar + &
            matmul(transpose(curl), curl_gradient_bar)
        curl_total_bar = curl_bar + matmul(curl_gradient_bar, transpose(gradient)) + &
            matmul(transpose(divergence), divergence_curl_bar)
        divergence_total_bar = divergence_bar + &
            matmul(divergence_curl_bar, transpose(curl))
        call scatter_local_bars( &
            scalar_map, hcurl_map, hdiv_map, l2_map, gradient_total_bar, &
            curl_total_bar, divergence_total_bar, local_gradient_bar, &
            local_curl_bar, local_divergence_bar)
        deallocate(gradient, curl, divergence, gradient_total_bar, curl_total_bar, &
            divergence_total_bar)
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine assemble_glued_feec_sequence_vjp

    subroutine assemble_global_maps( &
            local_gradient, local_curl, local_divergence, scalar_map, hcurl_map, &
            hdiv_map, l2_map, scalar_count, hcurl_count, hdiv_count, l2_count, &
            gradient, curl, divergence, status)
        real(dp), intent(in) :: local_gradient(:, :, :), local_curl(:, :, :)
        real(dp), intent(in) :: local_divergence(:, :, :)
        integer, intent(in) :: scalar_map(:, :), hcurl_map(:, :), hdiv_map(:, :)
        integer, intent(in) :: l2_map(:, :)
        integer, intent(in) :: scalar_count, hcurl_count, hdiv_count, l2_count
        real(dp), intent(out) :: gradient(:, :), curl(:, :), divergence(:, :)
        type(fortsparse_status_t), intent(out) :: status
        integer :: cell, row, column, global_row, global_column

        gradient = 0.0_dp
        curl = 0.0_dp
        divergence = 0.0_dp
        if (.not. valid_inputs( &
            local_gradient, local_curl, local_divergence, scalar_map, hcurl_map, &
            hdiv_map, l2_map, scalar_count, hcurl_count, hdiv_count, l2_count, &
            gradient, curl, divergence)) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "glued FEEC maps have incompatible arrays")
            return
        end if

        do cell = 1, size(local_gradient, 3)
            do row = 1, size(local_gradient, 1)
                global_row = abs(hcurl_map(row, cell))
                do column = 1, size(local_gradient, 2)
                    global_column = abs(scalar_map(column, cell))
                    gradient(global_row, global_column) = &
                        gradient(global_row, global_column) + &
                        signed_product(hcurl_map(row, cell), scalar_map(column, cell)) * &
                        local_gradient(row, column, cell)
                end do
            end do
            do row = 1, size(local_curl, 1)
                global_row = abs(hdiv_map(row, cell))
                do column = 1, size(local_curl, 2)
                    global_column = abs(hcurl_map(column, cell))
                    curl(global_row, global_column) = &
                        curl(global_row, global_column) + &
                        signed_product(hdiv_map(row, cell), hcurl_map(column, cell)) * &
                        local_curl(row, column, cell)
                end do
            end do
            do row = 1, size(local_divergence, 1)
                global_row = abs(l2_map(row, cell))
                do column = 1, size(local_divergence, 2)
                    global_column = abs(hdiv_map(column, cell))
                    divergence(global_row, global_column) = &
                        divergence(global_row, global_column) + &
                        signed_product(l2_map(row, cell), hdiv_map(column, cell)) * &
                        local_divergence(row, column, cell)
                end do
            end do
        end do
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine assemble_global_maps

    subroutine scatter_local_bars( &
            scalar_map, hcurl_map, hdiv_map, l2_map, gradient_bar, curl_bar, &
            divergence_bar, local_gradient_bar, local_curl_bar, local_divergence_bar)
        integer, intent(in) :: scalar_map(:, :), hcurl_map(:, :), hdiv_map(:, :)
        integer, intent(in) :: l2_map(:, :)
        real(dp), intent(in) :: gradient_bar(:, :), curl_bar(:, :), divergence_bar(:, :)
        real(dp), intent(out) :: local_gradient_bar(:, :, :), local_curl_bar(:, :, :)
        real(dp), intent(out) :: local_divergence_bar(:, :, :)
        integer :: cell, row, column

        local_gradient_bar = 0.0_dp
        local_curl_bar = 0.0_dp
        local_divergence_bar = 0.0_dp
        do cell = 1, size(scalar_map, 2)
            do row = 1, size(hcurl_map, 1)
                do column = 1, size(scalar_map, 1)
                    local_gradient_bar(row, column, cell) = &
                        signed_product(hcurl_map(row, cell), scalar_map(column, cell)) * &
                        gradient_bar(abs(hcurl_map(row, cell)), abs(scalar_map(column, cell)))
                end do
            end do
            do row = 1, size(hdiv_map, 1)
                do column = 1, size(hcurl_map, 1)
                    local_curl_bar(row, column, cell) = &
                        signed_product(hdiv_map(row, cell), hcurl_map(column, cell)) * &
                        curl_bar(abs(hdiv_map(row, cell)), abs(hcurl_map(column, cell)))
                end do
            end do
            do row = 1, size(l2_map, 1)
                do column = 1, size(hdiv_map, 1)
                    local_divergence_bar(row, column, cell) = &
                        signed_product(l2_map(row, cell), hdiv_map(column, cell)) * &
                        divergence_bar(abs(l2_map(row, cell)), abs(hdiv_map(column, cell)))
                end do
            end do
        end do
    end subroutine scatter_local_bars

    pure real(dp) function signed_product(first, second) result(value)
        integer, intent(in) :: first, second

        value = real(sign(1, first)*sign(1, second), dp)
    end function signed_product

    logical function valid_local_direction( &
            local_gradient, local_curl, local_divergence, local_gradient_dot, &
            local_curl_dot, local_divergence_dot) result(valid)
        real(dp), intent(in) :: local_gradient(:, :, :), local_curl(:, :, :)
        real(dp), intent(in) :: local_divergence(:, :, :)
        real(dp), intent(in) :: local_gradient_dot(:, :, :), local_curl_dot(:, :, :)
        real(dp), intent(in) :: local_divergence_dot(:, :, :)

        valid = all(shape(local_gradient_dot) == shape(local_gradient)) .and. &
            all(shape(local_curl_dot) == shape(local_curl)) .and. &
            all(shape(local_divergence_dot) == shape(local_divergence)) .and. &
            all(ieee_is_finite(local_gradient_dot)) .and. &
            all(ieee_is_finite(local_curl_dot)) .and. &
            all(ieee_is_finite(local_divergence_dot))
    end function valid_local_direction

    logical function valid_output_bars( &
            gradient, curl, divergence, gradient_bar, curl_bar, divergence_bar, &
            curl_gradient_bar, divergence_curl_bar, local_gradient_bar, &
            local_curl_bar, local_divergence_bar) result(valid)
        real(dp), intent(in) :: gradient(:, :), curl(:, :), divergence(:, :)
        real(dp), intent(in) :: gradient_bar(:, :), curl_bar(:, :), divergence_bar(:, :)
        real(dp), intent(in) :: curl_gradient_bar(:, :), divergence_curl_bar(:, :)
        real(dp), intent(in) :: local_gradient_bar(:, :, :), local_curl_bar(:, :, :)
        real(dp), intent(in) :: local_divergence_bar(:, :, :)

        valid = all(shape(gradient_bar) == shape(gradient)) .and. &
            all(shape(curl_bar) == shape(curl)) .and. &
            all(shape(divergence_bar) == shape(divergence)) .and. &
            all(shape(curl_gradient_bar) == [size(curl, 1), size(gradient, 2)]) .and. &
            all(shape(divergence_curl_bar) == [size(divergence, 1), size(curl, 2)]) .and. &
            all(ieee_is_finite(gradient_bar)) .and. all(ieee_is_finite(curl_bar)) .and. &
            all(ieee_is_finite(divergence_bar)) .and. &
            all(ieee_is_finite(curl_gradient_bar)) .and. &
            all(ieee_is_finite(divergence_curl_bar)) .and. &
            size(local_gradient_bar, 1) > 0 .and. size(local_curl_bar, 1) > 0 .and. &
            size(local_divergence_bar, 1) > 0
    end function valid_output_bars

    logical function valid_inputs( &
            local_gradient, local_curl, local_divergence, scalar_map, hcurl_map, &
            hdiv_map, l2_map, scalar_count, hcurl_count, hdiv_count, l2_count, &
            gradient, curl, divergence) result(valid)
        real(dp), intent(in) :: local_gradient(:, :, :), local_curl(:, :, :)
        real(dp), intent(in) :: local_divergence(:, :, :)
        integer, intent(in) :: scalar_map(:, :), hcurl_map(:, :), hdiv_map(:, :)
        integer, intent(in) :: l2_map(:, :)
        integer, intent(in) :: scalar_count, hcurl_count, hdiv_count, l2_count
        real(dp), intent(in) :: gradient(:, :), curl(:, :), divergence(:, :)
        integer :: cells

        cells = size(local_gradient, 3)
        valid = scalar_count > 0 .and. hcurl_count > 0 .and. hdiv_count > 0 .and. &
            l2_count > 0 .and. cells > 0 .and. &
            size(local_gradient, 1) > 0 .and. size(local_gradient, 2) > 0 .and. &
            all(shape(local_curl) == [size(local_curl, 1), size(local_gradient, 1), cells]) .and. &
            size(local_curl, 1) > 0 .and. &
            all(shape(local_divergence) == [size(local_divergence, 1), size(local_curl, 1), cells]) .and. &
            size(local_divergence, 1) > 0 .and. &
            all(shape(scalar_map) == [size(local_gradient, 2), cells]) .and. &
            all(shape(hcurl_map) == [size(local_gradient, 1), cells]) .and. &
            all(shape(hdiv_map) == [size(local_curl, 1), cells]) .and. &
            all(shape(l2_map) == [size(local_divergence, 1), cells]) .and. &
            all(shape(gradient) == [hcurl_count, scalar_count]) .and. &
            all(shape(curl) == [hdiv_count, hcurl_count]) .and. &
            all(shape(divergence) == [l2_count, hdiv_count])
        if (.not. valid) return
        valid = all(abs(scalar_map) >= 1) .and. all(abs(scalar_map) <= scalar_count) .and. &
            all(abs(hcurl_map) >= 1) .and. all(abs(hcurl_map) <= hcurl_count) .and. &
            all(abs(hdiv_map) >= 1) .and. all(abs(hdiv_map) <= hdiv_count) .and. &
            all(abs(l2_map) >= 1) .and. all(abs(l2_map) <= l2_count) .and. &
            all(ieee_is_finite(local_gradient)) .and. all(ieee_is_finite(local_curl)) .and. &
            all(ieee_is_finite(local_divergence))
    end function valid_inputs

end module fortfem_glued_feec_sequence
