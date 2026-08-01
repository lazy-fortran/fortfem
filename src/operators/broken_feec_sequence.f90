module fortfem_broken_feec_sequence
    !! Block-diagonal de Rham maps for broken cell-local spaces.
    !!
    !! The local operators are caller-owned reference or physical-cell maps.
    !! This module only places independent cell blocks into broken global
    !! numbering and exposes the two exact-sequence compositions.  No
    !! inter-cell continuity, metric map, trace law, or constitutive model is
    !! selected here, so the same contract applies to DG, HDG, cut/XFEM, and
    !! IGA patch-local operators.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortfem_kinds, only: dp
    use fortsparse, only: fortsparse_status_t, status_set, &
        FORTSPARSE_INVALID_MATRIX, FORTSPARSE_OK
    implicit none
    private

    public :: assemble_broken_feec_sequence
    public :: assemble_broken_feec_sequence_jvp
    public :: assemble_broken_feec_sequence_vjp

contains

    subroutine assemble_broken_feec_sequence( &
            local_gradient, local_curl, local_divergence, gradient, curl, &
            divergence, curl_gradient, divergence_curl, status)
        !! Embed local de Rham maps in cell-major broken global numbering.
        real(dp), intent(in) :: local_gradient(:, :, :), local_curl(:, :, :)
        real(dp), intent(in) :: local_divergence(:, :, :)
        real(dp), intent(out) :: gradient(:, :), curl(:, :), divergence(:, :)
        real(dp), intent(out) :: curl_gradient(:, :), divergence_curl(:, :)
        type(fortsparse_status_t), intent(out) :: status

        integer :: cell, scalar_count, hcurl_count, hdiv_count, l2_count

        gradient = 0.0_dp
        curl = 0.0_dp
        divergence = 0.0_dp
        curl_gradient = 0.0_dp
        divergence_curl = 0.0_dp
        call validate_inputs( &
            local_gradient, local_curl, local_divergence, gradient, curl, &
            divergence, curl_gradient, divergence_curl, status)
        if (status%code /= FORTSPARSE_OK) return

        scalar_count = size(local_gradient, 2)
        hcurl_count = size(local_gradient, 1)
        hdiv_count = size(local_curl, 1)
        l2_count = size(local_divergence, 1)
        do cell = 1, size(local_gradient, 3)
            call copy_block(local_gradient(:, :, cell), gradient, &
                hcurl_count, scalar_count, cell)
            call copy_block(local_curl(:, :, cell), curl, hdiv_count, &
                hcurl_count, cell)
            call copy_block(local_divergence(:, :, cell), divergence, l2_count, &
                hdiv_count, cell)
        end do
        curl_gradient = matmul(curl, gradient)
        divergence_curl = matmul(divergence, curl)
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine assemble_broken_feec_sequence

    subroutine assemble_broken_feec_sequence_jvp( &
            local_gradient, local_curl, local_divergence, local_gradient_dot, &
            local_curl_dot, local_divergence_dot, gradient_dot, curl_dot, &
            divergence_dot, curl_gradient_dot, divergence_curl_dot, status)
        !! Apply the fixed-topology JVP of broken de Rham assembly.
        real(dp), intent(in) :: local_gradient(:, :, :), local_curl(:, :, :)
        real(dp), intent(in) :: local_divergence(:, :, :)
        real(dp), intent(in) :: local_gradient_dot(:, :, :), local_curl_dot(:, :, :)
        real(dp), intent(in) :: local_divergence_dot(:, :, :)
        real(dp), intent(out) :: gradient_dot(:, :), curl_dot(:, :)
        real(dp), intent(out) :: divergence_dot(:, :)
        real(dp), intent(out) :: curl_gradient_dot(:, :), divergence_curl_dot(:, :)
        type(fortsparse_status_t), intent(out) :: status

        integer :: cell, scalar_count, hcurl_count, hdiv_count, l2_count
        real(dp), allocatable :: gradient(:, :), curl(:, :), divergence(:, :)

        gradient_dot = 0.0_dp
        curl_dot = 0.0_dp
        divergence_dot = 0.0_dp
        curl_gradient_dot = 0.0_dp
        divergence_curl_dot = 0.0_dp
        call validate_inputs( &
            local_gradient, local_curl, local_divergence, gradient_dot, curl_dot, &
            divergence_dot, curl_gradient_dot, divergence_curl_dot, status)
        if (status%code /= FORTSPARSE_OK) return
        if (.not. valid_direction( &
            local_gradient_dot, local_curl_dot, local_divergence_dot, &
            local_gradient, local_curl, local_divergence)) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "broken FEEC JVP has incompatible local increments")
            return
        end if

        scalar_count = size(local_gradient, 2)
        hcurl_count = size(local_gradient, 1)
        hdiv_count = size(local_curl, 1)
        l2_count = size(local_divergence, 1)
        allocate(gradient(hcurl_count*size(local_gradient, 3), &
            scalar_count*size(local_gradient, 3)))
        allocate(curl(hdiv_count*size(local_gradient, 3), &
            hcurl_count*size(local_gradient, 3)))
        allocate(divergence(l2_count*size(local_gradient, 3), &
            hdiv_count*size(local_gradient, 3)))
        gradient = 0.0_dp
        curl = 0.0_dp
        divergence = 0.0_dp
        do cell = 1, size(local_gradient, 3)
            call copy_block(local_gradient(:, :, cell), gradient, &
                hcurl_count, scalar_count, cell)
            call copy_block(local_curl(:, :, cell), curl, hdiv_count, &
                hcurl_count, cell)
            call copy_block(local_divergence(:, :, cell), divergence, l2_count, &
                hdiv_count, cell)
            call copy_block(local_gradient_dot(:, :, cell), gradient_dot, &
                hcurl_count, scalar_count, cell)
            call copy_block(local_curl_dot(:, :, cell), curl_dot, hdiv_count, &
                hcurl_count, cell)
            call copy_block(local_divergence_dot(:, :, cell), divergence_dot, &
                l2_count, hdiv_count, cell)
        end do
        curl_gradient_dot = matmul(curl_dot, gradient) + &
            matmul(curl, gradient_dot)
        divergence_curl_dot = matmul(divergence_dot, curl) + &
            matmul(divergence, curl_dot)
        deallocate(gradient, curl, divergence)
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine assemble_broken_feec_sequence_jvp

    subroutine assemble_broken_feec_sequence_vjp( &
            local_gradient, local_curl, local_divergence, gradient_bar, curl_bar, &
            divergence_bar, curl_gradient_bar, divergence_curl_bar, &
            local_gradient_bar, local_curl_bar, local_divergence_bar, status)
        !! Apply the real reverse product of broken maps and compositions.
        real(dp), intent(in) :: local_gradient(:, :, :), local_curl(:, :, :)
        real(dp), intent(in) :: local_divergence(:, :, :)
        real(dp), intent(in) :: gradient_bar(:, :), curl_bar(:, :)
        real(dp), intent(in) :: divergence_bar(:, :)
        real(dp), intent(in) :: curl_gradient_bar(:, :), divergence_curl_bar(:, :)
        real(dp), intent(out) :: local_gradient_bar(:, :, :), local_curl_bar(:, :, :)
        real(dp), intent(out) :: local_divergence_bar(:, :, :)
        type(fortsparse_status_t), intent(out) :: status

        integer :: cell, scalar_count, hcurl_count, hdiv_count, l2_count
        integer :: scalar_offset, hcurl_offset, hdiv_offset, l2_offset
        real(dp), allocatable :: gradient_local_bar(:, :), curl_local_bar(:, :)
        real(dp), allocatable :: divergence_local_bar(:, :)

        local_gradient_bar = 0.0_dp
        local_curl_bar = 0.0_dp
        local_divergence_bar = 0.0_dp
        call validate_inputs( &
            local_gradient, local_curl, local_divergence, gradient_bar, curl_bar, &
            divergence_bar, curl_gradient_bar, divergence_curl_bar, status)
        if (status%code /= FORTSPARSE_OK) return
        if (.not. valid_bars( &
            gradient_bar, curl_bar, divergence_bar, curl_gradient_bar, &
            divergence_curl_bar, local_gradient, local_curl, local_divergence, &
            local_gradient_bar, local_curl_bar, local_divergence_bar)) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "broken FEEC VJP has incompatible cotangents")
            return
        end if

        scalar_count = size(local_gradient, 2)
        hcurl_count = size(local_gradient, 1)
        hdiv_count = size(local_curl, 1)
        l2_count = size(local_divergence, 1)
        allocate(gradient_local_bar(hcurl_count, scalar_count))
        allocate(curl_local_bar(hdiv_count, hcurl_count))
        allocate(divergence_local_bar(l2_count, hdiv_count))
        do cell = 1, size(local_gradient, 3)
            scalar_offset = (cell - 1)*scalar_count
            hcurl_offset = (cell - 1)*hcurl_count
            hdiv_offset = (cell - 1)*hdiv_count
            l2_offset = (cell - 1)*l2_count
            gradient_local_bar = gradient_bar( &
                hcurl_offset + 1:hcurl_offset + hcurl_count, &
                scalar_offset + 1:scalar_offset + scalar_count)
            curl_local_bar = curl_bar( &
                hdiv_offset + 1:hdiv_offset + hdiv_count, &
                hcurl_offset + 1:hcurl_offset + hcurl_count)
            divergence_local_bar = divergence_bar( &
                l2_offset + 1:l2_offset + l2_count, &
                hdiv_offset + 1:hdiv_offset + hdiv_count)
            gradient_local_bar = gradient_local_bar + matmul(transpose( &
                local_curl(:, :, cell)), curl_gradient_bar( &
                hdiv_offset + 1:hdiv_offset + hdiv_count, &
                scalar_offset + 1:scalar_offset + scalar_count))
            curl_local_bar = curl_local_bar + matmul(curl_gradient_bar( &
                hdiv_offset + 1:hdiv_offset + hdiv_count, &
                scalar_offset + 1:scalar_offset + scalar_count), transpose( &
                local_gradient(:, :, cell)))
            curl_local_bar = curl_local_bar + matmul(transpose( &
                local_divergence(:, :, cell)), divergence_curl_bar( &
                l2_offset + 1:l2_offset + l2_count, &
                hcurl_offset + 1:hcurl_offset + hcurl_count))
            divergence_local_bar = divergence_local_bar + matmul( &
                divergence_curl_bar(l2_offset + 1:l2_offset + l2_count, &
                hcurl_offset + 1:hcurl_offset + hcurl_count), transpose( &
                local_curl(:, :, cell)))
            local_gradient_bar(:, :, cell) = gradient_local_bar
            local_curl_bar(:, :, cell) = curl_local_bar
            local_divergence_bar(:, :, cell) = divergence_local_bar
        end do
        deallocate(gradient_local_bar, curl_local_bar, divergence_local_bar)
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine assemble_broken_feec_sequence_vjp

    subroutine copy_block(local_block, global_block, row_count, column_count, cell)
        real(dp), intent(in) :: local_block(:, :)
        real(dp), intent(inout) :: global_block(:, :)
        integer, intent(in) :: row_count, column_count, cell
        integer :: row_offset, column_offset

        row_offset = (cell - 1)*row_count
        column_offset = (cell - 1)*column_count
        global_block(row_offset + 1:row_offset + row_count, &
            column_offset + 1:column_offset + column_count) = local_block
    end subroutine copy_block

    subroutine validate_inputs( &
            local_gradient, local_curl, local_divergence, gradient, curl, &
            divergence, curl_gradient, divergence_curl, status)
        real(dp), intent(in) :: local_gradient(:, :, :), local_curl(:, :, :)
        real(dp), intent(in) :: local_divergence(:, :, :)
        real(dp), intent(in) :: gradient(:, :), curl(:, :), divergence(:, :)
        real(dp), intent(in) :: curl_gradient(:, :), divergence_curl(:, :)
        type(fortsparse_status_t), intent(out) :: status
        integer :: cells, scalar_count, hcurl_count, hdiv_count, l2_count

        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            "broken FEEC maps have incompatible dimensions")
        cells = size(local_gradient, 3)
        scalar_count = size(local_gradient, 2)
        hcurl_count = size(local_gradient, 1)
        hdiv_count = size(local_curl, 1)
        l2_count = size(local_divergence, 1)
        if (cells < 1 .or. scalar_count < 1 .or. hcurl_count < 1 .or. &
            hdiv_count < 1 .or. l2_count < 1 .or. &
            size(local_curl, 2) /= hcurl_count .or. &
            size(local_curl, 3) /= cells .or. &
            size(local_divergence, 2) /= hdiv_count .or. &
            size(local_divergence, 3) /= cells .or. &
            size(gradient, 1) /= hcurl_count*cells .or. &
            size(gradient, 2) /= scalar_count*cells .or. &
            size(curl, 1) /= hdiv_count*cells .or. &
            size(curl, 2) /= hcurl_count*cells .or. &
            size(divergence, 1) /= l2_count*cells .or. &
            size(divergence, 2) /= hdiv_count*cells .or. &
            size(curl_gradient, 1) /= hdiv_count*cells .or. &
            size(curl_gradient, 2) /= scalar_count*cells .or. &
            size(divergence_curl, 1) /= l2_count*cells .or. &
            size(divergence_curl, 2) /= hcurl_count*cells) return
        if (any(.not. ieee_is_finite(local_gradient)) .or. &
            any(.not. ieee_is_finite(local_curl)) .or. &
            any(.not. ieee_is_finite(local_divergence)) .or. &
            any(.not. ieee_is_finite(gradient)) .or. &
            any(.not. ieee_is_finite(curl)) .or. &
            any(.not. ieee_is_finite(divergence)) .or. &
            any(.not. ieee_is_finite(curl_gradient)) .or. &
            any(.not. ieee_is_finite(divergence_curl))) return
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine validate_inputs

    pure logical function valid_direction( &
            local_gradient_dot, local_curl_dot, local_divergence_dot, &
            local_gradient, local_curl, local_divergence) result(valid)
        real(dp), intent(in) :: local_gradient_dot(:, :, :), local_curl_dot(:, :, :)
        real(dp), intent(in) :: local_divergence_dot(:, :, :)
        real(dp), intent(in) :: local_gradient(:, :, :), local_curl(:, :, :)
        real(dp), intent(in) :: local_divergence(:, :, :)

        valid = all(shape(local_gradient_dot) == shape(local_gradient)) .and. &
            all(shape(local_curl_dot) == shape(local_curl)) .and. &
            all(shape(local_divergence_dot) == shape(local_divergence)) .and. &
            all(ieee_is_finite(local_gradient_dot)) .and. &
            all(ieee_is_finite(local_curl_dot)) .and. &
            all(ieee_is_finite(local_divergence_dot))
    end function valid_direction

    pure logical function valid_bars( &
            gradient_bar, curl_bar, divergence_bar, curl_gradient_bar, &
            divergence_curl_bar, local_gradient, local_curl, local_divergence, &
            local_gradient_bar, local_curl_bar, local_divergence_bar) result(valid)
        real(dp), intent(in) :: gradient_bar(:, :), curl_bar(:, :)
        real(dp), intent(in) :: divergence_bar(:, :)
        real(dp), intent(in) :: curl_gradient_bar(:, :), divergence_curl_bar(:, :)
        real(dp), intent(in) :: local_gradient(:, :, :), local_curl(:, :, :)
        real(dp), intent(in) :: local_divergence(:, :, :)
        real(dp), intent(in) :: local_gradient_bar(:, :, :), local_curl_bar(:, :, :)
        real(dp), intent(in) :: local_divergence_bar(:, :, :)
        integer :: cells, scalar_count, hcurl_count, hdiv_count, l2_count

        cells = size(local_gradient, 3)
        scalar_count = size(local_gradient, 2)
        hcurl_count = size(local_gradient, 1)
        hdiv_count = size(local_curl, 1)
        l2_count = size(local_divergence, 1)
        valid = size(gradient_bar, 1) == hcurl_count*cells .and. &
            size(gradient_bar, 2) == scalar_count*cells .and. &
            size(curl_bar, 1) == hdiv_count*cells .and. &
            size(curl_bar, 2) == hcurl_count*cells .and. &
            size(divergence_bar, 1) == l2_count*cells .and. &
            size(divergence_bar, 2) == hdiv_count*cells .and. &
            size(curl_gradient_bar, 1) == hdiv_count*cells .and. &
            size(curl_gradient_bar, 2) == scalar_count*cells .and. &
            size(divergence_curl_bar, 1) == l2_count*cells .and. &
            size(divergence_curl_bar, 2) == hcurl_count*cells .and. &
            all(shape(local_gradient_bar) == shape(local_gradient)) .and. &
            all(shape(local_curl_bar) == shape(local_curl)) .and. &
            all(shape(local_divergence_bar) == shape(local_divergence)) .and. &
            all(ieee_is_finite(gradient_bar)) .and. all(ieee_is_finite(curl_bar)) .and. &
            all(ieee_is_finite(divergence_bar)) .and. &
            all(ieee_is_finite(curl_gradient_bar)) .and. &
            all(ieee_is_finite(divergence_curl_bar))
    end function valid_bars

end module fortfem_broken_feec_sequence
