module fortfem_enriched_feec_sequence
    !! Matrix-level FEEC composition for fixed-topology enriched spaces.
    !!
    !! The four maps are caller-owned coefficient maps from an enriched
    !! space into its base space.  They may be shifted XFEM/XIGA maps, a
    !! Piola/IGA coefficient map, or an independently generated projection;
    !! no geometry, numbering, or physical PDE is selected here.  For
    !! example, with ``S`` mapping enriched H1 coefficients to base H1
    !! coefficients and ``V`` doing the same for H(curl), the enriched
    !! gradient is ``V**T * G * S``.  The H(div) and L2 maps are composed in
    !! the same way.  Non-commuting enriched maps are deliberately reported
    !! through the two composition defects rather than silently corrected.
    !!
    !! Every value composition has a product-rule JVP and a real Frobenius
    !! VJP.  The maps are dense in this neutral reference layer; a caller can
    !! supply CSC or generated actions through the same algebra without
    !! imposing a sparse storage format here.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortfem_kinds, only: dp
    use fortsparse, only: fortsparse_status_t, status_set, &
        FORTSPARSE_INVALID_MATRIX, FORTSPARSE_OK
    implicit none
    private

    public :: assemble_enriched_feec_sequence
    public :: assemble_enriched_feec_sequence_jvp
    public :: assemble_enriched_feec_sequence_vjp

contains

    subroutine assemble_enriched_feec_sequence( &
            gradient, curl, divergence, scalar_map, hcurl_map, hdiv_map, l2_map, &
            enriched_gradient, enriched_curl, enriched_divergence, &
            curl_gradient, divergence_curl, status)
        !! Compose base incidence maps with four enriched coefficient maps.
        real(dp), intent(in) :: gradient(:, :), curl(:, :), divergence(:, :)
        real(dp), intent(in) :: scalar_map(:, :), hcurl_map(:, :), hdiv_map(:, :)
        real(dp), intent(in) :: l2_map(:, :)
        real(dp), intent(out) :: enriched_gradient(:, :), enriched_curl(:, :)
        real(dp), intent(out) :: enriched_divergence(:, :)
        real(dp), intent(out) :: curl_gradient(:, :), divergence_curl(:, :)
        type(fortsparse_status_t), intent(out) :: status

        enriched_gradient = 0.0_dp
        enriched_curl = 0.0_dp
        enriched_divergence = 0.0_dp
        curl_gradient = 0.0_dp
        divergence_curl = 0.0_dp
        call validate_value_inputs( &
            gradient, curl, divergence, scalar_map, hcurl_map, hdiv_map, l2_map, &
            enriched_gradient, enriched_curl, enriched_divergence, &
            curl_gradient, divergence_curl, status)
        if (status%code /= FORTSPARSE_OK) return

        enriched_gradient = matmul(transpose(hcurl_map), matmul(gradient, scalar_map))
        enriched_curl = matmul(transpose(hdiv_map), matmul(curl, hcurl_map))
        enriched_divergence = matmul( &
            transpose(l2_map), matmul(divergence, hdiv_map))
        curl_gradient = matmul(enriched_curl, enriched_gradient)
        divergence_curl = matmul(enriched_divergence, enriched_curl)
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine assemble_enriched_feec_sequence

    subroutine assemble_enriched_feec_sequence_jvp( &
            gradient, curl, divergence, scalar_map, hcurl_map, hdiv_map, l2_map, &
            gradient_dot, curl_dot, divergence_dot, scalar_map_dot, &
            hcurl_map_dot, hdiv_map_dot, l2_map_dot, enriched_gradient_dot, &
            enriched_curl_dot, enriched_divergence_dot, curl_gradient_dot, &
            divergence_curl_dot, status)
        !! Apply the exact product rule to the enriched maps and defects.
        real(dp), intent(in) :: gradient(:, :), curl(:, :), divergence(:, :)
        real(dp), intent(in) :: scalar_map(:, :), hcurl_map(:, :), hdiv_map(:, :)
        real(dp), intent(in) :: l2_map(:, :)
        real(dp), intent(in) :: gradient_dot(:, :), curl_dot(:, :), divergence_dot(:, :)
        real(dp), intent(in) :: scalar_map_dot(:, :), hcurl_map_dot(:, :)
        real(dp), intent(in) :: hdiv_map_dot(:, :), l2_map_dot(:, :)
        real(dp), intent(out) :: enriched_gradient_dot(:, :), enriched_curl_dot(:, :)
        real(dp), intent(out) :: enriched_divergence_dot(:, :)
        real(dp), intent(out) :: curl_gradient_dot(:, :), divergence_curl_dot(:, :)
        type(fortsparse_status_t), intent(out) :: status
        real(dp), allocatable :: enriched_gradient(:, :), enriched_curl(:, :)
        real(dp), allocatable :: enriched_divergence(:, :)

        enriched_gradient_dot = 0.0_dp
        enriched_curl_dot = 0.0_dp
        enriched_divergence_dot = 0.0_dp
        curl_gradient_dot = 0.0_dp
        divergence_curl_dot = 0.0_dp
        call validate_jvp_inputs( &
            gradient, curl, divergence, scalar_map, hcurl_map, hdiv_map, l2_map, &
            gradient_dot, curl_dot, divergence_dot, scalar_map_dot, hcurl_map_dot, &
            hdiv_map_dot, l2_map_dot, enriched_gradient_dot, enriched_curl_dot, &
            enriched_divergence_dot, curl_gradient_dot, divergence_curl_dot, status)
        if (status%code /= FORTSPARSE_OK) return

        allocate(enriched_gradient(size(hcurl_map, 2), size(scalar_map, 2)))
        allocate(enriched_curl(size(hdiv_map, 2), size(hcurl_map, 2)))
        allocate(enriched_divergence(size(l2_map, 2), size(hdiv_map, 2)))
        enriched_gradient = matmul(transpose(hcurl_map), matmul(gradient, scalar_map))
        enriched_curl = matmul(transpose(hdiv_map), matmul(curl, hcurl_map))
        enriched_divergence = matmul( &
            transpose(l2_map), matmul(divergence, hdiv_map))

        enriched_gradient_dot = matmul(transpose(hcurl_map_dot), &
            matmul(gradient, scalar_map)) + matmul(transpose(hcurl_map), &
            matmul(gradient_dot, scalar_map)) + matmul(transpose(hcurl_map), &
            matmul(gradient, scalar_map_dot))
        enriched_curl_dot = matmul(transpose(hdiv_map_dot), matmul(curl, hcurl_map)) + &
            matmul(transpose(hdiv_map), matmul(curl_dot, hcurl_map)) + &
            matmul(transpose(hdiv_map), matmul(curl, hcurl_map_dot))
        enriched_divergence_dot = matmul(transpose(l2_map_dot), &
            matmul(divergence, hdiv_map)) + matmul(transpose(l2_map), &
            matmul(divergence_dot, hdiv_map)) + matmul(transpose(l2_map), &
            matmul(divergence, hdiv_map_dot))
        curl_gradient_dot = matmul(enriched_curl_dot, enriched_gradient) + &
            matmul(enriched_curl, enriched_gradient_dot)
        divergence_curl_dot = matmul(enriched_divergence_dot, enriched_curl) + &
            matmul(enriched_divergence, enriched_curl_dot)
        deallocate(enriched_gradient, enriched_curl, enriched_divergence)
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine assemble_enriched_feec_sequence_jvp

    subroutine assemble_enriched_feec_sequence_vjp( &
            gradient, curl, divergence, scalar_map, hcurl_map, hdiv_map, l2_map, &
            enriched_gradient_bar, enriched_curl_bar, enriched_divergence_bar, &
            curl_gradient_bar, divergence_curl_bar, gradient_bar, curl_bar, &
            divergence_bar, scalar_map_bar, hcurl_map_bar, hdiv_map_bar, l2_map_bar, &
            status)
        !! Apply the real reverse product of enriched maps and defects.
        real(dp), intent(in) :: gradient(:, :), curl(:, :), divergence(:, :)
        real(dp), intent(in) :: scalar_map(:, :), hcurl_map(:, :), hdiv_map(:, :)
        real(dp), intent(in) :: l2_map(:, :)
        real(dp), intent(in) :: enriched_gradient_bar(:, :), enriched_curl_bar(:, :)
        real(dp), intent(in) :: enriched_divergence_bar(:, :)
        real(dp), intent(in) :: curl_gradient_bar(:, :), divergence_curl_bar(:, :)
        real(dp), intent(out) :: gradient_bar(:, :), curl_bar(:, :), divergence_bar(:, :)
        real(dp), intent(out) :: scalar_map_bar(:, :), hcurl_map_bar(:, :)
        real(dp), intent(out) :: hdiv_map_bar(:, :), l2_map_bar(:, :)
        type(fortsparse_status_t), intent(out) :: status
        real(dp), allocatable :: enriched_gradient(:, :), enriched_curl(:, :)
        real(dp), allocatable :: enriched_divergence(:, :)
        real(dp), allocatable :: gradient_total_bar(:, :), curl_total_bar(:, :)
        real(dp), allocatable :: divergence_total_bar(:, :)

        gradient_bar = 0.0_dp
        curl_bar = 0.0_dp
        divergence_bar = 0.0_dp
        scalar_map_bar = 0.0_dp
        hcurl_map_bar = 0.0_dp
        hdiv_map_bar = 0.0_dp
        l2_map_bar = 0.0_dp
        call validate_vjp_inputs( &
            gradient, curl, divergence, scalar_map, hcurl_map, hdiv_map, l2_map, &
            enriched_gradient_bar, enriched_curl_bar, enriched_divergence_bar, &
            curl_gradient_bar, divergence_curl_bar, gradient_bar, curl_bar, &
            divergence_bar, scalar_map_bar, hcurl_map_bar, hdiv_map_bar, l2_map_bar, &
            status)
        if (status%code /= FORTSPARSE_OK) return

        allocate(enriched_gradient(size(hcurl_map, 2), size(scalar_map, 2)))
        allocate(enriched_curl(size(hdiv_map, 2), size(hcurl_map, 2)))
        allocate(enriched_divergence(size(l2_map, 2), size(hdiv_map, 2)))
        allocate(gradient_total_bar(size(hcurl_map, 2), size(scalar_map, 2)))
        allocate(curl_total_bar(size(hdiv_map, 2), size(hcurl_map, 2)))
        allocate(divergence_total_bar(size(l2_map, 2), size(hdiv_map, 2)))
        enriched_gradient = matmul(transpose(hcurl_map), matmul(gradient, scalar_map))
        enriched_curl = matmul(transpose(hdiv_map), matmul(curl, hcurl_map))
        enriched_divergence = matmul( &
            transpose(l2_map), matmul(divergence, hdiv_map))
        gradient_total_bar = enriched_gradient_bar + matmul( &
            transpose(enriched_curl), curl_gradient_bar)
        curl_total_bar = enriched_curl_bar + matmul(curl_gradient_bar, &
            transpose(enriched_gradient)) + matmul(transpose(enriched_divergence), &
            divergence_curl_bar)
        divergence_total_bar = enriched_divergence_bar + matmul( &
            divergence_curl_bar, transpose(enriched_curl))

        gradient_bar = matmul(hcurl_map, matmul(gradient_total_bar, &
            transpose(scalar_map)))
        curl_bar = matmul(hdiv_map, matmul(curl_total_bar, transpose(hcurl_map)))
        divergence_bar = matmul(l2_map, matmul(divergence_total_bar, &
            transpose(hdiv_map)))
        scalar_map_bar = matmul(transpose(gradient), matmul(hcurl_map, &
            gradient_total_bar))
        hcurl_map_bar = matmul(gradient, matmul(scalar_map, &
            transpose(gradient_total_bar))) + matmul(transpose(curl), &
            matmul(hdiv_map, curl_total_bar))
        hdiv_map_bar = matmul(curl, matmul(hcurl_map, transpose(curl_total_bar))) + &
            matmul(transpose(divergence), matmul(l2_map, divergence_total_bar))
        l2_map_bar = matmul(divergence, matmul(hdiv_map, &
            transpose(divergence_total_bar)))
        deallocate(enriched_gradient, enriched_curl, enriched_divergence, &
            gradient_total_bar, curl_total_bar, divergence_total_bar)
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine assemble_enriched_feec_sequence_vjp

    subroutine validate_value_inputs( &
            gradient, curl, divergence, scalar_map, hcurl_map, hdiv_map, l2_map, &
            enriched_gradient, enriched_curl, enriched_divergence, &
            curl_gradient, divergence_curl, status)
        real(dp), intent(in) :: gradient(:, :), curl(:, :), divergence(:, :)
        real(dp), intent(in) :: scalar_map(:, :), hcurl_map(:, :), hdiv_map(:, :)
        real(dp), intent(in) :: l2_map(:, :)
        real(dp), intent(in) :: enriched_gradient(:, :), enriched_curl(:, :)
        real(dp), intent(in) :: enriched_divergence(:, :), curl_gradient(:, :)
        real(dp), intent(in) :: divergence_curl(:, :)
        type(fortsparse_status_t), intent(out) :: status
        integer :: n0, n1, n2, n3, e0, e1, e2, e3

        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            "enriched FEEC sequence has incompatible dimensions")
        n1 = size(gradient, 1)
        n0 = size(gradient, 2)
        n2 = size(curl, 1)
        n3 = size(divergence, 1)
        e0 = size(scalar_map, 2)
        e1 = size(hcurl_map, 2)
        e2 = size(hdiv_map, 2)
        e3 = size(l2_map, 2)
        if (min(n0, n1, n2, n3, e0, e1, e2, e3) < 1) return
        if (size(curl, 2) /= n1 .or. size(divergence, 2) /= n2 .or. &
            size(scalar_map, 1) /= n0 .or. size(hcurl_map, 1) /= n1 .or. &
            size(hdiv_map, 1) /= n2 .or. size(l2_map, 1) /= n3 .or. &
            any(shape(enriched_gradient) /= [e1, e0]) .or. &
            any(shape(enriched_curl) /= [e2, e1]) .or. &
            any(shape(enriched_divergence) /= [e3, e2]) .or. &
            any(shape(curl_gradient) /= [e2, e0]) .or. &
            any(shape(divergence_curl) /= [e3, e1])) return
        if (any(.not. ieee_is_finite(gradient)) .or. any(.not. ieee_is_finite(curl)) .or. &
            any(.not. ieee_is_finite(divergence)) .or. &
            any(.not. ieee_is_finite(scalar_map)) .or. &
            any(.not. ieee_is_finite(hcurl_map)) .or. &
            any(.not. ieee_is_finite(hdiv_map)) .or. &
            any(.not. ieee_is_finite(l2_map))) return
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine validate_value_inputs

    subroutine validate_jvp_inputs( &
            gradient, curl, divergence, scalar_map, hcurl_map, hdiv_map, l2_map, &
            gradient_dot, curl_dot, divergence_dot, scalar_map_dot, hcurl_map_dot, &
            hdiv_map_dot, l2_map_dot, enriched_gradient_dot, enriched_curl_dot, &
            enriched_divergence_dot, curl_gradient_dot, divergence_curl_dot, status)
        real(dp), intent(in) :: gradient(:, :), curl(:, :), divergence(:, :)
        real(dp), intent(in) :: scalar_map(:, :), hcurl_map(:, :), hdiv_map(:, :)
        real(dp), intent(in) :: l2_map(:, :), gradient_dot(:, :), curl_dot(:, :)
        real(dp), intent(in) :: divergence_dot(:, :), scalar_map_dot(:, :)
        real(dp), intent(in) :: hcurl_map_dot(:, :), hdiv_map_dot(:, :), l2_map_dot(:, :)
        real(dp), intent(in) :: enriched_gradient_dot(:, :), enriched_curl_dot(:, :)
        real(dp), intent(in) :: enriched_divergence_dot(:, :), curl_gradient_dot(:, :)
        real(dp), intent(in) :: divergence_curl_dot(:, :)
        type(fortsparse_status_t), intent(out) :: status

        call validate_value_inputs( &
            gradient, curl, divergence, scalar_map, hcurl_map, hdiv_map, l2_map, &
            enriched_gradient_dot, enriched_curl_dot, enriched_divergence_dot, &
            curl_gradient_dot, divergence_curl_dot, status)
        if (status%code /= FORTSPARSE_OK) return
        if (.not. same_shape_finite(gradient, gradient_dot) .or. &
            .not. same_shape_finite(curl, curl_dot) .or. &
            .not. same_shape_finite(divergence, divergence_dot) .or. &
            .not. same_shape_finite(scalar_map, scalar_map_dot) .or. &
            .not. same_shape_finite(hcurl_map, hcurl_map_dot) .or. &
            .not. same_shape_finite(hdiv_map, hdiv_map_dot) .or. &
            .not. same_shape_finite(l2_map, l2_map_dot)) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "enriched FEEC sequence JVP has incompatible increments")
        end if
    end subroutine validate_jvp_inputs

    subroutine validate_vjp_inputs( &
            gradient, curl, divergence, scalar_map, hcurl_map, hdiv_map, l2_map, &
            enriched_gradient_bar, enriched_curl_bar, enriched_divergence_bar, &
            curl_gradient_bar, divergence_curl_bar, gradient_bar, curl_bar, &
            divergence_bar, scalar_map_bar, hcurl_map_bar, hdiv_map_bar, l2_map_bar, &
            status)
        real(dp), intent(in) :: gradient(:, :), curl(:, :), divergence(:, :)
        real(dp), intent(in) :: scalar_map(:, :), hcurl_map(:, :), hdiv_map(:, :)
        real(dp), intent(in) :: l2_map(:, :), enriched_gradient_bar(:, :)
        real(dp), intent(in) :: enriched_curl_bar(:, :), enriched_divergence_bar(:, :)
        real(dp), intent(in) :: curl_gradient_bar(:, :), divergence_curl_bar(:, :)
        real(dp), intent(in) :: gradient_bar(:, :), curl_bar(:, :), divergence_bar(:, :)
        real(dp), intent(in) :: scalar_map_bar(:, :), hcurl_map_bar(:, :)
        real(dp), intent(in) :: hdiv_map_bar(:, :), l2_map_bar(:, :)
        type(fortsparse_status_t), intent(out) :: status
        integer :: e0, e1, e2, e3, n0, n1, n2, n3

        call validate_value_inputs( &
            gradient, curl, divergence, scalar_map, hcurl_map, hdiv_map, l2_map, &
            enriched_gradient_bar, enriched_curl_bar, enriched_divergence_bar, &
            curl_gradient_bar, divergence_curl_bar, status)
        if (status%code /= FORTSPARSE_OK) return
        n1 = size(gradient, 1); n0 = size(gradient, 2)
        n2 = size(curl, 1); n3 = size(divergence, 1)
        e0 = size(scalar_map, 2); e1 = size(hcurl_map, 2)
        e2 = size(hdiv_map, 2); e3 = size(l2_map, 2)
        if (any(shape(gradient_bar) /= [n1, n0]) .or. &
            any(shape(curl_bar) /= [n2, n1]) .or. &
            any(shape(divergence_bar) /= [n3, n2]) .or. &
            any(shape(scalar_map_bar) /= [n0, e0]) .or. &
            any(shape(hcurl_map_bar) /= [n1, e1]) .or. &
            any(shape(hdiv_map_bar) /= [n2, e2]) .or. &
            any(shape(l2_map_bar) /= [n3, e3]) .or. &
            any(.not. ieee_is_finite(enriched_gradient_bar)) .or. &
            any(.not. ieee_is_finite(enriched_curl_bar)) .or. &
            any(.not. ieee_is_finite(enriched_divergence_bar)) .or. &
            any(.not. ieee_is_finite(curl_gradient_bar)) .or. &
            any(.not. ieee_is_finite(divergence_curl_bar))) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "enriched FEEC sequence VJP has incompatible cotangents")
        end if
    end subroutine validate_vjp_inputs

    pure logical function same_shape_finite(reference, direction) result(valid)
        real(dp), intent(in) :: reference(:, :), direction(:, :)

        valid = all(shape(reference) == shape(direction)) .and. &
            all(ieee_is_finite(direction))
    end function same_shape_finite

end module fortfem_enriched_feec_sequence
