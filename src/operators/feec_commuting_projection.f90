module fortfem_feec_commuting_projection
    !! Differentiable commuting-diagram diagnostics for FEEC maps.
    !!
    !! A discrete differential and a projection commute when the discrete
    !! differential applied after projection equals projection applied after
    !! the continuous differential.  The three returned defects are neutral
    !! matrix products and therefore apply to simplicial, cut, IGA,
    !! multipatch, and periodic constructions without choosing a geometry or
    !! constitutive law.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortfem_kinds, only: dp
    use fortsparse, only: fortsparse_status_t, status_set, &
        FORTSPARSE_INVALID_MATRIX, FORTSPARSE_OK
    implicit none
    private

    public :: assemble_feec_commuting_projection
    public :: assemble_feec_commuting_projection_jvp
    public :: assemble_feec_commuting_projection_vjp

contains

    subroutine assemble_feec_commuting_projection( &
            gradient_discrete, curl_discrete, divergence_discrete, &
            gradient_continuous, curl_continuous, divergence_continuous, &
            projection_scalar, projection_hcurl, projection_hdiv, projection_l2, &
            gradient_defect, curl_defect, divergence_defect, status)
        !! Return discrete-minus-continuous commuting-diagram defects.
        real(dp), intent(in) :: gradient_discrete(:, :), curl_discrete(:, :)
        real(dp), intent(in) :: divergence_discrete(:, :)
        real(dp), intent(in) :: gradient_continuous(:, :), curl_continuous(:, :)
        real(dp), intent(in) :: divergence_continuous(:, :)
        real(dp), intent(in) :: projection_scalar(:, :), projection_hcurl(:, :)
        real(dp), intent(in) :: projection_hdiv(:, :), projection_l2(:, :)
        real(dp), intent(out) :: gradient_defect(:, :), curl_defect(:, :)
        real(dp), intent(out) :: divergence_defect(:, :)
        type(fortsparse_status_t), intent(out) :: status

        gradient_defect = 0.0_dp
        curl_defect = 0.0_dp
        divergence_defect = 0.0_dp
        call validate_commuting_inputs( &
            gradient_discrete, curl_discrete, divergence_discrete, &
            gradient_continuous, curl_continuous, divergence_continuous, &
            projection_scalar, projection_hcurl, projection_hdiv, projection_l2, &
            gradient_defect, curl_defect, divergence_defect, status)
        if (status%code /= FORTSPARSE_OK) return

        gradient_defect = matmul(gradient_discrete, projection_scalar) - &
            matmul(projection_hcurl, gradient_continuous)
        curl_defect = matmul(curl_discrete, projection_hcurl) - &
            matmul(projection_hdiv, curl_continuous)
        divergence_defect = matmul(divergence_discrete, projection_hdiv) - &
            matmul(projection_l2, divergence_continuous)
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine assemble_feec_commuting_projection

    subroutine assemble_feec_commuting_projection_jvp( &
            gradient_discrete, curl_discrete, divergence_discrete, &
            gradient_continuous, curl_continuous, divergence_continuous, &
            projection_scalar, projection_hcurl, projection_hdiv, projection_l2, &
            gradient_discrete_dot, curl_discrete_dot, divergence_discrete_dot, &
            gradient_continuous_dot, curl_continuous_dot, divergence_continuous_dot, &
            projection_scalar_dot, projection_hcurl_dot, projection_hdiv_dot, &
            projection_l2_dot, gradient_defect_dot, curl_defect_dot, &
            divergence_defect_dot, status)
        !! Apply the product-rule JVP to all commuting defects.
        real(dp), intent(in) :: gradient_discrete(:, :), curl_discrete(:, :)
        real(dp), intent(in) :: divergence_discrete(:, :)
        real(dp), intent(in) :: gradient_continuous(:, :), curl_continuous(:, :)
        real(dp), intent(in) :: divergence_continuous(:, :)
        real(dp), intent(in) :: projection_scalar(:, :), projection_hcurl(:, :)
        real(dp), intent(in) :: projection_hdiv(:, :), projection_l2(:, :)
        real(dp), intent(in) :: gradient_discrete_dot(:, :), curl_discrete_dot(:, :)
        real(dp), intent(in) :: divergence_discrete_dot(:, :)
        real(dp), intent(in) :: gradient_continuous_dot(:, :)
        real(dp), intent(in) :: curl_continuous_dot(:, :), divergence_continuous_dot(:, :)
        real(dp), intent(in) :: projection_scalar_dot(:, :)
        real(dp), intent(in) :: projection_hcurl_dot(:, :), projection_hdiv_dot(:, :)
        real(dp), intent(in) :: projection_l2_dot(:, :)
        real(dp), intent(out) :: gradient_defect_dot(:, :), curl_defect_dot(:, :)
        real(dp), intent(out) :: divergence_defect_dot(:, :)
        type(fortsparse_status_t), intent(out) :: status

        gradient_defect_dot = 0.0_dp
        curl_defect_dot = 0.0_dp
        divergence_defect_dot = 0.0_dp
        call validate_commuting_inputs( &
            gradient_discrete, curl_discrete, divergence_discrete, &
            gradient_continuous, curl_continuous, divergence_continuous, &
            projection_scalar, projection_hcurl, projection_hdiv, projection_l2, &
            gradient_defect_dot, curl_defect_dot, divergence_defect_dot, status)
        if (status%code /= FORTSPARSE_OK) return
        if (.not. valid_commuting_direction( &
            gradient_discrete_dot, curl_discrete_dot, divergence_discrete_dot, &
            gradient_continuous_dot, curl_continuous_dot, divergence_continuous_dot, &
            projection_scalar_dot, projection_hcurl_dot, projection_hdiv_dot, &
            projection_l2_dot, gradient_discrete, curl_discrete, divergence_discrete, &
            gradient_continuous, curl_continuous, divergence_continuous, &
            projection_scalar, projection_hcurl, projection_hdiv, projection_l2)) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "commuting projection JVP has incompatible increments")
            return
        end if

        gradient_defect_dot = matmul(gradient_discrete_dot, projection_scalar) + &
            matmul(gradient_discrete, projection_scalar_dot) - &
            matmul(projection_hcurl_dot, gradient_continuous) - &
            matmul(projection_hcurl, gradient_continuous_dot)
        curl_defect_dot = matmul(curl_discrete_dot, projection_hcurl) + &
            matmul(curl_discrete, projection_hcurl_dot) - &
            matmul(projection_hdiv_dot, curl_continuous) - &
            matmul(projection_hdiv, curl_continuous_dot)
        divergence_defect_dot = matmul(divergence_discrete_dot, projection_hdiv) + &
            matmul(divergence_discrete, projection_hdiv_dot) - &
            matmul(projection_l2_dot, divergence_continuous) - &
            matmul(projection_l2, divergence_continuous_dot)
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine assemble_feec_commuting_projection_jvp

    subroutine assemble_feec_commuting_projection_vjp( &
            gradient_discrete, curl_discrete, divergence_discrete, &
            gradient_continuous, curl_continuous, divergence_continuous, &
            projection_scalar, projection_hcurl, projection_hdiv, projection_l2, &
            gradient_defect_bar, curl_defect_bar, divergence_defect_bar, &
            gradient_discrete_bar, curl_discrete_bar, divergence_discrete_bar, &
            gradient_continuous_bar, curl_continuous_bar, divergence_continuous_bar, &
            projection_scalar_bar, projection_hcurl_bar, projection_hdiv_bar, &
            projection_l2_bar, status)
        !! Apply the real reverse product to every commuting-defect input.
        real(dp), intent(in) :: gradient_discrete(:, :), curl_discrete(:, :)
        real(dp), intent(in) :: divergence_discrete(:, :)
        real(dp), intent(in) :: gradient_continuous(:, :), curl_continuous(:, :)
        real(dp), intent(in) :: divergence_continuous(:, :)
        real(dp), intent(in) :: projection_scalar(:, :), projection_hcurl(:, :)
        real(dp), intent(in) :: projection_hdiv(:, :), projection_l2(:, :)
        real(dp), intent(in) :: gradient_defect_bar(:, :), curl_defect_bar(:, :)
        real(dp), intent(in) :: divergence_defect_bar(:, :)
        real(dp), intent(out) :: gradient_discrete_bar(:, :), curl_discrete_bar(:, :)
        real(dp), intent(out) :: divergence_discrete_bar(:, :)
        real(dp), intent(out) :: gradient_continuous_bar(:, :), curl_continuous_bar(:, :)
        real(dp), intent(out) :: divergence_continuous_bar(:, :)
        real(dp), intent(out) :: projection_scalar_bar(:, :)
        real(dp), intent(out) :: projection_hcurl_bar(:, :), projection_hdiv_bar(:, :)
        real(dp), intent(out) :: projection_l2_bar(:, :)
        type(fortsparse_status_t), intent(out) :: status

        gradient_discrete_bar = 0.0_dp
        curl_discrete_bar = 0.0_dp
        divergence_discrete_bar = 0.0_dp
        gradient_continuous_bar = 0.0_dp
        curl_continuous_bar = 0.0_dp
        divergence_continuous_bar = 0.0_dp
        projection_scalar_bar = 0.0_dp
        projection_hcurl_bar = 0.0_dp
        projection_hdiv_bar = 0.0_dp
        projection_l2_bar = 0.0_dp
        call validate_commuting_inputs( &
            gradient_discrete, curl_discrete, divergence_discrete, &
            gradient_continuous, curl_continuous, divergence_continuous, &
            projection_scalar, projection_hcurl, projection_hdiv, projection_l2, &
            gradient_defect_bar, curl_defect_bar, divergence_defect_bar, status)
        if (status%code /= FORTSPARSE_OK) return
        if (.not. valid_commuting_cotangents( &
            gradient_defect_bar, curl_defect_bar, divergence_defect_bar, &
            gradient_discrete, curl_discrete, divergence_discrete, &
            gradient_continuous, curl_continuous, divergence_continuous, &
            projection_scalar, projection_hcurl, projection_hdiv, projection_l2, &
            gradient_discrete_bar, curl_discrete_bar, divergence_discrete_bar, &
            gradient_continuous_bar, curl_continuous_bar, divergence_continuous_bar, &
            projection_scalar_bar, projection_hcurl_bar, projection_hdiv_bar, &
            projection_l2_bar)) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "commuting projection VJP has incompatible cotangents")
            return
        end if

        gradient_discrete_bar = matmul(gradient_defect_bar, transpose(projection_scalar))
        projection_scalar_bar = matmul(transpose(gradient_discrete), gradient_defect_bar)
        gradient_continuous_bar = -matmul(transpose(projection_hcurl), gradient_defect_bar)
        projection_hcurl_bar = -matmul(gradient_defect_bar, &
            transpose(gradient_continuous)) + matmul(transpose(curl_discrete), &
            curl_defect_bar)
        curl_discrete_bar = matmul(curl_defect_bar, transpose(projection_hcurl))
        curl_continuous_bar = -matmul(transpose(projection_hdiv), curl_defect_bar)
        projection_hdiv_bar = -matmul(curl_defect_bar, transpose(curl_continuous)) + &
            matmul(transpose(divergence_discrete), divergence_defect_bar)
        divergence_discrete_bar = matmul(divergence_defect_bar, transpose(projection_hdiv))
        divergence_continuous_bar = -matmul(transpose(projection_l2), divergence_defect_bar)
        projection_l2_bar = -matmul(divergence_defect_bar, &
            transpose(divergence_continuous))
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine assemble_feec_commuting_projection_vjp

    pure logical function valid_commuting_direction( &
            gd_dot, cd_dot, dd_dot, gc_dot, cc_dot, dc_dot, p0_dot, p1_dot, &
            p2_dot, p3_dot, gd, cd, dd, gc, cc, dc, p0, p1, p2, p3) result(valid)
        real(dp), intent(in) :: gd_dot(:, :), cd_dot(:, :), dd_dot(:, :)
        real(dp), intent(in) :: gc_dot(:, :), cc_dot(:, :), dc_dot(:, :)
        real(dp), intent(in) :: p0_dot(:, :), p1_dot(:, :), p2_dot(:, :), p3_dot(:, :)
        real(dp), intent(in) :: gd(:, :), cd(:, :), dd(:, :), gc(:, :), cc(:, :), dc(:, :)
        real(dp), intent(in) :: p0(:, :), p1(:, :), p2(:, :), p3(:, :)

        valid = same_shape(gd_dot, gd) .and. same_shape(cd_dot, cd) .and. &
            same_shape(dd_dot, dd) .and. same_shape(gc_dot, gc) .and. &
            same_shape(cc_dot, cc) .and. same_shape(dc_dot, dc) .and. &
            same_shape(p0_dot, p0) .and. same_shape(p1_dot, p1) .and. &
            same_shape(p2_dot, p2) .and. same_shape(p3_dot, p3) .and. &
            all(ieee_is_finite(gd_dot)) .and. all(ieee_is_finite(cd_dot)) .and. &
            all(ieee_is_finite(dd_dot)) .and. all(ieee_is_finite(gc_dot)) .and. &
            all(ieee_is_finite(cc_dot)) .and. all(ieee_is_finite(dc_dot)) .and. &
            all(ieee_is_finite(p0_dot)) .and. all(ieee_is_finite(p1_dot)) .and. &
            all(ieee_is_finite(p2_dot)) .and. all(ieee_is_finite(p3_dot))
    end function valid_commuting_direction

    pure logical function valid_commuting_cotangents( &
            gd_bar, cd_bar, dd_bar, gd, cd, dd, gc, cc, dc, p0, p1, p2, p3, &
            gd_out, cd_out, dd_out, gc_out, cc_out, dc_out, p0_out, p1_out, &
            p2_out, p3_out) result(valid)
        real(dp), intent(in) :: gd_bar(:, :), cd_bar(:, :), dd_bar(:, :)
        real(dp), intent(in) :: gd(:, :), cd(:, :), dd(:, :), gc(:, :), cc(:, :), dc(:, :)
        real(dp), intent(in) :: p0(:, :), p1(:, :), p2(:, :), p3(:, :)
        real(dp), intent(in) :: gd_out(:, :), cd_out(:, :), dd_out(:, :)
        real(dp), intent(in) :: gc_out(:, :), cc_out(:, :), dc_out(:, :)
        real(dp), intent(in) :: p0_out(:, :), p1_out(:, :), p2_out(:, :), p3_out(:, :)

        valid = size(gd_bar, 1) == size(gd, 1) .and. &
            size(gd_bar, 2) == size(gc, 2) .and. &
            size(cd_bar, 1) == size(cd, 1) .and. &
            size(cd_bar, 2) == size(cc, 2) .and. &
            size(dd_bar, 1) == size(dd, 1) .and. &
            size(dd_bar, 2) == size(dc, 2) .and. &
            same_shape(gc_out, gc) .and. &
            same_shape(cc_out, cc) .and. same_shape(dc_out, dc) .and. &
            same_shape(p0_out, p0) .and. same_shape(p1_out, p1) .and. &
            same_shape(p2_out, p2) .and. same_shape(p3_out, p3) .and. &
            same_shape(gd_out, gd) .and. same_shape(cd_out, cd) .and. &
            same_shape(dd_out, dd) .and. all(ieee_is_finite(gd_bar)) .and. &
            all(ieee_is_finite(cd_bar)) .and. all(ieee_is_finite(dd_bar))
    end function valid_commuting_cotangents

    pure logical function same_shape(first, second) result(same)
        real(dp), intent(in) :: first(:, :), second(:, :)

        same = all(shape(first) == shape(second))
    end function same_shape

    subroutine validate_commuting_inputs( &
            gd, cd, dd, gc, cc, dc, p0, p1, p2, p3, defect_g, defect_c, defect_d, status)
        real(dp), intent(in) :: gd(:, :), cd(:, :), dd(:, :), gc(:, :), cc(:, :), dc(:, :)
        real(dp), intent(in) :: p0(:, :), p1(:, :), p2(:, :), p3(:, :)
        real(dp), intent(in) :: defect_g(:, :), defect_c(:, :), defect_d(:, :)
        type(fortsparse_status_t), intent(out) :: status

        integer :: n0d, n1d, n2d, n3d, n0c, n1c, n2c, n3c

        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            "commuting projection maps have incompatible dimensions")
        n1d = size(gd, 1)
        n0d = size(gd, 2)
        n2d = size(cd, 1)
        n3d = size(dd, 1)
        n1c = size(gc, 1)
        n0c = size(gc, 2)
        n2c = size(cc, 1)
        n3c = size(dc, 1)
        if (n0d < 1 .or. n1d < 1 .or. n2d < 1 .or. n3d < 1 .or. &
            n0c < 1 .or. n1c < 1 .or. n2c < 1 .or. n3c < 1) return
        if (size(cd, 2) /= n1d .or. size(dd, 2) /= n2d .or. &
            size(cc, 2) /= n1c .or. size(dc, 2) /= n2c) return
        if (size(p0, 1) /= n0d .or. size(p0, 2) /= n0c .or. &
            size(p1, 1) /= n1d .or. size(p1, 2) /= n1c .or. &
            size(p2, 1) /= n2d .or. size(p2, 2) /= n2c .or. &
            size(p3, 1) /= n3d .or. size(p3, 2) /= n3c) return
        if (size(defect_g, 1) /= n1d .or. size(defect_g, 2) /= n0c .or. &
            size(defect_c, 1) /= n2d .or. size(defect_c, 2) /= n1c .or. &
            size(defect_d, 1) /= n3d .or. size(defect_d, 2) /= n2c) return
        if (.not. all_finite_commuting_inputs( &
            gd, cd, dd, gc, cc, dc, p0, p1, p2, p3)) return
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine validate_commuting_inputs

    pure logical function all_finite_commuting_inputs( &
            gd, cd, dd, gc, cc, dc, p0, p1, p2, p3) result(valid)
        real(dp), intent(in) :: gd(:, :), cd(:, :), dd(:, :), gc(:, :), cc(:, :), dc(:, :)
        real(dp), intent(in) :: p0(:, :), p1(:, :), p2(:, :), p3(:, :)

        valid = all(ieee_is_finite(gd)) .and. all(ieee_is_finite(cd)) .and. &
            all(ieee_is_finite(dd)) .and. all(ieee_is_finite(gc)) .and. &
            all(ieee_is_finite(cc)) .and. all(ieee_is_finite(dc)) .and. &
            all(ieee_is_finite(p0)) .and. all(ieee_is_finite(p1)) .and. &
            all(ieee_is_finite(p2)) .and. all(ieee_is_finite(p3))
    end function all_finite_commuting_inputs

end module fortfem_feec_commuting_projection
