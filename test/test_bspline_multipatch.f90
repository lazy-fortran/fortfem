program test_bspline_multipatch
    use check, only: check_condition, check_summary
    use fortfem_api, only: &
        BSPLINE_FACE_X_MAX, BSPLINE_FACE_X_MIN, BSPLINE_FACE_Y_MAX, &
        build_bspline_feec_2d_interface_dofs, &
        build_bspline_feec_2d_operators
    use fortfem_kinds, only: dp
    implicit none

    integer, allocatable :: h1_left(:), h1_right(:)
    integer, allocatable :: hcurl_left(:), hcurl_right(:), hcurl_sign(:)
    real(dp), allocatable :: coefficients_left(:), coefficients_right(:)
    real(dp), allocatable :: curl_left(:, :), curl_right(:, :)
    real(dp), allocatable :: gradient_left(:, :), gradient_right(:, :)
    real(dp), allocatable :: trace_gradient_left(:), trace_gradient_right(:)
    real(dp), parameter :: knots_long(8) = [ &
        0.0_dp, 0.0_dp, 0.0_dp, 0.3_dp, 0.7_dp, 1.0_dp, 1.0_dp, 1.0_dp]
    real(dp), parameter :: knots_short(7) = [ &
        0.0_dp, 0.0_dp, 0.0_dp, 0.4_dp, 1.0_dp, 1.0_dp, 1.0_dp]
    real(dp) :: trace_values(5)
    integer :: status

    call build_bspline_feec_2d_interface_dofs( &
        5, 4, 5, 4, BSPLINE_FACE_X_MAX, BSPLINE_FACE_X_MIN, .false., &
        h1_left, h1_right, hcurl_left, hcurl_right, hcurl_sign, status)
    call check_condition(status == 0 .and. &
        all(h1_left == [5, 10, 15, 20]) .and. &
        all(h1_right == [1, 6, 11, 16]), &
        "Multipatch H1 trace map identifies an aligned vertical interface")
    call check_condition( &
        all(hcurl_left == [21, 26, 31]) .and. &
        all(hcurl_right == [17, 22, 27]) .and. all(hcurl_sign == 1), &
        "Multipatch Hcurl map identifies aligned tangential one-forms")

    call build_bspline_feec_2d_interface_dofs( &
        5, 4, 4, 5, BSPLINE_FACE_Y_MAX, BSPLINE_FACE_X_MIN, .true., &
        h1_left, h1_right, hcurl_left, hcurl_right, hcurl_sign, status)
    call check_condition(status == 0 .and. &
        all(h1_left == [16, 17, 18, 19, 20]) .and. &
        all(h1_right == [17, 13, 9, 5, 1]), &
        "Multipatch H1 trace map handles a rotated, reversed interface")
    call check_condition( &
        all(hcurl_left == [13, 14, 15, 16]) .and. &
        all(hcurl_right == [28, 24, 20, 16]) .and. all(hcurl_sign == -1), &
        "Reversed interfaces flip covariant tangential trace orientation")

    call build_bspline_feec_2d_operators( &
        knots_long, knots_short, 2, 2, gradient_left, curl_left, status)
    call build_bspline_feec_2d_operators( &
        knots_short, knots_long, 2, 2, gradient_right, curl_right, status)
    allocate(coefficients_left(20), coefficients_right(20))
    coefficients_left = 0.0_dp
    coefficients_right = 0.0_dp
    trace_values = [0.0_dp, 0.15_dp, 0.5_dp, 0.85_dp, 1.0_dp]
    coefficients_left(h1_left) = trace_values
    coefficients_right(h1_right) = trace_values
    trace_gradient_left = matmul(gradient_left, coefficients_left)
    trace_gradient_right = matmul(gradient_right, coefficients_right)
    call check_condition(maxval(abs( &
        trace_gradient_right(hcurl_right) - &
        real(hcurl_sign, dp)*trace_gradient_left(hcurl_left))) < 2.0e-14_dp, &
        "Multipatch orientation makes the discrete trace gradient commute")

    call build_bspline_feec_2d_interface_dofs( &
        5, 4, 5, 5, BSPLINE_FACE_X_MAX, BSPLINE_FACE_X_MIN, .false., &
        h1_left, h1_right, hcurl_left, hcurl_right, hcurl_sign, status)
    call check_condition(status /= 0, &
        "Multipatch trace map rejects incompatible interface dimensions")

    call check_summary("Isogeometric multipatch trace topology")
end program test_bspline_multipatch
