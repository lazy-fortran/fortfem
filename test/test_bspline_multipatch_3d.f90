program test_bspline_multipatch_3d
    use check, only: check_condition, check_summary
    use fortfem_api, only: &
        BSPLINE_FACE_X_MIN, BSPLINE_FACE_Z_MAX, &
        build_bspline_feec_3d_interface_dofs, &
        build_bspline_feec_3d_operators
    use fortfem_kinds, only: dp
    implicit none

    real(dp), parameter :: knots_long(8) = [ &
        0.0_dp, 0.0_dp, 0.0_dp, 0.3_dp, 0.7_dp, 1.0_dp, 1.0_dp, 1.0_dp]
    real(dp), parameter :: knots_short(7) = [ &
        0.0_dp, 0.0_dp, 0.0_dp, 0.4_dp, 1.0_dp, 1.0_dp, 1.0_dp]
    integer, allocatable :: h1_left(:), h1_right(:)
    integer, allocatable :: hcurl_left(:), hcurl_right(:), hcurl_sign(:)
    integer, allocatable :: hdiv_left(:), hdiv_right(:), hdiv_sign(:)
    real(dp), allocatable :: coefficients_left(:), coefficients_right(:)
    real(dp), allocatable :: curl_left(:, :), curl_right(:, :)
    real(dp), allocatable :: divergence_left(:, :), divergence_right(:, :)
    real(dp), allocatable :: field_left(:), field_right(:)
    real(dp), allocatable :: gradient_left(:, :), gradient_right(:, :)
    real(dp) :: trace_values(20), edge_values(31)
    integer :: entry, status

    call build_bspline_feec_3d_interface_dofs( &
        5, 4, 4, 4, 5, 4, BSPLINE_FACE_Z_MAX, BSPLINE_FACE_X_MIN, .false., &
        .true., .false., h1_left, h1_right, hcurl_left, hcurl_right, &
        hcurl_sign, hdiv_left, hdiv_right, hdiv_sign, status)
    call check_condition(status == 0 .and. size(h1_left) == 20 .and. &
        size(hcurl_left) == 31 .and. size(hdiv_left) == 12, &
        "3D multipatch trace maps have tensor-product face dimensions")

    call build_bspline_feec_3d_operators( &
        knots_long, knots_short, knots_short, 2, 2, 2, gradient_left, &
        curl_left, divergence_left, status)
    call build_bspline_feec_3d_operators( &
        knots_short, knots_long, knots_short, 2, 2, 2, gradient_right, &
        curl_right, divergence_right, status)
    allocate(coefficients_left(80), coefficients_right(80))
    coefficients_left = 0.0_dp
    coefficients_right = 0.0_dp
    trace_values = [(real(entry, dp), entry = 1, 20)]
    coefficients_left(h1_left) = trace_values
    coefficients_right(h1_right) = trace_values
    field_left = matmul(gradient_left, coefficients_left)
    field_right = matmul(gradient_right, coefficients_right)
    call check_condition(maxval(abs( &
        field_right(hcurl_right) - &
        real(hcurl_sign, dp)*field_left(hcurl_left))) < 3.0e-13_dp, &
        "3D multipatch edge signs make the gradient trace commute")

    deallocate(coefficients_left, coefficients_right)
    allocate( &
        coefficients_left(size(curl_left, 2)), &
        coefficients_right(size(curl_right, 2)), source=0.0_dp)
    edge_values = [(real(entry, dp)/7.0_dp, entry = 1, 31)]
    coefficients_left(hcurl_left) = edge_values
    coefficients_right(hcurl_right) = real(hcurl_sign, dp)*edge_values
    field_left = matmul(curl_left, coefficients_left)
    field_right = matmul(curl_right, coefficients_right)
    call check_condition(maxval(abs( &
        field_right(hdiv_right) - &
        real(hdiv_sign, dp)*field_left(hdiv_left))) < 3.0e-13_dp, &
        "3D multipatch face signs make the curl trace commute")

    call build_bspline_feec_3d_interface_dofs( &
        5, 4, 4, 4, 4, 5, BSPLINE_FACE_Z_MAX, BSPLINE_FACE_X_MIN, .true., &
        .false., .true., h1_left, h1_right, hcurl_left, hcurl_right, &
        hcurl_sign, hdiv_left, hdiv_right, hdiv_sign, status)
    call build_bspline_feec_3d_operators( &
        knots_short, knots_short, knots_long, 2, 2, 2, gradient_right, &
        curl_right, divergence_right, status)
    deallocate(coefficients_left, coefficients_right)
    allocate(coefficients_left(80), coefficients_right(80))
    coefficients_left = 0.0_dp
    coefficients_right = 0.0_dp
    coefficients_left(h1_left) = trace_values
    coefficients_right(h1_right) = trace_values
    field_left = matmul(gradient_left, coefficients_left)
    field_right = matmul(gradient_right, coefficients_right)
    call check_condition(maxval(abs( &
        field_right(hcurl_right) - &
        real(hcurl_sign, dp)*field_left(hcurl_left))) < 3.0e-13_dp, &
        "Swapped and reversed face map preserves gradient traces")
    deallocate(coefficients_left, coefficients_right)
    allocate( &
        coefficients_left(size(curl_left, 2)), &
        coefficients_right(size(curl_right, 2)), source=0.0_dp)
    coefficients_left(hcurl_left) = edge_values
    coefficients_right(hcurl_right) = real(hcurl_sign, dp)*edge_values
    field_left = matmul(curl_left, coefficients_left)
    field_right = matmul(curl_right, coefficients_right)
    call check_condition(maxval(abs( &
        field_right(hdiv_right) - &
        real(hdiv_sign, dp)*field_left(hdiv_left))) < 3.0e-13_dp, &
        "Swapped and reversed face map preserves curl traces")

    call build_bspline_feec_3d_interface_dofs( &
        5, 4, 4, 4, 4, 4, BSPLINE_FACE_Z_MAX, BSPLINE_FACE_X_MIN, .false., &
        .false., .false., h1_left, h1_right, hcurl_left, hcurl_right, &
        hcurl_sign, hdiv_left, hdiv_right, hdiv_sign, status)
    call check_condition(status /= 0, &
        "3D multipatch trace map rejects incompatible face dimensions")

    call check_summary("Isogeometric 3D multipatch trace complex")
end program test_bspline_multipatch_3d
