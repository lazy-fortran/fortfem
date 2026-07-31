program test_fci_bilinear_interpolation_map_ad
    use check, only: check_condition, check_summary
    use fortfem_api, only: &
        build_fci_bilinear_interpolation_map_2d, &
        build_fci_bilinear_interpolation_map_2d_jvp, &
        build_fci_bilinear_interpolation_map_2d_vjp
    use fortfem_kinds, only: dp
    use fortsparse, only: fortsparse_status_t
    implicit none

    real(dp), parameter :: source_x(3) = [0.0_dp, 1.0_dp, 2.0_dp]
    real(dp), parameter :: source_y(2) = [0.0_dp, 2.0_dp]
    real(dp), parameter :: target_x(2) = [0.25_dp, 1.5_dp]
    real(dp), parameter :: target_y(2) = [0.5_dp, 1.0_dp]
    real(dp), parameter :: source_x_dot(3) = [0.01_dp, -0.02_dp, 0.03_dp]
    real(dp), parameter :: source_y_dot(2) = [0.04_dp, -0.03_dp]
    real(dp), parameter :: target_x_dot(2) = [0.05_dp, -0.04_dp]
    real(dp), parameter :: target_y_dot(2) = [0.02_dp, -0.01_dp]
    real(dp), parameter :: map_bar(2, 6) = reshape([ &
        0.7_dp, -0.3_dp, 0.2_dp, 0.1_dp, 0.6_dp, -0.4_dp, &
        -0.1_dp, 0.6_dp, 0.4_dp, -0.2_dp, 0.3_dp, 0.5_dp], [2, 6])
    real(dp), parameter :: endpoint_x(1) = [1.0_dp]
    real(dp), parameter :: endpoint_y(1) = [0.5_dp]
    real(dp) :: interpolation_map(2, 6), map_dot(2, 6)
    real(dp) :: map_plus(2, 6), map_minus(2, 6), map_fd(2, 6)
    real(dp) :: source_x_bar(3), source_y_bar(2), target_x_bar(2)
    real(dp) :: target_y_bar(2), left, right, epsilon
    type(fortsparse_status_t) :: status

    epsilon = 1.0e-6_dp
    call build_fci_bilinear_interpolation_map_2d( &
        source_x, source_y, target_x, target_y, interpolation_map, status)
    call check_condition(status%code == 0, &
        "FCI bilinear AD primal map accepts fixed-topology points")
    call build_fci_bilinear_interpolation_map_2d_jvp( &
        source_x, source_y, target_x, target_y, source_x_dot, source_y_dot, &
        target_x_dot, target_y_dot, map_dot, status)
    call check_condition(status%code == 0, &
        "FCI bilinear JVP accepts fixed-topology coordinate tangents")

    call build_fci_bilinear_interpolation_map_2d( &
        source_x + epsilon*source_x_dot, source_y + epsilon*source_y_dot, &
        target_x + epsilon*target_x_dot, target_y + epsilon*target_y_dot, &
        map_plus, status)
    call build_fci_bilinear_interpolation_map_2d( &
        source_x - epsilon*source_x_dot, source_y - epsilon*source_y_dot, &
        target_x - epsilon*target_x_dot, target_y - epsilon*target_y_dot, &
        map_minus, status)
    map_fd = (map_plus - map_minus)/(2.0_dp*epsilon)
    call check_condition(maxval(abs(map_dot - map_fd)) < 1.0e-8_dp, &
        "FCI bilinear JVP matches the independent central-difference oracle")

    call build_fci_bilinear_interpolation_map_2d_vjp( &
        source_x, source_y, target_x, target_y, map_bar, source_x_bar, &
        source_y_bar, target_x_bar, target_y_bar, status)
    left = sum(map_bar*map_dot)
    right = dot_product(source_x_bar, source_x_dot) + &
        dot_product(source_y_bar, source_y_dot) + &
        dot_product(target_x_bar, target_x_dot) + &
        dot_product(target_y_bar, target_y_dot)
    call check_condition(status%code == 0 .and. &
        abs(left - right) < 1.0e-12_dp, &
        "FCI bilinear VJP satisfies the real dot-product identity")

    call build_fci_bilinear_interpolation_map_2d_jvp( &
        source_x, source_y, endpoint_x, endpoint_y, source_x_dot, source_y_dot, &
        [0.0_dp], [0.0_dp], map_dot(1:1, :), status)
    call check_condition(status%code /= 0, &
        "FCI bilinear JVP rejects a grid-line topology event")
    call build_fci_bilinear_interpolation_map_2d_vjp( &
        source_x, source_y, endpoint_x, endpoint_y, map_bar(1:1, :), &
        source_x_bar, source_y_bar, target_x_bar(1:1), target_y_bar(1:1), &
        status)
    call check_condition(status%code /= 0, &
        "FCI bilinear VJP rejects a grid-line topology event")
    call check_summary("FCI bilinear interpolation derivatives")
end program test_fci_bilinear_interpolation_map_ad
