program test_fci_bilinear_maps_2d
    use check, only: check_condition, check_summary
    use fortfem_fci_interpolation_map, only: &
        build_fci_bilinear_interpolation_maps_2d, &
        build_fci_bilinear_interpolation_maps_2d_jvp, &
        build_fci_bilinear_interpolation_maps_2d_vjp
    use fortfem_kinds, only: dp
    use fortsparse, only: fortsparse_status_t
    implicit none

    real(dp), parameter :: source_x(3) = [0.0_dp, 1.0_dp, 2.0_dp]
    real(dp), parameter :: source_y(2) = [0.0_dp, 2.0_dp]
    real(dp), parameter :: forward_x(2, 2) = reshape([ &
        0.25_dp, 1.50_dp, 0.50_dp, 1.25_dp], [2, 2])
    real(dp), parameter :: forward_y(2, 2) = reshape([ &
        0.50_dp, 1.00_dp, 1.25_dp, 0.75_dp], [2, 2])
    real(dp), parameter :: backward_x(2, 2) = reshape([ &
        0.75_dp, 1.25_dp, 0.25_dp, 1.75_dp], [2, 2])
    real(dp), parameter :: backward_y(2, 2) = reshape([ &
        0.25_dp, 1.50_dp, 1.25_dp, 0.50_dp], [2, 2])
    real(dp), parameter :: source_x_dot(3) = [0.01_dp, -0.02_dp, 0.03_dp]
    real(dp), parameter :: source_y_dot(2) = [0.04_dp, -0.03_dp]
    real(dp), parameter :: forward_x_dot(2, 2) = reshape([ &
        0.05_dp, -0.04_dp, 0.03_dp, -0.02_dp], [2, 2])
    real(dp), parameter :: forward_y_dot(2, 2) = reshape([ &
        0.02_dp, -0.01_dp, -0.03_dp, 0.04_dp], [2, 2])
    real(dp), parameter :: backward_x_dot(2, 2) = reshape([ &
        -0.03_dp, 0.04_dp, -0.05_dp, 0.06_dp], [2, 2])
    real(dp), parameter :: backward_y_dot(2, 2) = reshape([ &
        0.06_dp, -0.05_dp, 0.04_dp, -0.03_dp], [2, 2])
    real(dp), parameter :: forward_map_bar(2, 6, 2) = reshape([ &
        0.7_dp, -0.3_dp, 0.2_dp, 0.1_dp, 0.6_dp, -0.4_dp, &
        -0.1_dp, 0.6_dp, 0.4_dp, -0.2_dp, 0.3_dp, 0.5_dp, &
        0.2_dp, -0.5_dp, 0.4_dp, 0.8_dp, -0.6_dp, 0.1_dp, &
        -0.3_dp, 0.9_dp, -0.2_dp, 0.5_dp, 0.7_dp, -0.8_dp], [2, 6, 2])
    real(dp), parameter :: backward_map_bar(2, 6, 2) = reshape([ &
        -0.2_dp, 0.4_dp, -0.6_dp, 0.8_dp, -1.0_dp, 1.2_dp, &
        0.3_dp, -0.5_dp, 0.7_dp, -0.9_dp, 1.1_dp, -1.3_dp, &
        0.4_dp, -0.6_dp, 0.8_dp, -1.0_dp, 1.2_dp, -1.4_dp, &
        -0.5_dp, 0.7_dp, -0.9_dp, 1.1_dp, -1.3_dp, 1.5_dp], [2, 6, 2])
    real(dp), parameter :: epsilon = 1.0e-6_dp
    real(dp) :: forward_map(2, 6, 2), backward_map(2, 6, 2)
    real(dp) :: forward_map_dot(2, 6, 2), backward_map_dot(2, 6, 2)
    real(dp) :: forward_map_plus(2, 6, 2), backward_map_plus(2, 6, 2)
    real(dp) :: forward_map_minus(2, 6, 2), backward_map_minus(2, 6, 2)
    real(dp) :: source_x_bar(3), source_y_bar(2)
    real(dp) :: forward_x_bar(2, 2), forward_y_bar(2, 2)
    real(dp) :: backward_x_bar(2, 2), backward_y_bar(2, 2)
    real(dp) :: left_pairing, right_pairing
    real(dp), parameter :: topology_x(2, 2) = reshape([ &
        1.0_dp, 1.5_dp, 0.5_dp, 1.25_dp], [2, 2])
    type(fortsparse_status_t) :: status

    call build_fci_bilinear_interpolation_maps_2d( &
        source_x, source_y, forward_x, forward_y, backward_x, backward_y, &
        forward_map, backward_map, status)
    call check_condition(status%code == 0, &
        "FCI map adapter accepts traced forward/backward endpoints")
    call check_condition(maxval(abs(sum(forward_map, dim=2) - 1.0_dp)) < &
        1.0e-14_dp .and. maxval(abs(sum(backward_map, dim=2) - 1.0_dp)) < &
        1.0e-14_dp, "FCI map adapter preserves partition of unity")

    call build_fci_bilinear_interpolation_maps_2d_jvp( &
        source_x, source_y, forward_x, forward_y, backward_x, backward_y, &
        source_x_dot, source_y_dot, forward_x_dot, forward_y_dot, &
        backward_x_dot, backward_y_dot, forward_map_dot, backward_map_dot, &
        status)
    call build_fci_bilinear_interpolation_maps_2d( &
        source_x + epsilon*source_x_dot, source_y + epsilon*source_y_dot, &
        forward_x + epsilon*forward_x_dot, forward_y + epsilon*forward_y_dot, &
        backward_x + epsilon*backward_x_dot, backward_y + epsilon*backward_y_dot, &
        forward_map_plus, backward_map_plus, status)
    call build_fci_bilinear_interpolation_maps_2d( &
        source_x - epsilon*source_x_dot, source_y - epsilon*source_y_dot, &
        forward_x - epsilon*forward_x_dot, forward_y - epsilon*forward_y_dot, &
        backward_x - epsilon*backward_x_dot, backward_y - epsilon*backward_y_dot, &
        forward_map_minus, backward_map_minus, status)
    call check_condition(status%code == 0 .and. &
        maxval(abs(forward_map_dot - (forward_map_plus - forward_map_minus)/ &
        (2.0_dp*epsilon))) < 1.0e-8_dp .and. &
        maxval(abs(backward_map_dot - (backward_map_plus - backward_map_minus)/ &
        (2.0_dp*epsilon))) < 1.0e-8_dp, &
        "FCI map adapter JVP matches the central-difference oracle")

    call build_fci_bilinear_interpolation_maps_2d_vjp( &
        source_x, source_y, forward_x, forward_y, backward_x, backward_y, &
        forward_map_bar, backward_map_bar, source_x_bar, source_y_bar, &
        forward_x_bar, forward_y_bar, backward_x_bar, backward_y_bar, status)
    left_pairing = sum(forward_map_bar*forward_map_dot) + &
        sum(backward_map_bar*backward_map_dot)
    right_pairing = dot_product(source_x_bar, source_x_dot) + &
        dot_product(source_y_bar, source_y_dot) + &
        sum(forward_x_bar*forward_x_dot) + sum(forward_y_bar*forward_y_dot) + &
        sum(backward_x_bar*backward_x_dot) + &
        sum(backward_y_bar*backward_y_dot)
    call check_condition(status%code == 0 .and. &
        abs(left_pairing - right_pairing) < 1.0e-12_dp, &
        "FCI map adapter VJP satisfies the full dot-product identity")

    call build_fci_bilinear_interpolation_maps_2d_jvp( &
        source_x, source_y, topology_x, forward_y, backward_x, backward_y, &
        source_x_dot, source_y_dot, forward_x_dot, forward_y_dot, &
        backward_x_dot, backward_y_dot, forward_map_dot, backward_map_dot, &
        status)
    call check_condition(status%code /= 0, &
        "FCI map adapter JVP rejects a grid-line topology event")
    call check_summary("FCI bilinear map adapter")
end program test_fci_bilinear_maps_2d
