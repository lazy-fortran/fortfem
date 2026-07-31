program test_fci_interpolation_map_ad
    use check, only: check_condition, check_summary
    use fortfem_api, only: &
        build_fci_linear_interpolation_map_1d_jvp, &
        build_fci_linear_interpolation_map_1d_vjp
    use fortfem_kinds, only: dp
    use fortsparse, only: fortsparse_status_t
    implicit none

    real(dp), parameter :: source_coordinates(3) = [0.0_dp, 1.0_dp, 2.0_dp]
    real(dp), parameter :: target_coordinates(2) = [0.25_dp, 1.5_dp]
    real(dp), parameter :: source_dot(3) = [0.1_dp, -0.2_dp, 0.3_dp]
    real(dp), parameter :: target_dot(2) = [0.05_dp, -0.04_dp]
    real(dp), parameter :: expected_map_dot(2, 3) = reshape([ &
        -0.025_dp, 0.0_dp, 0.025_dp, &
        0.09_dp, 0.0_dp, -0.09_dp], [2, 3])
    real(dp), parameter :: map_bar(2, 3) = reshape([ &
        0.7_dp, -0.3_dp, 0.2_dp, &
        -0.1_dp, 0.6_dp, 0.4_dp], [2, 3])
    real(dp) :: map_dot(2, 3), source_bar(3), target_bar(2)
    real(dp), parameter :: endpoint_target(1) = [0.0_dp]
    real(dp) :: endpoint_map_dot(1, 3), endpoint_source_bar(3)
    real(dp) :: endpoint_target_bar(1)
    type(fortsparse_status_t) :: status
    real(dp) :: left, right

    call build_fci_linear_interpolation_map_1d_jvp( &
        source_coordinates, target_coordinates, source_dot, target_dot, &
        map_dot, status)
    call check_condition(status%code == 0, &
        "FCI interpolation JVP accepts a fixed stencil topology")
    call check_condition(maxval(abs(map_dot - expected_map_dot)) < 1.0e-14_dp, &
        "FCI interpolation JVP matches the quotient-rule oracle")

    call build_fci_linear_interpolation_map_1d_vjp( &
        source_coordinates, target_coordinates, map_bar, source_bar, &
        target_bar, status)
    left = sum(map_bar*map_dot)
    right = dot_product(source_bar, source_dot) + &
        dot_product(target_bar, target_dot)
    call check_condition(status%code == 0 .and. &
        abs(left - right) < 1.0e-14_dp, &
        "FCI interpolation VJP satisfies the real dot-product identity")

    call build_fci_linear_interpolation_map_1d_jvp( &
        source_coordinates, endpoint_target, source_dot, [0.0_dp], &
        endpoint_map_dot, status)
    call check_condition(status%code /= 0, &
        "FCI interpolation JVP rejects a stencil-topology boundary")
    call build_fci_linear_interpolation_map_1d_vjp( &
        source_coordinates, endpoint_target, endpoint_map_dot, &
        endpoint_source_bar, endpoint_target_bar, status)
    call check_condition(status%code /= 0, &
        "FCI interpolation VJP rejects a stencil-topology boundary")
    call check_summary("FCI linear interpolation derivatives")
end program test_fci_interpolation_map_ad
