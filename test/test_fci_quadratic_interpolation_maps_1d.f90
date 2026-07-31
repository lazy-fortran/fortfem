program test_fci_quadratic_interpolation_maps_1d
    use check, only: check_condition, check_summary
    use fortfem_api, only: &
        build_fci_quadratic_interpolation_maps_1d, &
        build_fci_quadratic_interpolation_maps_1d_jvp, &
        build_fci_quadratic_interpolation_maps_1d_vjp
    use fortfem_kinds, only: dp
    use fortsparse, only: fortsparse_status_t
    implicit none

    real(dp), parameter :: source_coordinates(5) = [ &
        0.0_dp, 1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp]
    real(dp), parameter :: target_coordinates(3, 2) = reshape([ &
        0.25_dp, 1.75_dp, 3.2_dp, 0.4_dp, 1.6_dp, 3.4_dp], [3, 2])
    integer, parameter :: stencil_indices(3, 3, 2) = reshape([ &
        1, 2, 3, 2, 3, 4, 3, 4, 5, &
        1, 2, 3, 2, 3, 4, 3, 4, 5], [3, 3, 2])
    real(dp), parameter :: source_coordinates_dot(5) = [ &
        0.01_dp, -0.02_dp, 0.015_dp, -0.01_dp, 0.02_dp]
    real(dp), parameter :: target_coordinates_dot(3, 2) = reshape([ &
        0.03_dp, -0.02_dp, 0.01_dp, -0.01_dp, 0.02_dp, -0.03_dp], [3, 2])
    real(dp), parameter :: step = 2.0e-7_dp
    real(dp) :: interpolation_maps(3, 5, 2), interpolation_maps_dot(3, 5, 2)
    real(dp) :: interpolation_maps_plus(3, 5, 2), interpolation_maps_minus(3, 5, 2)
    real(dp) :: interpolation_maps_bar(3, 5, 2)
    real(dp) :: source_coordinates_bar(5), target_coordinates_bar(3, 2)
    real(dp) :: expected(3, 5), weights(3), values(3), lhs, rhs
    real(dp) :: source_coordinates_bad(5)
    type(fortsparse_status_t) :: status
    integer :: segment, row, node

    call build_fci_quadratic_interpolation_maps_1d( &
        source_coordinates, target_coordinates, stencil_indices, &
        interpolation_maps, status)
    call check_condition(status%code == 0, &
        "FCI batched quadratic interpolation accepts fixed segment stencils")

    do segment = 1, 2
        expected = 0.0_dp
        do row = 1, 3
            call quadratic_weights( &
                source_coordinates(stencil_indices(:, row, segment)), &
                target_coordinates(row, segment), weights)
            do node = 1, 3
                expected(row, stencil_indices(node, row, segment)) = weights(node)
            end do
        end do
        call check_condition(maxval(abs(interpolation_maps(:, :, segment) - &
            expected)) < 1.0e-14_dp, &
            "FCI batched quadratic map matches the independent Lagrange oracle")
        call check_condition(maxval(abs(sum(interpolation_maps(:, :, segment), &
            dim=2) - 1.0_dp)) < 1.0e-14_dp, &
            "FCI batched quadratic weights preserve partition of unity")
        values = matmul(interpolation_maps(:, :, segment), source_coordinates**2)
        call check_condition(maxval(abs(values - target_coordinates(:, segment)**2)) < &
            1.0e-14_dp, "FCI batched quadratic map reproduces a quadratic field")
    end do

    call build_fci_quadratic_interpolation_maps_1d_jvp( &
        source_coordinates, target_coordinates, stencil_indices, &
        source_coordinates_dot, target_coordinates_dot, interpolation_maps_dot, &
        status)
    call build_fci_quadratic_interpolation_maps_1d( &
        source_coordinates + step*source_coordinates_dot, &
        target_coordinates + step*target_coordinates_dot, stencil_indices, &
        interpolation_maps_plus, status)
    call build_fci_quadratic_interpolation_maps_1d( &
        source_coordinates - step*source_coordinates_dot, &
        target_coordinates - step*target_coordinates_dot, stencil_indices, &
        interpolation_maps_minus, status)
    call check_condition(maxval(abs(interpolation_maps_dot - &
        (interpolation_maps_plus - interpolation_maps_minus)/(2.0_dp*step))) < &
        3.0e-9_dp, "FCI batched quadratic JVP matches central differences")

    interpolation_maps_bar = reshape([ &
        0.4_dp, -0.3_dp, 0.2_dp, -0.5_dp, 0.6_dp, -0.7_dp, &
        0.8_dp, -0.9_dp, 1.0_dp, 1.1_dp, -1.2_dp, 1.3_dp, &
        -1.4_dp, 1.5_dp, -1.6_dp, 0.2_dp, -0.4_dp, 0.6_dp, &
        -0.8_dp, 1.0_dp, -1.2_dp, 1.4_dp, -1.6_dp, 1.8_dp, &
        -2.0_dp, 2.2_dp, -2.4_dp, 2.6_dp, -2.8_dp, 3.0_dp], [3, 5, 2])
    call build_fci_quadratic_interpolation_maps_1d_vjp( &
        source_coordinates, target_coordinates, stencil_indices, &
        interpolation_maps_bar, source_coordinates_bar, target_coordinates_bar, &
        status)
    lhs = sum(interpolation_maps_bar*interpolation_maps_dot)
    rhs = dot_product(source_coordinates_bar, source_coordinates_dot) + &
        sum(target_coordinates_bar*target_coordinates_dot)
    call check_condition(abs(lhs - rhs) < &
        3.0e-13_dp*max(1.0_dp, abs(lhs), abs(rhs)), &
        "FCI batched quadratic VJP satisfies the real dot-product identity")

    source_coordinates_bad = source_coordinates
    source_coordinates_bad(3) = source_coordinates_bad(2)
    call build_fci_quadratic_interpolation_maps_1d( &
        source_coordinates_bad, target_coordinates, stencil_indices, &
        interpolation_maps, status)
    call check_condition(status%code /= 0, &
        "FCI batched quadratic map rejects a degenerate local stencil")
    call check_summary("FCI batched quadratic interpolation maps")

contains

    pure subroutine quadratic_weights(nodes, target, weights)
        real(dp), intent(in) :: nodes(3), target
        real(dp), intent(out) :: weights(3)

        weights(1) = (target - nodes(2))*(target - nodes(3))/ &
            ((nodes(1) - nodes(2))*(nodes(1) - nodes(3)))
        weights(2) = (target - nodes(1))*(target - nodes(3))/ &
            ((nodes(2) - nodes(1))*(nodes(2) - nodes(3)))
        weights(3) = (target - nodes(1))*(target - nodes(2))/ &
            ((nodes(3) - nodes(1))*(nodes(3) - nodes(2)))
    end subroutine quadratic_weights

end program test_fci_quadratic_interpolation_maps_1d
