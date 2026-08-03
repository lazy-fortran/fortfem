program test_fci_quadratic_interpolation_map
    use check, only: check_condition, check_summary
    use fortfem_feec, only: build_fci_quadratic_interpolation_map_1d, &
        build_fci_quadratic_interpolation_map_1d_jvp, &
        build_fci_quadratic_interpolation_map_1d_vjp
    use fortfem_kinds, only: dp
    use fortsparse, only: fortsparse_status_t
    implicit none

    real(dp), parameter :: source_coordinates(5) = [ &
        0.0_dp, 1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp]
    real(dp), parameter :: target_coordinates(3) = [0.25_dp, 1.75_dp, 3.2_dp]
    integer, parameter :: stencil_indices(3, 3) = reshape([ &
        1, 2, 3, 2, 3, 4, 3, 4, 5], [3, 3])
    real(dp) :: interpolation_map(3, 5), interpolation_map_dot(3, 5)
    real(dp) :: interpolation_map_bar(3, 5)
    real(dp) :: interpolation_map_minus(3, 5), interpolation_map_plus(3, 5)
    real(dp) :: source_coordinates_dot(5), source_coordinates_bar(5)
    real(dp) :: target_coordinates_dot(3), target_coordinates_bar(3)
    real(dp) :: quadratic_values(3), lhs, rhs
    real(dp) :: expected(3, 5), weights(3)
    integer, parameter :: bad_stencils(3, 3) = reshape([ &
        1, 2, 3, 2, 2, 4, 3, 4, 5], [3, 3])
    type(fortsparse_status_t) :: status
    integer :: row, node

    call build_fci_quadratic_interpolation_map_1d( &
        source_coordinates, target_coordinates, stencil_indices, &
        interpolation_map, status)
    call check_condition(status%code == 0, &
        "FCI quadratic interpolation accepts fixed stencils")
    expected = 0.0_dp
    do row = 1, size(target_coordinates)
        call quadratic_weights( &
            source_coordinates(stencil_indices(:, row)), &
            target_coordinates(row), weights)
        do node = 1, 3
            expected(row, stencil_indices(node, row)) = weights(node)
        end do
    end do
    call check_condition(maxval(abs(interpolation_map - expected)) < 1.0e-14_dp, &
        "FCI quadratic interpolation matches the independent Lagrange oracle")
    call check_condition(maxval(abs(sum(interpolation_map, dim=2) - 1.0_dp)) < &
        1.0e-14_dp, "FCI quadratic weights preserve partition of unity")
    quadratic_values = matmul(interpolation_map, source_coordinates**2)
    call check_condition(maxval(abs(quadratic_values - &
        target_coordinates**2)) < 1.0e-14_dp, &
        "FCI quadratic weights reproduce a quadratic field")

    source_coordinates_dot = [0.01_dp, -0.02_dp, 0.015_dp, -0.01_dp, 0.02_dp]
    target_coordinates_dot = [0.03_dp, -0.02_dp, 0.01_dp]
    call build_fci_quadratic_interpolation_map_1d_jvp( &
        source_coordinates, target_coordinates, stencil_indices, &
        source_coordinates_dot, target_coordinates_dot, interpolation_map_dot, &
        status)
    call build_fci_quadratic_interpolation_map_1d( &
        source_coordinates + 2.0e-7_dp*source_coordinates_dot, &
        target_coordinates + 2.0e-7_dp*target_coordinates_dot, stencil_indices, &
        interpolation_map_plus, status)
    call build_fci_quadratic_interpolation_map_1d( &
        source_coordinates - 2.0e-7_dp*source_coordinates_dot, &
        target_coordinates - 2.0e-7_dp*target_coordinates_dot, stencil_indices, &
        interpolation_map_minus, status)
    call check_condition(maxval(abs(interpolation_map_dot - &
        (interpolation_map_plus - interpolation_map_minus)/(4.0e-7_dp))) < &
        3.0e-9_dp, "FCI quadratic JVP matches central differences")

    interpolation_map_bar = reshape([ &
        0.4_dp, -0.3_dp, 0.2_dp, -0.5_dp, 0.6_dp, -0.7_dp, &
        0.8_dp, -0.9_dp, 1.0_dp, 1.1_dp, -1.2_dp, 1.3_dp, &
        -1.4_dp, 1.5_dp, -1.6_dp], [3, 5])
    call build_fci_quadratic_interpolation_map_1d_vjp( &
        source_coordinates, target_coordinates, stencil_indices, &
        interpolation_map_bar, source_coordinates_bar, target_coordinates_bar, &
        status)
    lhs = sum(interpolation_map_bar*interpolation_map_dot)
    rhs = dot_product(source_coordinates_bar, source_coordinates_dot) + &
        dot_product(target_coordinates_bar, target_coordinates_dot)
    call check_condition(abs(lhs - rhs) < &
        3.0e-13_dp*max(1.0_dp, abs(lhs), abs(rhs)), &
        "FCI quadratic VJP satisfies the real dot-product identity")

    call build_fci_quadratic_interpolation_map_1d( &
        source_coordinates, target_coordinates, bad_stencils, interpolation_map, &
        status)
    call check_condition(status%code /= 0, &
        "FCI quadratic interpolation rejects repeated stencil nodes")
    call check_summary("FCI quadratic interpolation map")

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

end program test_fci_quadratic_interpolation_map
