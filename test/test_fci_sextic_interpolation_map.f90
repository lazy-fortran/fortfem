program test_fci_sextic_interpolation_map
    use check, only: check_condition, check_summary
    use fortfem_api, only: &
        build_fci_sextic_interpolation_map_1d, &
        build_fci_sextic_interpolation_map_1d_jvp, &
        build_fci_sextic_interpolation_map_1d_vjp
    use fortfem_kinds, only: dp
    use fortsparse, only: fortsparse_status_t
    implicit none

    real(dp), parameter :: source_coordinates(9) = [ &
        0.0_dp, 0.7_dp, 1.6_dp, 2.8_dp, 4.1_dp, 5.7_dp, 7.4_dp, 9.2_dp, &
        11.3_dp]
    real(dp), parameter :: target_coordinates(3) = [0.3_dp, 2.2_dp, 6.3_dp]
    integer, parameter :: stencil_indices(7, 3) = reshape([ &
        1, 2, 3, 4, 5, 6, 7, 2, 3, 4, 5, 6, 7, 8, &
        3, 4, 5, 6, 7, 8, 9], [7, 3])
    integer, parameter :: bad_stencils(7, 3) = reshape([ &
        1, 2, 3, 4, 5, 6, 7, 2, 3, 3, 5, 6, 7, 8, &
        3, 4, 5, 6, 7, 8, 9], [7, 3])
    real(dp), parameter :: step = 2.0e-7_dp
    real(dp) :: interpolation_map(3, 9), interpolation_map_dot(3, 9)
    real(dp) :: interpolation_map_bar(3, 9)
    real(dp) :: interpolation_map_plus(3, 9), interpolation_map_minus(3, 9)
    real(dp) :: source_coordinates_dot(9), source_coordinates_bar(9)
    real(dp) :: target_coordinates_dot(3), target_coordinates_bar(3)
    real(dp) :: expected(3, 9), weights(7), sextic_values(3), lhs, rhs
    real(dp) :: repeated_source_coordinates(9)
    type(fortsparse_status_t) :: status
    integer :: row, node

    call build_fci_sextic_interpolation_map_1d( &
        source_coordinates, target_coordinates, stencil_indices, &
        interpolation_map, status)
    call check_condition(status%code == 0, &
        "FCI sextic interpolation accepts fixed stencils")
    expected = 0.0_dp
    do row = 1, size(target_coordinates)
        call sextic_weights( &
            source_coordinates(stencil_indices(:, row)), &
            target_coordinates(row), weights)
        do node = 1, 7
            expected(row, stencil_indices(node, row)) = weights(node)
        end do
    end do
    call check_condition(maxval(abs(interpolation_map - expected)) < 5.0e-13_dp, &
        "FCI sextic interpolation matches the independent Lagrange oracle")
    call check_condition(maxval(abs(sum(interpolation_map, dim=2) - 1.0_dp)) < &
        5.0e-13_dp, "FCI sextic weights preserve partition of unity")
    sextic_values = matmul(interpolation_map, source_coordinates**6)
    call check_condition(maxval(abs(sextic_values - target_coordinates**6)) < &
        5.0e-9_dp, "FCI sextic weights reproduce a sextic field")

    source_coordinates_dot = [0.01_dp, -0.02_dp, 0.015_dp, -0.01_dp, 0.02_dp, &
        -0.015_dp, 0.012_dp, -0.009_dp, 0.011_dp]
    target_coordinates_dot = [0.03_dp, -0.02_dp, 0.01_dp]
    call build_fci_sextic_interpolation_map_1d_jvp( &
        source_coordinates, target_coordinates, stencil_indices, &
        source_coordinates_dot, target_coordinates_dot, interpolation_map_dot, &
        status)
    call build_fci_sextic_interpolation_map_1d( &
        source_coordinates + step*source_coordinates_dot, &
        target_coordinates + step*target_coordinates_dot, stencil_indices, &
        interpolation_map_plus, status)
    call build_fci_sextic_interpolation_map_1d( &
        source_coordinates - step*source_coordinates_dot, &
        target_coordinates - step*target_coordinates_dot, stencil_indices, &
        interpolation_map_minus, status)
    call check_condition(maxval(abs(interpolation_map_dot - &
        (interpolation_map_plus - interpolation_map_minus)/(2.0_dp*step))) < &
        5.0e-6_dp, "FCI sextic JVP matches central differences")

    interpolation_map_bar = reshape([ &
        0.4_dp, -0.3_dp, 0.2_dp, -0.5_dp, 0.6_dp, -0.7_dp, 0.8_dp, -0.9_dp, &
        1.0_dp, 1.1_dp, -1.2_dp, 1.3_dp, -1.4_dp, 1.5_dp, -1.6_dp, 1.7_dp, &
        -1.8_dp, 1.9_dp, -2.0_dp, 2.1_dp, -2.2_dp, 2.3_dp, -2.4_dp, 2.5_dp, &
        -2.6_dp, 2.7_dp, -2.8_dp], [3, 9])
    call build_fci_sextic_interpolation_map_1d_vjp( &
        source_coordinates, target_coordinates, stencil_indices, &
        interpolation_map_bar, source_coordinates_bar, target_coordinates_bar, &
        status)
    lhs = sum(interpolation_map_bar*interpolation_map_dot)
    rhs = dot_product(source_coordinates_bar, source_coordinates_dot) + &
        dot_product(target_coordinates_bar, target_coordinates_dot)
    call check_condition(abs(lhs - rhs) < &
        5.0e-10_dp*max(1.0_dp, abs(lhs), abs(rhs)), &
        "FCI sextic VJP satisfies the real dot-product identity")

    call build_fci_sextic_interpolation_map_1d( &
        source_coordinates, target_coordinates, bad_stencils, interpolation_map, &
        status)
    call check_condition(status%code /= 0, &
        "FCI sextic interpolation rejects repeated stencil nodes")
    repeated_source_coordinates = source_coordinates
    repeated_source_coordinates(1) = repeated_source_coordinates(3)
    call build_fci_sextic_interpolation_map_1d( &
        repeated_source_coordinates, target_coordinates, stencil_indices, &
        interpolation_map, status)
    call check_condition(status%code /= 0, &
        "FCI sextic interpolation rejects repeated source coordinates")
    call check_summary("FCI sextic interpolation map")

contains

    pure subroutine sextic_weights(nodes, target, weights)
        real(dp), intent(in) :: nodes(7), target
        real(dp), intent(out) :: weights(7)
        integer :: node, other

        do node = 1, 7
            weights(node) = 1.0_dp
            do other = 1, 7
                if (other == node) cycle
                weights(node) = weights(node)*(target - nodes(other))/ &
                    (nodes(node) - nodes(other))
            end do
        end do
    end subroutine sextic_weights

end program test_fci_sextic_interpolation_map
