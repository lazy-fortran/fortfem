program test_fci_cubic_interpolation_map
    use check, only: check_condition, check_summary
    use fortfem_feec, only: build_fci_cubic_interpolation_map_1d, &
        build_fci_cubic_interpolation_map_1d_jvp, &
        build_fci_cubic_interpolation_map_1d_vjp
    use fortfem_kinds, only: dp
    use fortsparse, only: fortsparse_status_t
    implicit none

    real(dp), parameter :: source_coordinates(6) = [ &
        0.0_dp, 1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp]
    real(dp), parameter :: target_coordinates(3) = [0.3_dp, 1.7_dp, 3.4_dp]
    integer, parameter :: stencil_indices(4, 3) = reshape([ &
        1, 2, 3, 4, 2, 3, 4, 5, 3, 4, 5, 6], [4, 3])
    real(dp), parameter :: step = 2.0e-7_dp
    real(dp) :: interpolation_map(3, 6), interpolation_map_dot(3, 6)
    real(dp) :: interpolation_map_bar(3, 6)
    real(dp) :: interpolation_map_plus(3, 6), interpolation_map_minus(3, 6)
    real(dp) :: source_coordinates_dot(6), source_coordinates_bar(6)
    real(dp) :: repeated_source_coordinates(6)
    real(dp) :: target_coordinates_dot(3), target_coordinates_bar(3)
    real(dp) :: expected(3, 6), weights(4), cubic_values(3), lhs, rhs
    integer, parameter :: bad_stencils(4, 3) = reshape([ &
        1, 2, 3, 4, 2, 3, 3, 5, 3, 4, 5, 6], [4, 3])
    type(fortsparse_status_t) :: status
    integer :: row, node
    logical :: all_passed

    all_passed = .true.

    call build_fci_cubic_interpolation_map_1d( &
        source_coordinates, target_coordinates, stencil_indices, &
        interpolation_map, status)
    call record_condition(status%code == 0, &
        "FCI cubic interpolation accepts fixed stencils")
    expected = 0.0_dp
    do row = 1, size(target_coordinates)
        call cubic_weights( &
            source_coordinates(stencil_indices(:, row)), &
            target_coordinates(row), weights)
        do node = 1, 4
            expected(row, stencil_indices(node, row)) = weights(node)
        end do
    end do
    call record_condition(maxval(abs(interpolation_map - expected)) < 1.0e-14_dp, &
        "FCI cubic interpolation matches the independent Lagrange oracle")
    call record_condition(maxval(abs(sum(interpolation_map, dim=2) - 1.0_dp)) < &
        1.0e-14_dp, "FCI cubic weights preserve partition of unity")
    cubic_values = matmul(interpolation_map, source_coordinates**3)
    call record_condition(maxval(abs(cubic_values - target_coordinates**3)) < &
        1.0e-14_dp, "FCI cubic weights reproduce a cubic field")

    source_coordinates_dot = [0.01_dp, -0.02_dp, 0.015_dp, -0.01_dp, 0.02_dp, &
        -0.015_dp]
    target_coordinates_dot = [0.03_dp, -0.02_dp, 0.01_dp]
    call build_fci_cubic_interpolation_map_1d_jvp( &
        source_coordinates, target_coordinates, stencil_indices, &
        source_coordinates_dot, target_coordinates_dot, interpolation_map_dot, &
        status)
    call build_fci_cubic_interpolation_map_1d( &
        source_coordinates + step*source_coordinates_dot, &
        target_coordinates + step*target_coordinates_dot, stencil_indices, &
        interpolation_map_plus, status)
    call build_fci_cubic_interpolation_map_1d( &
        source_coordinates - step*source_coordinates_dot, &
        target_coordinates - step*target_coordinates_dot, stencil_indices, &
        interpolation_map_minus, status)
    call record_condition(maxval(abs(interpolation_map_dot - &
        (interpolation_map_plus - interpolation_map_minus)/(2.0_dp*step))) < &
        4.0e-9_dp, "FCI cubic JVP matches central differences")

    interpolation_map_bar = reshape([ &
        0.4_dp, -0.3_dp, 0.2_dp, -0.5_dp, 0.6_dp, -0.7_dp, &
        0.8_dp, -0.9_dp, 1.0_dp, 1.1_dp, -1.2_dp, 1.3_dp, &
        -1.4_dp, 1.5_dp, -1.6_dp, 1.7_dp, -1.8_dp, 1.9_dp], [3, 6])
    call build_fci_cubic_interpolation_map_1d_vjp( &
        source_coordinates, target_coordinates, stencil_indices, &
        interpolation_map_bar, source_coordinates_bar, target_coordinates_bar, &
        status)
    lhs = sum(interpolation_map_bar*interpolation_map_dot)
    rhs = dot_product(source_coordinates_bar, source_coordinates_dot) + &
        dot_product(target_coordinates_bar, target_coordinates_dot)
    call record_condition(abs(lhs - rhs) < &
        4.0e-13_dp*max(1.0_dp, abs(lhs), abs(rhs)), &
        "FCI cubic VJP satisfies the real dot-product identity")

    call build_fci_cubic_interpolation_map_1d( &
        source_coordinates, target_coordinates, bad_stencils, interpolation_map, &
        status)
    call record_condition(status%code /= 0, &
        "FCI cubic interpolation rejects repeated stencil nodes")
    repeated_source_coordinates = source_coordinates
    repeated_source_coordinates(1) = repeated_source_coordinates(3)
    call build_fci_cubic_interpolation_map_1d( &
        repeated_source_coordinates, target_coordinates, stencil_indices, &
        interpolation_map, status)
    call record_condition(status%code /= 0, &
        "FCI cubic interpolation rejects repeated source coordinates")
    call check_summary("FCI cubic interpolation map")
    if (.not. all_passed) error stop 1

contains

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

    pure subroutine cubic_weights(nodes, target, weights)
        real(dp), intent(in) :: nodes(4), target
        real(dp), intent(out) :: weights(4)
        integer :: node, other

        do node = 1, 4
            weights(node) = 1.0_dp
            do other = 1, 4
                if (other == node) cycle
                weights(node) = weights(node)*(target - nodes(other))/ &
                    (nodes(node) - nodes(other))
            end do
        end do
    end subroutine cubic_weights

end program test_fci_cubic_interpolation_map
