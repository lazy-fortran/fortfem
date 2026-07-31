program test_fci_quadratic_interpolation_map
    use check, only: check_condition, check_summary
    use fortfem_api, only: build_fci_quadratic_interpolation_map_1d
    use fortfem_kinds, only: dp
    use fortsparse, only: fortsparse_status_t
    implicit none

    real(dp), parameter :: source_coordinates(5) = [ &
        0.0_dp, 1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp]
    real(dp), parameter :: target_coordinates(3) = [0.25_dp, 1.75_dp, 3.2_dp]
    integer, parameter :: stencil_indices(3, 3) = reshape([ &
        1, 2, 3, 2, 3, 4, 3, 4, 5], [3, 3])
    real(dp) :: interpolation_map(3, 5), quadratic_values(3)
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
