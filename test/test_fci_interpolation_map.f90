program test_fci_interpolation_map
    use check, only: check_condition, check_summary
    use fortfem_api, only: build_fci_linear_interpolation_map_1d
    use fortfem_kinds, only: dp
    use fortsparse, only: fortsparse_status_t
    implicit none

    real(dp), parameter :: source_coordinates(3) = [0.0_dp, 1.0_dp, 2.0_dp]
    real(dp), parameter :: target_coordinates(5) = [ &
        0.0_dp, 0.25_dp, 1.0_dp, 1.5_dp, 2.0_dp]
    real(dp), parameter :: source_values(3) = [3.0_dp, 5.0_dp, 7.0_dp]
    real(dp), parameter :: expected_values(5) = [ &
        3.0_dp, 3.5_dp, 5.0_dp, 6.0_dp, 7.0_dp]
    real(dp) :: interpolation_map(5, 3), interpolated_values(5)
    real(dp) :: bad_map(4, 3)
    real(dp), parameter :: outside_target(1) = [-0.1_dp]
    real(dp), parameter :: repeated_source(3) = [0.0_dp, 1.0_dp, 1.0_dp]
    type(fortsparse_status_t) :: status
    integer :: row

    call build_fci_linear_interpolation_map_1d( &
        source_coordinates, target_coordinates, interpolation_map, status)
    call check_condition(status%code == 0, &
        "FCI linear interpolation accepts in-range monotone coordinates")
    interpolated_values = matmul(interpolation_map, source_values)
    call check_condition(maxval(abs(interpolated_values - expected_values)) < &
        1.0e-14_dp, "FCI interpolation reproduces an affine field")
    do row = 1, size(target_coordinates)
        call check_condition(abs(sum(interpolation_map(row, :)) - 1.0_dp) < &
            1.0e-14_dp, "FCI interpolation weights form a partition of unity")
    end do

    call build_fci_linear_interpolation_map_1d( &
        source_coordinates, outside_target, interpolation_map(1:1, :), status)
    call check_condition(status%code /= 0, &
        "FCI interpolation rejects targets outside the source interval")
    call build_fci_linear_interpolation_map_1d( &
        repeated_source, target_coordinates(1:1), interpolation_map(1:1, :), &
        status)
    call check_condition(status%code /= 0, &
        "FCI interpolation rejects repeated source coordinates")
    call build_fci_linear_interpolation_map_1d( &
        source_coordinates, target_coordinates, bad_map, status)
    call check_condition(status%code /= 0, &
        "FCI interpolation rejects an incompatible output shape")
    call check_summary("FCI linear interpolation map")
end program test_fci_interpolation_map
