program test_fci_bilinear_interpolation_map
    use check, only: check_condition, check_summary
    use fortfem_feec, only: build_fci_bilinear_interpolation_map_2d
    use fortfem_kinds, only: dp
    use fortsparse, only: fortsparse_status_t
    implicit none

    real(dp), parameter :: source_x(3) = [0.0_dp, 1.0_dp, 2.0_dp]
    real(dp), parameter :: source_y(2) = [0.0_dp, 2.0_dp]
    real(dp), parameter :: target_x(2) = [0.25_dp, 1.5_dp]
    real(dp), parameter :: target_y(2) = [0.5_dp, 1.0_dp]
    real(dp), parameter :: source_values(6) = [ &
        5.0_dp, 7.0_dp, 9.0_dp, 2.0_dp, 4.0_dp, 6.0_dp]
    real(dp), parameter :: expected_values(2) = [4.75_dp, 6.5_dp]
    real(dp) :: interpolation_map(2, 6), interpolated_values(2)
    real(dp) :: bad_map(2, 5)
    real(dp), parameter :: outside_x(1) = [-0.1_dp]
    real(dp), parameter :: one_target_y(1) = [0.5_dp]
    type(fortsparse_status_t) :: status
    integer :: row

    call build_fci_bilinear_interpolation_map_2d( &
        source_x, source_y, target_x, target_y, interpolation_map, status)
    call check_condition(status%code == 0, &
        "FCI bilinear interpolation accepts in-range Cartesian coordinates")
    interpolated_values = matmul(interpolation_map, source_values)
    call check_condition(maxval(abs(interpolated_values - expected_values)) < &
        1.0e-14_dp, "FCI bilinear interpolation reproduces an affine field")
    do row = 1, size(target_x)
        call check_condition(abs(sum(interpolation_map(row, :)) - 1.0_dp) < &
            1.0e-14_dp, "FCI bilinear weights form a partition of unity")
    end do

    call build_fci_bilinear_interpolation_map_2d( &
        source_x, source_y, outside_x, one_target_y, interpolation_map(1:1, :), &
        status)
    call check_condition(status%code /= 0, &
        "FCI bilinear interpolation rejects targets outside the source box")
    call build_fci_bilinear_interpolation_map_2d( &
        source_x, source_y, target_x, target_y, bad_map, status)
    call check_condition(status%code /= 0, &
        "FCI bilinear interpolation rejects an incompatible output shape")
    call check_summary("FCI bilinear interpolation map")
end program test_fci_bilinear_interpolation_map
