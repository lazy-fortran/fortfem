program test_fci_support_volumes
    use check, only: check_condition, check_summary
    use fortfem_api, only: compute_fci_staggered_flux_box_volumes
    use fortfem_kinds, only: dp
    use fortsparse, only: fortsparse_status_t
    implicit none

    real(dp), parameter :: forward_expansion(3) = [0.2_dp, 0.5_dp, 0.7_dp]
    real(dp), parameter :: backward_expansion(3) = [0.3_dp, 0.1_dp, 0.4_dp]
    real(dp), parameter :: base_cell_area(3) = [2.0_dp, 3.0_dp, 4.0_dp]
    real(dp), parameter :: toroidal_field(3) = [5.0_dp, 2.0_dp, 1.5_dp]
    real(dp), parameter :: expected_volumes(3) = [5.0_dp, 3.6_dp, 6.6_dp]
    real(dp) :: staggered_volumes(3)
    real(dp), parameter :: bad_field(2) = [1.0_dp, 2.0_dp]
    type(fortsparse_status_t) :: status

    call compute_fci_staggered_flux_box_volumes( &
        forward_expansion, backward_expansion, base_cell_area, toroidal_field, &
        staggered_volumes, status)
    call check_condition(status%code == 0, &
        "FCI support volumes accept positive geometric factors")
    call check_condition(maxval(abs(staggered_volumes - expected_volumes)) < &
        1.0e-14_dp, "FCI support volumes match the flux-expansion oracle")

    call compute_fci_staggered_flux_box_volumes( &
        forward_expansion, bad_field, base_cell_area, toroidal_field, &
        staggered_volumes, status)
    call check_condition(status%code /= 0, &
        "FCI support volumes reject incompatible array sizes")
    call check_summary("FCI support volumes")
end program test_fci_support_volumes
