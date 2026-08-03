program test_fci_support_volumes
    use check, only: check_condition, check_summary
    use fortfem_feec, only: &
        compute_fci_staggered_flux_box_volumes, &
        compute_fci_staggered_flux_box_volumes_jvp, &
        compute_fci_staggered_flux_box_volumes_vjp
    use fortfem_kinds, only: dp
    use fortsparse, only: fortsparse_status_t
    implicit none

    real(dp), parameter :: forward_expansion(3) = [0.2_dp, 0.5_dp, 0.7_dp]
    real(dp), parameter :: backward_expansion(3) = [0.3_dp, 0.1_dp, 0.4_dp]
    real(dp), parameter :: base_cell_area(3) = [2.0_dp, 3.0_dp, 4.0_dp]
    real(dp), parameter :: toroidal_field(3) = [5.0_dp, 2.0_dp, 1.5_dp]
    real(dp), parameter :: expected_volumes(3) = [5.0_dp, 3.6_dp, 6.6_dp]
    real(dp), parameter :: forward_dot(3) = [0.01_dp, -0.02_dp, 0.03_dp]
    real(dp), parameter :: backward_dot(3) = [-0.04_dp, 0.05_dp, 0.02_dp]
    real(dp), parameter :: area_dot(3) = [0.1_dp, -0.2_dp, 0.3_dp]
    real(dp), parameter :: field_dot(3) = [0.2_dp, 0.1_dp, -0.1_dp]
    real(dp), parameter :: expected_volume_dot(3) = [0.15_dp, 0.12_dp, 0.355_dp]
    real(dp), parameter :: volume_bar(3) = [0.4_dp, -0.3_dp, 0.8_dp]
    real(dp) :: staggered_volumes(3)
    real(dp) :: staggered_volume_dot(3)
    real(dp) :: forward_bar(3), backward_bar(3), area_bar(3), field_bar(3)
    real(dp) :: vjp_left, vjp_right
    real(dp), parameter :: bad_field(2) = [1.0_dp, 2.0_dp]
    type(fortsparse_status_t) :: status

    call compute_fci_staggered_flux_box_volumes( &
        forward_expansion, backward_expansion, base_cell_area, toroidal_field, &
        staggered_volumes, status)
    call check_condition(status%code == 0, &
        "FCI support volumes accept positive geometric factors")
    call check_condition(maxval(abs(staggered_volumes - expected_volumes)) < &
        1.0e-14_dp, "FCI support volumes match the flux-expansion oracle")

    call compute_fci_staggered_flux_box_volumes_jvp( &
        forward_expansion, backward_expansion, base_cell_area, toroidal_field, &
        forward_dot, backward_dot, area_dot, field_dot, staggered_volume_dot, &
        status)
    call check_condition(status%code == 0 .and. &
        maxval(abs(staggered_volume_dot - expected_volume_dot)) < 1.0e-14_dp, &
        "FCI support-volume JVP matches the product-rule oracle")

    call compute_fci_staggered_flux_box_volumes_vjp( &
        forward_expansion, backward_expansion, base_cell_area, toroidal_field, &
        volume_bar, forward_bar, backward_bar, area_bar, field_bar, status)
    vjp_left = dot_product(volume_bar, staggered_volume_dot)
    vjp_right = dot_product(forward_bar, forward_dot) + &
        dot_product(backward_bar, backward_dot) + dot_product(area_bar, area_dot) + &
        dot_product(field_bar, field_dot)
    call check_condition(status%code == 0 .and. &
        abs(vjp_left - vjp_right) < 1.0e-14_dp, &
        "FCI support-volume VJP satisfies the real dot-product identity")

    call compute_fci_staggered_flux_box_volumes( &
        forward_expansion, bad_field, base_cell_area, toroidal_field, &
        staggered_volumes, status)
    call check_condition(status%code /= 0, &
        "FCI support volumes reject incompatible array sizes")
    call check_summary("FCI support volumes")
end program test_fci_support_volumes
