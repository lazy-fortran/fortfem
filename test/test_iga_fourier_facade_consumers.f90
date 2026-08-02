program test_iga_fourier_facade_consumers
    !! Compile and behavior gate for the IGA/Fourier canonical facades.
    !!
    !! The example clients intentionally import these namespaces directly;
    !! this test keeps the migration honest without depending on the legacy
    !! compatibility umbrella.
    use check, only: check_condition, check_summary
    use fortfem_core, only: cartesian_to_toroidal, &
        evaluate_axis_regular_radial_basis, evaluate_nested_surface_geometry, &
        toroidal_point_to_cartesian
    use fortfem_feec, only: &
        apply_fci_parallel_diffusion, &
        apply_fci_parallel_gradient, &
        assemble_bspline_h1_operator_csc, assemble_tensor_diffusion_matrix, &
        compute_fci_curved_polygon_cell_areas_2d, &
        compute_fci_curved_quadrilateral_cell_areas_2d, &
        compute_fci_cubic_curved_polygon_cell_areas_2d, &
        compute_fci_polygon_cell_areas_2d, &
        compute_fci_quadrilateral_cell_areas_2d, &
        compute_fci_quartic_curved_polygon_cell_areas_2d, &
        compute_fci_quintic_curved_polygon_cell_areas_2d, &
        compute_fci_sextic_curved_polygon_cell_areas_2d, &
        evaluate_bspline_basis, evaluate_field_aligned_flux, &
        trace_fci_field_line_rk4
    use fortfem_fourier, only: &
        scalar_reluctivity_curvilinear_fourier_coefficients, toroidal_p
    use fortfem_kinds, only: dp
    use fortfem_time, only: &
        advance_bspline_jorek_poloidal_flux_midpoint_steps
    use fortsparse, only: fortsparse_status_t
    implicit none

    real(dp), parameter :: direction(3) = [0.0_dp, 0.0_dp, 1.0_dp]
    real(dp), parameter :: gradient(3) = [1.0_dp, -2.0_dp, 3.0_dp]
    real(dp), parameter :: expected_flux(3) = [1.0_dp, -2.0_dp, 30.0_dp]
    real(dp), parameter :: metric(3, 3) = reshape([ &
        1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, &
        0.0_dp, 0.0_dp, 4.0_dp], [3, 3])
    real(dp) :: flux(3), curl_weight, mass_tensor(2, 2)
    real(dp) :: harmonic, eta, theta, phi, point(3)
    integer :: status
    type(fortsparse_status_t) :: sparse_status

    call evaluate_field_aligned_flux(10.0_dp, 1.0_dp, direction, gradient, &
        flux, sparse_status)
    call check_condition(sparse_status%code == 0 .and. &
        maxval(abs(flux - expected_flux)) < 1.0e-14_dp, &
        "FEEC facade evaluates field-aligned flux")

    call scalar_reluctivity_curvilinear_fourier_coefficients( &
        metric, 1.0_dp, curl_weight, mass_tensor, status)
    call check_condition(status == 0 .and. abs(curl_weight - 2.0_dp) < &
        1.0e-14_dp, "Fourier facade evaluates curvilinear metric weights")

    harmonic = toroidal_p(1, 0, cosh(0.7_dp))
    call check_condition(harmonic == harmonic, &
        "Fourier facade exposes toroidal harmonics")

    call toroidal_point_to_cartesian(2.0_dp, 0.7_dp, 0.2_dp, 0.4_dp, point)
    call cartesian_to_toroidal(point, 2.0_dp, eta, theta, phi)
    call check_condition(maxval(abs([eta, theta, phi] - &
        [0.7_dp, 0.2_dp, 0.4_dp])) < 1.0e-12_dp, &
        "core facade preserves toroidal coordinate round trip")

    call check_summary("IGA/Fourier facade consumers")
end program test_iga_fourier_facade_consumers
