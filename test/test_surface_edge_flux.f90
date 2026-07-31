program test_surface_edge_flux
    use check, only: check_condition, check_summary
    use fortfem_api, only: assemble_surface_edge_flux, &
        assemble_surface_edge_flux_jvp, assemble_surface_edge_flux_vjp
    use fortfem_kinds, only: dp
    use fortsparse, only: fortsparse_status_t
    implicit none

    real(dp), parameter :: edge_conormal(2, 2, 3) = reshape([ &
        1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 2.0_dp, 0.0_dp, &
        0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp], [2, 2, 3])
    real(dp), parameter :: edge_weights(2, 2) = reshape([ &
        2.0_dp, 1.0_dp, 3.0_dp, 4.0_dp], [2, 2])
    real(dp), parameter :: surface_current(2, 2, 3) = reshape([ &
        3.0_dp, 4.0_dp, 5.0_dp, 1.0_dp, -2.0_dp, 2.0_dp, &
        2.0_dp, -1.0_dp, 3.0_dp, -2.0_dp, 1.0_dp, 4.0_dp], [2, 2, 3])
    real(dp), parameter :: expected_flux(2) = [-2.0_dp, 12.0_dp]
    real(dp), parameter :: edge_conormal_dot(2, 2, 3) = reshape([ &
        0.2_dp, -0.1_dp, 0.3_dp, -0.4_dp, 0.5_dp, -0.6_dp, &
        0.7_dp, -0.8_dp, 0.9_dp, -1.0_dp, 1.1_dp, -1.2_dp], [2, 2, 3])
    real(dp), parameter :: edge_weights_dot(2, 2) = reshape([ &
        0.1_dp, -0.2_dp, 0.3_dp, -0.4_dp], [2, 2])
    real(dp), parameter :: surface_current_dot(2, 2, 3) = reshape([ &
        -0.1_dp, 0.2_dp, -0.3_dp, 0.4_dp, -0.5_dp, 0.6_dp, &
        -0.7_dp, 0.8_dp, -0.9_dp, 1.0_dp, -1.1_dp, 1.2_dp], [2, 2, 3])
    real(dp), parameter :: edge_flux_bar(2) = [0.8_dp, -0.7_dp]
    real(dp), parameter :: eps = 1.0e-6_dp
    real(dp) :: edge_flux(2), edge_flux_dot(2), edge_flux_plus(2)
    real(dp) :: edge_flux_minus(2), expected_dot(2)
    real(dp) :: edge_conormal_bar(2, 2, 3), edge_weights_bar(2, 2)
    real(dp) :: surface_current_bar(2, 2, 3), lhs, rhs
    real(dp) :: reversed_flux(2)
    real(dp), parameter :: bad_weights(2, 2) = reshape([ &
        2.0_dp, -1.0_dp, 3.0_dp, 4.0_dp], [2, 2])
    type(fortsparse_status_t) :: status
    integer :: quadrature, edge

    call assemble_surface_edge_flux( &
        edge_conormal, edge_weights, surface_current, edge_flux, status)
    call check_condition(status%code == 0, &
        "surface edge flux accepts positive oriented edge quadrature")
    call check_condition(maxval(abs(edge_flux - expected_flux)) < 1.0e-14_dp, &
        "surface edge flux matches the independent contraction oracle")

    call assemble_surface_edge_flux( &
        -edge_conormal, edge_weights, surface_current, reversed_flux, status)
    call check_condition(status%code == 0 .and. &
        maxval(abs(reversed_flux + edge_flux)) < 1.0e-14_dp, &
        "reversing each edge conormal reverses its oriented flux")

    expected_dot = 0.0_dp
    do quadrature = 1, 2
        do edge = 1, 2
            expected_dot(edge) = expected_dot(edge) + &
                edge_weights_dot(quadrature, edge)*dot_product( &
                edge_conormal(quadrature, edge, :), &
                surface_current(quadrature, edge, :)) + &
                edge_weights(quadrature, edge)*dot_product( &
                edge_conormal_dot(quadrature, edge, :), &
                surface_current(quadrature, edge, :)) + &
                edge_weights(quadrature, edge)*dot_product( &
                edge_conormal(quadrature, edge, :), &
                surface_current_dot(quadrature, edge, :))
        end do
    end do
    call assemble_surface_edge_flux_jvp( &
        edge_conormal, edge_weights, surface_current, edge_conormal_dot, &
        edge_weights_dot, surface_current_dot, edge_flux_dot, status)
    call assemble_surface_edge_flux( &
        edge_conormal + eps*edge_conormal_dot, &
        edge_weights + eps*edge_weights_dot, &
        surface_current + eps*surface_current_dot, edge_flux_plus, status)
    call assemble_surface_edge_flux( &
        edge_conormal - eps*edge_conormal_dot, &
        edge_weights - eps*edge_weights_dot, &
        surface_current - eps*surface_current_dot, edge_flux_minus, status)
    call check_condition(status%code == 0 .and. &
        maxval(abs(edge_flux_dot - expected_dot)) < 1.0e-14_dp .and. &
        maxval(abs((edge_flux_plus - edge_flux_minus)/(2.0_dp*eps) - &
        edge_flux_dot)) < 1.0e-8_dp, &
        "surface edge flux JVP matches product rule and finite differences")

    call assemble_surface_edge_flux_vjp( &
        edge_conormal, edge_weights, surface_current, edge_flux_bar, &
        edge_conormal_bar, edge_weights_bar, surface_current_bar, status)
    lhs = dot_product(edge_flux_bar, edge_flux_dot)
    rhs = sum(edge_conormal_bar*edge_conormal_dot) + &
        sum(edge_weights_bar*edge_weights_dot) + &
        sum(surface_current_bar*surface_current_dot)
    call check_condition(status%code == 0 .and. abs(lhs - rhs) < 1.0e-12_dp, &
        "surface edge flux VJP satisfies the real dot-product identity")

    call assemble_surface_edge_flux( &
        edge_conormal, bad_weights, surface_current, edge_flux, status)
    call check_condition(status%code /= 0, &
        "surface edge flux rejects non-positive edge quadrature")

    call check_summary("surface edge flux")
end program test_surface_edge_flux
