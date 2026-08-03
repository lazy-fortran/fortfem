program test_surface_edge_balance
    use check, only: check_condition, check_summary
    use fortfem_feec, only: assemble_surface_edge_flux_balance, &
        assemble_surface_edge_flux_balance_jvp, &
        assemble_surface_edge_flux_balance_vjp
    use fortfem_kinds, only: dp
    use fortsparse, only: fortsparse_status_t
    implicit none

    integer, parameter :: edge_boundary(3, 2) = reshape([ &
        -1, 1, 0, 0, -1, 1], [3, 2])
    integer, parameter :: closed_boundary(3, 3) = reshape([ &
        -1, 1, 0, 0, -1, 1, 1, 0, -1], [3, 3])
    integer, parameter :: invalid_boundary(2, 1) = reshape([1, 0], [2, 1])
    real(dp), parameter :: edge_flux(2) = [2.0_dp, -3.0_dp]
    real(dp), parameter :: edge_flux_dot(2) = [0.25_dp, -0.4_dp]
    real(dp), parameter :: expected_balance(3) = [-2.0_dp, 5.0_dp, -3.0_dp]
    real(dp), parameter :: balance_bar(3) = [0.2_dp, -0.4_dp, 0.7_dp]
    real(dp) :: balance(3), global_balance
    real(dp) :: balance_dot(3), global_balance_dot
    real(dp) :: edge_flux_bar(2), lhs, rhs
    real(dp) :: closed_flux(3) = [1.0_dp, -2.0_dp, 0.5_dp]
    real(dp) :: closed_balance(3), closed_global
    type(fortsparse_status_t) :: status

    call assemble_surface_edge_flux_balance( &
        edge_boundary, edge_flux, balance, global_balance, status)
    call check_condition(status%code == 0 .and. &
        maxval(abs(balance - expected_balance)) < 1.0e-14_dp .and. &
        abs(global_balance) < 1.0e-14_dp, &
        "open surface edge balance matches the independent incidence oracle")

    call assemble_surface_edge_flux_balance_jvp( &
        edge_boundary, edge_flux, edge_flux_dot, balance_dot, &
        global_balance_dot, status)
    call check_condition(status%code == 0 .and. &
        abs(global_balance_dot) < 1.0e-14_dp .and. &
        maxval(abs(balance_dot - [-0.25_dp, 0.65_dp, -0.4_dp])) < &
        1.0e-14_dp, "surface edge-balance JVP preserves global conservation")

    call assemble_surface_edge_flux_balance_vjp( &
        edge_boundary, edge_flux, balance_bar, 0.3_dp, edge_flux_bar, status)
    lhs = dot_product(balance_bar, balance) + 0.3_dp*global_balance
    rhs = dot_product(edge_flux_bar, edge_flux)
    call check_condition(status%code == 0 .and. abs(lhs - rhs) < 1.0e-14_dp, &
        "surface edge-balance VJP satisfies the real dot-product identity")

    call assemble_surface_edge_flux_balance( &
        closed_boundary, closed_flux, closed_balance, closed_global, status)
    call check_condition(status%code == 0 .and. abs(closed_global) < 1.0e-14_dp, &
        "closed surface edge cycle has zero global divergence")

    call assemble_surface_edge_flux_balance( &
        invalid_boundary, edge_flux(1:1), balance(1:2), global_balance, status)
    call check_condition(status%code /= 0, &
        "surface edge balance rejects a non-conservative incidence column")

    call check_summary("surface edge flux balance")
end program test_surface_edge_balance
