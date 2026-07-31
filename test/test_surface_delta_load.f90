program test_surface_delta_load
    use check, only: check_condition, check_summary
    use fortfem_api, only: assemble_surface_delta_load, &
        assemble_surface_delta_load_jvp, assemble_surface_delta_load_vjp
    use fortfem_kinds, only: dp
    use fortsparse, only: fortsparse_status_t
    implicit none

    real(dp), parameter :: trace_basis(3, 2) = reshape([ &
        1.0_dp, 0.0_dp, 0.5_dp, 0.0_dp, 1.0_dp, 0.5_dp], [3, 2])
    real(dp), parameter :: surface_weights(3) = [2.0_dp, 1.0_dp, 2.0_dp]
    real(dp), parameter :: surface_source(3) = [3.0_dp, 4.0_dp, 5.0_dp]
    real(dp), parameter :: expected_load(2) = [11.0_dp, 9.0_dp]
    real(dp), parameter :: trace_basis_dot(3, 2) = reshape([ &
        0.1_dp, -0.2_dp, 0.3_dp, -0.4_dp, 0.5_dp, -0.6_dp], [3, 2])
    real(dp), parameter :: surface_weights_dot(3) = [0.2_dp, -0.3_dp, 0.4_dp]
    real(dp), parameter :: surface_source_dot(3) = [-0.1_dp, 0.2_dp, 0.3_dp]
    real(dp), parameter :: load_bar(2) = [0.7_dp, -0.4_dp]
    real(dp), parameter :: eps = 1.0e-6_dp
    real(dp) :: load(2), load_bad(1), load_dot(2), load_plus(2), load_minus(2)
    real(dp) :: trace_basis_bar(3, 2), surface_weights_bar(3)
    real(dp) :: surface_source_bar(3), lhs, rhs, oracle_load_dot(2)
    real(dp), parameter :: bad_weights(2) = [1.0_dp, 2.0_dp]
    type(fortsparse_status_t) :: status

    call assemble_surface_delta_load( &
        trace_basis, surface_weights, surface_source, load, status)
    call check_condition(status%code == 0, &
        "surface delta load accepts a positive trace quadrature")
    call check_condition(maxval(abs(load - expected_load)) < 1.0e-14_dp, &
        "surface delta load matches the explicit trace-transpose oracle")

    oracle_load_dot = 0.0_dp
    oracle_load_dot = matmul(transpose(trace_basis_dot), &
        surface_weights*surface_source) + &
        matmul(transpose(trace_basis), surface_weights_dot*surface_source) + &
        matmul(transpose(trace_basis), surface_weights*surface_source_dot)
    call assemble_surface_delta_load_jvp( &
        trace_basis, surface_weights, surface_source, trace_basis_dot, &
        surface_weights_dot, surface_source_dot, load_dot, status)
    call assemble_surface_delta_load( &
        trace_basis + eps*trace_basis_dot, &
        surface_weights + eps*surface_weights_dot, &
        surface_source + eps*surface_source_dot, load_plus, status)
    call assemble_surface_delta_load( &
        trace_basis - eps*trace_basis_dot, &
        surface_weights - eps*surface_weights_dot, &
        surface_source - eps*surface_source_dot, load_minus, status)
    call check_condition(status%code == 0 .and. &
        maxval(abs(load_dot - oracle_load_dot)) < 1.0e-14_dp .and. &
        maxval(abs((load_plus - load_minus)/(2.0_dp*eps) - load_dot)) < 1.0e-8_dp, &
        "surface delta load JVP matches the independent product-rule oracle")

    call assemble_surface_delta_load_vjp( &
        trace_basis, surface_weights, surface_source, load_bar, &
        trace_basis_bar, surface_weights_bar, surface_source_bar, status)
    lhs = dot_product(load_bar, load_dot)
    rhs = sum(trace_basis_bar*trace_basis_dot) + &
        dot_product(surface_weights_bar, surface_weights_dot) + &
        dot_product(surface_source_bar, surface_source_dot)
    call check_condition(status%code == 0 .and. abs(lhs - rhs) < 1.0e-12_dp, &
        "surface delta load VJP satisfies the real dot-product identity")

    call assemble_surface_delta_load( &
        trace_basis, bad_weights, surface_source, load, status)
    call check_condition(status%code /= 0, &
        "surface delta load rejects incompatible quadrature sizes")
    call assemble_surface_delta_load( &
        trace_basis, surface_weights, surface_source, load_bad, status)
    call check_condition(status%code /= 0, &
        "surface delta load rejects incompatible output size")
    call check_summary("surface delta load")
end program test_surface_delta_load
