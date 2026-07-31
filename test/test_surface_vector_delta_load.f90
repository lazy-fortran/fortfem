program test_surface_vector_delta_load
    use check, only: check_condition, check_summary
    use fortfem_api, only: assemble_surface_vector_delta_load, &
        assemble_surface_vector_delta_load_jvp, &
        assemble_surface_vector_delta_load_vjp
    use fortfem_kinds, only: dp
    use fortsparse, only: fortsparse_status_t
    implicit none

    real(dp), parameter :: tangential_trace_basis(2, 2, 3) = reshape([ &
        1.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, &
        0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp, &
        0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp], [2, 2, 3])
    real(dp), parameter :: surface_weights(2) = [2.0_dp, 1.0_dp]
    real(dp), parameter :: surface_current(2, 3) = reshape([ &
        3.0_dp, 4.0_dp, 5.0_dp, 1.0_dp, -2.0_dp, 2.0_dp], [2, 3])
    real(dp), parameter :: expected_load(2) = [7.0_dp, 1.0_dp]
    real(dp), parameter :: tangential_trace_basis_dot(2, 2, 3) = reshape([ &
        0.2_dp, -0.1_dp, 0.3_dp, -0.4_dp, 0.5_dp, -0.6_dp, &
        0.7_dp, -0.8_dp, 0.9_dp, -1.0_dp, 1.1_dp, -1.2_dp], [2, 2, 3])
    real(dp), parameter :: surface_weights_dot(2) = [0.3_dp, -0.2_dp]
    real(dp), parameter :: surface_current_dot(2, 3) = reshape([ &
        -0.1_dp, 0.2_dp, -0.3_dp, 0.4_dp, -0.5_dp, 0.6_dp], [2, 3])
    real(dp), parameter :: load_bar(2) = [0.8_dp, -0.7_dp]
    real(dp), parameter :: eps = 1.0e-6_dp
    real(dp) :: load(2), load_bad(1), load_dot(2), load_plus(2), load_minus(2)
    real(dp) :: trace_basis_bar(2, 2, 3), surface_weights_bar(2)
    real(dp) :: surface_current_bar(2, 3), lhs, rhs, oracle_load_dot(2)
    real(dp), parameter :: bad_weights(1) = [1.0_dp]
    type(fortsparse_status_t) :: status
    integer :: quadrature, dof

    call assemble_surface_vector_delta_load( &
        tangential_trace_basis, surface_weights, surface_current, load, status)
    call check_condition(status%code == 0, &
        "surface vector delta load accepts tangential trace data")
    call check_condition(maxval(abs(load - expected_load)) < 1.0e-14_dp, &
        "surface vector delta load matches the current-sheet pairing oracle")

    oracle_load_dot = 0.0_dp
    do quadrature = 1, 2
        do dof = 1, 2
            oracle_load_dot(dof) = oracle_load_dot(dof) + &
                surface_weights_dot(quadrature)*dot_product( &
                    tangential_trace_basis(quadrature, dof, :), &
                    surface_current(quadrature, :)) + &
                surface_weights(quadrature)*dot_product( &
                    tangential_trace_basis_dot(quadrature, dof, :), &
                    surface_current(quadrature, :)) + &
                surface_weights(quadrature)*dot_product( &
                    tangential_trace_basis(quadrature, dof, :), &
                    surface_current_dot(quadrature, :))
        end do
    end do
    call assemble_surface_vector_delta_load_jvp( &
        tangential_trace_basis, surface_weights, surface_current, &
        tangential_trace_basis_dot, surface_weights_dot, surface_current_dot, &
        load_dot, status)
    call assemble_surface_vector_delta_load( &
        tangential_trace_basis + eps*tangential_trace_basis_dot, &
        surface_weights + eps*surface_weights_dot, &
        surface_current + eps*surface_current_dot, load_plus, status)
    call assemble_surface_vector_delta_load( &
        tangential_trace_basis - eps*tangential_trace_basis_dot, &
        surface_weights - eps*surface_weights_dot, &
        surface_current - eps*surface_current_dot, load_minus, status)
    call check_condition(status%code == 0 .and. &
        maxval(abs(load_dot - oracle_load_dot)) < 1.0e-14_dp .and. &
        maxval(abs((load_plus - load_minus)/(2.0_dp*eps) - load_dot)) < 1.0e-8_dp, &
        "surface vector delta load JVP matches the independent product-rule oracle")

    call assemble_surface_vector_delta_load_vjp( &
        tangential_trace_basis, surface_weights, surface_current, load_bar, &
        trace_basis_bar, surface_weights_bar, surface_current_bar, status)
    lhs = dot_product(load_bar, load_dot)
    rhs = sum(trace_basis_bar*tangential_trace_basis_dot) + &
        dot_product(surface_weights_bar, surface_weights_dot) + &
        sum(surface_current_bar*surface_current_dot)
    call check_condition(status%code == 0 .and. abs(lhs - rhs) < 1.0e-12_dp, &
        "surface vector delta load VJP satisfies the real dot-product identity")

    call assemble_surface_vector_delta_load( &
        tangential_trace_basis, bad_weights, surface_current, load, status)
    call check_condition(status%code /= 0, &
        "surface vector delta load rejects incompatible quadrature sizes")
    call assemble_surface_vector_delta_load( &
        tangential_trace_basis, surface_weights, surface_current, load_bad, status)
    call check_condition(status%code /= 0, &
        "surface vector delta load rejects incompatible output size")
    call check_summary("surface vector delta load")
end program test_surface_vector_delta_load
