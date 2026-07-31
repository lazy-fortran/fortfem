program test_surface_current_balance
    use check, only: check_condition, check_summary
    use fortfem_api, only: assemble_surface_current_junction_balance, &
        assemble_surface_current_junction_balance_jvp, &
        assemble_surface_current_junction_balance_vjp
    use fortfem_kinds, only: dp
    use fortsparse, only: fortsparse_status_t
    implicit none

    integer, parameter :: junction_incidence(3, 2) = reshape([ &
        -1, 1, 0, 0, -1, 1], [3, 2])
    integer, parameter :: invalid_incidence(2, 1) = reshape([1, 0], [2, 1])
    real(dp), parameter :: manifold_current(3, 2) = reshape([ &
        1.0_dp, 2.0_dp, 3.0_dp, -2.0_dp, 1.0_dp, 4.0_dp], [3, 2])
    real(dp), parameter :: manifold_current_dot(3, 2) = reshape([ &
        0.1_dp, -0.2_dp, 0.3_dp, 0.4_dp, 0.5_dp, -0.1_dp], [3, 2])
    real(dp), parameter :: expected_balance(3, 3) = reshape([ &
        -1.0_dp, -2.0_dp, -3.0_dp, 3.0_dp, 1.0_dp, -1.0_dp, &
        -2.0_dp, 1.0_dp, 4.0_dp], [3, 3])
    real(dp), parameter :: junction_balance_bar(3, 3) = reshape([ &
        0.2_dp, -0.3_dp, 0.5_dp, -0.4_dp, 0.1_dp, 0.7_dp, &
        0.6_dp, -0.2_dp, -0.8_dp], [3, 3])
    real(dp), parameter :: global_balance_bar(3) = [0.3_dp, -0.6_dp, 0.4_dp]
    real(dp) :: junction_balance(3, 3), global_balance(3)
    real(dp) :: junction_balance_dot(3, 3), global_balance_dot(3)
    real(dp) :: manifold_current_bar(3, 2), lhs, rhs
    real(dp) :: closed_current(3, 1), closed_balance(3, 0), closed_global(3)
    integer, allocatable :: closed_incidence(:, :)
    type(fortsparse_status_t) :: status

    call assemble_surface_current_junction_balance( &
        junction_incidence, manifold_current, junction_balance, &
        global_balance, status)
    call check_condition(status%code == 0 .and. &
        maxval(abs(junction_balance - expected_balance)) < 1.0e-14_dp .and. &
        maxval(abs(global_balance)) < 1.0e-14_dp, &
        "open manifold current ledger matches the independent balance oracle")

    call assemble_surface_current_junction_balance_jvp( &
        junction_incidence, manifold_current, manifold_current_dot, &
        junction_balance_dot, global_balance_dot, status)
    call check_condition(status%code == 0 .and. &
        maxval(abs(global_balance_dot)) < 1.0e-14_dp, &
        "fixed-topology current ledger JVP preserves global conservation")

    call assemble_surface_current_junction_balance_vjp( &
        junction_incidence, manifold_current, junction_balance_bar, &
        global_balance_bar, manifold_current_bar, status)
    lhs = sum(junction_balance_bar*junction_balance_dot) + &
        dot_product(global_balance_bar, global_balance_dot)
    rhs = sum(manifold_current_bar*manifold_current_dot)
    call check_condition(status%code == 0 .and. abs(lhs - rhs) < 1.0e-14_dp, &
        "current ledger VJP satisfies the real dot-product identity")

    closed_current(:, 1) = [2.0_dp, -1.0_dp, 0.5_dp]
    allocate(closed_incidence(0, 1))
    call assemble_surface_current_junction_balance( &
        closed_incidence, closed_current, closed_balance, closed_global, status)
    call check_condition(status%code == 0 .and. &
        maxval(abs(closed_global)) < 1.0e-14_dp, &
        "boundaryless manifold has zero junction current balance")

    call assemble_surface_current_junction_balance( &
        invalid_incidence, closed_current, junction_balance(:, 1:2), &
        global_balance, status)
    call check_condition(status%code /= 0, &
        "current ledger rejects non-conservative incidence columns")

    call check_summary("surface current junction balance")
end program test_surface_current_balance
