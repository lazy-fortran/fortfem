program test_surface_current_constraints
    use check, only: check_condition, check_summary
    use fortfem_boundary, only: &
        assemble_surface_current_loop_constraints, &
        assemble_surface_current_loop_constraints_jvp, &
        assemble_surface_current_loop_constraints_vjp
    use fortfem_kinds, only: dp
    use fortsparse, only: fortsparse_status_t
    implicit none

    integer, parameter :: loop_basis(2, 3) = reshape([ &
        1, 0, -1, 1, 0, -1], [2, 3])
    real(dp), parameter :: manifold_current(3, 3) = reshape([ &
        1.0_dp, 2.0_dp, -1.0_dp, -2.0_dp, 0.5_dp, 3.0_dp, &
        4.0_dp, -1.0_dp, 2.0_dp], [3, 3])
    real(dp), parameter :: prescribed_current(3, 2) = reshape([ &
        0.5_dp, 1.0_dp, -2.0_dp, 1.0_dp, -0.5_dp, 0.25_dp], [3, 2])
    real(dp), parameter :: manifold_current_dot(3, 3) = reshape([ &
        0.1_dp, -0.2_dp, 0.3_dp, 0.4_dp, 0.5_dp, -0.1_dp, &
        -0.3_dp, 0.2_dp, 0.6_dp], [3, 3])
    real(dp), parameter :: prescribed_current_dot(3, 2) = reshape([ &
        -0.2_dp, 0.1_dp, 0.4_dp, 0.3_dp, -0.5_dp, 0.2_dp], [3, 2])
    real(dp), parameter :: expected_residual(3, 2) = reshape([ &
        2.5_dp, 0.5_dp, -2.0_dp, -7.0_dp, 2.0_dp, 0.75_dp], [3, 2])
    real(dp), parameter :: residual_bar(3, 2) = reshape([ &
        0.2_dp, -0.3_dp, 0.5_dp, -0.4_dp, 0.1_dp, 0.7_dp], [3, 2])
    real(dp) :: residual(3, 2), residual_dot(3, 2)
    real(dp) :: manifold_current_bar(3, 3), prescribed_current_bar(3, 2)
    real(dp) :: lhs, rhs
    type(fortsparse_status_t) :: status

    call assemble_surface_current_loop_constraints( &
        loop_basis, manifold_current, prescribed_current, residual, status)
    call check_condition(status%code == 0 .and. &
        maxval(abs(residual - expected_residual)) < 1.0e-14_dp, &
        "loop-current residual matches the independent cycle oracle")

    call assemble_surface_current_loop_constraints_jvp( &
        loop_basis, manifold_current, prescribed_current, &
        manifold_current_dot, prescribed_current_dot, residual_dot, status)
    call check_condition(status%code == 0 .and. &
        maxval(abs(residual_dot - matmul(manifold_current_dot, &
        transpose(real(loop_basis, dp))) + &
        prescribed_current_dot)) < 1.0e-14_dp, &
        "fixed-topology loop-current JVP matches the linear oracle")

    call assemble_surface_current_loop_constraints_vjp( &
        loop_basis, manifold_current, prescribed_current, residual_bar, &
        manifold_current_bar, prescribed_current_bar, status)
    lhs = sum(residual_bar*residual_dot)
    rhs = sum(manifold_current_bar*manifold_current_dot) + &
        sum(prescribed_current_bar*prescribed_current_dot)
    call check_condition(status%code == 0 .and. abs(lhs - rhs) < 1.0e-14_dp, &
        "loop-current VJP satisfies the real dot-product identity")

    call assemble_surface_current_loop_constraints( &
        reshape([1], [1, 1]), manifold_current(:, 1:2), &
        prescribed_current(:, 1:1), residual(:, 1:1), status)
    call check_condition(status%code /= 0, &
        "loop-current constraints reject a basis with incompatible shape")

    call check_summary("surface current loop constraints")
end program test_surface_current_constraints
