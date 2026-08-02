program test_complex_boundary_trace_residual
    use, intrinsic :: iso_fortran_env, only: real64
    use check, only: check_condition, check_summary
    use fortfem_api, only: &
        assemble_complex_boundary_trace_residual, &
        assemble_complex_boundary_trace_residual_jvp, &
        assemble_complex_boundary_trace_residual_vjp
    use fortsparse, only: fortsparse_status_t
    implicit none

    integer, parameter :: dp = real64
    integer, parameter :: state_count = 4, normal_count = 2, tangential_count = 3
    complex(dp) :: normal_trace(normal_count, state_count)
    complex(dp) :: tangential_trace(tangential_count, state_count)
    complex(dp) :: state(state_count), normal_target(normal_count)
    complex(dp) :: tangential_target(tangential_count)
    real(dp) :: normal_weights(normal_count), tangential_weights(tangential_count)
    complex(dp) :: normal_trace_dot(normal_count, state_count)
    complex(dp) :: tangential_trace_dot(tangential_count, state_count)
    complex(dp) :: state_dot(state_count), normal_target_dot(normal_count)
    complex(dp) :: tangential_target_dot(tangential_count)
    real(dp) :: normal_weights_dot(normal_count), tangential_weights_dot(tangential_count)
    complex(dp) :: normal_residual(normal_count), tangential_residual(tangential_count)
    complex(dp) :: normal_residual_dot(normal_count)
    complex(dp) :: tangential_residual_dot(tangential_count)
    complex(dp) :: normal_residual_plus(normal_count)
    complex(dp) :: tangential_residual_plus(tangential_count)
    complex(dp) :: normal_residual_bar(normal_count)
    complex(dp) :: tangential_residual_bar(tangential_count)
    complex(dp) :: normal_trace_bar(normal_count, state_count)
    complex(dp) :: tangential_trace_bar(tangential_count, state_count)
    complex(dp) :: state_bar(state_count)
    complex(dp) :: normal_target_bar(normal_count)
    complex(dp) :: tangential_target_bar(tangential_count)
    real(dp) :: normal_weights_bar(normal_count)
    real(dp) :: tangential_weights_bar(tangential_count)
    real(dp) :: epsilon, finite_difference_error, lhs, rhs
    type(fortsparse_status_t) :: status
    logical :: all_passed

    all_passed = .true.
    call build_system( &
        normal_trace, tangential_trace, normal_weights, tangential_weights, state, &
        normal_target, tangential_target)
    call assemble_complex_boundary_trace_residual( &
        normal_trace, tangential_trace, normal_weights, tangential_weights, state, &
        normal_target, tangential_target, normal_residual, tangential_residual, status)
    call record_condition(status%code == 0, &
        "complex boundary trace residual assembles")
    call record_condition(maxval(abs(normal_residual - normal_weights*( &
        matmul(normal_trace, state) - normal_target))) < 1.0e-14_dp, &
        "normal trace residual matches the independent oracle")
    call record_condition(maxval(abs(tangential_residual - tangential_weights*( &
        matmul(tangential_trace, state) - tangential_target))) < 1.0e-14_dp, &
        "tangential trace residual matches the independent oracle")

    normal_trace_dot = reshape([ &
        cmplx(0.03_dp, -0.02_dp, dp), cmplx(-0.02_dp, 0.01_dp, dp), &
        cmplx(0.01_dp, 0.04_dp, dp), cmplx(0.04_dp, -0.01_dp, dp), &
        cmplx(-0.01_dp, 0.02_dp, dp), cmplx(0.02_dp, -0.03_dp, dp), &
        cmplx(0.05_dp, 0.02_dp, dp), cmplx(-0.03_dp, 0.01_dp, dp)], &
        shape(normal_trace_dot))
    tangential_trace_dot = reshape([ &
        cmplx(-0.02_dp, 0.01_dp, dp), cmplx(0.01_dp, 0.03_dp, dp), &
        cmplx(0.03_dp, -0.02_dp, dp), cmplx(-0.01_dp, 0.04_dp, dp), &
        cmplx(0.04_dp, -0.03_dp, dp), cmplx(-0.03_dp, 0.02_dp, dp), &
        cmplx(0.02_dp, 0.05_dp, dp), cmplx(0.05_dp, -0.01_dp, dp), &
        cmplx(-0.04_dp, 0.02_dp, dp), cmplx(0.03_dp, -0.01_dp, dp), &
        cmplx(0.01_dp, 0.02_dp, dp), cmplx(-0.02_dp, 0.04_dp, dp)], &
        shape(tangential_trace_dot))
    state_dot = [ &
        cmplx(0.02_dp, -0.03_dp, dp), cmplx(-0.03_dp, 0.01_dp, dp), &
        cmplx(0.01_dp, 0.04_dp, dp), cmplx(0.04_dp, -0.02_dp, dp)]
    normal_target_dot = [ &
        cmplx(-0.01_dp, 0.02_dp, dp), cmplx(0.03_dp, -0.01_dp, dp)]
    tangential_target_dot = [ &
        cmplx(-0.02_dp, 0.04_dp, dp), cmplx(0.01_dp, -0.03_dp, dp), &
        cmplx(0.03_dp, 0.01_dp, dp)]
    normal_weights_dot = [0.04_dp, -0.03_dp]
    tangential_weights_dot = [0.02_dp, -0.01_dp, 0.03_dp]
    call assemble_complex_boundary_trace_residual_jvp( &
        normal_trace, tangential_trace, normal_weights, tangential_weights, state, &
        normal_target, tangential_target, normal_trace_dot, tangential_trace_dot, &
        normal_weights_dot, tangential_weights_dot, state_dot, normal_target_dot, &
        tangential_target_dot, normal_residual_dot, tangential_residual_dot, status)
    call record_condition(status%code == 0, &
        "complex boundary trace residual JVP assembles")

    epsilon = 1.0e-7_dp
    call assemble_complex_boundary_trace_residual( &
        normal_trace + epsilon*normal_trace_dot, &
        tangential_trace + epsilon*tangential_trace_dot, &
        normal_weights + epsilon*normal_weights_dot, &
        tangential_weights + epsilon*tangential_weights_dot, state + epsilon*state_dot, &
        normal_target + epsilon*normal_target_dot, &
        tangential_target + epsilon*tangential_target_dot, normal_residual_plus, &
        tangential_residual_plus, status)
    finite_difference_error = max( &
        maxval(abs(normal_residual_dot - &
        (normal_residual_plus - normal_residual)/epsilon)), &
        maxval(abs(tangential_residual_dot - &
        (tangential_residual_plus - tangential_residual)/epsilon)))
    call record_condition(finite_difference_error < 2.0e-8_dp, &
        "complex boundary trace JVP matches a forward difference")

    normal_residual_bar = [ &
        cmplx(0.2_dp, -0.3_dp, dp), cmplx(-0.3_dp, 0.4_dp, dp)]
    tangential_residual_bar = [ &
        cmplx(0.4_dp, 0.1_dp, dp), cmplx(-0.1_dp, 0.5_dp, dp), &
        cmplx(0.5_dp, -0.2_dp, dp)]
    call assemble_complex_boundary_trace_residual_vjp( &
        normal_trace, tangential_trace, normal_weights, tangential_weights, state, &
        normal_target, tangential_target, normal_residual_bar, &
        tangential_residual_bar, normal_trace_bar, tangential_trace_bar, &
        normal_weights_bar, tangential_weights_bar, state_bar, normal_target_bar, &
        tangential_target_bar, status)
    call record_condition(status%code == 0, &
        "complex boundary trace residual VJP assembles")
    lhs = real(sum(conjg(normal_residual_bar)*normal_residual_dot) + &
        sum(conjg(tangential_residual_bar)*tangential_residual_dot), dp)
    rhs = real(sum(conjg(normal_trace_bar)*normal_trace_dot) + &
        sum(conjg(tangential_trace_bar)*tangential_trace_dot) + &
        dot_product(normal_weights_bar, normal_weights_dot) + &
        dot_product(tangential_weights_bar, tangential_weights_dot) + &
        sum(conjg(state_bar)*state_dot) + &
        sum(conjg(normal_target_bar)*normal_target_dot) + &
        sum(conjg(tangential_target_bar)*tangential_target_dot), dp)
    call record_condition(abs(lhs - rhs) < 2.0e-12_dp, &
        "complex boundary trace VJP satisfies the real adjoint identity")

    call check_summary("complex boundary trace residual")
    if (.not. all_passed) error stop 1

contains

    subroutine build_system( &
            normal_trace, tangential_trace, normal_weights, tangential_weights, &
            state, normal_target, tangential_target)
        complex(dp), intent(out) :: normal_trace(:, :), tangential_trace(:, :)
        real(dp), intent(out) :: normal_weights(:), tangential_weights(:)
        complex(dp), intent(out) :: state(:), normal_target(:), tangential_target(:)

        normal_trace = reshape([ &
            cmplx(1.0_dp, 0.2_dp, dp), cmplx(-0.2_dp, 0.1_dp, dp), &
            cmplx(0.3_dp, -0.3_dp, dp), cmplx(0.4_dp, 0.2_dp, dp), &
            cmplx(-0.5_dp, 0.4_dp, dp), cmplx(0.6_dp, -0.1_dp, dp), &
            cmplx(0.7_dp, 0.3_dp, dp), cmplx(0.1_dp, -0.2_dp, dp)], &
            shape(normal_trace))
        tangential_trace = reshape([ &
            cmplx(0.4_dp, -0.1_dp, dp), cmplx(-0.1_dp, 0.2_dp, dp), &
            cmplx(0.2_dp, 0.3_dp, dp), cmplx(0.8_dp, -0.2_dp, dp), &
            cmplx(-0.6_dp, 0.1_dp, dp), cmplx(0.3_dp, 0.4_dp, dp), &
            cmplx(0.5_dp, -0.3_dp, dp), cmplx(-0.7_dp, 0.2_dp, dp), &
            cmplx(0.2_dp, 0.1_dp, dp), cmplx(0.6_dp, -0.4_dp, dp), &
            cmplx(-0.3_dp, 0.2_dp, dp), cmplx(0.1_dp, 0.5_dp, dp)], &
            shape(tangential_trace))
        normal_weights = [1.5_dp, 2.0_dp]
        tangential_weights = [0.8_dp, 1.3_dp, 1.7_dp]
        state = [ &
            cmplx(0.6_dp, -0.4_dp, dp), cmplx(-0.4_dp, 0.2_dp, dp), &
            cmplx(0.2_dp, 0.5_dp, dp), cmplx(0.9_dp, -0.1_dp, dp)]
        normal_target = [ &
            cmplx(0.1_dp, 0.2_dp, dp), cmplx(-0.3_dp, 0.4_dp, dp)]
        tangential_target = [ &
            cmplx(0.5_dp, -0.2_dp, dp), cmplx(-0.1_dp, 0.3_dp, dp), &
            cmplx(0.2_dp, -0.4_dp, dp)]
    end subroutine build_system

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_complex_boundary_trace_residual
