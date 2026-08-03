program test_coupled_field_residual
    use, intrinsic :: iso_fortran_env, only: real64
    use check, only: check_condition, check_summary
    use fortfem_feec, only: &
        assemble_coupled_field_residual, &
        assemble_coupled_field_residual_jvp, &
        assemble_coupled_field_residual_vjp
    use fortsparse, only: fortsparse_status_t
    implicit none

    integer, parameter :: dp = real64
    integer, parameter :: field_count = 3, state_count = 4
    integer, parameter :: constraint_count = 2
    real(dp) :: field_operator(field_count, state_count)
    real(dp) :: constraint_operator(constraint_count, state_count)
    real(dp) :: state(state_count), source(field_count)
    real(dp) :: constraint_target(constraint_count)
    real(dp) :: field_operator_dot(field_count, state_count)
    real(dp) :: constraint_operator_dot(constraint_count, state_count)
    real(dp) :: state_dot(state_count), source_dot(field_count)
    real(dp) :: constraint_target_dot(constraint_count)
    real(dp) :: field_residual(field_count), constraint_residual(constraint_count)
    real(dp) :: field_residual_dot(field_count)
    real(dp) :: constraint_residual_dot(constraint_count)
    real(dp) :: field_residual_plus(field_count)
    real(dp) :: constraint_residual_plus(constraint_count)
    real(dp) :: field_residual_bar(field_count)
    real(dp) :: constraint_residual_bar(constraint_count)
    real(dp) :: field_operator_bar(field_count, state_count)
    real(dp) :: constraint_operator_bar(constraint_count, state_count)
    real(dp) :: state_bar(state_count), source_bar(field_count)
    real(dp) :: constraint_target_bar(constraint_count)
    real(dp) :: epsilon, finite_difference_error, lhs, rhs
    type(fortsparse_status_t) :: status
    logical :: all_passed

    all_passed = .true.
    call build_system( &
        field_operator, constraint_operator, state, source, constraint_target)
    call assemble_coupled_field_residual( &
        field_operator, state, source, constraint_operator, constraint_target, &
        field_residual, constraint_residual, status)
    call record_condition(status%code == 0, "coupled field residual assembles")
    call record_condition(maxval(abs(field_residual - &
        (matmul(field_operator, state) - source))) < 1.0e-14_dp, &
        "field residual matches its matrix oracle")
    call record_condition(maxval(abs(constraint_residual - &
        (matmul(constraint_operator, state) - constraint_target))) < 1.0e-14_dp, &
        "constraint residual matches its matrix oracle")

    field_operator_dot = reshape([ &
        0.03_dp, -0.02_dp, 0.01_dp, 0.04_dp, -0.01_dp, 0.02_dp, &
        0.05_dp, -0.03_dp, 0.02_dp, 0.01_dp, -0.04_dp, 0.03_dp], &
        shape(field_operator_dot))
    constraint_operator_dot = reshape([ &
        -0.02_dp, 0.01_dp, 0.03_dp, -0.01_dp, 0.04_dp, -0.03_dp, &
        0.02_dp, 0.05_dp], shape(constraint_operator_dot))
    state_dot = [0.02_dp, -0.03_dp, 0.01_dp, 0.04_dp]
    source_dot = [-0.01_dp, 0.03_dp, -0.02_dp]
    constraint_target_dot = [0.05_dp, -0.04_dp]
    call assemble_coupled_field_residual_jvp( &
        field_operator, state, source, constraint_operator, constraint_target, &
        field_operator_dot, state_dot, source_dot, constraint_operator_dot, &
        constraint_target_dot, field_residual_dot, constraint_residual_dot, status)
    call record_condition(status%code == 0, "coupled field residual JVP assembles")

    epsilon = 1.0e-7_dp
    call assemble_coupled_field_residual( &
        field_operator + epsilon*field_operator_dot, state + epsilon*state_dot, &
        source + epsilon*source_dot, &
        constraint_operator + epsilon*constraint_operator_dot, &
        constraint_target + epsilon*constraint_target_dot, &
        field_residual_plus, constraint_residual_plus, status)
    finite_difference_error = max( &
        maxval(abs(field_residual_dot - &
        (field_residual_plus - field_residual)/epsilon)), &
        maxval(abs(constraint_residual_dot - &
        (constraint_residual_plus - constraint_residual)/epsilon)))
    call record_condition(finite_difference_error < 2.0e-8_dp, &
        "coupled field residual JVP matches a forward difference")

    field_residual_bar = [0.2_dp, -0.3_dp, 0.4_dp]
    constraint_residual_bar = [-0.1_dp, 0.5_dp]
    call assemble_coupled_field_residual_vjp( &
        field_operator, state, source, constraint_operator, constraint_target, &
        field_residual_bar, constraint_residual_bar, field_operator_bar, state_bar, &
        source_bar, constraint_operator_bar, constraint_target_bar, status)
    call record_condition(status%code == 0, "coupled field residual VJP assembles")
    lhs = dot_product(field_residual_bar, field_residual_dot) + &
        dot_product(constraint_residual_bar, constraint_residual_dot)
    rhs = sum(field_operator_bar*field_operator_dot) + &
        sum(constraint_operator_bar*constraint_operator_dot) + &
        dot_product(state_bar, state_dot) + dot_product(source_bar, source_dot) + &
        dot_product(constraint_target_bar, constraint_target_dot)
    call record_condition(abs(lhs - rhs) < 2.0e-13_dp, &
        "coupled field residual VJP satisfies the adjoint identity")

    call assemble_coupled_field_residual( &
        field_operator(:, :field_count - 1), state, source, constraint_operator, &
        constraint_target, field_residual, constraint_residual, status)
    call record_condition(status%code /= 0, &
        "incompatible field dimensions are rejected")

    call check_summary("coupled field residual")
    if (.not. all_passed) error stop 1

contains

    subroutine build_system( &
            field_operator, constraint_operator, state, source, constraint_target)
        real(dp), intent(out) :: field_operator(:, :), constraint_operator(:, :)
        real(dp), intent(out) :: state(:), source(:), constraint_target(:)

        field_operator = reshape([ &
            1.0_dp, -0.2_dp, 0.3_dp, 0.4_dp, -0.5_dp, 0.6_dp, &
            0.7_dp, 0.1_dp, -0.8_dp, 0.2_dp, 0.9_dp, -0.3_dp], &
            shape(field_operator))
        constraint_operator = reshape([ &
            0.4_dp, -0.1_dp, 0.2_dp, 0.8_dp, -0.6_dp, 0.3_dp, &
            0.5_dp, -0.7_dp], shape(constraint_operator))
        state = [0.6_dp, -0.4_dp, 0.2_dp, 0.9_dp]
        source = [0.1_dp, -0.3_dp, 0.5_dp]
        constraint_target = [-0.2_dp, 0.4_dp]
    end subroutine build_system

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_coupled_field_residual
