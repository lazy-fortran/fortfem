program test_block_2x2_residual
    use, intrinsic :: iso_fortran_env, only: real64
    use check, only: check_condition, check_summary
    use fortfem_api, only: &
        assemble_block_2x2_residual, &
        assemble_block_2x2_residual_jvp, &
        assemble_block_2x2_residual_vjp
    use fortsparse, only: fortsparse_status_t
    implicit none

    integer, parameter :: dp = real64
    integer, parameter :: first_rows = 4, second_rows = 3
    integer, parameter :: first_state = 3, second_state = 2
    real(dp) :: aa(first_rows, first_state), ab(first_rows, second_state)
    real(dp) :: ba(second_rows, first_state), bb(second_rows, second_state)
    real(dp) :: x(first_state), y(second_state)
    real(dp) :: rhs_first(first_rows), rhs_second(second_rows)
    real(dp) :: aa_dot(first_rows, first_state), ab_dot(first_rows, second_state)
    real(dp) :: ba_dot(second_rows, first_state), bb_dot(second_rows, second_state)
    real(dp) :: x_dot(first_state), y_dot(second_state)
    real(dp) :: rhs_first_dot(first_rows), rhs_second_dot(second_rows)
    real(dp) :: residual_first(first_rows), residual_second(second_rows)
    real(dp) :: residual_first_dot(first_rows), residual_second_dot(second_rows)
    real(dp) :: residual_first_plus(first_rows), residual_second_plus(second_rows)
    real(dp) :: residual_first_bar(first_rows), residual_second_bar(second_rows)
    real(dp) :: aa_bar(first_rows, first_state), ab_bar(first_rows, second_state)
    real(dp) :: ba_bar(second_rows, first_state), bb_bar(second_rows, second_state)
    real(dp) :: x_bar(first_state), y_bar(second_state)
    real(dp) :: rhs_first_bar(first_rows), rhs_second_bar(second_rows)
    real(dp) :: epsilon, finite_difference_error, lhs, rhs
    type(fortsparse_status_t) :: status
    logical :: all_passed

    all_passed = .true.
    call build_system(aa, ab, ba, bb, x, y, rhs_first, rhs_second)
    call assemble_block_2x2_residual(aa, ab, ba, bb, x, y, rhs_first, rhs_second, &
        residual_first, residual_second, status)
    call record_condition(status%code == 0, "2x2 block residual assembles")
    call record_condition(maxval(abs(residual_first - &
        (matmul(aa, x) + matmul(ab, y) - rhs_first))) < 1.0e-14_dp, &
        "first block residual matches the matrix oracle")
    call record_condition(maxval(abs(residual_second - &
        (matmul(ba, x) + matmul(bb, y) - rhs_second))) < 1.0e-14_dp, &
        "second block residual matches the matrix oracle")

    aa_dot = reshape([ &
        0.03_dp, -0.02_dp, 0.01_dp, 0.04_dp, -0.01_dp, 0.02_dp, &
        0.05_dp, -0.03_dp, 0.02_dp, 0.01_dp, -0.04_dp, 0.03_dp], &
        shape(aa_dot))
    ab_dot = reshape([0.02_dp, -0.01_dp, 0.04_dp, -0.03_dp, 0.01_dp, 0.05_dp, &
        -0.02_dp, 0.03_dp], shape(ab_dot))
    ba_dot = reshape([-0.01_dp, 0.03_dp, 0.02_dp, 0.04_dp, -0.02_dp, 0.01_dp, &
        0.05_dp, -0.04_dp, 0.03_dp], shape(ba_dot))
    bb_dot = reshape([0.04_dp, -0.03_dp, 0.01_dp, 0.02_dp, -0.05_dp, 0.03_dp], &
        shape(bb_dot))
    x_dot = [0.02_dp, -0.03_dp, 0.01_dp]
    y_dot = [0.04_dp, -0.02_dp]
    rhs_first_dot = [-0.01_dp, 0.03_dp, -0.02_dp, 0.04_dp]
    rhs_second_dot = [0.05_dp, -0.04_dp, 0.01_dp]
    call assemble_block_2x2_residual_jvp(aa, ab, ba, bb, x, y, rhs_first, rhs_second, &
        aa_dot, ab_dot, ba_dot, bb_dot, x_dot, y_dot, rhs_first_dot, rhs_second_dot, &
        residual_first_dot, residual_second_dot, status)
    call record_condition(status%code == 0, "2x2 block residual JVP assembles")

    epsilon = 1.0e-7_dp
    call assemble_block_2x2_residual(aa + epsilon*aa_dot, ab + epsilon*ab_dot, &
        ba + epsilon*ba_dot, bb + epsilon*bb_dot, x + epsilon*x_dot, y + epsilon*y_dot, &
        rhs_first + epsilon*rhs_first_dot, rhs_second + epsilon*rhs_second_dot, &
        residual_first_plus, residual_second_plus, status)
    finite_difference_error = max( &
        maxval(abs(residual_first_dot - &
        (residual_first_plus - residual_first)/epsilon)), &
        maxval(abs(residual_second_dot - &
        (residual_second_plus - residual_second)/epsilon)))
    call record_condition(finite_difference_error < 3.0e-8_dp, &
        "2x2 block residual JVP matches a forward difference")

    residual_first_bar = [0.2_dp, -0.3_dp, 0.4_dp, -0.1_dp]
    residual_second_bar = [-0.1_dp, 0.5_dp, 0.3_dp]
    call assemble_block_2x2_residual_vjp(aa, ab, ba, bb, x, y, rhs_first, rhs_second, &
        residual_first_bar, residual_second_bar, aa_bar, ab_bar, ba_bar, bb_bar, &
        x_bar, y_bar, rhs_first_bar, rhs_second_bar, status)
    call record_condition(status%code == 0, "2x2 block residual VJP assembles")
    lhs = dot_product(residual_first_bar, residual_first_dot) + &
        dot_product(residual_second_bar, residual_second_dot)
    rhs = sum(aa_bar*aa_dot) + sum(ab_bar*ab_dot) + sum(ba_bar*ba_dot) + &
        sum(bb_bar*bb_dot) + dot_product(x_bar, x_dot) + dot_product(y_bar, y_dot) + &
        dot_product(rhs_first_bar, rhs_first_dot) + dot_product(rhs_second_bar, rhs_second_dot)
    call record_condition(abs(lhs - rhs) < 3.0e-13_dp, &
        "2x2 block residual VJP satisfies the adjoint identity")

    call assemble_block_2x2_residual(aa(:first_rows - 1, :), ab, ba, bb, x, y, &
        rhs_first, rhs_second, residual_first, residual_second, status)
    call record_condition(status%code /= 0, "incompatible 2x2 block dimensions are rejected")

    call check_summary("2x2 block residual")
    if (.not. all_passed) error stop 1

contains

    subroutine build_system(aa, ab, ba, bb, x, y, rhs_first, rhs_second)
        real(dp), intent(out) :: aa(:, :), ab(:, :), ba(:, :), bb(:, :)
        real(dp), intent(out) :: x(:), y(:), rhs_first(:), rhs_second(:)

        aa = reshape([ &
            1.0_dp, -0.2_dp, 0.3_dp, 0.4_dp, -0.5_dp, 0.6_dp, &
            0.7_dp, 0.1_dp, -0.8_dp, 0.2_dp, 0.9_dp, -0.3_dp], shape(aa))
        ab = reshape([0.4_dp, -0.1_dp, 0.2_dp, 0.8_dp, -0.6_dp, 0.3_dp, &
            0.5_dp, -0.7_dp], shape(ab))
        ba = reshape([0.2_dp, -0.3_dp, 0.5_dp, 0.7_dp, -0.4_dp, 0.1_dp, &
            0.6_dp, 0.8_dp, -0.2_dp], shape(ba))
        bb = reshape([-0.5_dp, 0.2_dp, 0.4_dp, -0.7_dp, 0.3_dp, 0.9_dp], shape(bb))
        x = [0.6_dp, -0.4_dp, 0.2_dp]
        y = [0.9_dp, -0.1_dp]
        rhs_first = [0.1_dp, -0.3_dp, 0.5_dp, -0.2_dp]
        rhs_second = [-0.2_dp, 0.4_dp, 0.3_dp]
    end subroutine build_system

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_block_2x2_residual
