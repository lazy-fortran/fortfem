program test_maxwell_fem_bem_state_ad
    use check, only: check_condition, check_summary
    use fortfem_boundary, only: &
        solve_maxwell_fem_bem_linear_state, &
        solve_maxwell_fem_bem_linear_state_jvp, &
        solve_maxwell_fem_bem_linear_state_vjp
    use fortfem_kinds, only: dp
    implicit none

    integer, parameter :: n = 4
    real(dp), parameter :: step = 2.0e-5_dp
    complex(dp) :: matrix(n, n), matrix_dot(n, n)
    complex(dp) :: rhs(n), rhs_dot(n), state_bar(n)
    complex(dp), allocatable :: state(:), state_dot(:)
    complex(dp), allocatable :: state_plus(:), state_minus(:)
    complex(dp), allocatable :: matrix_bar(:, :), rhs_bar(:)
    real(dp) :: jvp_error, lhs, rhs_identity
    integer :: row, column, status, status_plus, status_minus

    matrix = cmplx(0.0_dp, 0.0_dp, dp)
    matrix_dot = cmplx(0.0_dp, 0.0_dp, dp)
    do row = 1, n
        matrix(row, row) = cmplx(2.5_dp + 0.2_dp*real(row, dp), &
            0.1_dp*real(row, dp), dp)
        matrix_dot(row, row) = cmplx(0.03_dp*real(row, dp), &
            -0.02_dp*real(row, dp), dp)
        do column = 1, n
            if (column == row) cycle
            matrix(row, column) = cmplx(0.04_dp*real(row + column, dp), &
                -0.01_dp*real(row - column, dp), dp)
            matrix_dot(row, column) = cmplx(0.005_dp*real(row - column, dp), &
                0.004_dp*real(row + column, dp), dp)
        end do
    end do
    rhs = [cmplx(0.7_dp, -0.2_dp, dp), cmplx(-0.4_dp, 0.5_dp, dp), &
        cmplx(0.3_dp, 0.15_dp, dp), cmplx(-0.2_dp, -0.35_dp, dp)]
    rhs_dot = [cmplx(-0.02_dp, 0.03_dp, dp), cmplx(0.04_dp, -0.01_dp, dp), &
        cmplx(0.015_dp, 0.025_dp, dp), cmplx(-0.01_dp, 0.02_dp, dp)]
    state_bar = [cmplx(0.2_dp, -0.1_dp, dp), cmplx(-0.3_dp, 0.25_dp, dp), &
        cmplx(0.15_dp, 0.35_dp, dp), cmplx(-0.1_dp, 0.05_dp, dp)]

    call solve_maxwell_fem_bem_linear_state(matrix, rhs, state, status)
    call check_condition(status == 0, "mixed FEM-BEM linear state solves")
    if (status /= 0) error stop 1
    call solve_maxwell_fem_bem_linear_state_jvp( &
        matrix, rhs, matrix_dot, rhs_dot, state_dot, status)
    call check_condition(status == 0, "mixed FEM-BEM linear state JVP solves")
    if (status /= 0) error stop 1
    call solve_maxwell_fem_bem_linear_state( &
        matrix + step*matrix_dot, rhs + step*rhs_dot, state_plus, status_plus)
    call solve_maxwell_fem_bem_linear_state( &
        matrix - step*matrix_dot, rhs - step*rhs_dot, state_minus, status_minus)
    jvp_error = maxval(abs( &
        state_dot - (state_plus - state_minus)/(2.0_dp*step)))
    call check_condition( &
        status_plus == 0 .and. status_minus == 0 .and. jvp_error < 2.0e-10_dp, &
        "mixed FEM-BEM linear state JVP matches central differences")

    call solve_maxwell_fem_bem_linear_state_vjp( &
        matrix, rhs, state_bar, matrix_bar, rhs_bar, status)
    call check_condition(status == 0, "mixed FEM-BEM linear state VJP solves")
    if (status /= 0) error stop 1
    lhs = real(sum(conjg(state_bar)*state_dot), dp)
    rhs_identity = real(sum(conjg(matrix_bar)*matrix_dot) + &
        sum(conjg(rhs_bar)*rhs_dot), dp)
    call check_condition( &
        abs(lhs - rhs_identity) < 3.0e-13_dp*max(1.0_dp, abs(lhs), abs(rhs_identity)), &
        "mixed FEM-BEM linear state products obey the real-complex adjoint identity")
    call check_summary("Maxwell FEM-BEM linear state differentiation")
end program test_maxwell_fem_bem_state_ad
