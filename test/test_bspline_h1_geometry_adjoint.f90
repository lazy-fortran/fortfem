program test_bspline_h1_geometry_adjoint
    use, intrinsic :: iso_fortran_env, only: error_unit
    use fortfem_api, only: assemble_bspline_h1_operator_csc, &
        assemble_bspline_h1_operator_csc_jvp, &
        assemble_bspline_h1_operator_csc_vjp, sparse_direct_factor_csc, &
        sparse_direct_factor_t, sparse_direct_factor_transpose_csc, &
        sparse_direct_free, sparse_direct_solve_factored, &
        sparse_direct_solve_factored_jvp, sparse_direct_solve_factored_vjp
    use fortfem_kinds, only: dp
    use fortsparse, only: csc_t, fortsparse_status_t
    implicit none

    real(dp), parameter :: knots_x(7) = [ &
        0.0_dp, 0.0_dp, 0.0_dp, 0.5_dp, 1.0_dp, 1.0_dp, 1.0_dp]
    real(dp), parameter :: knots_y(7) = knots_x
    real(dp), parameter :: step = 1.0e-6_dp
    real(dp), allocatable :: control_points(:, :, :)
    real(dp), allocatable :: control_points_bar(:, :, :)
    real(dp), allocatable :: control_points_dot(:, :, :)
    real(dp), allocatable :: control_points_minus(:, :, :)
    real(dp), allocatable :: control_points_plus(:, :, :)
    real(dp), allocatable :: weights(:, :), weights_bar(:, :)
    real(dp), allocatable :: weights_dot(:, :), weights_minus(:, :)
    real(dp), allocatable :: weights_plus(:, :)
    real(dp), allocatable :: rhs(:), rhs_bar(:), rhs_dot(:)
    real(dp), allocatable :: solution(:), solution_bar(:), solution_dot(:)
    real(dp), allocatable :: solution_minus(:), solution_plus(:)
    real(dp), allocatable :: matrix_values_bar(:)
    real(dp), allocatable :: x_points(:), y_points(:)
    real(dp) :: adjoint_geometry, adjoint_state, central_difference
    type(csc_t) :: matrix, matrix_bar, matrix_dot, matrix_minus, matrix_plus
    type(fortsparse_status_t) :: sparse_status
    type(sparse_direct_factor_t) :: factor, factor_minus, factor_plus
    type(sparse_direct_factor_t) :: transpose_factor
    integer :: failures, i, ix, iy, status

    failures = 0
    x_points = greville_abscissae(knots_x, 2)
    y_points = greville_abscissae(knots_y, 2)
    allocate (control_points(2, size(x_points), size(y_points)))
    allocate (weights(size(x_points), size(y_points)))
    do iy = 1, size(y_points)
        do ix = 1, size(x_points)
            control_points(:, ix, iy) = [ &
                1.0_dp + 1.2_dp*x_points(ix), &
                0.8_dp*y_points(iy) + 0.1_dp*x_points(ix)]
            weights(ix, iy) = 1.0_dp + &
                0.05_dp*sin(real(3*ix + 5*iy, dp))
        end do
    end do
    allocate (control_points_dot, mold=control_points)
    allocate (weights_dot, mold=weights)
    control_points_dot = reshape([ &
        (0.1_dp*sin(real(7*i + 2, dp)), i=1, size(control_points_dot))], &
        shape(control_points_dot))
    weights_dot = reshape([ &
        (0.03_dp*cos(real(5*i - 1, dp)), i=1, size(weights_dot))], &
        shape(weights_dot))

    call assemble_bspline_h1_operator_csc( &
        knots_x, knots_y, 2, 2, control_points, weights, 4, matrix, &
        sparse_status, stiffness_coefficient=1.0_dp, mass_coefficient=0.5_dp)
    call check(sparse_status%code == 0, "primal assembly succeeds")
    call assemble_bspline_h1_operator_csc_jvp( &
        knots_x, knots_y, 2, 2, control_points, weights, &
        control_points_dot, weights_dot, 4, matrix_dot, sparse_status, &
        stiffness_coefficient=1.0_dp, mass_coefficient=0.5_dp)
    call check(sparse_status%code == 0, "assembly JVP succeeds")

    allocate (rhs(matrix%nrow), rhs_dot(matrix%nrow))
    allocate (solution(matrix%nrow), solution_dot(matrix%nrow))
    allocate (solution_bar(matrix%nrow), rhs_bar(matrix%nrow))
    allocate (matrix_values_bar(matrix%nnz))
    rhs = [(1.0_dp + 0.1_dp*cos(real(2*i, dp)), i=1, matrix%nrow)]
    rhs_dot = 0.0_dp
    solution_bar = [(sin(real(3*i + 1, dp)), i=1, matrix%nrow)]
    call sparse_direct_factor_csc( &
        factor, matrix%nrow, matrix%col_ptr, matrix%row_idx, matrix%val, status)
    call check(status == 0, "primal sparse factorization succeeds")
    call sparse_direct_solve_factored(factor, rhs, solution, status)
    call check(status == 0, "primal state solve succeeds")
    call sparse_direct_solve_factored_jvp( &
        factor, matrix%nrow, matrix_dot%col_ptr, matrix_dot%row_idx, &
        matrix_dot%val, solution, rhs_dot, solution_dot, status)
    call check(status == 0, "state JVP succeeds")

    call sparse_direct_factor_transpose_csc( &
        transpose_factor, matrix%nrow, matrix%col_ptr, matrix%row_idx, &
        matrix%val, status)
    call check(status == 0, "adjoint factorization succeeds")
    call sparse_direct_solve_factored_vjp( &
        transpose_factor, matrix%nrow, matrix%col_ptr, matrix%row_idx, &
        solution, solution_bar, rhs_bar, matrix_values_bar, status)
    call check(status == 0, "state VJP succeeds")
    matrix_bar = matrix
    matrix_bar%val = matrix_values_bar
    allocate (control_points_bar, mold=control_points)
    allocate (weights_bar, mold=weights)
    call assemble_bspline_h1_operator_csc_vjp( &
        knots_x, knots_y, 2, 2, control_points, weights, 4, matrix_bar, &
        control_points_bar, weights_bar, sparse_status, &
        stiffness_coefficient=1.0_dp, mass_coefficient=0.5_dp)
    call check(sparse_status%code == 0, "assembly VJP succeeds")

    adjoint_state = dot_product(solution_bar, solution_dot)
    adjoint_geometry = sum(control_points_bar*control_points_dot) + &
        sum(weights_bar*weights_dot)
    call check(abs(adjoint_state - adjoint_geometry) < 2.0e-10_dp, &
        "geometry-to-state JVP and adjoint satisfy the chain-rule identity")

    allocate (control_points_plus, mold=control_points)
    allocate (control_points_minus, mold=control_points)
    allocate (weights_plus, mold=weights)
    allocate (weights_minus, mold=weights)
    control_points_plus = control_points + step*control_points_dot
    control_points_minus = control_points - step*control_points_dot
    weights_plus = weights + step*weights_dot
    weights_minus = weights - step*weights_dot
    call assemble_bspline_h1_operator_csc( &
        knots_x, knots_y, 2, 2, control_points_plus, weights_plus, 4, &
        matrix_plus, sparse_status, stiffness_coefficient=1.0_dp, &
        mass_coefficient=0.5_dp)
    call assemble_bspline_h1_operator_csc( &
        knots_x, knots_y, 2, 2, control_points_minus, weights_minus, 4, &
        matrix_minus, sparse_status, stiffness_coefficient=1.0_dp, &
        mass_coefficient=0.5_dp)
    allocate (solution_plus(matrix%nrow), solution_minus(matrix%nrow))
    call sparse_direct_factor_csc( &
        factor_plus, matrix_plus%nrow, matrix_plus%col_ptr, &
        matrix_plus%row_idx, matrix_plus%val, status)
    call sparse_direct_solve_factored( &
        factor_plus, rhs, solution_plus, status)
    call sparse_direct_factor_csc( &
        factor_minus, matrix_minus%nrow, matrix_minus%col_ptr, &
        matrix_minus%row_idx, matrix_minus%val, status)
    call sparse_direct_solve_factored( &
        factor_minus, rhs, solution_minus, status)
    central_difference = dot_product( &
        solution_bar, solution_plus - solution_minus)/(2.0_dp*step)
    call check(abs(adjoint_geometry - central_difference) < 8.0e-8_dp, &
        "geometry-to-state adjoint matches an independent central difference")

    call sparse_direct_free(factor)
    call sparse_direct_free(factor_plus)
    call sparse_direct_free(factor_minus)
    call sparse_direct_free(transpose_factor)
    if (failures > 0) then
        write (error_unit, "(i0,a)") failures, " test(s) failed"
        stop 1
    end if
    write (*, "(a)") "PASS"

contains

    pure function greville_abscissae(knots, degree) result(points)
        real(dp), intent(in) :: knots(:)
        integer, intent(in) :: degree
        real(dp), allocatable :: points(:)
        integer :: basis

        allocate (points(size(knots) - degree - 1))
        do basis = 1, size(points)
            points(basis) = &
                sum(knots(basis + 1:basis + degree))/real(degree, dp)
        end do
    end function greville_abscissae

    subroutine check(condition, label)
        logical, intent(in) :: condition
        character(*), intent(in) :: label

        if (condition) return
        failures = failures + 1
        write (error_unit, "(a)") "FAIL: "//label
    end subroutine check

end program test_bspline_h1_geometry_adjoint
