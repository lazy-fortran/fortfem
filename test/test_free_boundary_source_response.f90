program test_free_boundary_source_response
    use, intrinsic :: iso_fortran_env, only: real64
    use check, only: check_condition, check_summary
    use fortfem_boundary, only: &
        assemble_free_boundary_source_response, &
        assemble_free_boundary_source_response_jvp, &
        assemble_free_boundary_source_response_vjp
    use fortsparse, only: fortsparse_status_t, FORTSPARSE_OK
    implicit none

    integer, parameter :: dp = real64
    integer, parameter :: nq = 2, nc = 2, ns = 2
    real(dp), parameter :: epsilon = 1.0e-7_dp
    real(dp) :: basis(nq, nc, ns), basis_dot(nq, nc, ns)
    real(dp) :: coefficients(ns), coefficients_dot(ns)
    real(dp) :: target(nq, nc), target_dot(nq, nc), weights(nq), weights_dot(nq)
    real(dp) :: sheet(nq, nc), sheet_dot(nq, nc)
    real(dp) :: basis_plus(nq, nc, ns), basis_minus(nq, nc, ns)
    real(dp) :: coefficients_plus(ns), coefficients_minus(ns)
    real(dp) :: target_plus(nq, nc), target_minus(nq, nc)
    real(dp) :: weights_plus(nq), weights_minus(nq)
    real(dp) :: sheet_plus(nq, nc), sheet_minus(nq, nc)
    real(dp), allocatable :: trace(:,:), residual(:,:), trace_dot(:,:), residual_dot(:,:)
    real(dp), allocatable :: trace_plus(:,:), trace_minus(:,:)
    real(dp), allocatable :: residual_plus(:,:), residual_minus(:,:)
    real(dp), allocatable :: basis_bar(:,:,:), coefficients_bar(:), target_bar(:,:)
    real(dp), allocatable :: weights_bar(:), sheet_bar(:,:)
    real(dp) :: trace_bar(nq, nc), residual_bar(nq, nc)
    real(dp) :: lhs, rhs, finite_difference_error
    type(fortsparse_status_t) :: status

    basis = reshape([ &
        1.0_dp, -0.2_dp, 0.4_dp, 0.7_dp, &
        -0.3_dp, 0.5_dp, 0.8_dp, -0.6_dp], shape(basis))
    basis_dot = reshape([ &
        0.02_dp, -0.01_dp, 0.03_dp, 0.04_dp, &
        -0.05_dp, 0.01_dp, 0.02_dp, -0.03_dp], shape(basis_dot))
    coefficients = [0.6_dp, -0.4_dp]
    coefficients_dot = [-0.03_dp, 0.02_dp]
    target = reshape([0.1_dp, -0.2_dp, 0.3_dp, 0.4_dp], shape(target))
    target_dot = reshape([-0.02_dp, 0.01_dp, 0.04_dp, -0.03_dp], shape(target_dot))
    weights = [1.5_dp, 0.75_dp]
    weights_dot = [0.03_dp, -0.02_dp]
    sheet = reshape([0.05_dp, -0.1_dp, 0.02_dp, 0.07_dp], shape(sheet))
    sheet_dot = reshape([-0.01_dp, 0.03_dp, 0.02_dp, -0.04_dp], shape(sheet_dot))

    call assemble_free_boundary_source_response( &
        basis, coefficients, target, weights, trace, residual, status, sheet)
    call check_condition(status%code == FORTSPARSE_OK .and. &
        maxval(abs(trace - expected_trace(basis, coefficients))) < 1.0e-14_dp .and. &
        maxval(abs(residual - expected_residual(trace, target, weights, sheet))) < &
        1.0e-14_dp, "source-to-free-boundary composition matches independent oracle")

    call assemble_free_boundary_source_response_jvp( &
        basis, coefficients, target, weights, basis_dot, coefficients_dot, &
        target_dot, weights_dot, trace_dot, residual_dot, status, sheet, sheet_dot)
    basis_plus = basis + epsilon*basis_dot
    basis_minus = basis - epsilon*basis_dot
    coefficients_plus = coefficients + epsilon*coefficients_dot
    coefficients_minus = coefficients - epsilon*coefficients_dot
    target_plus = target + epsilon*target_dot
    target_minus = target - epsilon*target_dot
    weights_plus = weights + epsilon*weights_dot
    weights_minus = weights - epsilon*weights_dot
    sheet_plus = sheet + epsilon*sheet_dot
    sheet_minus = sheet - epsilon*sheet_dot
    call assemble_free_boundary_source_response( &
        basis_plus, coefficients_plus, target_plus, weights_plus, &
        trace_plus, residual_plus, status, sheet_plus)
    call assemble_free_boundary_source_response( &
        basis_minus, coefficients_minus, target_minus, weights_minus, &
        trace_minus, residual_minus, status, sheet_minus)
    finite_difference_error = max( &
        maxval(abs(trace_dot - (trace_plus - trace_minus)/(2.0_dp*epsilon))), &
        maxval(abs(residual_dot - (residual_plus - residual_minus)/(2.0_dp*epsilon))))
    call check_condition(status%code == FORTSPARSE_OK .and. &
        finite_difference_error < 2.0e-8_dp, &
        "source-to-free-boundary JVP matches central differences")

    trace_bar = reshape([0.7_dp, -0.1_dp, 0.2_dp, 0.6_dp], shape(trace_bar))
    residual_bar = reshape([-0.4_dp, 0.3_dp, 0.5_dp, -0.2_dp], shape(residual_bar))
    call assemble_free_boundary_source_response_vjp( &
        basis, coefficients, target, weights, trace_bar, residual_bar, &
        basis_bar, coefficients_bar, target_bar, weights_bar, status, sheet, sheet_bar)
    lhs = sum(trace_bar*trace_dot) + sum(residual_bar*residual_dot)
    rhs = sum(basis_bar*basis_dot) + sum(coefficients_bar*coefficients_dot) + &
        sum(target_bar*target_dot) + sum(weights_bar*weights_dot) + &
        sum(sheet_bar*sheet_dot)
    call check_condition(status%code == FORTSPARSE_OK .and. abs(lhs - rhs) < 2.0e-12_dp, &
        "source-to-free-boundary VJP satisfies the transpose oracle")

    call check_summary("free-boundary source response")

contains

    pure function expected_trace(basis, coefficients) result(trace)
        real(dp), intent(in) :: basis(:, :, :), coefficients(:)
        real(dp) :: trace(size(basis, 1), size(basis, 2))
        integer :: q, c, s

        trace = 0.0_dp
        do q = 1, size(basis, 1)
            do c = 1, size(basis, 2)
                do s = 1, size(coefficients)
                    trace(q, c) = trace(q, c) + basis(q, c, s)*coefficients(s)
                end do
            end do
        end do
    end function expected_trace

    pure function expected_residual(trace, target, weights, sheet) result(residual)
        real(dp), intent(in) :: trace(:, :), target(:, :), weights(:), sheet(:, :)
        real(dp) :: residual(size(trace, 1), size(trace, 2))
        integer :: q

        do q = 1, size(weights)
            residual(q, :) = weights(q)*(trace(q, :) - target(q, :) - sheet(q, :))
        end do
    end function expected_residual

end program test_free_boundary_source_response
