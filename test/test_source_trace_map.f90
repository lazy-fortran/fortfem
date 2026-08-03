program test_source_trace_map
    use, intrinsic :: iso_fortran_env, only: real64
    use check, only: check_condition, check_summary
    use fortfem_boundary, only: &
        evaluate_source_trace_map, evaluate_source_trace_map_jvp, &
        evaluate_source_trace_map_vjp, &
        evaluate_weighted_source_trace_reciprocity_defect
    use fortsparse, only: fortsparse_status_t, FORTSPARSE_OK
    implicit none

    integer, parameter :: dp = real64, n = 2
    real(dp), parameter :: epsilon = 1.0e-7_dp
    complex(dp) :: source_matrix(n, n), source_matrix_dot(n, n)
    complex(dp) :: source_coefficients(n), source_coefficients_dot(n)
    complex(dp) :: source_matrix_plus(n, n), source_matrix_minus(n, n)
    complex(dp) :: source_coefficients_plus(n), source_coefficients_minus(n)
    complex(dp) :: trace_bar(n)
    complex(dp), allocatable :: trace(:), trace_dot_expected(:)
    complex(dp), allocatable :: trace_plus(:), trace_minus(:)
    complex(dp), allocatable :: source_matrix_bar(:, :), source_coefficients_bar(:)
    complex(dp) :: source_one(n), source_two(n), target_one(n), target_two(n)
    real(dp) :: weights(n), reciprocity_defect, expected_defect
    real(dp) :: lhs, rhs, finite_difference_error
    complex(dp) :: first_pairing, second_pairing
    type(fortsparse_status_t) :: status

    source_matrix = reshape([ &
        cmplx(1.2_dp, 0.1_dp, dp), cmplx(-0.3_dp, 0.4_dp, dp), &
        cmplx(0.5_dp, -0.2_dp, dp), cmplx(0.8_dp, 0.3_dp, dp)], shape(source_matrix))
    source_matrix_dot = reshape([ &
        cmplx(0.03_dp, -0.01_dp, dp), cmplx(-0.02_dp, 0.04_dp, dp), &
        cmplx(0.01_dp, 0.02_dp, dp), cmplx(-0.05_dp, 0.03_dp, dp)], &
        shape(source_matrix_dot))
    source_coefficients = [cmplx(0.4_dp, -0.2_dp, dp), cmplx(-0.6_dp, 0.3_dp, dp)]
    source_coefficients_dot = [cmplx(-0.02_dp, 0.03_dp, dp), &
        cmplx(0.05_dp, -0.01_dp, dp)]

    call evaluate_source_trace_map( &
        source_matrix, source_coefficients, trace, status)
    call check_condition(status%code == FORTSPARSE_OK .and. &
        maxval(abs(trace - matmul(source_matrix, source_coefficients))) < 1.0e-14_dp, &
        "source-trace value matches matrix-vector oracle")

    call evaluate_source_trace_map_jvp( &
        source_matrix, source_coefficients, source_matrix_dot, &
        source_coefficients_dot, trace_dot_expected, status)
    source_matrix_plus = source_matrix + epsilon*source_matrix_dot
    source_matrix_minus = source_matrix - epsilon*source_matrix_dot
    source_coefficients_plus = source_coefficients + epsilon*source_coefficients_dot
    source_coefficients_minus = source_coefficients - epsilon*source_coefficients_dot
    call evaluate_source_trace_map( &
        source_matrix_plus, source_coefficients_plus, trace_plus, status)
    call evaluate_source_trace_map( &
        source_matrix_minus, source_coefficients_minus, trace_minus, status)
    finite_difference_error = maxval(abs(trace_dot_expected - &
        (trace_plus - trace_minus)/(2.0_dp*epsilon)))
    call check_condition(status%code == FORTSPARSE_OK .and. &
        finite_difference_error < 2.0e-8_dp, &
        "source-trace JVP matches central differences")

    trace_bar = [cmplx(0.7_dp, -0.2_dp, dp), cmplx(-0.1_dp, 0.6_dp, dp)]
    call evaluate_source_trace_map_vjp( &
        source_matrix, source_coefficients, trace_bar, source_matrix_bar, &
        source_coefficients_bar, status)
    lhs = real(sum(conjg(trace_bar)*trace_dot_expected), dp)
    rhs = real(sum(conjg(source_matrix_bar)*source_matrix_dot) + &
        sum(conjg(source_coefficients_bar)*source_coefficients_dot), dp)
    call check_condition(status%code == FORTSPARSE_OK .and. &
        abs(lhs - rhs) < 2.0e-12_dp, &
        "source-trace VJP satisfies the real complex adjoint identity")

    source_one = [cmplx(0.2_dp, 0.1_dp, dp), cmplx(-0.4_dp, 0.3_dp, dp)]
    source_two = [cmplx(-0.3_dp, 0.5_dp, dp), cmplx(0.6_dp, -0.2_dp, dp)]
    target_one = [cmplx(0.8_dp, -0.1_dp, dp), cmplx(-0.2_dp, 0.7_dp, dp)]
    target_two = [cmplx(-0.5_dp, 0.4_dp, dp), cmplx(0.3_dp, 0.2_dp, dp)]
    weights = [1.5_dp, 0.75_dp]
    first_pairing = sum(target_one*weights*matmul(source_matrix, source_two))
    second_pairing = sum(target_two*weights*matmul(source_matrix, source_one))
    expected_defect = abs(first_pairing - second_pairing)/ &
        max(1.0_dp, abs(first_pairing), abs(second_pairing))
    call evaluate_weighted_source_trace_reciprocity_defect( &
        source_matrix, source_one, target_one, source_two, target_two, weights, &
        reciprocity_defect, status)
    call check_condition(status%code == FORTSPARSE_OK .and. &
        abs(reciprocity_defect - expected_defect) < 1.0e-14_dp, &
        "weighted reciprocity defect matches transpose-work oracle")

    call check_summary("source-trace map")
end program test_source_trace_map
