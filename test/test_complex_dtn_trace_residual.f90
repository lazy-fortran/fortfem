program test_complex_dtn_trace_residual
    use, intrinsic :: ieee_arithmetic, only: ieee_quiet_nan, ieee_value
    use check, only: check_condition, check_summary
    use fortfem_boundary, only: &
        assemble_complex_dtn_trace_residual, &
        assemble_complex_dtn_trace_residual_jvp, &
        assemble_complex_dtn_trace_residual_vjp
    use fortfem_kinds, only: dp
    use fortsparse, only: FORTSPARSE_OK, fortsparse_status_t
    implicit none

    integer, parameter :: source_count = 2, trace_count = 3
    real(dp), parameter :: difference_step = 1.0e-7_dp
    complex(dp) :: fem_trace(trace_count), fem_trace_bar(trace_count)
    complex(dp) :: fem_trace_dot(trace_count)
    complex(dp) :: map(trace_count, source_count)
    complex(dp) :: map_bar(trace_count, source_count)
    complex(dp) :: map_dot(trace_count, source_count)
    complex(dp) :: map_dot_invalid(trace_count, source_count)
    complex(dp) :: map_invalid(trace_count, source_count)
    complex(dp) :: residual(trace_count), residual_bar(trace_count)
    complex(dp) :: residual_bar_invalid(trace_count)
    complex(dp) :: residual_dot(trace_count), residual_expected(trace_count)
    complex(dp) :: residual_minus(trace_count), residual_plus(trace_count)
    complex(dp) :: short_residual(trace_count - 1)
    complex(dp) :: source(source_count), source_bar(source_count)
    complex(dp) :: source_dot(source_count)
    real(dp) :: forward_pairing, reverse_pairing
    real(dp) :: weights(trace_count), weights_bar(trace_count)
    real(dp) :: weights_dot(trace_count), weights_invalid(trace_count)
    type(fortsparse_status_t) :: status, status_minus, status_plus
    integer :: column, row
    logical :: all_passed

    all_passed = .true.
    map = reshape([ &
        cmplx(0.8_dp, -0.2_dp, dp), cmplx(-0.3_dp, 0.4_dp, dp), &
        cmplx(0.5_dp, 0.1_dp, dp), cmplx(-0.6_dp, -0.3_dp, dp), &
        cmplx(0.2_dp, 0.7_dp, dp), cmplx(0.9_dp, -0.5_dp, dp)], &
        shape(map))
    map_dot = reshape([ &
        cmplx(0.03_dp, 0.01_dp, dp), cmplx(-0.02_dp, 0.04_dp, dp), &
        cmplx(0.01_dp, -0.03_dp, dp), cmplx(-0.04_dp, 0.02_dp, dp), &
        cmplx(0.05_dp, -0.01_dp, dp), cmplx(-0.01_dp, -0.02_dp, dp)], &
        shape(map_dot))
    source = [cmplx(0.4_dp, -0.3_dp, dp), cmplx(-0.7_dp, 0.2_dp, dp)]
    source_dot = [cmplx(-0.02_dp, 0.05_dp, dp), &
        cmplx(0.04_dp, -0.01_dp, dp)]
    fem_trace = [cmplx(0.6_dp, 0.2_dp, dp), &
        cmplx(-0.1_dp, 0.8_dp, dp), cmplx(0.3_dp, -0.5_dp, dp)]
    fem_trace_dot = [cmplx(0.01_dp, -0.04_dp, dp), &
        cmplx(0.03_dp, 0.02_dp, dp), cmplx(-0.05_dp, 0.01_dp, dp)]
    weights = [0.7_dp, 1.4_dp, 2.1_dp]
    weights_dot = [0.04_dp, -0.03_dp, 0.02_dp]

    call assemble_complex_dtn_trace_residual( &
        map, source, fem_trace, weights, residual, status)
    residual_expected = fem_trace
    do row = 1, trace_count
        do column = 1, source_count
            residual_expected(row) = residual_expected(row) - &
                map(row, column)*source(column)
        end do
        residual_expected(row) = weights(row)*residual_expected(row)
    end do
    call record_condition(status%code == FORTSPARSE_OK .and. &
        maxval(abs(residual - residual_expected)) < 1.0e-14_dp, &
        "complex FEM/DtN trace residual matches the manual loop oracle")

    call assemble_complex_dtn_trace_residual_jvp( &
        map, source, fem_trace, weights, map_dot, source_dot, fem_trace_dot, &
        weights_dot, residual_dot, status)
    call assemble_complex_dtn_trace_residual( &
        map + difference_step*map_dot, source + difference_step*source_dot, &
        fem_trace + difference_step*fem_trace_dot, &
        weights + difference_step*weights_dot, residual_plus, status_plus)
    call assemble_complex_dtn_trace_residual( &
        map - difference_step*map_dot, source - difference_step*source_dot, &
        fem_trace - difference_step*fem_trace_dot, &
        weights - difference_step*weights_dot, residual_minus, status_minus)
    call record_condition(status%code == FORTSPARSE_OK .and. &
        status_plus%code == FORTSPARSE_OK .and. &
        status_minus%code == FORTSPARSE_OK, &
        "complex FEM/DtN trace JVP accepts valid perturbations")
    call record_condition(maxval(abs(residual_dot - &
        (residual_plus - residual_minus)/(2.0_dp*difference_step))) < &
        2.0e-9_dp, "complex FEM/DtN trace JVP matches full reassembly")

    residual_bar = [cmplx(-0.2_dp, 0.5_dp, dp), &
        cmplx(0.7_dp, -0.1_dp, dp), cmplx(-0.4_dp, -0.3_dp, dp)]
    call assemble_complex_dtn_trace_residual_vjp( &
        map, source, fem_trace, weights, residual_bar, map_bar, source_bar, &
        fem_trace_bar, weights_bar, status)
    forward_pairing = real(sum(conjg(residual_bar)*residual_dot), dp)
    reverse_pairing = real(sum(conjg(map_bar)*map_dot) + &
        sum(conjg(source_bar)*source_dot) + &
        sum(conjg(fem_trace_bar)*fem_trace_dot), dp) + &
        dot_product(weights_bar, weights_dot)
    call record_condition(status%code == FORTSPARSE_OK .and. &
        abs(forward_pairing - reverse_pairing) < 2.0e-12_dp, &
        "complex FEM/DtN trace VJP satisfies the real adjoint identity")

    call assemble_complex_dtn_trace_residual( &
        map, source, fem_trace, weights, short_residual, status)
    call record_condition(status%code /= FORTSPARSE_OK .and. &
        all(short_residual == cmplx(0.0_dp, 0.0_dp, dp)), &
        "complex FEM/DtN trace residual rejects incompatible shapes")

    map_invalid = map
    map_invalid(1, 1) = cmplx( &
        ieee_value(0.0_dp, ieee_quiet_nan), 0.0_dp, dp)
    call assemble_complex_dtn_trace_residual( &
        map_invalid, source, fem_trace, weights, residual, status)
    call record_condition(status%code /= FORTSPARSE_OK .and. &
        all(residual == cmplx(0.0_dp, 0.0_dp, dp)), &
        "complex FEM/DtN trace residual rejects non-finite maps")

    weights_invalid = weights
    weights_invalid(2) = 0.0_dp
    call assemble_complex_dtn_trace_residual( &
        map, source, fem_trace, weights_invalid, residual, status)
    call record_condition(status%code /= FORTSPARSE_OK, &
        "complex FEM/DtN trace residual rejects nonpositive work weights")

    map_dot_invalid = map_dot
    map_dot_invalid(2, 1) = cmplx( &
        0.0_dp, ieee_value(0.0_dp, ieee_quiet_nan), dp)
    call assemble_complex_dtn_trace_residual_jvp( &
        map, source, fem_trace, weights, map_dot_invalid, source_dot, &
        fem_trace_dot, weights_dot, residual_dot, status)
    call record_condition(status%code /= FORTSPARSE_OK .and. &
        all(residual_dot == cmplx(0.0_dp, 0.0_dp, dp)), &
        "complex FEM/DtN trace JVP rejects non-finite directions")

    residual_bar_invalid = residual_bar
    residual_bar_invalid(3) = cmplx( &
        ieee_value(0.0_dp, ieee_quiet_nan), 0.0_dp, dp)
    call assemble_complex_dtn_trace_residual_vjp( &
        map, source, fem_trace, weights, residual_bar_invalid, map_bar, &
        source_bar, fem_trace_bar, weights_bar, status)
    call record_condition(status%code /= FORTSPARSE_OK .and. &
        all(map_bar == cmplx(0.0_dp, 0.0_dp, dp)) .and. &
        all(source_bar == cmplx(0.0_dp, 0.0_dp, dp)) .and. &
        all(fem_trace_bar == cmplx(0.0_dp, 0.0_dp, dp)) .and. &
        all(weights_bar == 0.0_dp), &
        "complex FEM/DtN trace VJP rejects non-finite adjoints")

    call check_summary("complex FEM/DtN boundary trace residual")
    if (.not. all_passed) error stop 1

contains

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_complex_dtn_trace_residual
