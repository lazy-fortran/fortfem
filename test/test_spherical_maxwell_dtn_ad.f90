program test_spherical_maxwell_dtn_ad
    use check, only: check_condition, check_summary
    use fortfem_boundary, only: apply_spherical_maxwell_dtn, &
        apply_spherical_maxwell_dtn_jvp, apply_spherical_maxwell_dtn_vjp
    use fortfem_kinds, only: dp
    implicit none

    integer, parameter :: count = 3
    real(dp), parameter :: step = 2.0e-7_dp
    integer :: degrees(count), status, status_plus, status_minus
    complex(dp) :: te_trace(count), tm_trace(count)
    complex(dp) :: te_trace_dot(count), tm_trace_dot(count)
    complex(dp) :: te_trace_bar(count), tm_trace_bar(count)
    complex(dp) :: te_result(count), tm_result(count)
    complex(dp) :: te_result_dot(count), tm_result_dot(count)
    complex(dp) :: te_result_bar(count), tm_result_bar(count)
    complex(dp) :: te_plus(count), tm_plus(count), te_minus(count), tm_minus(count)
    real(dp) :: wavenumber, radius, wavenumber_dot, radius_dot
    real(dp) :: wavenumber_bar, radius_bar, lhs, rhs, relative_error
    integer :: coefficient

    degrees = [1, 2, 4]
    wavenumber = 2.1_dp
    radius = 1.4_dp
    wavenumber_dot = 0.13_dp
    radius_dot = -0.09_dp
    do coefficient = 1, count
        te_trace(coefficient) = cmplx( &
            0.2_dp*coefficient, -0.07_dp*coefficient, dp)
        tm_trace(coefficient) = cmplx( &
            -0.11_dp*coefficient, 0.05_dp*coefficient, dp)
        te_trace_dot(coefficient) = cmplx( &
            -0.03_dp*coefficient, 0.04_dp*coefficient, dp)
        tm_trace_dot(coefficient) = cmplx( &
            0.02_dp*coefficient, -0.01_dp*coefficient, dp)
        te_result_bar(coefficient) = cmplx( &
            0.06_dp*coefficient, 0.02_dp*coefficient, dp)
        tm_result_bar(coefficient) = cmplx( &
            -0.04_dp*coefficient, 0.03_dp*coefficient, dp)
    end do

    call apply_spherical_maxwell_dtn( &
        degrees, te_trace, tm_trace, wavenumber, radius, te_result, tm_result, &
        status)
    call apply_spherical_maxwell_dtn_jvp( &
        degrees, te_trace, tm_trace, wavenumber, radius, te_trace_dot, &
        tm_trace_dot, wavenumber_dot, radius_dot, te_result_dot, &
        tm_result_dot, status)
    call apply_spherical_maxwell_dtn( &
        degrees, te_trace + step*te_trace_dot, tm_trace + step*tm_trace_dot, &
        wavenumber + step*wavenumber_dot, radius + step*radius_dot, te_plus, &
        tm_plus, status_plus)
    call apply_spherical_maxwell_dtn( &
        degrees, te_trace - step*te_trace_dot, tm_trace - step*tm_trace_dot, &
        wavenumber - step*wavenumber_dot, radius - step*radius_dot, te_minus, &
        tm_minus, status_minus)
    relative_error = max( &
        maxval(abs(te_result_dot - (te_plus - te_minus)/(2.0_dp*step))), &
        maxval(abs(tm_result_dot - (tm_plus - tm_minus)/(2.0_dp*step))))/ &
        max(1.0_dp, maxval(abs(te_result_dot)), maxval(abs(tm_result_dot)))
    call check_condition( &
        status == 0 .and. status_plus == 0 .and. status_minus == 0, &
        "Spherical Maxwell DtN JVP accepts valid modes")
    call check_condition( &
        relative_error < 2.0e-8_dp, &
        "Spherical Maxwell DtN JVP matches a complete central difference")

    call apply_spherical_maxwell_dtn_vjp( &
        degrees, te_trace, tm_trace, wavenumber, radius, te_result_bar, &
        tm_result_bar, te_trace_bar, tm_trace_bar, wavenumber_bar, radius_bar, &
        status)
    lhs = real(sum(conjg(te_result_bar)*te_result_dot) + &
        sum(conjg(tm_result_bar)*tm_result_dot), dp)
    rhs = real(sum(conjg(te_trace_bar)*te_trace_dot) + &
        sum(conjg(tm_trace_bar)*tm_trace_dot), dp) + &
        wavenumber_bar*wavenumber_dot + radius_bar*radius_dot
    call check_condition(status == 0, "Spherical Maxwell DtN VJP succeeds")
    call check_condition( &
        abs(lhs - rhs) < 2.0e-11_dp*max(1.0_dp, abs(lhs), abs(rhs)), &
        "Spherical Maxwell DtN products satisfy the complex adjoint identity")

    call check_summary("Spherical Maxwell DtN derivatives")
end program test_spherical_maxwell_dtn_ad
