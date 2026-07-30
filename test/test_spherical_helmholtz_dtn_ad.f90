program test_spherical_helmholtz_dtn_ad
    use check, only: check_condition, check_summary
    use fortfem_api, only: apply_spherical_helmholtz_dtn, &
        apply_spherical_helmholtz_dtn_jvp, apply_spherical_helmholtz_dtn_vjp
    use fortfem_kinds, only: dp
    implicit none

    integer, parameter :: count = 4
    real(dp), parameter :: step = 2.0e-7_dp
    integer :: degrees(count), status, status_plus, status_minus
    complex(dp) :: trace(count), trace_dot(count), trace_bar(count)
    complex(dp) :: result(count), result_dot(count), result_bar(count)
    complex(dp) :: plus(count), minus(count)
    real(dp) :: wavenumber, radius, wavenumber_dot, radius_dot
    real(dp) :: wavenumber_bar, radius_bar, lhs, rhs, relative_error
    integer :: coefficient

    degrees = [0, 1, 2, 4]
    wavenumber = 2.1_dp
    radius = 1.4_dp
    wavenumber_dot = 0.13_dp
    radius_dot = -0.09_dp
    do coefficient = 1, count
        trace(coefficient) = cmplx( &
            0.2_dp*coefficient, -0.07_dp*coefficient, dp)
        trace_dot(coefficient) = cmplx( &
            -0.03_dp*coefficient, 0.04_dp*coefficient, dp)
        result_bar(coefficient) = cmplx( &
            0.06_dp*coefficient, 0.02_dp*coefficient, dp)
    end do

    call apply_spherical_helmholtz_dtn( &
        degrees, trace, wavenumber, radius, result, status)
    call apply_spherical_helmholtz_dtn_jvp( &
        degrees, trace, wavenumber, radius, trace_dot, wavenumber_dot, &
        radius_dot, result_dot, status)
    call apply_spherical_helmholtz_dtn( &
        degrees, trace + step*trace_dot, wavenumber + step*wavenumber_dot, &
        radius + step*radius_dot, plus, status_plus)
    call apply_spherical_helmholtz_dtn( &
        degrees, trace - step*trace_dot, wavenumber - step*wavenumber_dot, &
        radius - step*radius_dot, minus, status_minus)
    relative_error = maxval(abs( &
        result_dot - (plus - minus)/(2.0_dp*step)))/ &
        max(1.0_dp, maxval(abs(result_dot)))
    call check_condition( &
        status == 0 .and. status_plus == 0 .and. status_minus == 0, &
        "Spherical Helmholtz DtN JVP accepts valid modes")
    call check_condition( &
        relative_error < 2.0e-8_dp, &
        "Spherical Helmholtz DtN JVP matches a complete central difference")

    call apply_spherical_helmholtz_dtn_vjp( &
        degrees, trace, wavenumber, radius, result_bar, trace_bar, &
        wavenumber_bar, radius_bar, status)
    lhs = real(sum(conjg(result_bar)*result_dot), dp)
    rhs = real(sum(conjg(trace_bar)*trace_dot), dp) + &
        wavenumber_bar*wavenumber_dot + radius_bar*radius_dot
    call check_condition(status == 0, "Spherical Helmholtz DtN VJP succeeds")
    call check_condition( &
        abs(lhs - rhs) < 2.0e-11_dp*max(1.0_dp, abs(lhs), abs(rhs)), &
        "Spherical Helmholtz DtN products satisfy the complex adjoint identity")

    call check_summary("Spherical Helmholtz DtN derivatives")
end program test_spherical_helmholtz_dtn_ad
