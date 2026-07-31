program test_circular_helmholtz_dtn_ad
    use check, only: check_condition, check_summary
    use fortfem_api, only: &
        apply_circular_helmholtz_dtn, apply_circular_helmholtz_dtn_jvp, &
        apply_circular_helmholtz_dtn_vjp, circular_helmholtz_dtn_eigenvalue, &
        circular_helmholtz_dtn_eigenvalue_jvp, &
        circular_helmholtz_dtn_eigenvalue_vjp
    use fortfem_kinds, only: dp
    implicit none

    integer, parameter :: point_count = 8, maximum_mode = 2
    real(dp), parameter :: step = 2.0e-6_dp
    real(dp) :: wavenumber, radius, wavenumber_dot, radius_dot
    real(dp) :: wavenumber_bar, radius_bar, discarded_bar
    real(dp) :: discarded, discarded_dot, discarded_plus, discarded_minus
    real(dp) :: lhs, rhs, relative_error
    complex(dp) :: trace(point_count), trace_dot(point_count)
    complex(dp) :: normal_derivative(point_count), normal_derivative_dot(point_count)
    complex(dp) :: normal_derivative_plus(point_count), normal_derivative_minus(point_count)
    complex(dp) :: normal_derivative_bar(point_count), trace_bar(point_count)
    complex(dp) :: eigenvalue, eigenvalue_dot, eigenvalue_plus, eigenvalue_minus
    complex(dp) :: eigenvalue_bar
    integer :: index, status

    wavenumber = 1.7_dp
    radius = 0.9_dp
    wavenumber_dot = 0.13_dp
    radius_dot = -0.08_dp
    do index = 1, point_count
        trace(index) = cmplx( &
            0.2_dp*cos(0.4_dp*real(index, dp)), &
            -0.15_dp*sin(0.3_dp*real(index, dp)), dp)
        trace_dot(index) = cmplx( &
            -0.07_dp*sin(0.2_dp*real(index, dp)), &
            0.04_dp*cos(0.6_dp*real(index, dp)), dp)
        normal_derivative_bar(index) = cmplx( &
            0.11_dp*cos(0.5_dp*real(index, dp)), &
            -0.08_dp*sin(0.7_dp*real(index, dp)), dp)
    end do

    call circular_helmholtz_dtn_eigenvalue( &
        -3, wavenumber, radius, eigenvalue, status)
    call circular_helmholtz_dtn_eigenvalue_jvp( &
        -3, wavenumber, radius, wavenumber_dot, radius_dot, eigenvalue_dot, &
        status)
    call circular_helmholtz_dtn_eigenvalue( &
        -3, wavenumber + step*wavenumber_dot, radius + step*radius_dot, &
        eigenvalue_plus, status)
    call circular_helmholtz_dtn_eigenvalue( &
        -3, wavenumber - step*wavenumber_dot, radius - step*radius_dot, &
        eigenvalue_minus, status)
    relative_error = abs(eigenvalue_dot - &
        (eigenvalue_plus - eigenvalue_minus)/(2.0_dp*step))/ &
        max(1.0_dp, abs(eigenvalue_dot))
    call check_condition(status == 0 .and. relative_error < 3.0e-7_dp, &
        "Circular DtN eigenvalue JVP matches re-evaluation")

    eigenvalue_bar = cmplx(0.37_dp, -0.22_dp, dp)
    call circular_helmholtz_dtn_eigenvalue_vjp( &
        -3, wavenumber, radius, eigenvalue_bar, wavenumber_bar, radius_bar, &
        status)
    lhs = real(conjg(eigenvalue_bar)*eigenvalue_dot, dp)
    rhs = wavenumber_bar*wavenumber_dot + radius_bar*radius_dot
    call check_condition(status == 0 .and. &
        abs(lhs - rhs) < 2.0e-8_dp*max(1.0_dp, abs(lhs), abs(rhs)), &
        "Circular DtN eigenvalue VJP satisfies the adjoint identity")

    call apply_circular_helmholtz_dtn( &
        trace, wavenumber, radius, maximum_mode, normal_derivative, discarded, &
        status)
    call apply_circular_helmholtz_dtn_jvp( &
        trace, wavenumber, radius, maximum_mode, trace_dot, wavenumber_dot, &
        radius_dot, normal_derivative_dot, discarded_dot, status)
    call apply_circular_helmholtz_dtn( &
        trace + step*trace_dot, wavenumber + step*wavenumber_dot, &
        radius + step*radius_dot, maximum_mode, normal_derivative_plus, &
        discarded_plus, status)
    call apply_circular_helmholtz_dtn( &
        trace - step*trace_dot, wavenumber - step*wavenumber_dot, &
        radius - step*radius_dot, maximum_mode, normal_derivative_minus, &
        discarded_minus, status)
    relative_error = max( &
        maxval(abs(normal_derivative_dot - &
        (normal_derivative_plus - normal_derivative_minus)/(2.0_dp*step))), &
        abs(discarded_dot - (discarded_plus - discarded_minus)/(2.0_dp*step)))/ &
        max(1.0_dp, maxval(abs(normal_derivative_dot)), abs(discarded_dot))
    call check_condition(status == 0 .and. relative_error < 2.0e-6_dp, &
        "Circular DtN apply JVP matches re-evaluation")

    discarded_bar = 0.23_dp
    call apply_circular_helmholtz_dtn_vjp( &
        trace, wavenumber, radius, maximum_mode, normal_derivative, &
        normal_derivative_bar, discarded, discarded_bar, trace_bar, &
        wavenumber_bar, radius_bar, status)
    lhs = real(sum(conjg(normal_derivative_bar)*normal_derivative_dot), dp) + &
        discarded_bar*discarded_dot
    rhs = real(sum(conjg(trace_bar)*trace_dot), dp) + &
        wavenumber_bar*wavenumber_dot + radius_bar*radius_dot
    call check_condition(status == 0 .and. &
        abs(lhs - rhs) < 2.0e-7_dp*max(1.0_dp, abs(lhs), abs(rhs)), &
        "Circular DtN apply VJP satisfies the adjoint identity")

    call check_summary("Differentiable circular Helmholtz DtN")
end program test_circular_helmholtz_dtn_ad
