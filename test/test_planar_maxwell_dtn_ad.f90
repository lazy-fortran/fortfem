program test_planar_maxwell_dtn_ad
    use fortfem_boundary, only: apply_planar_maxwell_dtn, &
        apply_planar_maxwell_dtn_jvp, apply_planar_maxwell_dtn_vjp
    use fortfem_kinds, only: dp
    use check, only: check_condition, check_summary
    implicit none

    integer, parameter :: nx = 3, ny = 3
    complex(dp) :: trace(2, nx, ny), trace_dot(2, nx, ny)
    complex(dp) :: derivative(2, nx, ny), derivative_dot(2, nx, ny)
    complex(dp) :: derivative_bar(2, nx, ny), trace_bar(2, nx, ny)
    complex(dp) :: plus(2, nx, ny), minus(2, nx, ny)
    real(dp), parameter :: step = 2.0e-7_dp
    real(dp) :: wave_number, length_x, length_y
    real(dp) :: wave_number_dot, length_x_dot, length_y_dot
    real(dp) :: wave_number_bar, length_x_bar, length_y_bar
    real(dp) :: lhs, rhs, relative_error
    integer :: component, i, j, status, status_plus, status_minus

    wave_number = 2.3_dp
    length_x = 4.7_dp
    length_y = 5.2_dp
    wave_number_dot = 0.17_dp
    length_x_dot = -0.11_dp
    length_y_dot = 0.08_dp
    do j = 1, ny
        do i = 1, nx
            do component = 1, 2
                trace(component, i, j) = cmplx( &
                    0.13_dp*component + 0.07_dp*i - 0.04_dp*j, &
                    -0.09_dp*component + 0.03_dp*i + 0.05_dp*j, dp)
                trace_dot(component, i, j) = cmplx( &
                    -0.05_dp*component + 0.02_dp*i + 0.01_dp*j, &
                    0.04_dp*component - 0.03_dp*i + 0.02_dp*j, dp)
                derivative_bar(component, i, j) = cmplx( &
                    0.06_dp*component - 0.01_dp*i + 0.03_dp*j, &
                    -0.02_dp*component + 0.04_dp*i - 0.01_dp*j, dp)
            end do
        end do
    end do

    call apply_planar_maxwell_dtn( &
        trace, wave_number, length_x, length_y, derivative, status)
    call apply_planar_maxwell_dtn_jvp( &
        trace, wave_number, length_x, length_y, trace_dot, wave_number_dot, &
        length_x_dot, length_y_dot, derivative_dot, status)
    call apply_planar_maxwell_dtn( &
        trace + step*trace_dot, wave_number + step*wave_number_dot, &
        length_x + step*length_x_dot, length_y + step*length_y_dot, &
        plus, status_plus)
    call apply_planar_maxwell_dtn( &
        trace - step*trace_dot, wave_number - step*wave_number_dot, &
        length_x - step*length_x_dot, length_y - step*length_y_dot, &
        minus, status_minus)
    relative_error = maxval(abs(derivative_dot - (plus - minus)/(2.0_dp*step)))/ &
        max(1.0_dp, maxval(abs(derivative_dot)))
    call check_condition( &
        status == 0 .and. status_plus == 0 .and. status_minus == 0, &
        "Planar Maxwell DtN JVP accepts valid inputs")
    call check_condition( &
        relative_error < 2.0e-8_dp, &
        "Planar Maxwell DtN JVP matches a complete central difference")

    call apply_planar_maxwell_dtn_vjp( &
        trace, wave_number, length_x, length_y, derivative_bar, trace_bar, &
        wave_number_bar, length_x_bar, length_y_bar, status)
    lhs = real(sum(conjg(derivative_bar)*derivative_dot), dp)
    rhs = real(sum(conjg(trace_bar)*trace_dot), dp) + &
        wave_number_bar*wave_number_dot + length_x_bar*length_x_dot + &
        length_y_bar*length_y_dot
    call check_condition(status == 0, "Planar Maxwell DtN VJP succeeds")
    call check_condition( &
        abs(lhs - rhs) < 2.0e-10_dp*max(1.0_dp, abs(lhs), abs(rhs)), &
        "Planar Maxwell DtN products satisfy the complex adjoint identity")

    call check_summary("Planar Maxwell DtN derivatives")
end program test_planar_maxwell_dtn_ad
