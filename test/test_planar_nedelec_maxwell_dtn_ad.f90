program test_planar_nedelec_maxwell_dtn_ad
    use check, only: check_condition, check_summary
    use fortfem_api, only: pullback_planar_maxwell_dtn_form, &
        pullback_planar_maxwell_dtn_form_jvp, &
        pullback_planar_maxwell_dtn_form_vjp
    use fortfem_kinds, only: dp
    implicit none

    integer, parameter :: nx = 3, ny = 3, trace_count = 2*nx*ny
    integer, parameter :: dof_count = 4
    real(dp), parameter :: step = 2.0e-7_dp
    complex(dp) :: sampling(trace_count, dof_count)
    complex(dp) :: sampling_dot(trace_count, dof_count)
    complex(dp) :: sampling_bar(trace_count, dof_count)
    complex(dp) :: form(dof_count, dof_count)
    complex(dp) :: form_dot(dof_count, dof_count)
    complex(dp) :: form_bar(dof_count, dof_count)
    complex(dp) :: plus(dof_count, dof_count), minus(dof_count, dof_count)
    real(dp) :: wave_number, length_x, length_y
    real(dp) :: wave_number_dot, length_x_dot, length_y_dot
    real(dp) :: wave_number_bar, length_x_bar, length_y_bar
    real(dp) :: lhs, rhs, relative_error
    integer :: column, row, status, status_plus, status_minus

    wave_number = 2.3_dp
    length_x = 4.7_dp
    length_y = 5.2_dp
    wave_number_dot = 0.17_dp
    length_x_dot = -0.11_dp
    length_y_dot = 0.08_dp
    do column = 1, dof_count
        do row = 1, trace_count
            sampling(row, column) = cmplx( &
                0.02_dp*row - 0.03_dp*column, &
                -0.01_dp*row + 0.015_dp*column, dp)
            sampling_dot(row, column) = cmplx( &
                -0.007_dp*row + 0.011_dp*column, &
                0.005_dp*row + 0.003_dp*column, dp)
        end do
        do row = 1, dof_count
            form_bar(row, column) = cmplx( &
                0.04_dp*row - 0.03_dp*column, &
                -0.02_dp*row + 0.01_dp*column, dp)
        end do
    end do

    call pullback_planar_maxwell_dtn_form( &
        sampling, nx, ny, wave_number, length_x, length_y, form, status)
    call pullback_planar_maxwell_dtn_form_jvp( &
        sampling, nx, ny, wave_number, length_x, length_y, sampling_dot, &
        wave_number_dot, length_x_dot, length_y_dot, form_dot, status)
    call pullback_planar_maxwell_dtn_form( &
        sampling + step*sampling_dot, nx, ny, &
        wave_number + step*wave_number_dot, length_x + step*length_x_dot, &
        length_y + step*length_y_dot, plus, status_plus)
    call pullback_planar_maxwell_dtn_form( &
        sampling - step*sampling_dot, nx, ny, &
        wave_number - step*wave_number_dot, length_x - step*length_x_dot, &
        length_y - step*length_y_dot, minus, status_minus)
    relative_error = maxval(abs(form_dot - (plus - minus)/(2.0_dp*step)))/ &
        max(1.0_dp, maxval(abs(form_dot)))
    call check_condition( &
        status == 0 .and. status_plus == 0 .and. status_minus == 0, &
        "Nedelec-DtN pullback JVP accepts valid inputs")
    call check_condition( &
        relative_error < 3.0e-8_dp, &
        "Nedelec-DtN pullback JVP matches a complete central difference")

    call pullback_planar_maxwell_dtn_form_vjp( &
        sampling, nx, ny, wave_number, length_x, length_y, form_bar, &
        sampling_bar, wave_number_bar, length_x_bar, length_y_bar, status)
    lhs = real(sum(conjg(form_bar)*form_dot), dp)
    rhs = real(sum(conjg(sampling_bar)*sampling_dot), dp) + &
        wave_number_bar*wave_number_dot + length_x_bar*length_x_dot + &
        length_y_bar*length_y_dot
    call check_condition(status == 0, "Nedelec-DtN pullback VJP succeeds")
    call check_condition( &
        abs(lhs - rhs) < 3.0e-10_dp*max(1.0_dp, abs(lhs), abs(rhs)), &
        "Nedelec-DtN pullback products satisfy the complex adjoint identity")

    call check_summary("Planar Nedelec-Maxwell DtN derivatives")
end program test_planar_nedelec_maxwell_dtn_ad
