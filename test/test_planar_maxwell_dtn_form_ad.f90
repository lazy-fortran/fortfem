program test_planar_maxwell_dtn_form_ad
    use check, only: check_condition, check_summary
    use fortfem_boundary, only: assemble_planar_maxwell_dtn_form, &
        assemble_planar_maxwell_dtn_form_jvp, &
        assemble_planar_maxwell_dtn_form_vjp
    use fortfem_kinds, only: dp
    implicit none

    integer, parameter :: nx = 3, ny = 3, n = 2*nx*ny
    real(dp), parameter :: step = 2.0e-7_dp
    complex(dp), allocatable :: form(:, :), form_dot(:, :)
    complex(dp), allocatable :: plus(:, :), minus(:, :)
    complex(dp) :: form_bar(n, n)
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
    do column = 1, n
        do row = 1, n
            form_bar(row, column) = cmplx( &
                0.003_dp*row - 0.002_dp*column, &
                -0.001_dp*row + 0.004_dp*column, dp)
        end do
    end do

    call assemble_planar_maxwell_dtn_form( &
        nx, ny, wave_number, length_x, length_y, form, status)
    call assemble_planar_maxwell_dtn_form_jvp( &
        nx, ny, wave_number, length_x, length_y, wave_number_dot, &
        length_x_dot, length_y_dot, form_dot, status)
    call assemble_planar_maxwell_dtn_form( &
        nx, ny, wave_number + step*wave_number_dot, &
        length_x + step*length_x_dot, length_y + step*length_y_dot, &
        plus, status_plus)
    call assemble_planar_maxwell_dtn_form( &
        nx, ny, wave_number - step*wave_number_dot, &
        length_x - step*length_x_dot, length_y - step*length_y_dot, &
        minus, status_minus)
    relative_error = maxval(abs(form_dot - (plus - minus)/(2.0_dp*step)))/ &
        max(1.0_dp, maxval(abs(form_dot)))
    call check_condition( &
        status == 0 .and. status_plus == 0 .and. status_minus == 0, &
        "Planar Maxwell DtN form JVP accepts valid inputs")
    call check_condition( &
        relative_error < 2.0e-8_dp, &
        "Planar Maxwell DtN form JVP matches a complete central difference")

    call assemble_planar_maxwell_dtn_form_vjp( &
        nx, ny, wave_number, length_x, length_y, form_bar, wave_number_bar, &
        length_x_bar, length_y_bar, status)
    lhs = real(sum(conjg(form_bar)*form_dot), dp)
    rhs = wave_number_bar*wave_number_dot + length_x_bar*length_x_dot + &
        length_y_bar*length_y_dot
    call check_condition(status == 0, "Planar Maxwell DtN form VJP succeeds")
    call check_condition( &
        abs(lhs - rhs) < 2.0e-10_dp*max(1.0_dp, abs(lhs), abs(rhs)), &
        "Planar Maxwell DtN form products satisfy the complex adjoint identity")

    call check_summary("Planar Maxwell DtN form derivatives")
end program test_planar_maxwell_dtn_form_ad
