program test_planar_maxwell_dtn_system_ad
    use check, only: check_condition, check_summary
    use fortfem_api, only: solve_planar_maxwell_dtn_system, &
        solve_planar_maxwell_dtn_system_jvp, &
        solve_planar_maxwell_dtn_system_vjp
    use fortfem_kinds, only: dp
    implicit none

    integer, parameter :: nx = 3, ny = 3, trace_count = 2*nx*ny
    integer, parameter :: dof_count = 4
    real(dp), parameter :: step = 2.0e-7_dp
    complex(dp) :: volume_matrix(dof_count, dof_count)
    complex(dp) :: volume_matrix_dot(dof_count, dof_count)
    complex(dp) :: volume_matrix_bar(dof_count, dof_count)
    complex(dp) :: sampling(trace_count, dof_count)
    complex(dp) :: sampling_dot(trace_count, dof_count)
    complex(dp) :: sampling_bar(trace_count, dof_count)
    complex(dp) :: load(dof_count), load_dot(dof_count), load_bar(dof_count)
    complex(dp) :: state(dof_count), state_dot(dof_count), state_bar(dof_count)
    complex(dp) :: plus(dof_count), minus(dof_count)
    real(dp) :: wave_number, length_x, length_y
    real(dp) :: wave_number_dot, length_x_dot, length_y_dot
    real(dp) :: wave_number_bar, length_x_bar, length_y_bar
    real(dp) :: lhs, rhs, relative_error
    integer :: column, row, status, status_plus, status_minus

    volume_matrix = cmplx(0.0_dp, 0.0_dp, dp)
    volume_matrix_dot = cmplx(0.0_dp, 0.0_dp, dp)
    do column = 1, dof_count
        do row = 1, dof_count
            if (row == column) then
                volume_matrix(row, column) = cmplx( &
                    8.0_dp + 0.4_dp*row, 0.2_dp*row, dp)
            else
                volume_matrix(row, column) = cmplx( &
                    0.03_dp*(row + column), -0.01_dp*(row - column), dp)
            end if
            volume_matrix_dot(row, column) = cmplx( &
                0.004_dp*row - 0.003_dp*column, &
                -0.002_dp*row + 0.001_dp*column, dp)
        end do
    end do
    do column = 1, dof_count
        do row = 1, trace_count
            sampling(row, column) = cmplx( &
                0.01_dp*row - 0.015_dp*column, &
                -0.004_dp*row + 0.007_dp*column, dp)
            sampling_dot(row, column) = cmplx( &
                -0.003_dp*row + 0.005_dp*column, &
                0.002_dp*row + 0.001_dp*column, dp)
        end do
        load(column) = cmplx(0.3_dp*column, -0.1_dp*column, dp)
        load_dot(column) = cmplx(-0.04_dp*column, 0.02_dp*column, dp)
        state_bar(column) = cmplx(0.05_dp*column, -0.03_dp*column, dp)
    end do
    wave_number = 2.3_dp
    length_x = 4.7_dp
    length_y = 5.2_dp
    wave_number_dot = 0.17_dp
    length_x_dot = -0.11_dp
    length_y_dot = 0.08_dp

    call solve_planar_maxwell_dtn_system( &
        volume_matrix, sampling, nx, ny, wave_number, length_x, length_y, &
        load, state, status)
    call solve_planar_maxwell_dtn_system_jvp( &
        volume_matrix, sampling, nx, ny, wave_number, length_x, length_y, &
        load, volume_matrix_dot, sampling_dot, wave_number_dot, length_x_dot, &
        length_y_dot, load_dot, state_dot, status)
    call solve_planar_maxwell_dtn_system( &
        volume_matrix + step*volume_matrix_dot, &
        sampling + step*sampling_dot, nx, ny, &
        wave_number + step*wave_number_dot, length_x + step*length_x_dot, &
        length_y + step*length_y_dot, load + step*load_dot, plus, status_plus)
    call solve_planar_maxwell_dtn_system( &
        volume_matrix - step*volume_matrix_dot, &
        sampling - step*sampling_dot, nx, ny, &
        wave_number - step*wave_number_dot, length_x - step*length_x_dot, &
        length_y - step*length_y_dot, load - step*load_dot, minus, status_minus)
    relative_error = maxval(abs( &
        state_dot - (plus - minus)/(2.0_dp*step)))/ &
        max(1.0_dp, maxval(abs(state_dot)))
    call check_condition( &
        status == 0 .and. status_plus == 0 .and. status_minus == 0, &
        "Maxwell FEM-DtN state JVP accepts valid inputs")
    call check_condition( &
        relative_error < 3.0e-9_dp, &
        "Maxwell FEM-DtN state JVP matches reassemble-refactor difference")

    call solve_planar_maxwell_dtn_system_vjp( &
        volume_matrix, sampling, nx, ny, wave_number, length_x, length_y, &
        load, state, state_bar, volume_matrix_bar, sampling_bar, &
        wave_number_bar, length_x_bar, length_y_bar, load_bar, status)
    lhs = real(sum(conjg(state_bar)*state_dot), dp)
    rhs = real(sum(conjg(volume_matrix_bar)*volume_matrix_dot) + &
        sum(conjg(sampling_bar)*sampling_dot) + &
        sum(conjg(load_bar)*load_dot), dp) + &
        wave_number_bar*wave_number_dot + length_x_bar*length_x_dot + &
        length_y_bar*length_y_dot
    call check_condition(status == 0, "Maxwell FEM-DtN state VJP succeeds")
    call check_condition( &
        abs(lhs - rhs) < 3.0e-11_dp*max(1.0_dp, abs(lhs), abs(rhs)), &
        "Maxwell FEM-DtN state products satisfy the complete adjoint identity")

    call check_summary("Planar Maxwell FEM-DtN state derivatives")
end program test_planar_maxwell_dtn_system_ad
