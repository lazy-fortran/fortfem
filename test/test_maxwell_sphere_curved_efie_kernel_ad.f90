program test_maxwell_sphere_curved_efie_kernel_ad
    use check, only: check_condition, check_summary
    use fortfem_api, only: &
        assemble_maxwell_sphere_curved_efie_rwg_3d, &
        assemble_maxwell_sphere_curved_efie_imaginary_rwg_3d, &
        assemble_maxwell_sphere_efie_wave_number_jvp, &
        assemble_maxwell_sphere_efie_wave_number_vjp, &
        assemble_maxwell_sphere_efie_imaginary_decay_jvp, &
        assemble_maxwell_sphere_efie_imaginary_decay_vjp, &
        generate_sphere_surface_mesh
    use fortfem_kinds, only: dp
    implicit none

    real(dp), parameter :: radius = 1.0_dp
    real(dp), parameter :: wave_number = 0.4_dp, decay_rate = 0.6_dp
    real(dp), parameter :: impedance = 1.7_dp, step = 2.0e-6_dp
    real(dp) :: wave_number_dot, wave_number_bar, decay_rate_dot, decay_rate_bar
    real(dp) :: jvp_error, lhs, rhs
    complex(dp), allocatable :: matrix(:, :), matrix_dot(:, :)
    complex(dp), allocatable :: matrix_minus(:, :), matrix_plus(:, :)
    complex(dp), allocatable :: matrix_bar(:, :)
    real(dp), allocatable :: vertices(:, :)
    integer, allocatable :: triangles(:, :)
    integer :: i, j, status, status_minus, status_plus
    logical :: all_passed

    all_passed = .true.
    wave_number_dot = -0.031_dp
    decay_rate_dot = -0.027_dp
    call generate_sphere_surface_mesh(radius, 0, vertices, triangles)

    call assemble_maxwell_sphere_efie_wave_number_jvp( &
        vertices, triangles, radius, wave_number, impedance, 2, 1.0e-5_dp, 1, &
        wave_number_dot, matrix, matrix_dot, status)
    call assemble_maxwell_sphere_curved_efie_rwg_3d( &
        vertices, triangles, radius, wave_number + step*wave_number_dot, &
        impedance, 2, 1.0e-5_dp, 1, matrix_plus, status_plus)
    call assemble_maxwell_sphere_curved_efie_rwg_3d( &
        vertices, triangles, radius, wave_number - step*wave_number_dot, &
        impedance, 2, 1.0e-5_dp, 1, matrix_minus, status_minus)
    if (status == 0) then
        if (status_plus == 0 .and. status_minus == 0) then
            jvp_error = maxval(abs(matrix_dot - &
                (matrix_plus - matrix_minus)/(2.0_dp*step)))
        else
            jvp_error = huge(1.0_dp)
        end if
    else
        jvp_error = huge(1.0_dp)
    end if
    call record_condition(status == 0, &
        "spherical propagating EFIE wave-number JVP succeeds")
    call record_condition(status_plus == 0 .and. status_minus == 0, &
        "perturbed spherical propagating EFIE assembles")
    call record_condition(jvp_error < 2.0e-6_dp, &
        "spherical propagating EFIE wave-number JVP matches reassembly")

    call fill_matrix_bar(matrix_bar, matrix)
    call assemble_maxwell_sphere_efie_wave_number_vjp( &
        vertices, triangles, radius, wave_number, impedance, 2, 1.0e-5_dp, 1, &
        matrix_bar, wave_number_bar, status)
    lhs = real(sum(conjg(matrix_bar)*matrix_dot), dp)
    rhs = wave_number_bar*wave_number_dot
    call record_condition(status == 0, &
        "spherical propagating EFIE wave-number VJP succeeds")
    call record_condition(abs(lhs - rhs) < &
        2.0e-7_dp*max(1.0_dp, abs(lhs), abs(rhs)), &
        "spherical propagating EFIE wave-number products obey the adjoint identity")
    deallocate(matrix_bar)

    call assemble_maxwell_sphere_efie_imaginary_decay_jvp( &
        vertices, triangles, radius, decay_rate, impedance, 2, 1.0e-5_dp, 1, &
        decay_rate_dot, matrix, matrix_dot, status)
    call assemble_maxwell_sphere_curved_efie_imaginary_rwg_3d( &
        vertices, triangles, radius, decay_rate + step*decay_rate_dot, &
        impedance, 2, 1.0e-5_dp, 1, matrix_plus, status_plus)
    call assemble_maxwell_sphere_curved_efie_imaginary_rwg_3d( &
        vertices, triangles, radius, decay_rate - step*decay_rate_dot, &
        impedance, 2, 1.0e-5_dp, 1, matrix_minus, status_minus)
    if (status == 0) then
        if (status_plus == 0 .and. status_minus == 0) then
            jvp_error = maxval(abs(matrix_dot - &
                (matrix_plus - matrix_minus)/(2.0_dp*step)))
        else
            jvp_error = huge(1.0_dp)
        end if
    else
        jvp_error = huge(1.0_dp)
    end if
    call record_condition(status == 0, &
        "spherical imaginary EFIE decay-rate JVP succeeds")
    call record_condition(status_plus == 0 .and. status_minus == 0, &
        "perturbed spherical imaginary EFIE assembles")
    call record_condition(jvp_error < 2.0e-6_dp, &
        "spherical imaginary EFIE decay-rate JVP matches reassembly")

    call fill_matrix_bar(matrix_bar, matrix)
    call assemble_maxwell_sphere_efie_imaginary_decay_vjp( &
        vertices, triangles, radius, decay_rate, impedance, 2, 1.0e-5_dp, 1, &
        matrix_bar, decay_rate_bar, status)
    lhs = real(sum(conjg(matrix_bar)*matrix_dot), dp)
    rhs = decay_rate_bar*decay_rate_dot
    call record_condition(status == 0, &
        "spherical imaginary EFIE decay-rate VJP succeeds")
    call record_condition(abs(lhs - rhs) < &
        2.0e-7_dp*max(1.0_dp, abs(lhs), abs(rhs)), &
        "spherical imaginary EFIE decay-rate products obey the adjoint identity")

    call check_summary("Differentiable spherical EFIE kernel parameters")
    if (.not. all_passed) error stop 1

contains

    subroutine fill_matrix_bar(matrix_bar_local, reference_matrix)
        complex(dp), allocatable, intent(out) :: matrix_bar_local(:, :)
        complex(dp), intent(in) :: reference_matrix(:, :)

        allocate(matrix_bar_local(size(reference_matrix, 1), &
            size(reference_matrix, 2)))
        do j = 1, size(matrix_bar_local, 2)
            do i = 1, size(matrix_bar_local, 1)
                matrix_bar_local(i, j) = cmplx( &
                    sin(real(i + 2*j, dp)), cos(real(2*i - j, dp)), dp)
            end do
        end do
    end subroutine fill_matrix_bar

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_maxwell_sphere_curved_efie_kernel_ad
