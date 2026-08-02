program test_maxwell_sphere_curved_efie_impedance_ad
    use check, only: check_condition, check_summary
    use fortfem_api, only: &
        assemble_maxwell_sphere_curved_efie_rwg_3d, &
        assemble_maxwell_sphere_curved_efie_imaginary_rwg_3d, &
        assemble_maxwell_sphere_efie_propagating_impedance_jvp, &
        assemble_maxwell_sphere_efie_propagating_impedance_vjp, &
        assemble_maxwell_sphere_efie_imaginary_impedance_jvp, &
        assemble_maxwell_sphere_efie_imaginary_impedance_vjp, &
        generate_sphere_surface_mesh
    use fortfem_kinds, only: dp
    implicit none

    real(dp), parameter :: radius = 1.0_dp
    real(dp), parameter :: wave_number = 0.4_dp, decay_rate = 0.6_dp
    real(dp), parameter :: impedance = 1.7_dp, step = 2.0e-6_dp
    real(dp) :: impedance_dot, impedance_bar
    real(dp) :: jvp_error, lhs, rhs
    complex(dp), allocatable :: matrix(:, :), matrix_dot(:, :)
    complex(dp), allocatable :: matrix_minus(:, :), matrix_plus(:, :)
    complex(dp), allocatable :: matrix_bar(:, :)
    real(dp), allocatable :: vertices(:, :)
    integer, allocatable :: triangles(:, :)
    integer :: i, j, status, status_minus, status_plus
    logical :: all_passed

    abstract interface
        subroutine impedance_vjp_interface( &
                vertices, triangles, radius, wave_or_decay, impedance, &
                quadrature_degree, tolerance, max_depth, matrix_bar, &
                impedance_bar, status)
            import dp
            real(dp), intent(in) :: vertices(:, :), radius, wave_or_decay
            real(dp), intent(in) :: impedance, tolerance
            integer, intent(in) :: triangles(:, :), quadrature_degree, max_depth
            complex(dp), intent(in) :: matrix_bar(:, :)
            real(dp), intent(out) :: impedance_bar
            integer, intent(out) :: status
        end subroutine impedance_vjp_interface
    end interface

    all_passed = .true.
    impedance_dot = -0.061_dp
    call generate_sphere_surface_mesh(radius, 0, vertices, triangles)

    call assemble_maxwell_sphere_efie_propagating_impedance_jvp( &
        vertices, triangles, radius, wave_number, impedance, 2, 1.0e-5_dp, 1, &
        impedance_dot, matrix, matrix_dot, status)
    call assemble_maxwell_sphere_curved_efie_rwg_3d( &
        vertices, triangles, radius, wave_number, impedance + &
        step*impedance_dot, 2, 1.0e-5_dp, 1, matrix_plus, status_plus)
    call assemble_maxwell_sphere_curved_efie_rwg_3d( &
        vertices, triangles, radius, wave_number, impedance - &
        step*impedance_dot, 2, 1.0e-5_dp, 1, matrix_minus, status_minus)
    if (status == 0 .and. status_plus == 0 .and. status_minus == 0) then
        jvp_error = maxval(abs(matrix_dot - &
            (matrix_plus - matrix_minus)/(2.0_dp*step)))
    else
        jvp_error = huge(1.0_dp)
    end if
    call record_condition(status == 0, &
        "spherical propagating EFIE impedance JVP succeeds")
    call record_condition(status_plus == 0 .and. status_minus == 0, &
        "perturbed spherical propagating EFIE assembles")
    call record_condition(jvp_error < 2.0e-8_dp, &
        "spherical propagating EFIE impedance JVP matches reassembly")
    call check_impedance_vjp( &
        assemble_maxwell_sphere_efie_propagating_impedance_vjp, matrix, &
        matrix_dot, "spherical propagating EFIE")

    call assemble_maxwell_sphere_efie_imaginary_impedance_jvp( &
        vertices, triangles, radius, decay_rate, impedance, 2, 1.0e-5_dp, 1, &
        impedance_dot, matrix, matrix_dot, status)
    call assemble_maxwell_sphere_curved_efie_imaginary_rwg_3d( &
        vertices, triangles, radius, decay_rate, impedance + &
        step*impedance_dot, 2, 1.0e-5_dp, 1, matrix_plus, status_plus)
    call assemble_maxwell_sphere_curved_efie_imaginary_rwg_3d( &
        vertices, triangles, radius, decay_rate, impedance - &
        step*impedance_dot, 2, 1.0e-5_dp, 1, matrix_minus, status_minus)
    if (status == 0 .and. status_plus == 0 .and. status_minus == 0) then
        jvp_error = maxval(abs(matrix_dot - &
            (matrix_plus - matrix_minus)/(2.0_dp*step)))
    else
        jvp_error = huge(1.0_dp)
    end if
    call record_condition(status == 0, &
        "spherical imaginary EFIE impedance JVP succeeds")
    call record_condition(status_plus == 0 .and. status_minus == 0, &
        "perturbed spherical imaginary EFIE assembles")
    call record_condition(jvp_error < 2.0e-8_dp, &
        "spherical imaginary EFIE impedance JVP matches reassembly")
    call check_impedance_vjp( &
        assemble_maxwell_sphere_efie_imaginary_impedance_vjp, matrix, &
        matrix_dot, "spherical imaginary EFIE")

    call check_summary("Differentiable spherical EFIE impedance")
    if (.not. all_passed) error stop 1

contains

    subroutine check_impedance_vjp(vjp, matrix, matrix_dot, label)
        procedure(impedance_vjp_interface) :: vjp
        complex(dp), intent(in) :: matrix(:, :), matrix_dot(:, :)
        character(*), intent(in) :: label

        real(dp) :: lhs_local, rhs_local, impedance_bar_local
        integer :: local_status

        allocate(matrix_bar(size(matrix, 1), size(matrix, 2)))
        do j = 1, size(matrix_bar, 2)
            do i = 1, size(matrix_bar, 1)
                matrix_bar(i, j) = cmplx( &
                    sin(real(i + 2*j, dp)), cos(real(2*i - j, dp)), dp)
            end do
        end do
        call vjp(vertices, triangles, radius, merge(wave_number, decay_rate, &
            index(label, "imaginary") == 0), impedance, 2, 1.0e-5_dp, 1, &
            matrix_bar, impedance_bar_local, local_status)
        lhs_local = real(sum(conjg(matrix_bar)*matrix_dot), dp)
        rhs_local = impedance_bar_local*impedance_dot
        call record_condition(local_status == 0, trim(label)// &
            " impedance VJP succeeds")
        call record_condition(abs(lhs_local - rhs_local) < &
            2.0e-9_dp*max(1.0_dp, abs(lhs_local), abs(rhs_local)), trim(label)// &
            " impedance products obey the adjoint identity")
        deallocate(matrix_bar)
    end subroutine check_impedance_vjp

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_maxwell_sphere_curved_efie_impedance_ad
