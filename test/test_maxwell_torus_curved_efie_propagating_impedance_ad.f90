program test_maxwell_torus_curved_efie_propagating_impedance_ad
    use check, only: check_condition, check_summary
    use fortfem_api, only: &
        assemble_maxwell_torus_curved_efie_rwg_3d, &
        assemble_maxwell_torus_efie_propagating_impedance_jvp, &
        assemble_maxwell_torus_efie_propagating_impedance_vjp, &
        generate_torus_surface_mesh
    use fortfem_kinds, only: dp
    implicit none

    real(dp), parameter :: major_radius = 2.0_dp, minor_radius = 0.6_dp
    real(dp), parameter :: wave_number = 0.4_dp, impedance = 1.7_dp
    real(dp), parameter :: step = 2.0e-6_dp
    real(dp) :: impedance_dot, impedance_bar
    real(dp) :: jvp_error, lhs, rhs
    complex(dp), allocatable :: matrix(:, :), matrix_dot(:, :)
    complex(dp), allocatable :: matrix_minus(:, :), matrix_plus(:, :)
    complex(dp), allocatable :: matrix_bar(:, :)
    real(dp), allocatable :: vertices(:, :), parameters(:, :)
    integer, allocatable :: triangles(:, :)
    integer :: i, j, status, status_minus, status_plus
    logical :: all_passed

    all_passed = .true.
    impedance_dot = -0.061_dp
    call generate_torus_surface_mesh( &
        major_radius, minor_radius, 3, 4, vertices, triangles, parameters)

    call assemble_maxwell_torus_efie_propagating_impedance_jvp( &
        vertices, triangles, parameters, major_radius, minor_radius, &
        wave_number, impedance, 2, 1.0e-5_dp, 1, impedance_dot, matrix, &
        matrix_dot, status)
    call assemble_maxwell_torus_curved_efie_rwg_3d( &
        vertices, triangles, parameters, major_radius, minor_radius, &
        wave_number, impedance + step*impedance_dot, 2, 1.0e-5_dp, 1, &
        matrix_plus, status_plus)
    call assemble_maxwell_torus_curved_efie_rwg_3d( &
        vertices, triangles, parameters, major_radius, minor_radius, &
        wave_number, impedance - step*impedance_dot, 2, 1.0e-5_dp, 1, &
        matrix_minus, status_minus)
    if (status == 0 .and. status_plus == 0 .and. status_minus == 0) then
        jvp_error = maxval(abs(matrix_dot - &
            (matrix_plus - matrix_minus)/(2.0_dp*step)))
    else
        jvp_error = huge(1.0_dp)
    end if
    call record_condition(status == 0, &
        "toroidal propagating EFIE impedance JVP succeeds")
    call record_condition(status_plus == 0 .and. status_minus == 0, &
        "perturbed toroidal propagating EFIE assembles")
    call record_condition(jvp_error < 2.0e-8_dp, &
        "toroidal propagating EFIE impedance JVP matches reassembly")

    allocate(matrix_bar(size(matrix, 1), size(matrix, 2)))
    do j = 1, size(matrix_bar, 2)
        do i = 1, size(matrix_bar, 1)
            matrix_bar(i, j) = cmplx( &
                sin(real(i + 2*j, dp)), cos(real(2*i - j, dp)), dp)
        end do
    end do
    call assemble_maxwell_torus_efie_propagating_impedance_vjp( &
        vertices, triangles, parameters, major_radius, minor_radius, &
        wave_number, impedance, 2, 1.0e-5_dp, 1, matrix_bar, impedance_bar, &
        status)
    lhs = real(sum(conjg(matrix_bar)*matrix_dot), dp)
    rhs = impedance_bar*impedance_dot
    call record_condition(status == 0, &
        "toroidal propagating EFIE impedance VJP succeeds")
    call record_condition(abs(lhs - rhs) < &
        2.0e-9_dp*max(1.0_dp, abs(lhs), abs(rhs)), &
        "toroidal propagating EFIE impedance products obey the adjoint identity")

    call check_summary("Differentiable toroidal propagating EFIE impedance")
    if (.not. all_passed) error stop 1

contains

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_maxwell_torus_curved_efie_propagating_impedance_ad
