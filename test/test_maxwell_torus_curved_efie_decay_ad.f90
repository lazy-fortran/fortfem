program test_maxwell_torus_curved_efie_decay_ad
    use check, only: check_condition, check_summary
    use fortfem_api, only: &
        assemble_maxwell_torus_curved_efie_imaginary_rwg_3d, &
        assemble_maxwell_torus_efie_imaginary_decay_jvp, &
        assemble_maxwell_torus_efie_imaginary_decay_vjp, &
        generate_torus_surface_mesh
    use fortfem_kinds, only: dp
    implicit none

    real(dp), parameter :: major_radius = 2.0_dp, minor_radius = 0.6_dp
    real(dp), parameter :: decay_rate = 0.6_dp, impedance = 1.7_dp
    real(dp), parameter :: step = 2.0e-6_dp
    real(dp) :: decay_rate_dot, decay_rate_bar
    real(dp) :: jvp_error, lhs, rhs
    complex(dp), allocatable :: matrix(:, :), matrix_dot(:, :)
    complex(dp), allocatable :: matrix_minus(:, :), matrix_plus(:, :)
    complex(dp), allocatable :: matrix_bar(:, :)
    real(dp), allocatable :: vertices(:, :), parameters(:, :)
    integer, allocatable :: triangles(:, :)
    integer :: i, j, status, status_minus, status_plus
    logical :: all_passed

    all_passed = .true.
    decay_rate_dot = -0.047_dp
    call generate_torus_surface_mesh( &
        major_radius, minor_radius, 3, 4, vertices, triangles, parameters)

    call assemble_maxwell_torus_efie_imaginary_decay_jvp( &
        vertices, triangles, parameters, major_radius, minor_radius, &
        decay_rate, impedance, 2, 1.0e-5_dp, 1, decay_rate_dot, matrix, &
        matrix_dot, status)
    call assemble_maxwell_torus_curved_efie_imaginary_rwg_3d( &
        vertices, triangles, parameters, major_radius, minor_radius, &
        decay_rate + step*decay_rate_dot, impedance, 2, 1.0e-5_dp, 1, &
        matrix_plus, status_plus)
    call assemble_maxwell_torus_curved_efie_imaginary_rwg_3d( &
        vertices, triangles, parameters, major_radius, minor_radius, &
        decay_rate - step*decay_rate_dot, impedance, 2, 1.0e-5_dp, 1, &
        matrix_minus, status_minus)
    if (status == 0 .and. status_plus == 0 .and. status_minus == 0) then
        jvp_error = maxval(abs(matrix_dot - &
            (matrix_plus - matrix_minus)/(2.0_dp*step)))
    else
        jvp_error = huge(1.0_dp)
    end if
    call record_condition(status == 0, &
        "toroidal imaginary EFIE decay JVP succeeds")
    call record_condition(status_plus == 0 .and. status_minus == 0, &
        "perturbed toroidal imaginary EFIE decay assembles")
    call record_condition(jvp_error < 2.0e-6_dp, &
        "toroidal imaginary EFIE decay JVP matches reassembly")

    allocate(matrix_bar(size(matrix, 1), size(matrix, 2)))
    do j = 1, size(matrix_bar, 2)
        do i = 1, size(matrix_bar, 1)
            matrix_bar(i, j) = cmplx( &
                sin(real(i + 2*j, dp)), cos(real(2*i - j, dp)), dp)
        end do
    end do
    call assemble_maxwell_torus_efie_imaginary_decay_vjp( &
        vertices, triangles, parameters, major_radius, minor_radius, &
        decay_rate, impedance, 2, 1.0e-5_dp, 1, matrix_bar, decay_rate_bar, &
        status)
    lhs = real(sum(conjg(matrix_bar)*matrix_dot), dp)
    rhs = decay_rate_bar*decay_rate_dot
    call record_condition(status == 0, &
        "toroidal imaginary EFIE decay VJP succeeds")
    call record_condition(abs(lhs - rhs) < &
        2.0e-7_dp*max(1.0_dp, abs(lhs), abs(rhs)), &
        "toroidal imaginary EFIE decay products obey the adjoint identity")

    call check_summary("Differentiable toroidal imaginary decay EFIE")
    if (.not. all_passed) error stop 1

contains

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_maxwell_torus_curved_efie_decay_ad
