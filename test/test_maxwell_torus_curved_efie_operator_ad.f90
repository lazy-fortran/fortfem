program test_maxwell_torus_curved_efie_operator_ad
    use check, only: check_condition, check_summary
    use fortfem_core, only: generate_torus_surface_mesh
    use fortfem_kinds, only: dp
    use fortfem_maxwell_torus_curved_rwg, only: &
        assemble_maxwell_torus_curved_efie_rwg_3d, &
        assemble_maxwell_torus_curved_efie_rwg_3d_jvp, &
        assemble_maxwell_torus_curved_efie_rwg_3d_vjp
    implicit none

    real(dp), parameter :: impedance = 1.7_dp, impedance_dot = 0.053_dp
    real(dp), parameter :: major_radius = 2.0_dp, minor_radius = 0.6_dp
    real(dp), parameter :: step = 2.0e-6_dp, wave_number = 0.4_dp
    real(dp), parameter :: wave_number_dot = -0.031_dp
    complex(dp), allocatable :: matrix(:, :), matrix_bar(:, :)
    complex(dp), allocatable :: matrix_dot(:, :), matrix_minus(:, :)
    complex(dp), allocatable :: matrix_plus(:, :)
    real(dp), allocatable :: parameters(:, :), vertices(:, :)
    integer, allocatable :: triangles(:, :)
    real(dp) :: impedance_bar, lhs, primal_error, rhs, tangent_error
    real(dp) :: wave_number_bar
    integer :: i, j, status, status_minus, status_plus
    logical :: all_passed

    all_passed = .true.
    call generate_torus_surface_mesh( &
        major_radius, minor_radius, 3, 4, vertices, triangles, parameters)

    call assemble_maxwell_torus_curved_efie_rwg_3d_jvp( &
        vertices, triangles, parameters, major_radius, minor_radius, &
        wave_number, impedance, quadrature_degree=2, tolerance=1.0e-5_dp, &
        max_depth=1, wave_number_dot=wave_number_dot, &
        impedance_dot=impedance_dot, matrix=matrix, matrix_dot=matrix_dot, &
        status=status)
    call assemble_maxwell_torus_curved_efie_rwg_3d( &
        vertices, triangles, parameters, major_radius, minor_radius, &
        wave_number + step*wave_number_dot, &
        impedance + step*impedance_dot, 2, 1.0e-5_dp, 1, matrix_plus, &
        status_plus)
    call assemble_maxwell_torus_curved_efie_rwg_3d( &
        vertices, triangles, parameters, major_radius, minor_radius, &
        wave_number - step*wave_number_dot, &
        impedance - step*impedance_dot, 2, 1.0e-5_dp, 1, matrix_minus, &
        status_minus)
    if (status == 0 .and. status_plus == 0 .and. status_minus == 0) then
        primal_error = maxval(abs( &
            matrix - 0.5_dp*(matrix_plus + matrix_minus)))
        tangent_error = maxval(abs( &
            matrix_dot - (matrix_plus - matrix_minus)/(2.0_dp*step)))
    else
        primal_error = huge(1.0_dp)
        tangent_error = huge(1.0_dp)
    end if
    call record_condition(status == 0, &
        "assembled curved-torus EFIE JVP succeeds")
    call record_condition(status_plus == 0 .and. status_minus == 0, &
        "perturbed assembled curved-torus EFIE operators reassemble")
    call record_condition(primal_error < 2.0e-10_dp, &
        "assembled EFIE JVP returns the independently assembled primal")
    call record_condition(tangent_error < 2.0e-6_dp, &
        "assembled EFIE JVP matches central reassembly")

    allocate(matrix_bar(size(matrix, 1), size(matrix, 2)))
    do j = 1, size(matrix_bar, 2)
        do i = 1, size(matrix_bar, 1)
            matrix_bar(i, j) = cmplx( &
                sin(real(i + 2*j, dp)), cos(real(2*i - j, dp)), dp)
        end do
    end do
    call assemble_maxwell_torus_curved_efie_rwg_3d_vjp( &
        vertices, triangles, parameters, major_radius, minor_radius, &
        wave_number, impedance, 2, 1.0e-5_dp, 1, matrix_bar, &
        wave_number_bar, impedance_bar, status)
    lhs = real(sum(conjg(matrix_bar)*matrix_dot), dp)
    rhs = wave_number_bar*wave_number_dot + impedance_bar*impedance_dot
    call record_condition(status == 0, &
        "assembled curved-torus EFIE VJP succeeds")
    call record_condition(abs(lhs - rhs) < &
        2.0e-7_dp*max(1.0_dp, abs(lhs), abs(rhs)), &
        "assembled EFIE JVP and VJP obey the real complex adjoint identity")

    call check_summary("Differentiable assembled curved-torus EFIE operator")
    if (.not. all_passed) error stop 1

contains

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_maxwell_torus_curved_efie_operator_ad
