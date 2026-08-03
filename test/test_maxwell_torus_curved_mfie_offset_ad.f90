program test_maxwell_torus_curved_mfie_offset_ad
    use check, only: check_condition, check_summary
    use fortfem_boundary, only: &
        assemble_maxwell_torus_curved_mfie_offset_trace_rwg_rbc_3d, &
        assemble_maxwell_torus_mfie_offset_jvp, &
        assemble_maxwell_torus_mfie_offset_vjp
    use fortfem_core, only: generate_torus_surface_mesh
    use fortfem_kinds, only: dp
    implicit none

    real(dp), parameter :: major_radius = 2.0_dp, minor_radius = 0.6_dp
    real(dp), parameter :: wave_number = 0.4_dp, relative_offset = 0.12_dp
    real(dp), parameter :: relative_offset_dot = -0.037_dp, step = 2.0e-6_dp
    real(dp) :: wave_number_dot, wave_number_bar
    real(dp) :: relative_offset_bar, lhs, rhs, jvp_error
    complex(dp), allocatable :: matrix(:, :), matrix_dot(:, :)
    complex(dp), allocatable :: matrix_minus(:, :), matrix_plus(:, :)
    complex(dp), allocatable :: matrix_bar(:, :)
    real(dp), allocatable :: vertices(:, :), parameters(:, :)
    integer, allocatable :: triangles(:, :)
    integer :: i, j, status, status_minus, status_plus
    logical :: all_passed

    all_passed = .true.
    wave_number_dot = 0.023_dp
    call generate_torus_surface_mesh( &
        major_radius, minor_radius, 3, 3, vertices, triangles, parameters)

    call assemble_maxwell_torus_mfie_offset_jvp( &
        vertices, triangles, parameters, major_radius, minor_radius, &
        wave_number, 2, relative_offset, relative_offset_dot, matrix, &
        matrix_dot, status, wave_number_dot=wave_number_dot)
    call assemble_maxwell_torus_curved_mfie_offset_trace_rwg_rbc_3d( &
        vertices, triangles, parameters, major_radius, minor_radius, &
        wave_number + step*wave_number_dot, 2, &
        relative_offset + step*relative_offset_dot, &
        matrix_plus, status_plus)
    call assemble_maxwell_torus_curved_mfie_offset_trace_rwg_rbc_3d( &
        vertices, triangles, parameters, major_radius, minor_radius, &
        wave_number - step*wave_number_dot, 2, &
        relative_offset - step*relative_offset_dot, &
        matrix_minus, status_minus)
    jvp_error = maxval(abs(matrix_dot - &
        (matrix_plus - matrix_minus)/(2.0_dp*step)))
    call record_condition(status == 0, &
        "toroidal MFIE offset-trace JVP succeeds")
    call record_condition(status_plus == 0 .and. status_minus == 0, &
        "perturbed toroidal MFIE offset traces assemble")
    call record_condition(jvp_error < 5.0e-8_dp, &
        "toroidal MFIE offset-trace JVP matches reassembly")

    allocate(matrix_bar(size(matrix, 1), size(matrix, 2)))
    do j = 1, size(matrix_bar, 2)
        do i = 1, size(matrix_bar, 1)
            matrix_bar(i, j) = cmplx( &
                sin(real(i + 2*j, dp)), cos(real(2*i - j, dp)), dp)
        end do
    end do
    call assemble_maxwell_torus_mfie_offset_vjp( &
        vertices, triangles, parameters, major_radius, minor_radius, &
        wave_number, 2, relative_offset, matrix_bar, relative_offset_bar, &
        status, wave_number_bar=wave_number_bar)
    lhs = real(sum(conjg(matrix_bar)*matrix_dot), dp)
    rhs = relative_offset_bar*relative_offset_dot + wave_number_bar*wave_number_dot
    call record_condition(status == 0, &
        "toroidal MFIE offset-trace VJP succeeds")
    call record_condition(abs(lhs - rhs) < &
        5.0e-8_dp*max(1.0_dp, abs(lhs), abs(rhs)), &
        "toroidal MFIE offset-trace products obey the adjoint identity")

    call check_summary("Differentiable toroidal MFIE offset trace")
    if (.not. all_passed) error stop 1

contains

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_maxwell_torus_curved_mfie_offset_ad
