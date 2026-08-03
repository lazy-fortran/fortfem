program test_maxwell_torus_curved_potential_operator_ad
    use check, only: check_condition, check_summary
    use fortfem_core, only: generate_torus_surface_mesh
    use fortfem_kinds, only: dp
    use fortfem_maxwell_torus_curved_rwg, only: &
        assemble_maxwell_torus_curved_potential_operators_rwg_3d, &
        assemble_maxwell_torus_curved_potential_operators_rwg_3d_jvp, &
        assemble_maxwell_torus_curved_potential_operators_rwg_3d_vjp
    implicit none

    real(dp), parameter :: major_radius = 2.0_dp, minor_radius = 0.6_dp
    real(dp), parameter :: wave_number = 0.4_dp, step = 2.0e-6_dp
    real(dp), parameter :: wave_number_dot = -0.031_dp
    complex(dp), allocatable :: scalar(:, :), scalar_bar(:, :)
    complex(dp), allocatable :: scalar_dot(:, :), scalar_minus(:, :)
    complex(dp), allocatable :: scalar_plus(:, :), vector(:, :)
    complex(dp), allocatable :: vector_bar(:, :), vector_dot(:, :)
    complex(dp), allocatable :: vector_minus(:, :), vector_plus(:, :)
    real(dp), allocatable :: parameters(:, :), vertices(:, :)
    integer, allocatable :: triangles(:, :)
    real(dp) :: lhs, reassembly_error, rhs, tangent_error, wave_number_bar
    integer :: i, j, status, status_minus, status_plus
    logical :: all_passed

    all_passed = .true.
    call generate_torus_surface_mesh( &
        major_radius, minor_radius, 3, 4, vertices, triangles, parameters)

    call assemble_maxwell_torus_curved_potential_operators_rwg_3d_jvp( &
        vertices, triangles, parameters, major_radius, minor_radius, &
        wave_number, wave_number_dot, 2, 1.0e-5_dp, 1, .false., vector, &
        scalar, vector_dot, scalar_dot, status)
    call assemble_maxwell_torus_curved_potential_operators_rwg_3d( &
        vertices, triangles, parameters, major_radius, minor_radius, &
        wave_number + step*wave_number_dot, 2, 1.0e-5_dp, 1, &
        vector_plus, scalar_plus, status_plus)
    call assemble_maxwell_torus_curved_potential_operators_rwg_3d( &
        vertices, triangles, parameters, major_radius, minor_radius, &
        wave_number - step*wave_number_dot, 2, 1.0e-5_dp, 1, &
        vector_minus, scalar_minus, status_minus)
    reassembly_error = max( &
        maxval(abs(vector - 0.5_dp*(vector_plus + vector_minus))), &
        maxval(abs(scalar - 0.5_dp*(scalar_plus + scalar_minus))))
    tangent_error = max( &
        maxval(abs(vector_dot - (vector_plus - vector_minus)/(2.0_dp*step))), &
        maxval(abs(scalar_dot - (scalar_plus - scalar_minus)/(2.0_dp*step))))
    call record_condition(status == 0, &
        "curved-torus potential-operator block JVP succeeds")
    call record_condition(status_plus == 0 .and. status_minus == 0, &
        "perturbed curved-torus potential blocks reassemble")
    call record_condition(reassembly_error < 2.0e-10_dp, &
        "potential-operator JVP returns the independently assembled primal blocks")
    call record_condition(tangent_error < 2.0e-6_dp, &
        "potential-operator block JVP matches central reassembly")

    allocate( &
        vector_bar(size(vector, 1), size(vector, 2)), &
        scalar_bar(size(scalar, 1), size(scalar, 2)))
    do j = 1, size(vector_bar, 2)
        do i = 1, size(vector_bar, 1)
            vector_bar(i, j) = cmplx( &
                sin(real(i + 2*j, dp)), cos(real(2*i - j, dp)), dp)
            scalar_bar(i, j) = cmplx( &
                cos(real(3*i + j, dp)), sin(real(i - 3*j, dp)), dp)
        end do
    end do
    call assemble_maxwell_torus_curved_potential_operators_rwg_3d_vjp( &
        vertices, triangles, parameters, major_radius, minor_radius, &
        wave_number, 2, 1.0e-5_dp, 1, .false., vector_bar, scalar_bar, &
        wave_number_bar, status)
    lhs = real(sum(conjg(vector_bar)*vector_dot) + &
        sum(conjg(scalar_bar)*scalar_dot), dp)
    rhs = wave_number_bar*wave_number_dot
    call record_condition(status == 0, &
        "curved-torus potential-operator block VJP succeeds")
    call record_condition(abs(lhs - rhs) < &
        2.0e-7_dp*max(1.0_dp, abs(lhs), abs(rhs)), &
        "potential-operator JVP and VJP obey the real complex adjoint identity")

    call check_summary("Differentiable curved-torus potential-operator block")
    if (.not. all_passed) error stop 1

contains

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_maxwell_torus_curved_potential_operator_ad
