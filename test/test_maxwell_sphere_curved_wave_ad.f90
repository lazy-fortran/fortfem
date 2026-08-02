program test_maxwell_sphere_curved_wave_ad
    use check, only: check_condition, check_summary
    use fortfem_api, only: &
        evaluate_maxwell_sphere_curved_far_field_rwg_3d, &
        evaluate_maxwell_sphere_curved_far_field_rwg_3d_jvp, &
        evaluate_maxwell_sphere_curved_far_field_rwg_3d_vjp, &
        generate_sphere_surface_mesh
    use fortfem_kinds, only: dp
    use fortfem_maxwell_rwg_surface, only: build_maxwell_rwg_surface_space
    implicit none

    real(dp), parameter :: radius = 1.0_dp, wave_number = 0.45_dp
    real(dp), parameter :: impedance = 1.7_dp, step = 2.0e-6_dp
    real(dp) :: direction(3), direction_dot(3), direction_bar(3)
    real(dp) :: wave_number_dot, wave_number_bar, impedance_dot, impedance_bar
    real(dp) :: radius_dot, radius_bar, jvp_error, lhs, rhs
    complex(dp), allocatable :: coefficients(:), coefficients_dot(:)
    complex(dp), allocatable :: coefficients_bar(:)
    complex(dp) :: far_field(3), far_field_dot(3), far_field_bar(3)
    complex(dp) :: far_field_minus(3), far_field_plus(3)
    real(dp), allocatable :: vertices(:, :), vertices_dot(:, :), vertices_bar(:, :)
    integer, allocatable :: triangles(:, :), edge_triangles(:, :), edge_vertices(:, :)
    integer :: j, status, status_minus, status_plus
    logical :: all_passed

    all_passed = .true.
    direction = [0.0_dp, 0.0_dp, 1.0_dp]
    ! Keep the direction perturbation tangent to the unit sphere so the
    ! perturbed primal evaluations remain on the validated direction manifold.
    direction_dot = [0.013_dp, -0.009_dp, 0.0_dp]
    wave_number_dot = -0.031_dp
    impedance_dot = 0.023_dp
    radius_dot = -0.017_dp
    call generate_sphere_surface_mesh(radius, 0, vertices, triangles)
    call build_maxwell_rwg_surface_space( &
        vertices, triangles, edge_vertices, edge_triangles, status)
    call record_condition(status == 0, "spherical far-field topology builds")
    allocate( &
        coefficients(size(edge_vertices, 2)), &
        coefficients_dot(size(edge_vertices, 2)), &
        vertices_dot(size(vertices, 1), size(vertices, 2)))
    do j = 1, size(coefficients)
        coefficients(j) = cmplx(sin(real(2*j, dp)), cos(real(3*j, dp)), dp)
        coefficients_dot(j) = cmplx( &
            -0.017_dp*cos(real(j, dp)), 0.011_dp*sin(real(2*j, dp)), dp)
    end do
    do j = 1, size(vertices_dot, 2)
        vertices_dot(:, j) = [ &
            0.013_dp*sin(real(2*j, dp)), &
            -0.009_dp*cos(real(3*j, dp)), &
            0.011_dp*sin(real(5*j, dp))]
    end do

    call evaluate_maxwell_sphere_curved_far_field_rwg_3d_jvp( &
        vertices, triangles, radius, coefficients, direction, wave_number, &
        impedance, 3, vertices_dot, radius_dot, coefficients_dot, direction_dot, &
        wave_number_dot, impedance_dot, far_field, far_field_dot, status)
    call evaluate_maxwell_sphere_curved_far_field_rwg_3d( &
        vertices + step*vertices_dot, triangles, radius + step*radius_dot, &
        coefficients + step*coefficients_dot, direction + step*direction_dot, &
        wave_number + step*wave_number_dot, impedance + step*impedance_dot, 3, &
        far_field_plus, status_plus)
    call evaluate_maxwell_sphere_curved_far_field_rwg_3d( &
        vertices - step*vertices_dot, triangles, radius - step*radius_dot, &
        coefficients - step*coefficients_dot, direction - step*direction_dot, &
        wave_number - step*wave_number_dot, impedance - step*impedance_dot, 3, &
        far_field_minus, status_minus)
    if (status == 0) then
        if (status_plus == 0 .and. status_minus == 0) then
            jvp_error = maxval(abs(far_field_dot - &
                (far_field_plus - far_field_minus)/(2.0_dp*step)))
        else
            jvp_error = huge(1.0_dp)
        end if
    else
        jvp_error = huge(1.0_dp)
    end if
    call record_condition(status == 0, "spherical far-field JVP succeeds")
    call record_condition(status_plus == 0 .and. status_minus == 0, &
        "perturbed spherical far fields evaluate")
    call record_condition(jvp_error < 2.0e-6_dp, &
        "spherical far-field JVP matches reassembly")

    far_field_bar = [ &
        cmplx(0.31_dp, -0.27_dp, dp), cmplx(-0.19_dp, 0.23_dp, dp), &
        cmplx(0.11_dp, 0.17_dp, dp)]
    allocate(vertices_bar(size(vertices, 1), size(vertices, 2)))
    call evaluate_maxwell_sphere_curved_far_field_rwg_3d_vjp( &
        vertices, triangles, radius, coefficients, direction, wave_number, &
        impedance, 3, far_field_bar, vertices_bar, radius_bar, coefficients_bar, &
        direction_bar, wave_number_bar, impedance_bar, status)
    lhs = real(sum(conjg(far_field_bar)*far_field_dot), dp)
    rhs = sum(vertices_bar*vertices_dot) + radius_bar*radius_dot + &
        real(sum(conjg(coefficients_bar)*coefficients_dot), dp) + &
        sum(direction_bar*direction_dot) + wave_number_bar*wave_number_dot + &
        impedance_bar*impedance_dot
    call record_condition(status == 0, "spherical far-field VJP succeeds")
    call record_condition(abs(lhs - rhs) < &
        3.0e-8_dp*max(1.0_dp, abs(lhs), abs(rhs)), &
        "spherical far-field products obey the adjoint identity")

    call check_summary("Differentiable spherical Maxwell far field")
    if (.not. all_passed) error stop 1

contains

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_maxwell_sphere_curved_wave_ad
