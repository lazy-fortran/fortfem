program test_maxwell_magnetic_rwg_3d_ad
    use check, only: check_condition, check_summary
    use fortfem_boundary, only: &
        evaluate_maxwell_magnetic_field_rwg_3d, &
        evaluate_maxwell_magnetic_field_rwg_3d_jvp, &
        evaluate_maxwell_magnetic_field_rwg_3d_vjp
    use fortfem_feec, only: build_maxwell_rwg_surface_space
    use fortfem_kinds, only: dp
    implicit none

    real(dp), parameter :: observation(3) = [1.7_dp, -0.8_dp, 1.3_dp]
    real(dp), parameter :: wave_number = 1.4_dp
    real(dp), parameter :: step = 1.0e-6_dp
    real(dp) :: vertices(3, 4), vertices_dot(3, 4), vertices_bar(3, 4)
    integer :: triangles(3, 4)
    complex(dp) :: coefficients(6), coefficients_dot(6)
    complex(dp), allocatable :: coefficients_bar(:)
    complex(dp) :: field_dot(3), field_bar(3), field_plus(3), field_minus(3)
    real(dp) :: observation_dot(3), observation_bar(3)
    real(dp) :: wave_number_dot, wave_number_bar
    real(dp) :: derivative_error, adjoint_error, lhs, rhs
    integer :: basis, status
    logical :: all_passed

    all_passed = .true.
    vertices(:, 1) = [0.0_dp, 0.0_dp, 0.0_dp]
    vertices(:, 2) = [1.0_dp, 0.0_dp, 0.0_dp]
    vertices(:, 3) = [0.0_dp, 1.0_dp, 0.0_dp]
    vertices(:, 4) = [0.0_dp, 0.0_dp, 1.0_dp]
    triangles(:, 1) = [1, 3, 2]
    triangles(:, 2) = [1, 2, 4]
    triangles(:, 3) = [1, 4, 3]
    triangles(:, 4) = [2, 3, 4]
    do basis = 1, 4
        vertices_dot(:, basis) = [ &
            0.002_dp*basis, -0.0013_dp*basis, 0.0008_dp*basis]
    end do
    do basis = 1, 6
        coefficients(basis) = cmplx( &
            cos(0.31_dp*basis), sin(0.17_dp*basis), dp)
        coefficients_dot(basis) = cmplx( &
            0.11_dp*sin(0.13_dp*basis), -0.07_dp*cos(0.19_dp*basis), dp)
    end do
    observation_dot = [0.004_dp, -0.003_dp, 0.002_dp]
    wave_number_dot = 0.017_dp

    call evaluate_maxwell_magnetic_field_rwg_3d_jvp( &
        vertices, triangles, coefficients, observation, wave_number, 6, &
        vertices_dot, coefficients_dot, observation_dot, wave_number_dot, &
        field_dot, status)
    call record_condition(status == 0, &
        "flat RWG magnetic derivative accepts a valid direction")
    call evaluate_maxwell_magnetic_field_rwg_3d( &
        vertices + step*vertices_dot, triangles, coefficients + step*coefficients_dot, &
        observation + step*observation_dot, wave_number + step*wave_number_dot, 6, &
        field_plus, status)
    call evaluate_maxwell_magnetic_field_rwg_3d( &
        vertices - step*vertices_dot, triangles, coefficients - step*coefficients_dot, &
        observation - step*observation_dot, wave_number - step*wave_number_dot, 6, &
        field_minus, status)
    derivative_error = maxval(abs(field_dot - &
        (field_plus - field_minus)/(2.0_dp*step)))
    call record_condition(status == 0 .and. derivative_error < 2.0e-7_dp, &
        "flat RWG magnetic JVP matches central reassembly")

    field_bar = [ &
        cmplx(0.7_dp, -0.2_dp, dp), cmplx(-0.3_dp, 0.4_dp, dp), &
        cmplx(0.1_dp, 0.6_dp, dp)]
    call evaluate_maxwell_magnetic_field_rwg_3d_vjp( &
        vertices, triangles, coefficients, observation, wave_number, 6, field_bar, &
        vertices_bar, coefficients_bar, observation_bar, wave_number_bar, status)
    lhs = real(sum(conjg(field_bar)*field_dot), dp)
    rhs = sum(vertices_bar*vertices_dot) + &
        real(sum(conjg(coefficients_bar)*coefficients_dot), dp) + &
        dot_product(observation_bar, observation_dot) + wave_number_bar*wave_number_dot
    adjoint_error = abs(lhs - rhs)
    call record_condition(status == 0 .and. adjoint_error < &
        5.0e-8_dp*max(1.0_dp, abs(lhs), abs(rhs)), &
        "flat RWG magnetic VJP satisfies the real complex adjoint")

    call check_summary("flat RWG magnetic derivative")
    if (.not. all_passed) error stop 1

contains

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_maxwell_magnetic_rwg_3d_ad
