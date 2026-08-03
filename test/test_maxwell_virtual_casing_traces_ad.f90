program test_maxwell_virtual_casing_traces_ad
    use check, only: check_condition, check_summary
    use fortfem_boundary, only: &
        evaluate_maxwell_virtual_casing_rwg_3d, &
        evaluate_maxwell_virtual_casing_rwg_3d_jvp, &
        evaluate_maxwell_virtual_casing_rwg_3d_vjp
    use fortfem_kinds, only: dp
    implicit none

    real(dp), parameter :: wave_number = 1.4_dp
    real(dp), parameter :: step = 1.0e-6_dp
    real(dp) :: vertices(3, 4), vertices_dot(3, 4), vertices_bar(3, 4)
    real(dp) :: observations(3, 2), observations_dot(3, 2)
    real(dp) :: observations_bar(3, 2)
    real(dp) :: normals(3, 2), normals_dot(3, 2), normals_bar(3, 2)
    integer :: triangles(3, 4)
    complex(dp) :: coefficients(6), coefficients_dot(6)
    complex(dp), allocatable :: coefficients_bar(:)
    complex(dp) :: fields(3, 2), fields_dot(3, 2), fields_bar(3, 2)
    complex(dp) :: normal_fields(2), normal_fields_dot(2), normal_fields_bar(2)
    complex(dp) :: tangential_fields(3, 2), tangential_fields_dot(3, 2)
    complex(dp) :: tangential_fields_bar(3, 2)
    complex(dp) :: fields_plus(3, 2), fields_minus(3, 2)
    complex(dp) :: normal_fields_plus(2), normal_fields_minus(2)
    complex(dp) :: tangential_fields_plus(3, 2), tangential_fields_minus(3, 2)
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
    observations(:, 1) = [1.7_dp, -0.8_dp, 1.3_dp]
    observations(:, 2) = [-0.7_dp, 1.8_dp, 1.1_dp]
    observations_dot(:, 1) = [0.004_dp, -0.003_dp, 0.002_dp]
    observations_dot(:, 2) = [-0.002_dp, 0.005_dp, 0.003_dp]
    normals(:, 1) = [1.0_dp, 0.0_dp, 0.0_dp]
    normals(:, 2) = [0.0_dp, 1.0_dp, 0.0_dp]
    normals_dot(:, 1) = [0.0_dp, 0.013_dp, -0.009_dp]
    normals_dot(:, 2) = [0.017_dp, 0.0_dp, -0.011_dp]
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
    wave_number_dot = 0.017_dp

    call evaluate_maxwell_virtual_casing_rwg_3d( &
        vertices, triangles, coefficients, observations, normals, wave_number, 6, &
        fields, normal_fields, tangential_fields, status)
    call record_condition(status == 0, &
        "virtual-casing trace map accepts target normals")
    do basis = 1, 2
        call record_condition(maxval(abs(tangential_fields(:, basis) + &
            normal_fields(basis)*normals(:, basis) - fields(:, basis))) < &
            2.0e-13_dp, "virtual-casing traces reconstruct the full magnetic field")
    end do

    call evaluate_maxwell_virtual_casing_rwg_3d_jvp( &
        vertices, triangles, coefficients, observations, normals, wave_number, 6, &
        vertices_dot, coefficients_dot, observations_dot, normals_dot, &
        wave_number_dot, fields_dot, normal_fields_dot, tangential_fields_dot, status)
    call evaluate_maxwell_virtual_casing_rwg_3d( &
        vertices + step*vertices_dot, triangles, coefficients + step*coefficients_dot, &
        observations + step*observations_dot, normals + step*normals_dot, &
        wave_number + step*wave_number_dot, 6, fields_plus, normal_fields_plus, &
        tangential_fields_plus, status)
    call evaluate_maxwell_virtual_casing_rwg_3d( &
        vertices - step*vertices_dot, triangles, coefficients - step*coefficients_dot, &
        observations - step*observations_dot, normals - step*normals_dot, &
        wave_number - step*wave_number_dot, 6, fields_minus, normal_fields_minus, &
        tangential_fields_minus, status)
    derivative_error = max( &
        maxval(abs(fields_dot - (fields_plus - fields_minus)/(2.0_dp*step))), &
        maxval(abs(normal_fields_dot - &
        (normal_fields_plus - normal_fields_minus)/(2.0_dp*step))), &
        maxval(abs(tangential_fields_dot - &
        (tangential_fields_plus - tangential_fields_minus)/(2.0_dp*step))))
    call record_condition(status == 0 .and. derivative_error < 2.0e-7_dp, &
        "virtual-casing trace JVP matches central reassembly")

    fields_bar(:, 1) = [ &
        cmplx(0.7_dp, -0.2_dp, dp), cmplx(-0.3_dp, 0.4_dp, dp), &
        cmplx(0.1_dp, 0.6_dp, dp)]
    fields_bar(:, 2) = [ &
        cmplx(-0.2_dp, 0.5_dp, dp), cmplx(0.4_dp, 0.1_dp, dp), &
        cmplx(-0.6_dp, 0.3_dp, dp)]
    normal_fields_bar = [ &
        cmplx(0.3_dp, -0.5_dp, dp), cmplx(-0.4_dp, 0.2_dp, dp)]
    tangential_fields_bar(:, 1) = [ &
        cmplx(-0.1_dp, 0.3_dp, dp), cmplx(0.5_dp, -0.2_dp, dp), &
        cmplx(0.2_dp, 0.1_dp, dp)]
    tangential_fields_bar(:, 2) = [ &
        cmplx(0.4_dp, 0.2_dp, dp), cmplx(-0.2_dp, 0.6_dp, dp), &
        cmplx(0.1_dp, -0.3_dp, dp)]
    call evaluate_maxwell_virtual_casing_rwg_3d_vjp( &
        vertices, triangles, coefficients, observations, normals, wave_number, 6, &
        fields_bar, normal_fields_bar, tangential_fields_bar, vertices_bar, &
        coefficients_bar, observations_bar, normals_bar, wave_number_bar, status)
    lhs = real(sum(conjg(fields_bar)*fields_dot), dp) + &
        real(sum(conjg(normal_fields_bar)*normal_fields_dot), dp) + &
        real(sum(conjg(tangential_fields_bar)*tangential_fields_dot), dp)
    rhs = sum(vertices_bar*vertices_dot) + &
        real(sum(conjg(coefficients_bar)*coefficients_dot), dp) + &
        sum(observations_bar*observations_dot) + sum(normals_bar*normals_dot) + &
        wave_number_bar*wave_number_dot
    adjoint_error = abs(lhs - rhs)
    call record_condition(status == 0 .and. adjoint_error < &
        5.0e-8_dp*max(1.0_dp, abs(lhs), abs(rhs)), &
        "virtual-casing trace VJP satisfies the real complex adjoint")

    call check_summary("virtual-casing magnetic traces")
    if (.not. all_passed) error stop 1

contains

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_maxwell_virtual_casing_traces_ad
