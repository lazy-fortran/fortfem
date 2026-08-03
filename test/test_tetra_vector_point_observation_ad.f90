program test_tetra_vector_point_observation_ad
    use, intrinsic :: iso_fortran_env, only: error_unit
    use fortfem_feec, only: &
        evaluate_tetra_nedelec_interpolant_at_point, &
        evaluate_tetra_nedelec_interpolant_at_point_jvp, &
        evaluate_tetra_nedelec_interpolant_at_point_vjp, &
        evaluate_tetra_rt_interpolant_at_point, &
        evaluate_tetra_rt_interpolant_at_point_jvp, &
        evaluate_tetra_rt_interpolant_at_point_vjp, &
        initialize_tetra_nedelec_first_kind, initialize_tetra_rt, &
        tetra_nedelec_dof_count, tetra_nedelec_first_kind_t, &
        tetra_rt_dof_count, tetra_rt_t
    use fortfem_kinds, only: dp
    implicit none

    real(dp), parameter :: step = 2.0e-7_dp
    type(tetra_nedelec_first_kind_t) :: nedelec_basis
    type(tetra_rt_t) :: rt_basis
    real(dp), allocatable :: dofs(:), dofs_bar(:), dofs_dot(:)
    real(dp) :: point(3), point_bar(3), point_dot(3)
    real(dp) :: vertices(3, 4), vertices_bar(3, 4), vertices_dot(3, 4)
    integer :: failures, status

    failures = 0
    vertices = reshape([ &
        0.1_dp, -0.2_dp, 0.05_dp, 1.3_dp, 0.1_dp, -0.1_dp, &
        -0.15_dp, 1.2_dp, 0.2_dp, 0.2_dp, -0.1_dp, 1.4_dp], [3, 4])
    vertices_dot = reshape([ &
        0.01_dp, -0.02_dp, 0.015_dp, -0.03_dp, 0.025_dp, 0.01_dp, &
        0.02_dp, 0.015_dp, -0.025_dp, -0.01_dp, 0.03_dp, 0.02_dp], [3, 4])
    point = vertices(:, 1) + &
        matmul(vertices(:, 2:4) - spread(vertices(:, 1), 2, 3), &
        [0.21_dp, 0.18_dp, 0.27_dp])
    point_dot = [0.025_dp, -0.018_dp, 0.014_dp]

    call initialize_tetra_nedelec_first_kind(2, nedelec_basis, status)
    call initialize_arrays(tetra_nedelec_dof_count(nedelec_basis))
    call check_nedelec()
    call initialize_tetra_rt(1, rt_basis, status)
    call initialize_arrays(tetra_rt_dof_count(rt_basis))
    call check_rt()

    if (failures /= 0) then
        write(error_unit, "(i0,a)") failures, " test(s) failed"
        stop 1
    end if
    write(*, "(a)") "PASS"

contains

    subroutine initialize_arrays(count)
        integer, intent(in) :: count
        integer :: dof

        if (allocated(dofs)) deallocate(dofs, dofs_bar, dofs_dot)
        allocate(dofs(count), dofs_bar(count), dofs_dot(count))
        do dof = 1, count
            dofs(dof) = 0.035_dp*dof - 0.12_dp
            dofs_dot(dof) = -0.009_dp*dof + 0.028_dp
        end do
    end subroutine initialize_arrays

    subroutine check_nedelec()
        real(dp) :: curl_bar(3), curl_dot(3), curl_minus(3), curl_plus(3)
        real(dp) :: lhs, rhs, value_bar(3), value_dot(3)
        real(dp) :: value_minus(3), value_plus(3)

        value_bar = [0.6_dp, -0.35_dp, 0.22_dp]
        curl_bar = [0.17_dp, -0.21_dp, 0.13_dp]
        call evaluate_tetra_nedelec_interpolant_at_point_jvp( &
            vertices, nedelec_basis, dofs, point, vertices_dot, dofs_dot, &
            point_dot, value_dot, curl_dot, status)
        call evaluate_tetra_nedelec_interpolant_at_point( &
            vertices + step*vertices_dot, nedelec_basis, dofs + step*dofs_dot, &
            point + step*point_dot, value_plus, curl_plus, status)
        call evaluate_tetra_nedelec_interpolant_at_point( &
            vertices - step*vertices_dot, nedelec_basis, dofs - step*dofs_dot, &
            point - step*point_dot, value_minus, curl_minus, status)
        call check(maxval(abs(value_dot - &
            (value_plus - value_minus)/(2.0_dp*step))) < 3.0e-8_dp, &
            "physical-point Nedelec value JVP")
        call check(maxval(abs(curl_dot - &
            (curl_plus - curl_minus)/(2.0_dp*step))) < 3.0e-8_dp, &
            "physical-point Nedelec curl JVP")
        call evaluate_tetra_nedelec_interpolant_at_point_vjp( &
            vertices, nedelec_basis, dofs, point, value_bar, curl_bar, &
            vertices_bar, dofs_bar, point_bar, status)
        lhs = dot_product(value_bar, value_dot) + dot_product(curl_bar, curl_dot)
        rhs = sum(vertices_bar*vertices_dot) + dot_product(dofs_bar, dofs_dot) + &
            dot_product(point_bar, point_dot)
        call check(abs(lhs - rhs) < 3.0e-11_dp*max(1.0_dp, abs(lhs), abs(rhs)), &
            "physical-point Nedelec adjoint identity")
    end subroutine check_nedelec

    subroutine check_rt()
        real(dp) :: divergence_bar, divergence_dot, divergence_minus
        real(dp) :: divergence_plus, lhs, rhs, value_bar(3), value_dot(3)
        real(dp) :: value_minus(3), value_plus(3)

        value_bar = [0.6_dp, -0.35_dp, 0.22_dp]
        divergence_bar = -0.45_dp
        call evaluate_tetra_rt_interpolant_at_point_jvp( &
            vertices, rt_basis, dofs, point, vertices_dot, dofs_dot, &
            point_dot, value_dot, divergence_dot, status)
        call evaluate_tetra_rt_interpolant_at_point( &
            vertices + step*vertices_dot, rt_basis, dofs + step*dofs_dot, &
            point + step*point_dot, value_plus, divergence_plus, status)
        call evaluate_tetra_rt_interpolant_at_point( &
            vertices - step*vertices_dot, rt_basis, dofs - step*dofs_dot, &
            point - step*point_dot, value_minus, divergence_minus, status)
        call check(maxval(abs(value_dot - &
            (value_plus - value_minus)/(2.0_dp*step))) < 4.0e-8_dp, &
            "physical-point RT value JVP")
        call check(abs(divergence_dot - &
            (divergence_plus - divergence_minus)/(2.0_dp*step)) < 5.0e-8_dp, &
            "physical-point RT divergence JVP")
        call evaluate_tetra_rt_interpolant_at_point_vjp( &
            vertices, rt_basis, dofs, point, value_bar, divergence_bar, &
            vertices_bar, dofs_bar, point_bar, status)
        lhs = dot_product(value_bar, value_dot) + divergence_bar*divergence_dot
        rhs = sum(vertices_bar*vertices_dot) + dot_product(dofs_bar, dofs_dot) + &
            dot_product(point_bar, point_dot)
        call check(abs(lhs - rhs) < 3.0e-11_dp*max(1.0_dp, abs(lhs), abs(rhs)), &
            "physical-point RT adjoint identity")
    end subroutine check_rt

    subroutine check(condition, label)
        logical, intent(in) :: condition
        character(*), intent(in) :: label

        if (condition) return
        failures = failures + 1
        write(error_unit, "(a)") "FAIL: "//label
    end subroutine check

end program test_tetra_vector_point_observation_ad
