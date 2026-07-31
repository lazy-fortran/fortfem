program test_triangle_full_vector_point_observation_ad
    use, intrinsic :: iso_fortran_env, only: error_unit
    use fortfem_api, only: &
        evaluate_triangle_bdm_interpolant_at_point, &
        evaluate_triangle_bdm_interpolant_at_point_jvp, &
        evaluate_triangle_bdm_interpolant_at_point_vjp, &
        evaluate_triangle_nedelec_interpolant_at_point, &
        evaluate_triangle_nedelec_interpolant_at_point_jvp, &
        evaluate_triangle_nedelec_interpolant_at_point_vjp, &
        evaluate_triangle_nedelec_second_interpolant_at_point, &
        evaluate_triangle_nedelec_second_interpolant_at_point_jvp, &
        evaluate_triangle_nedelec_second_interpolant_at_point_vjp, &
        evaluate_triangle_rt_interpolant_at_point, &
        evaluate_triangle_rt_interpolant_at_point_jvp, &
        evaluate_triangle_rt_interpolant_at_point_vjp, &
        initialize_triangle_bdm, initialize_triangle_nedelec_first_kind, &
        initialize_triangle_nedelec_second_kind, &
        initialize_triangle_raviart_thomas, &
        triangle_bdm_basis_t, triangle_bdm_dof_count, &
        triangle_nedelec_dof_count, triangle_nedelec_first_kind_t, &
        triangle_nedelec_second_kind_dof_count, &
        triangle_nedelec_second_kind_t, triangle_rt_basis_t, &
        triangle_rt_dof_count
    use fortfem_kinds, only: dp
    implicit none

    real(dp), parameter :: step = 2.0e-6_dp
    type(triangle_bdm_basis_t) :: bdm_basis
    type(triangle_nedelec_first_kind_t) :: nedelec_first_basis
    type(triangle_nedelec_second_kind_t) :: nedelec_basis
    type(triangle_rt_basis_t) :: rt_basis
    real(dp), allocatable :: dofs(:), dofs_bar(:), dofs_dot(:)
    real(dp) :: point(2), point_bar(2), point_dot(2)
    real(dp) :: scalar, scalar_bar, scalar_dot, scalar_minus, scalar_plus
    real(dp) :: value(2), value_bar(2), value_dot(2)
    real(dp) :: value_minus(2), value_plus(2)
    real(dp) :: vertices(2, 3), vertices_bar(2, 3), vertices_dot(2, 3)
    integer :: failures, status

    failures = 0
    vertices = reshape([ &
        0.1_dp, -0.2_dp, 1.4_dp, 0.15_dp, -0.3_dp, 1.2_dp], [2, 3])
    vertices_dot = reshape([ &
        0.03_dp, -0.01_dp, -0.02_dp, 0.04_dp, 0.01_dp, -0.025_dp], [2, 3])
    point = [0.37_dp, 0.28_dp]
    point_dot = [-0.015_dp, 0.022_dp]
    value_bar = [0.7_dp, -0.4_dp]
    scalar_bar = 0.35_dp

    call initialize_triangle_bdm(2, bdm_basis, status)
    call initialize_arrays(triangle_bdm_dof_count(bdm_basis))
    call check_bdm()
    call initialize_triangle_raviart_thomas(2, rt_basis, status)
    call initialize_arrays(triangle_rt_dof_count(rt_basis))
    call check_rt()
    call initialize_triangle_nedelec_first_kind(2, nedelec_first_basis, status)
    call initialize_arrays(triangle_nedelec_dof_count(nedelec_first_basis))
    call check_nedelec_first()
    call initialize_triangle_nedelec_second_kind(2, nedelec_basis, status)
    call initialize_arrays( &
        triangle_nedelec_second_kind_dof_count(nedelec_basis))
    call check_nedelec()

    if (failures > 0) then
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
            dofs(dof) = 0.08_dp*dof - 0.17_dp
            dofs_dot(dof) = -0.015_dp*dof + 0.04_dp
        end do
    end subroutine initialize_arrays

    subroutine check_bdm()
        real(dp) :: lhs, rhs

        call evaluate_triangle_bdm_interpolant_at_point_jvp( &
            vertices, bdm_basis, dofs, point, vertices_dot, dofs_dot, &
            point_dot, value_dot, scalar_dot, status)
        call evaluate_triangle_bdm_interpolant_at_point( &
            vertices + step*vertices_dot, bdm_basis, dofs + step*dofs_dot, &
            point + step*point_dot, value_plus, scalar_plus, status)
        call evaluate_triangle_bdm_interpolant_at_point( &
            vertices - step*vertices_dot, bdm_basis, dofs - step*dofs_dot, &
            point - step*point_dot, value_minus, scalar_minus, status)
        call check_products("BDM")
        call evaluate_triangle_bdm_interpolant_at_point( &
            vertices, bdm_basis, dofs, point, value, scalar, status)
        call evaluate_triangle_bdm_interpolant_at_point_vjp( &
            vertices, bdm_basis, dofs, point, value_bar, scalar_bar, &
            vertices_bar, dofs_bar, point_bar, status)
        lhs = dot_product(value_bar, value_dot) + scalar_bar*scalar_dot
        rhs = sum(vertices_bar*vertices_dot) + dot_product(dofs_bar, dofs_dot) + &
            dot_product(point_bar, point_dot)
        call check(abs(lhs - rhs) < 3.0e-11_dp, &
            "BDM physical-point products obey adjoint identity")
    end subroutine check_bdm

    subroutine check_rt()
        real(dp) :: lhs, rhs

        call evaluate_triangle_rt_interpolant_at_point_jvp( &
            vertices, rt_basis, dofs, point, vertices_dot, dofs_dot, &
            point_dot, value_dot, scalar_dot, status)
        call evaluate_triangle_rt_interpolant_at_point( &
            vertices + step*vertices_dot, rt_basis, dofs + step*dofs_dot, &
            point + step*point_dot, value_plus, scalar_plus, status)
        call evaluate_triangle_rt_interpolant_at_point( &
            vertices - step*vertices_dot, rt_basis, dofs - step*dofs_dot, &
            point - step*point_dot, value_minus, scalar_minus, status)
        call check_products("Raviart-Thomas")
        call evaluate_triangle_rt_interpolant_at_point_vjp( &
            vertices, rt_basis, dofs, point, value_bar, scalar_bar, &
            vertices_bar, dofs_bar, point_bar, status)
        lhs = dot_product(value_bar, value_dot) + scalar_bar*scalar_dot
        rhs = sum(vertices_bar*vertices_dot) + dot_product(dofs_bar, dofs_dot) + &
            dot_product(point_bar, point_dot)
        call check(abs(lhs - rhs) < 3.0e-11_dp, &
            "Raviart-Thomas physical-point products obey adjoint identity")
    end subroutine check_rt

    subroutine check_nedelec_first()
        real(dp) :: lhs, rhs

        call evaluate_triangle_nedelec_interpolant_at_point_jvp( &
            vertices, nedelec_first_basis, dofs, point, vertices_dot, &
            dofs_dot, point_dot, value_dot, scalar_dot, status)
        call evaluate_triangle_nedelec_interpolant_at_point( &
            vertices + step*vertices_dot, nedelec_first_basis, &
            dofs + step*dofs_dot, point + step*point_dot, value_plus, &
            scalar_plus, status)
        call evaluate_triangle_nedelec_interpolant_at_point( &
            vertices - step*vertices_dot, nedelec_first_basis, &
            dofs - step*dofs_dot, point - step*point_dot, value_minus, &
            scalar_minus, status)
        call check_products("first-kind Nedelec")
        call evaluate_triangle_nedelec_interpolant_at_point_vjp( &
            vertices, nedelec_first_basis, dofs, point, value_bar, scalar_bar, &
            vertices_bar, dofs_bar, point_bar, status)
        lhs = dot_product(value_bar, value_dot) + scalar_bar*scalar_dot
        rhs = sum(vertices_bar*vertices_dot) + dot_product(dofs_bar, dofs_dot) + &
            dot_product(point_bar, point_dot)
        call check(abs(lhs - rhs) < 3.0e-11_dp, &
            "first-kind Nedelec physical-point products obey adjoint identity")
    end subroutine check_nedelec_first

    subroutine check_nedelec()
        real(dp) :: lhs, rhs

        call evaluate_triangle_nedelec_second_interpolant_at_point_jvp( &
            vertices, nedelec_basis, dofs, point, vertices_dot, dofs_dot, &
            point_dot, value_dot, scalar_dot, status)
        call evaluate_triangle_nedelec_second_interpolant_at_point( &
            vertices + step*vertices_dot, nedelec_basis, &
            dofs + step*dofs_dot, point + step*point_dot, value_plus, &
            scalar_plus, status)
        call evaluate_triangle_nedelec_second_interpolant_at_point( &
            vertices - step*vertices_dot, nedelec_basis, &
            dofs - step*dofs_dot, point - step*point_dot, value_minus, &
            scalar_minus, status)
        call check_products("second-kind Nedelec")
        call evaluate_triangle_nedelec_second_interpolant_at_point( &
            vertices, nedelec_basis, dofs, point, value, scalar, status)
        call evaluate_triangle_nedelec_second_interpolant_at_point_vjp( &
            vertices, nedelec_basis, dofs, point, value_bar, scalar_bar, &
            vertices_bar, dofs_bar, point_bar, status)
        lhs = dot_product(value_bar, value_dot) + scalar_bar*scalar_dot
        rhs = sum(vertices_bar*vertices_dot) + dot_product(dofs_bar, dofs_dot) + &
            dot_product(point_bar, point_dot)
        call check(abs(lhs - rhs) < 3.0e-11_dp, &
            "second-kind Nedelec physical-point products obey adjoint identity")
    end subroutine check_nedelec

    subroutine check_products(label)
        character(*), intent(in) :: label

        call check(status == 0, label//" physical-point products succeed")
        call check(maxval(abs( &
            value_dot - (value_plus - value_minus)/(2.0_dp*step))) < &
            3.0e-8_dp, label//" point value JVP matches full re-evaluation")
        call check(abs( &
            scalar_dot - (scalar_plus - scalar_minus)/(2.0_dp*step)) < &
            3.0e-8_dp, label//" point scalar JVP matches full re-evaluation")
    end subroutine check_products

    subroutine check(condition, label)
        logical, intent(in) :: condition
        character(*), intent(in) :: label

        if (condition) return
        failures = failures + 1
        write(error_unit, "(a)") "FAIL: "//label
    end subroutine check

end program test_triangle_full_vector_point_observation_ad
