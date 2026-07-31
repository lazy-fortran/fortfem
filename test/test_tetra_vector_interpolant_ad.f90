program test_tetra_vector_interpolant_ad
    use, intrinsic :: iso_fortran_env, only: error_unit
    use fortfem_api, only: evaluate_tetra_nedelec_interpolant, &
        evaluate_tetra_nedelec_interpolant_jvp, &
        evaluate_tetra_nedelec_interpolant_vjp, &
        evaluate_tetra_rt_interpolant, evaluate_tetra_rt_interpolant_jvp, &
        evaluate_tetra_rt_interpolant_vjp, &
        initialize_tetra_nedelec_first_kind, initialize_tetra_rt, &
        tetra_nedelec_dof_count, tetra_nedelec_first_kind_t, &
        tetra_rt_dof_count, tetra_rt_t
    use fortfem_kinds, only: dp
    implicit none

    real(dp), parameter :: step = 2.0e-7_dp
    type(tetra_nedelec_first_kind_t) :: nedelec_basis
    type(tetra_rt_t) :: rt_basis
    real(dp), allocatable :: dofs(:), dofs_bar(:), dofs_dot(:)
    real(dp) :: scalar, scalar_bar, scalar_dot, scalar_minus, scalar_plus
    real(dp) :: vector(3), vector_bar(3), vector_dot(3)
    real(dp) :: vector_minus(3), vector_plus(3)
    real(dp) :: vertices(3, 4), vertices_bar(3, 4), vertices_dot(3, 4)
    real(dp) :: reference_point(3)
    integer :: failures, status

    failures = 0
    vertices = reshape([ &
        0.1_dp, -0.2_dp, 0.05_dp, &
        1.3_dp, 0.1_dp, -0.1_dp, &
        -0.15_dp, 1.2_dp, 0.2_dp, &
        0.2_dp, -0.1_dp, 1.4_dp], [3, 4])
    vertices_dot = reshape([ &
        0.01_dp, -0.02_dp, 0.015_dp, &
        -0.03_dp, 0.025_dp, 0.01_dp, &
        0.02_dp, 0.015_dp, -0.025_dp, &
        -0.01_dp, 0.03_dp, 0.02_dp], [3, 4])
    reference_point = [0.21_dp, 0.18_dp, 0.27_dp]
    vector_bar = [0.6_dp, -0.35_dp, 0.22_dp]
    scalar_bar = -0.45_dp

    call initialize_tetra_nedelec_first_kind(2, nedelec_basis, status)
    call initialize_arrays(tetra_nedelec_dof_count(nedelec_basis))
    call check_nedelec()
    call initialize_tetra_rt(1, rt_basis, status)
    call initialize_arrays(tetra_rt_dof_count(rt_basis))
    call check_rt()

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
            dofs(dof) = 0.035_dp*dof - 0.12_dp
            dofs_dot(dof) = -0.009_dp*dof + 0.028_dp
        end do
    end subroutine initialize_arrays

    subroutine check_nedelec()
        real(dp) :: curl(3), curl_bar(3), curl_dot(3)
        real(dp) :: curl_minus(3), curl_plus(3), lhs, rhs

        curl_bar = [0.17_dp, -0.21_dp, 0.13_dp]
        call evaluate_tetra_nedelec_interpolant_jvp( &
            vertices, nedelec_basis, dofs, reference_point, vertices_dot, &
            dofs_dot, vector_dot, curl_dot, status)
        call evaluate_tetra_nedelec_interpolant( &
            vertices + step*vertices_dot, nedelec_basis, &
            dofs + step*dofs_dot, reference_point, vector_plus, curl_plus, &
            status)
        call evaluate_tetra_nedelec_interpolant( &
            vertices - step*vertices_dot, nedelec_basis, &
            dofs - step*dofs_dot, reference_point, vector_minus, curl_minus, &
            status)
        call check(status == 0, "tetrahedral Nedelec observation JVP succeeds")
        call check(maxval(abs( &
            vector_dot - (vector_plus - vector_minus)/(2.0_dp*step))) < &
            2.0e-8_dp, "tetrahedral Nedelec value JVP matches re-evaluation")
        call check(maxval(abs( &
            curl_dot - (curl_plus - curl_minus)/(2.0_dp*step))) < 2.0e-8_dp, &
            "tetrahedral Nedelec curl JVP matches re-evaluation")
        call evaluate_tetra_nedelec_interpolant( &
            vertices, nedelec_basis, dofs, reference_point, vector, curl, status)
        call evaluate_tetra_nedelec_interpolant_vjp( &
            vertices, nedelec_basis, dofs, reference_point, vector_bar, &
            curl_bar, vertices_bar, dofs_bar, status)
        lhs = dot_product(vector_bar, vector_dot) + &
            dot_product(curl_bar, curl_dot)
        rhs = sum(vertices_bar*vertices_dot) + dot_product(dofs_bar, dofs_dot)
        call check(status == 0, "tetrahedral Nedelec observation VJP succeeds")
        call check(abs(lhs - rhs) < 4.0e-12_dp*max(1.0_dp, abs(lhs), abs(rhs)), &
            "tetrahedral Nedelec observation products obey adjoint identity")
    end subroutine check_nedelec

    subroutine check_rt()
        real(dp) :: lhs, rhs

        call evaluate_tetra_rt_interpolant_jvp( &
            vertices, rt_basis, dofs, reference_point, vertices_dot, dofs_dot, &
            vector_dot, scalar_dot, status)
        call evaluate_tetra_rt_interpolant( &
            vertices + step*vertices_dot, rt_basis, dofs + step*dofs_dot, &
            reference_point, vector_plus, scalar_plus, status)
        call evaluate_tetra_rt_interpolant( &
            vertices - step*vertices_dot, rt_basis, dofs - step*dofs_dot, &
            reference_point, vector_minus, scalar_minus, status)
        call check(status == 0, "tetrahedral RT observation JVP succeeds")
        call check(maxval(abs( &
            vector_dot - (vector_plus - vector_minus)/(2.0_dp*step))) < &
            2.0e-8_dp, "tetrahedral RT value JVP matches re-evaluation")
        call check(abs( &
            scalar_dot - (scalar_plus - scalar_minus)/(2.0_dp*step)) < &
            2.0e-8_dp, "tetrahedral RT divergence JVP matches re-evaluation")
        call evaluate_tetra_rt_interpolant( &
            vertices, rt_basis, dofs, reference_point, vector, scalar, status)
        call evaluate_tetra_rt_interpolant_vjp( &
            vertices, rt_basis, dofs, reference_point, vector_bar, scalar_bar, &
            vertices_bar, dofs_bar, status)
        lhs = dot_product(vector_bar, vector_dot) + scalar_bar*scalar_dot
        rhs = sum(vertices_bar*vertices_dot) + dot_product(dofs_bar, dofs_dot)
        call check(status == 0, "tetrahedral RT observation VJP succeeds")
        call check(abs(lhs - rhs) < 4.0e-12_dp*max(1.0_dp, abs(lhs), abs(rhs)), &
            "tetrahedral RT observation products obey adjoint identity")
    end subroutine check_rt

    subroutine check(condition, label)
        logical, intent(in) :: condition
        character(*), intent(in) :: label

        if (condition) return
        failures = failures + 1
        write(error_unit, "(a)") "FAIL: "//label
    end subroutine check

end program test_tetra_vector_interpolant_ad
