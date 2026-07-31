program test_tetra_vector_basis_ad
    use, intrinsic :: iso_fortran_env, only: error_unit
    use fortfem_api, only: evaluate_tetra_nedelec_first_kind, &
        evaluate_tetra_nedelec_first_kind_jvp, &
        evaluate_tetra_nedelec_first_kind_vjp, evaluate_tetra_rt, &
        evaluate_tetra_rt_jvp, evaluate_tetra_rt_vjp, &
        initialize_tetra_nedelec_first_kind, initialize_tetra_rt, &
        tetra_nedelec_dof_count, tetra_nedelec_first_kind_t, &
        tetra_rt_dof_count, tetra_rt_t
    use fortfem_kinds, only: dp
    implicit none

    real(dp), parameter :: step = 2.0e-7_dp
    type(tetra_nedelec_first_kind_t) :: nedelec_basis
    type(tetra_rt_t) :: rt_basis
    integer :: failures, status

    failures = 0
    call initialize_tetra_nedelec_first_kind(2, nedelec_basis, status)
    call check(status == 0, "initialize Nedelec basis")
    call check_nedelec()
    call initialize_tetra_rt(1, rt_basis, status)
    call check(status == 0, "initialize RT basis")
    call check_rt()

    if (failures /= 0) then
        write(error_unit, "(i0,a)") failures, " test(s) failed"
        stop 1
    end if
    write(*, "(a)") "PASS"

contains

    subroutine check_nedelec()
        real(dp), allocatable :: curls(:, :), curls_bar(:, :), curls_dot(:, :)
        real(dp), allocatable :: curls_minus(:, :), curls_plus(:, :)
        real(dp), allocatable :: values(:, :), values_bar(:, :)
        real(dp), allocatable :: values_dot(:, :), values_minus(:, :)
        real(dp), allocatable :: values_plus(:, :)
        real(dp) :: lhs, point(3), point_bar(3), point_dot(3), rhs
        integer :: count, dof

        count = tetra_nedelec_dof_count(nedelec_basis)
        allocate(values(3, count), values_dot(3, count))
        allocate(values_minus(3, count), values_plus(3, count))
        allocate(values_bar(3, count))
        allocate(curls(3, count), curls_dot(3, count))
        allocate(curls_minus(3, count), curls_plus(3, count))
        allocate(curls_bar(3, count))
        point = [0.19_dp, 0.21_dp, 0.17_dp]
        point_dot = [0.13_dp, -0.07_dp, 0.09_dp]
        do dof = 1, count
            values_bar(:, dof) = [0.01_dp*dof, -0.006_dp*dof, 0.004_dp*dof]
            curls_bar(:, dof) = [-0.003_dp*dof, 0.005_dp*dof, 0.002_dp*dof]
        end do
        call evaluate_tetra_nedelec_first_kind_jvp( &
            nedelec_basis, point, point_dot, values_dot, curls_dot, status)
        call evaluate_tetra_nedelec_first_kind( &
            nedelec_basis, point + step*point_dot, values_plus, curls_plus, &
            status)
        call evaluate_tetra_nedelec_first_kind( &
            nedelec_basis, point - step*point_dot, values_minus, curls_minus, &
            status)
        call check(maxval(abs(values_dot - &
            (values_plus - values_minus)/(2.0_dp*step))) < 3.0e-8_dp, &
            "Nedelec basis JVP matches re-evaluation")
        call check(maxval(abs(curls_dot - &
            (curls_plus - curls_minus)/(2.0_dp*step))) < 3.0e-8_dp, &
            "Nedelec curl JVP matches re-evaluation")
        call evaluate_tetra_nedelec_first_kind( &
            nedelec_basis, point, values, curls, status)
        call evaluate_tetra_nedelec_first_kind_vjp( &
            nedelec_basis, point, values_bar, curls_bar, point_bar, status)
        lhs = sum(values_bar*values_dot) + sum(curls_bar*curls_dot)
        rhs = dot_product(point_bar, point_dot)
        call check(abs(lhs - rhs) < 2.0e-11_dp*max(1.0_dp, abs(lhs), abs(rhs)), &
            "Nedelec basis JVP and VJP obey adjoint identity")
    end subroutine check_nedelec

    subroutine check_rt()
        real(dp), allocatable :: divergences(:), divergences_bar(:)
        real(dp), allocatable :: divergences_dot(:), divergences_minus(:)
        real(dp), allocatable :: divergences_plus(:), values(:, :)
        real(dp), allocatable :: values_bar(:, :), values_dot(:, :)
        real(dp), allocatable :: values_minus(:, :), values_plus(:, :)
        real(dp) :: lhs, point(3), point_bar(3), point_dot(3), rhs
        integer :: count, dof

        count = tetra_rt_dof_count(rt_basis)
        allocate(values(3, count), values_dot(3, count))
        allocate(values_minus(3, count), values_plus(3, count))
        allocate(values_bar(3, count))
        allocate(divergences(count), divergences_dot(count))
        allocate(divergences_minus(count), divergences_plus(count))
        allocate(divergences_bar(count))
        point = [0.19_dp, 0.21_dp, 0.17_dp]
        point_dot = [0.13_dp, -0.07_dp, 0.09_dp]
        do dof = 1, count
            values_bar(:, dof) = [0.008_dp*dof, -0.005_dp*dof, 0.003_dp*dof]
            divergences_bar(dof) = -0.004_dp*dof
        end do
        call evaluate_tetra_rt_jvp( &
            rt_basis, point, point_dot, values_dot, divergences_dot, status)
        call evaluate_tetra_rt( &
            rt_basis, point + step*point_dot, values_plus, divergences_plus, &
            status)
        call evaluate_tetra_rt( &
            rt_basis, point - step*point_dot, values_minus, divergences_minus, &
            status)
        call check(maxval(abs(values_dot - &
            (values_plus - values_minus)/(2.0_dp*step))) < 3.0e-8_dp, &
            "RT basis JVP matches re-evaluation")
        call check_close(maxval(abs(divergences_dot - &
            (divergences_plus - divergences_minus)/(2.0_dp*step))), &
            5.0e-8_dp, "RT divergence JVP matches re-evaluation")
        call evaluate_tetra_rt(rt_basis, point, values, divergences, status)
        call evaluate_tetra_rt_vjp( &
            rt_basis, point, values_bar, divergences_bar, point_bar, status)
        lhs = sum(values_bar*values_dot) + &
            dot_product(divergences_bar, divergences_dot)
        rhs = dot_product(point_bar, point_dot)
        call check(abs(lhs - rhs) < 2.0e-11_dp*max(1.0_dp, abs(lhs), abs(rhs)), &
            "RT basis JVP and VJP obey adjoint identity")
    end subroutine check_rt

    subroutine check(condition, label)
        logical, intent(in) :: condition
        character(*), intent(in) :: label

        if (condition) return
        failures = failures + 1
        write(error_unit, "(a)") "FAIL: "//label
    end subroutine check

    subroutine check_close(error, tolerance, label)
        real(dp), intent(in) :: error, tolerance
        character(*), intent(in) :: label

        if (error < tolerance) return
        failures = failures + 1
        write(error_unit, "(a,es12.4)") "FAIL: "//label//" error=", error
    end subroutine check_close

end program test_tetra_vector_basis_ad
