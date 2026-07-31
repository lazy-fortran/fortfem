program test_tetra_lagrange_basis_ad
    use, intrinsic :: iso_fortran_env, only: error_unit
    use fortfem_api, only: evaluate_tetra_lagrange, &
        evaluate_tetra_lagrange_jvp, evaluate_tetra_lagrange_vjp, &
        initialize_tetra_lagrange, tetra_lagrange_dof_count, tetra_lagrange_t
    use fortfem_kinds, only: dp
    implicit none

    real(dp), parameter :: step = 2.0e-7_dp
    type(tetra_lagrange_t) :: basis
    real(dp), allocatable :: gradients(:, :), gradients_bar(:, :)
    real(dp), allocatable :: gradients_dot(:, :), gradients_minus(:, :)
    real(dp), allocatable :: gradients_plus(:, :), values(:), values_bar(:)
    real(dp), allocatable :: values_dot(:), values_minus(:), values_plus(:)
    real(dp) :: lhs, point(3), point_bar(3), point_dot(3), rhs
    integer :: count, dof, failures, status

    failures = 0
    call initialize_tetra_lagrange(4, basis, status)
    count = tetra_lagrange_dof_count(basis)
    allocate(values(count), values_bar(count), values_dot(count))
    allocate(values_minus(count), values_plus(count))
    allocate(gradients(3, count), gradients_bar(3, count))
    allocate(gradients_dot(3, count))
    allocate(gradients_minus(3, count), gradients_plus(3, count))
    point = [0.19_dp, 0.21_dp, 0.17_dp]
    point_dot = [0.13_dp, -0.07_dp, 0.09_dp]
    do dof = 1, count
        values_bar(dof) = 0.006_dp*dof - 0.03_dp
        gradients_bar(:, dof) = &
            [0.003_dp*dof, -0.002_dp*dof, 0.001_dp*dof]
    end do

    call evaluate_tetra_lagrange_jvp( &
        basis, point, point_dot, values_dot, gradients_dot, status)
    call evaluate_tetra_lagrange( &
        basis, point + step*point_dot, values_plus, gradients_plus, status)
    call evaluate_tetra_lagrange( &
        basis, point - step*point_dot, values_minus, gradients_minus, status)
    call check(maxval(abs(values_dot - &
        (values_plus - values_minus)/(2.0_dp*step))) < 5.0e-8_dp, &
        "tetrahedral Lagrange value JVP")
    call check(maxval(abs(gradients_dot - &
        (gradients_plus - gradients_minus)/(2.0_dp*step))) < 2.0e-7_dp, &
        "tetrahedral Lagrange gradient JVP")
    call evaluate_tetra_lagrange( &
        basis, point, values, gradients, status)
    call evaluate_tetra_lagrange_vjp( &
        basis, point, values_bar, gradients_bar, point_bar, status)
    lhs = dot_product(values_bar, values_dot) + sum(gradients_bar*gradients_dot)
    rhs = dot_product(point_bar, point_dot)
    call check(abs(lhs - rhs) < 3.0e-11_dp*max(1.0_dp, abs(lhs), abs(rhs)), &
        "tetrahedral Lagrange basis adjoint identity")

    if (failures /= 0) then
        write(error_unit, "(i0,a)") failures, " test(s) failed"
        stop 1
    end if
    write(*, "(a)") "PASS"

contains

    subroutine check(condition, label)
        logical, intent(in) :: condition
        character(*), intent(in) :: label

        if (condition) return
        failures = failures + 1
        write(error_unit, "(a)") "FAIL: "//label
    end subroutine check

end program test_tetra_lagrange_basis_ad
