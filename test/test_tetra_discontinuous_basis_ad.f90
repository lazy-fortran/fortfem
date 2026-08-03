program test_tetra_discontinuous_basis_ad
    use, intrinsic :: iso_fortran_env, only: error_unit
    use fortfem_tetra_discontinuous_arbitrary_order, only: &
        evaluate_tetra_discontinuous, &
        evaluate_tetra_discontinuous_jvp, evaluate_tetra_discontinuous_vjp, &
        initialize_tetra_discontinuous, tetra_discontinuous_dof_count, &
        tetra_discontinuous_t
    use fortfem_kinds, only: dp
    implicit none

    real(dp), parameter :: step = 2.0e-7_dp
    type(tetra_discontinuous_t) :: basis
    real(dp), allocatable :: values(:), values_bar(:), values_dot(:)
    real(dp), allocatable :: values_minus(:), values_plus(:)
    real(dp) :: lhs, point(3), point_bar(3), point_dot(3), rhs
    integer :: count, dof, failures, status

    failures = 0
    call initialize_tetra_discontinuous(6, basis, status)
    count = tetra_discontinuous_dof_count(basis)
    allocate(values(count), values_bar(count), values_dot(count))
    allocate(values_minus(count), values_plus(count))
    point = [0.19_dp, 0.21_dp, 0.17_dp]
    point_dot = [0.13_dp, -0.07_dp, 0.09_dp]
    do dof = 1, count
        values_bar(dof) = 0.003_dp*dof - 0.04_dp
    end do

    call evaluate_tetra_discontinuous_jvp( &
        basis, point, point_dot, values_dot, status)
    call evaluate_tetra_discontinuous( &
        basis, point(1) + step*point_dot(1), &
        point(2) + step*point_dot(2), point(3) + step*point_dot(3), &
        values_plus, status)
    call evaluate_tetra_discontinuous( &
        basis, point(1) - step*point_dot(1), &
        point(2) - step*point_dot(2), point(3) - step*point_dot(3), &
        values_minus, status)
    call check(maxval(abs(values_dot - &
        (values_plus - values_minus)/(2.0_dp*step))) < 2.0e-7_dp, &
        "modal discontinuous basis JVP")
    call evaluate_tetra_discontinuous( &
        basis, point(1), point(2), point(3), values, status)
    call evaluate_tetra_discontinuous_vjp( &
        basis, point, values_bar, point_bar, status)
    lhs = dot_product(values_bar, values_dot)
    rhs = dot_product(point_bar, point_dot)
    call check(abs(lhs - rhs) < 3.0e-11_dp*max(1.0_dp, abs(lhs), abs(rhs)), &
        "modal discontinuous basis adjoint identity")

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

end program test_tetra_discontinuous_basis_ad
