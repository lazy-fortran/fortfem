program test_tetra_lagrange_solution_ad
    use, intrinsic :: iso_fortran_env, only: error_unit
    use fortfem_api, only: evaluate_tetra_lagrange_solution, &
        evaluate_tetra_lagrange_solution_jvp, &
        evaluate_tetra_lagrange_solution_vjp, initialize_tetra_lagrange, &
        tetra_lagrange_dof_count, tetra_lagrange_t
    use fortfem_kinds, only: dp
    implicit none

    integer, parameter :: degree = 2
    real(dp), parameter :: step = 2.0e-7_dp
    type(tetra_lagrange_t) :: basis
    integer :: tetrahedra(4, 1)
    real(dp), allocatable :: solution(:), solution_bar(:), solution_dot(:)
    real(dp) :: gradient(3), gradient_bar(3), gradient_dot(3)
    real(dp) :: gradient_minus(3), gradient_plus(3)
    real(dp) :: lhs, reference_point(3), rhs, value, value_bar, value_dot
    real(dp) :: value_minus, value_plus
    real(dp) :: vertices(3, 4), vertices_bar(3, 4), vertices_dot(3, 4)
    integer :: dof, failures, status

    failures = 0
    tetrahedra(:, 1) = [1, 2, 3, 4]
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
    gradient_bar = [0.6_dp, -0.35_dp, 0.22_dp]
    value_bar = -0.45_dp

    call initialize_tetra_lagrange(degree, basis, status)
    allocate(solution(tetra_lagrange_dof_count(basis)))
    allocate(solution_bar(size(solution)), solution_dot(size(solution)))
    do dof = 1, size(solution)
        solution(dof) = 0.07_dp*dof - 0.19_dp
        solution_dot(dof) = -0.012_dp*dof + 0.035_dp
    end do

    call evaluate_tetra_lagrange_solution_jvp( &
        vertices, tetrahedra, degree, solution, 1, reference_point, &
        vertices_dot, solution_dot, value_dot, gradient_dot, status)
    call evaluate_tetra_lagrange_solution( &
        vertices + step*vertices_dot, tetrahedra, degree, &
        solution + step*solution_dot, 1, reference_point, value_plus, &
        gradient_plus, status)
    call evaluate_tetra_lagrange_solution( &
        vertices - step*vertices_dot, tetrahedra, degree, &
        solution - step*solution_dot, 1, reference_point, value_minus, &
        gradient_minus, status)
    call check(status == 0, "tetrahedral H1 observation JVP succeeds")
    call check(abs( &
        value_dot - (value_plus - value_minus)/(2.0_dp*step)) < 3.0e-9_dp, &
        "tetrahedral H1 value JVP matches full re-evaluation")
    call check(maxval(abs( &
        gradient_dot - (gradient_plus - gradient_minus)/(2.0_dp*step))) < &
        4.0e-9_dp, "tetrahedral H1 gradient JVP matches full re-evaluation")
    call evaluate_tetra_lagrange_solution( &
        vertices, tetrahedra, degree, solution, 1, reference_point, value, &
        gradient, status)
    call evaluate_tetra_lagrange_solution_vjp( &
        vertices, tetrahedra, degree, solution, 1, reference_point, &
        value_bar, gradient_bar, vertices_bar, solution_bar, status)
    lhs = value_bar*value_dot + dot_product(gradient_bar, gradient_dot)
    rhs = sum(vertices_bar*vertices_dot) + &
        dot_product(solution_bar, solution_dot)
    call check(status == 0, "tetrahedral H1 observation VJP succeeds")
    call check(abs(lhs - rhs) < 3.0e-13_dp*max(1.0_dp, abs(lhs), abs(rhs)), &
        "tetrahedral H1 observation products obey adjoint identity")

    if (failures > 0) then
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

end program test_tetra_lagrange_solution_ad
