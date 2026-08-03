program test_tetra_lagrange_state_ad
    use, intrinsic :: iso_fortran_env, only: error_unit
    use fortfem_feec, only: solve_tetra_lagrange_state, &
        solve_tetra_lagrange_state_jvp, solve_tetra_lagrange_state_vjp
    use fortfem_kinds, only: dp
    implicit none

    integer, parameter :: degree = 1, quadrature_degree = 4
    integer, parameter :: tetrahedra(4, 2) = reshape([ &
        1, 2, 3, 4, 2, 3, 4, 5], [4, 2])
    real(dp), parameter :: step = 2.0e-7_dp
    logical :: constrained(5)
    real(dp) :: vertices(3, 5), vertices_dot(3, 5), vertices_bar(3, 5)
    real(dp) :: constrained_values(5), constrained_values_bar(5)
    real(dp) :: constrained_values_dot(5)
    real(dp) :: load(5), load_bar(5), load_dot(5)
    real(dp) :: solution(5), solution_bar(5), solution_dot(5)
    real(dp) :: solution_minus(5), solution_plus(5)
    real(dp) :: stiffness_coefficient, stiffness_coefficient_bar
    real(dp) :: stiffness_coefficient_dot
    real(dp) :: mass_coefficient, mass_coefficient_bar, mass_coefficient_dot
    real(dp) :: lhs, rhs
    integer :: failures, status

    failures = 0
    vertices = reshape([ &
        0.0_dp, 0.0_dp, 0.0_dp, &
        1.0_dp, 0.0_dp, 0.0_dp, &
        0.0_dp, 1.0_dp, 0.0_dp, &
        0.0_dp, 0.0_dp, 1.0_dp, &
        1.0_dp, 1.0_dp, 1.0_dp], [3, 5])
    vertices_dot = reshape([ &
        0.01_dp, -0.02_dp, 0.015_dp, &
        -0.03_dp, 0.01_dp, 0.02_dp, &
        0.025_dp, -0.015_dp, 0.005_dp, &
        -0.01_dp, 0.03_dp, -0.02_dp, &
        0.02_dp, -0.01_dp, 0.025_dp], [3, 5])
    constrained = [.true., .false., .false., .false., .true.]
    constrained_values = [0.4_dp, 0.0_dp, 0.0_dp, 0.0_dp, -0.3_dp]
    constrained_values_dot = [-0.06_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.09_dp]
    load = [0.2_dp, 1.1_dp, -0.7_dp, 0.5_dp, 0.3_dp]
    load_dot = [-0.03_dp, 0.08_dp, 0.04_dp, -0.02_dp, 0.06_dp]
    solution_bar = [0.3_dp, -0.4_dp, 0.25_dp, -0.15_dp, 0.2_dp]
    stiffness_coefficient = 1.7_dp
    mass_coefficient = 0.8_dp
    stiffness_coefficient_dot = 0.13_dp
    mass_coefficient_dot = -0.09_dp

    call solve_tetra_lagrange_state( &
        vertices, tetrahedra, degree, quadrature_degree, &
        stiffness_coefficient, mass_coefficient, load, constrained, &
        constrained_values, solution, status)
    call check(status == 0, "tetra H1 constrained primal state succeeds")
    call solve_tetra_lagrange_state_jvp( &
        vertices, tetrahedra, degree, quadrature_degree, &
        stiffness_coefficient, mass_coefficient, load, constrained, &
        constrained_values, vertices_dot, stiffness_coefficient_dot, &
        mass_coefficient_dot, load_dot, constrained_values_dot, solution_dot, &
        status)
    call check(status == 0, "tetra H1 constrained JVP succeeds")
    call solve_tetra_lagrange_state( &
        vertices + step*vertices_dot, tetrahedra, degree, quadrature_degree, &
        stiffness_coefficient + step*stiffness_coefficient_dot, &
        mass_coefficient + step*mass_coefficient_dot, load + step*load_dot, &
        constrained, constrained_values + step*constrained_values_dot, &
        solution_plus, status)
    call solve_tetra_lagrange_state( &
        vertices - step*vertices_dot, tetrahedra, degree, quadrature_degree, &
        stiffness_coefficient - step*stiffness_coefficient_dot, &
        mass_coefficient - step*mass_coefficient_dot, load - step*load_dot, &
        constrained, constrained_values - step*constrained_values_dot, &
        solution_minus, status)
    call check(maxval(abs( &
        solution_dot - (solution_plus - solution_minus)/(2.0_dp*step))) < &
        2.0e-8_dp, &
        "tetra H1 JVP matches independent full reassembly")

    call solve_tetra_lagrange_state_vjp( &
        vertices, tetrahedra, degree, quadrature_degree, &
        stiffness_coefficient, mass_coefficient, load, constrained, &
        constrained_values, solution, solution_bar, vertices_bar, &
        stiffness_coefficient_bar, mass_coefficient_bar, load_bar, &
        constrained_values_bar, status)
    call check(status == 0, "tetra H1 constrained VJP succeeds")
    lhs = dot_product(solution_bar, solution_dot)
    rhs = sum(vertices_bar*vertices_dot) + &
        stiffness_coefficient_bar*stiffness_coefficient_dot + &
        mass_coefficient_bar*mass_coefficient_dot + &
        dot_product(load_bar, load_dot) + &
        dot_product(constrained_values_bar, constrained_values_dot)
    call check(abs(lhs - rhs) < 3.0e-10_dp*max(1.0_dp, abs(lhs), abs(rhs)), &
        "tetra H1 products accumulate the moving-mesh Dirichlet adjoint")

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

end program test_tetra_lagrange_state_ad
