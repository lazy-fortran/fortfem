program test_tetra_lagrange_sampled_state_ad
    use, intrinsic :: iso_fortran_env, only: error_unit
    use fortfem_api, only: solve_tetra_lagrange_sampled_state, &
        solve_tetra_lagrange_sampled_state_jvp, &
        solve_tetra_lagrange_sampled_state_vjp, tetra_duffy_quadrature
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
    real(dp), allocatable :: eta(:), source_gradients(:, :, :)
    real(dp), allocatable :: source_parameter_dot(:, :)
    real(dp), allocatable :: source_values(:, :), source_values_bar(:, :)
    real(dp), allocatable :: source_values_minus(:, :)
    real(dp), allocatable :: source_values_plus(:, :), weights(:), xi(:)
    real(dp), allocatable :: zeta(:)
    real(dp) :: solution(5), solution_bar(5), solution_dot(5)
    real(dp) :: solution_minus(5), solution_plus(5)
    real(dp) :: stiffness_coefficient, stiffness_coefficient_bar
    real(dp) :: stiffness_coefficient_dot
    real(dp) :: mass_coefficient, mass_coefficient_bar, mass_coefficient_dot
    real(dp) :: lhs, parameter, parameter_dot, rhs
    integer :: failures, status

    failures = 0
    vertices = reshape([ &
        0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, &
        0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, &
        1.0_dp, 1.0_dp, 1.0_dp], [3, 5])
    vertices_dot = reshape([ &
        0.01_dp, -0.02_dp, 0.015_dp, -0.03_dp, 0.01_dp, 0.02_dp, &
        0.025_dp, -0.015_dp, 0.005_dp, -0.01_dp, 0.03_dp, -0.02_dp, &
        0.02_dp, -0.01_dp, 0.025_dp], [3, 5])
    constrained = [.true., .false., .false., .false., .true.]
    constrained_values = [0.4_dp, 0.0_dp, 0.0_dp, 0.0_dp, -0.3_dp]
    constrained_values_dot = [-0.06_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.09_dp]
    solution_bar = [0.3_dp, -0.4_dp, 0.25_dp, -0.15_dp, 0.2_dp]
    stiffness_coefficient = 1.7_dp
    mass_coefficient = 0.8_dp
    stiffness_coefficient_dot = 0.13_dp
    mass_coefficient_dot = -0.09_dp
    parameter = 0.7_dp
    parameter_dot = -0.11_dp
    call tetra_duffy_quadrature( &
        quadrature_degree, xi, eta, zeta, weights, status)
    call check(status == 0, "sampled state quadrature succeeds")
    allocate(source_values(size(weights), size(tetrahedra, 2)))
    allocate(source_values_plus, mold=source_values)
    allocate(source_values_minus, mold=source_values)
    allocate(source_parameter_dot, mold=source_values)
    allocate(source_gradients(3, size(weights), size(tetrahedra, 2)))
    call sample_source(vertices, parameter, source_values)
    call sample_source( &
        vertices + step*vertices_dot, parameter + step*parameter_dot, &
        source_values_plus)
    call sample_source( &
        vertices - step*vertices_dot, parameter - step*parameter_dot, &
        source_values_minus)
    source_parameter_dot = parameter_dot
    source_gradients = spread(spread([1.0_dp, 2.0_dp, -0.5_dp], 2, &
        size(weights)), 3, size(tetrahedra, 2))

    call solve_tetra_lagrange_sampled_state( &
        vertices, tetrahedra, degree, quadrature_degree, &
        stiffness_coefficient, mass_coefficient, source_values, constrained, &
        constrained_values, solution, status)
    call check(status == 0, "sampled H1 primal state succeeds")
    call solve_tetra_lagrange_sampled_state_jvp( &
        vertices, tetrahedra, degree, quadrature_degree, &
        stiffness_coefficient, mass_coefficient, source_values, &
        source_gradients, constrained, constrained_values, vertices_dot, &
        stiffness_coefficient_dot, mass_coefficient_dot, &
        source_parameter_dot, constrained_values_dot, solution_dot, status)
    call check(status == 0, "sampled H1 state JVP succeeds")
    call solve_tetra_lagrange_sampled_state( &
        vertices + step*vertices_dot, tetrahedra, degree, quadrature_degree, &
        stiffness_coefficient + step*stiffness_coefficient_dot, &
        mass_coefficient + step*mass_coefficient_dot, source_values_plus, &
        constrained, constrained_values + step*constrained_values_dot, &
        solution_plus, status)
    call solve_tetra_lagrange_sampled_state( &
        vertices - step*vertices_dot, tetrahedra, degree, quadrature_degree, &
        stiffness_coefficient - step*stiffness_coefficient_dot, &
        mass_coefficient - step*mass_coefficient_dot, source_values_minus, &
        constrained, constrained_values - step*constrained_values_dot, &
        solution_minus, status)
    call check(maxval(abs( &
        solution_dot - (solution_plus - solution_minus)/(2.0_dp*step))) < &
        2.0e-8_dp, "sampled H1 JVP matches full Eulerian reassembly")

    allocate(source_values_bar, mold=source_values)
    call solve_tetra_lagrange_sampled_state_vjp( &
        vertices, tetrahedra, degree, quadrature_degree, &
        stiffness_coefficient, mass_coefficient, source_values, &
        source_gradients, constrained, constrained_values, solution, &
        solution_bar, vertices_bar, stiffness_coefficient_bar, &
        mass_coefficient_bar, source_values_bar, constrained_values_bar, &
        status)
    call check(status == 0, "sampled H1 state VJP succeeds")
    lhs = dot_product(solution_bar, solution_dot)
    rhs = sum(vertices_bar*vertices_dot) + &
        stiffness_coefficient_bar*stiffness_coefficient_dot + &
        mass_coefficient_bar*mass_coefficient_dot + &
        sum(source_values_bar*source_parameter_dot) + &
        dot_product(constrained_values_bar, constrained_values_dot)
    call check(abs(lhs - rhs) < 3.0e-10_dp*max(1.0_dp, abs(lhs), abs(rhs)), &
        "sampled H1 state products satisfy the Eulerian adjoint identity")

    if (failures > 0) then
        write(error_unit, "(i0,a)") failures, " test(s) failed"
        stop 1
    end if
    write(*, "(a)") "PASS"

contains

    subroutine sample_source(mesh_vertices, source_parameter, samples)
        real(dp), intent(in) :: mesh_vertices(:, :), source_parameter
        real(dp), intent(out) :: samples(:, :)

        real(dp) :: jacobian(3, 3), point(3), cell_vertices(3, 4)
        integer :: cell, node, quadrature_point

        do cell = 1, size(tetrahedra, 2)
            do node = 1, 4
                cell_vertices(:, node) = &
                    mesh_vertices(:, tetrahedra(node, cell))
            end do
            jacobian(:, 1) = cell_vertices(:, 2) - cell_vertices(:, 1)
            jacobian(:, 2) = cell_vertices(:, 3) - cell_vertices(:, 1)
            jacobian(:, 3) = cell_vertices(:, 4) - cell_vertices(:, 1)
            do quadrature_point = 1, size(weights)
                point = cell_vertices(:, 1) + matmul(jacobian, [ &
                    xi(quadrature_point), eta(quadrature_point), &
                    zeta(quadrature_point)])
                samples(quadrature_point, cell) = source_parameter + &
                    point(1) + 2.0_dp*point(2) - 0.5_dp*point(3)
            end do
        end do
    end subroutine sample_source

    subroutine check(condition, label)
        logical, intent(in) :: condition
        character(*), intent(in) :: label

        if (condition) return
        failures = failures + 1
        write(error_unit, "(a)") "FAIL: "//label
    end subroutine check

end program test_tetra_lagrange_sampled_state_ad
