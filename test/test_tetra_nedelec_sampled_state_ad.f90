program test_tetra_nedelec_sampled_state_ad
    use, intrinsic :: iso_fortran_env, only: error_unit
    use fortfem_feec, only: solve_tetra_nedelec_sampled_state, &
        solve_tetra_nedelec_sampled_state_jvp, &
        solve_tetra_nedelec_sampled_state_vjp, tetra_duffy_quadrature
    use fortfem_kinds, only: dp
    implicit none

    integer, parameter :: order = 1, quadrature_degree = 4, dof_count = 9
    integer, parameter :: tetrahedra(4, 2) = reshape([ &
        1, 2, 3, 4, 2, 3, 4, 5], [4, 2])
    real(dp), parameter :: step = 2.0e-7_dp
    logical :: constrained(dof_count)
    real(dp) :: vertices(3, 5), vertices_dot(3, 5), vertices_bar(3, 5)
    real(dp) :: constrained_values(dof_count), constrained_values_bar(dof_count)
    real(dp) :: constrained_values_dot(dof_count)
    real(dp), allocatable :: eta(:), source_gradients(:, :, :, :)
    real(dp), allocatable :: source_parameter_dot(:, :, :)
    real(dp), allocatable :: source_values(:, :, :), source_values_bar(:, :, :)
    real(dp), allocatable :: source_values_minus(:, :, :)
    real(dp), allocatable :: source_values_plus(:, :, :), weights(:), xi(:)
    real(dp), allocatable :: zeta(:)
    real(dp) :: state(dof_count), state_bar(dof_count), state_dot(dof_count)
    real(dp) :: state_minus(dof_count), state_plus(dof_count)
    real(dp) :: curl_coefficient, curl_coefficient_bar, curl_coefficient_dot
    real(dp) :: mass_coefficient, mass_coefficient_bar, mass_coefficient_dot
    real(dp) :: lhs, parameter, parameter_dot, rhs
    integer :: dof, failures, status

    failures = 0
    vertices = reshape([ &
        0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, &
        0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, &
        1.0_dp, 1.0_dp, 1.0_dp], [3, 5])
    vertices_dot = reshape([ &
        0.01_dp, -0.02_dp, 0.015_dp, -0.03_dp, 0.01_dp, 0.02_dp, &
        0.025_dp, -0.015_dp, 0.005_dp, -0.01_dp, 0.03_dp, -0.02_dp, &
        0.02_dp, -0.01_dp, 0.025_dp], [3, 5])
    constrained = .false.
    constrained([1, 9]) = .true.
    constrained_values = 0.0_dp
    constrained_values([1, 9]) = [0.2_dp, -0.15_dp]
    constrained_values_dot = 0.0_dp
    constrained_values_dot([1, 9]) = [-0.04_dp, 0.06_dp]
    do dof = 1, dof_count
        state_bar(dof) = 0.03_dp*dof - 0.002_dp*dof**2
    end do
    curl_coefficient = 1.7_dp
    mass_coefficient = 0.8_dp
    curl_coefficient_dot = 0.13_dp
    mass_coefficient_dot = -0.09_dp
    parameter = 0.7_dp
    parameter_dot = -0.11_dp
    call tetra_duffy_quadrature( &
        quadrature_degree, xi, eta, zeta, weights, status)
    call check(status == 0, "Maxwell state quadrature succeeds")
    allocate(source_values(3, size(weights), size(tetrahedra, 2)))
    allocate(source_values_plus, mold=source_values)
    allocate(source_values_minus, mold=source_values)
    allocate(source_values_bar, mold=source_values)
    allocate(source_parameter_dot, mold=source_values)
    allocate(source_gradients(3, 3, size(weights), size(tetrahedra, 2)))
    call sample_source(vertices, parameter, source_values)
    call sample_source( &
        vertices + step*vertices_dot, parameter + step*parameter_dot, &
        source_values_plus)
    call sample_source( &
        vertices - step*vertices_dot, parameter - step*parameter_dot, &
        source_values_minus)
    source_parameter_dot = spread(spread( &
        [parameter_dot, 2.0_dp*parameter_dot, -parameter_dot], 2, &
        size(weights)), 3, size(tetrahedra, 2))
    call fill_source_gradients(source_gradients)

    call solve_tetra_nedelec_sampled_state( &
        vertices, tetrahedra, order, quadrature_degree, curl_coefficient, &
        mass_coefficient, source_values, constrained, constrained_values, &
        state, status)
    call check(status == 0, "sampled Maxwell primal state succeeds")
    call solve_tetra_nedelec_sampled_state_jvp( &
        vertices, tetrahedra, order, quadrature_degree, curl_coefficient, &
        mass_coefficient, source_values, source_gradients, constrained, &
        constrained_values, vertices_dot, curl_coefficient_dot, &
        mass_coefficient_dot, source_parameter_dot, constrained_values_dot, &
        state_dot, status)
    call check(status == 0, "sampled Maxwell state JVP succeeds")
    call solve_tetra_nedelec_sampled_state( &
        vertices + step*vertices_dot, tetrahedra, order, quadrature_degree, &
        curl_coefficient + step*curl_coefficient_dot, &
        mass_coefficient + step*mass_coefficient_dot, source_values_plus, &
        constrained, constrained_values + step*constrained_values_dot, &
        state_plus, status)
    call solve_tetra_nedelec_sampled_state( &
        vertices - step*vertices_dot, tetrahedra, order, quadrature_degree, &
        curl_coefficient - step*curl_coefficient_dot, &
        mass_coefficient - step*mass_coefficient_dot, source_values_minus, &
        constrained, constrained_values - step*constrained_values_dot, &
        state_minus, status)
    call check(maxval(abs( &
        state_dot - (state_plus - state_minus)/(2.0_dp*step))) < 3.0e-8_dp, &
        "sampled Maxwell JVP matches full Eulerian reassembly")

    call solve_tetra_nedelec_sampled_state_vjp( &
        vertices, tetrahedra, order, quadrature_degree, curl_coefficient, &
        mass_coefficient, source_values, source_gradients, constrained, &
        constrained_values, state, state_bar, vertices_bar, &
        curl_coefficient_bar, mass_coefficient_bar, source_values_bar, &
        constrained_values_bar, status)
    call check(status == 0, "sampled Maxwell state VJP succeeds")
    lhs = dot_product(state_bar, state_dot)
    rhs = sum(vertices_bar*vertices_dot) + &
        curl_coefficient_bar*curl_coefficient_dot + &
        mass_coefficient_bar*mass_coefficient_dot + &
        sum(source_values_bar*source_parameter_dot) + &
        dot_product(constrained_values_bar, constrained_values_dot)
    call check(abs(lhs - rhs) < 5.0e-10_dp*max(1.0_dp, abs(lhs), abs(rhs)), &
        "sampled Maxwell products satisfy the Eulerian adjoint identity")

    if (failures > 0) then
        write(error_unit, "(i0,a)") failures, " test(s) failed"
        stop 1
    end if
    write(*, "(a)") "PASS"

contains

    subroutine sample_source(mesh_vertices, source_parameter, samples)
        real(dp), intent(in) :: mesh_vertices(:, :), source_parameter
        real(dp), intent(out) :: samples(:, :, :)

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
                samples(:, quadrature_point, cell) = [ &
                    source_parameter + point(1), &
                    2.0_dp*source_parameter - point(2) + point(3), &
                    -source_parameter + 0.5_dp*point(1) + 2.0_dp*point(3)]
            end do
        end do
    end subroutine sample_source

    subroutine fill_source_gradients(gradients)
        real(dp), intent(out) :: gradients(:, :, :, :)

        integer :: cell, quadrature_point

        do cell = 1, size(gradients, 4)
            do quadrature_point = 1, size(gradients, 3)
                gradients(:, :, quadrature_point, cell) = reshape([ &
                    1.0_dp, 0.0_dp, 0.5_dp, 0.0_dp, -1.0_dp, 0.0_dp, &
                    0.0_dp, 1.0_dp, 2.0_dp], [3, 3])
            end do
        end do
    end subroutine fill_source_gradients

    subroutine check(condition, label)
        logical, intent(in) :: condition
        character(*), intent(in) :: label

        if (condition) return
        failures = failures + 1
        write(error_unit, "(a)") "FAIL: "//label
    end subroutine check

end program test_tetra_nedelec_sampled_state_ad
