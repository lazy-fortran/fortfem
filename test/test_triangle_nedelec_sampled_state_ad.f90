program test_triangle_nedelec_sampled_state_ad
    use, intrinsic :: iso_fortran_env, only: error_unit
    use fortfem_api, only: &
        assemble_triangle_nedelec_vector_load_samples, &
        solve_triangle_nedelec_sampled_state, &
        solve_triangle_nedelec_sampled_state_jvp, &
        solve_triangle_nedelec_sampled_state_vjp, &
        triangle_duffy_quadrature
    use fortfem_kinds, only: dp
    use fortfem_mesh_2d, only: mesh_2d_t
    use fortsparse, only: fortsparse_status_t
    implicit none

    integer, parameter :: order = 1, quadrature_degree = 4
    real(dp), parameter :: step = 2.0e-7_dp
    type(mesh_2d_t) :: mesh, minus_mesh, plus_mesh
    real(dp) :: vertices_dot(2, 4), vertices_bar(2, 4)
    real(dp), allocatable :: state(:), state_bar(:), state_dot(:)
    real(dp), allocatable :: state_minus(:), state_plus(:)
    real(dp), allocatable :: source(:,:,:), source_bar(:,:,:)
    real(dp), allocatable :: source_dot(:,:,:), source_gradients(:,:,:,:)
    real(dp), allocatable :: source_minus(:,:,:), source_plus(:,:,:)
    real(dp), allocatable :: constrained_values(:)
    real(dp), allocatable :: constrained_values_bar(:)
    real(dp), allocatable :: constrained_values_dot(:), load(:)
    logical, allocatable :: constrained(:)
    real(dp), allocatable :: eta(:), weights(:), xi(:)
    type(fortsparse_status_t) :: sparse_status
    real(dp) :: curl_coefficient, curl_coefficient_bar, curl_coefficient_dot
    real(dp) :: mass_tensor(2, 2), mass_tensor_bar(2, 2)
    real(dp) :: mass_tensor_dot(2, 2)
    real(dp) :: lhs, point(2), point_dot(2), rhs
    integer :: cell, dof, failures, point_index, quadrature_status, status

    failures = 0
    call mesh%create_rectangular(2, 2, 0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp)
    call plus_mesh%create_rectangular(2, 2, 0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp)
    call minus_mesh%create_rectangular(2, 2, 0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp)
    vertices_dot = reshape([ &
        0.01_dp, -0.02_dp, -0.03_dp, 0.015_dp, &
        0.025_dp, -0.01_dp, 0.02_dp, 0.03_dp], [2, 4])
    plus_mesh%vertices = plus_mesh%vertices + step*vertices_dot
    minus_mesh%vertices = minus_mesh%vertices - step*vertices_dot
    curl_coefficient = 1.4_dp
    curl_coefficient_dot = -0.11_dp
    mass_tensor = reshape([1.3_dp, 0.12_dp, 0.12_dp, 0.9_dp], [2, 2])
    mass_tensor_dot = reshape([-0.08_dp, 0.03_dp, 0.03_dp, 0.06_dp], [2, 2])
    call triangle_duffy_quadrature( &
        quadrature_degree, xi, eta, weights, quadrature_status)
    call check(quadrature_status == 0, "triangle quadrature succeeds")
    allocate(source(2, size(weights), 2), source_dot(2, size(weights), 2))
    allocate(source_minus(2, size(weights), 2))
    allocate(source_plus(2, size(weights), 2))
    allocate(source_gradients(2, 2, size(weights), 2))
    do cell = 1, 2
        do point_index = 1, size(weights)
            source_gradients(:, :, point_index, cell) = reshape([ &
                0.2_dp*cell, -0.1_dp, 0.15_dp, 0.25_dp*cell], [2, 2])
            call physical_sample( &
                mesh%vertices, mesh%triangles, cell, point_index, point)
            call physical_sample( &
                vertices_dot, mesh%triangles, cell, point_index, point_dot)
            source(:, point_index, cell) = [0.4_dp, -0.3_dp*cell] + &
                matmul(source_gradients(:, :, point_index, cell), point)
            source_dot(:, point_index, cell) = [ &
                0.01_dp*cell, -0.012_dp*cell] + 0.0004_dp*point_index
            source_plus(:, point_index, cell) = source(:, point_index, cell) + &
                step*(source_dot(:, point_index, cell) + matmul( &
                source_gradients(:, :, point_index, cell), point_dot))
            source_minus(:, point_index, cell) = source(:, point_index, cell) - &
                step*(source_dot(:, point_index, cell) + matmul( &
                source_gradients(:, :, point_index, cell), point_dot))
        end do
    end do
    call assemble_triangle_nedelec_vector_load_samples( &
        mesh, order, quadrature_degree, source, load, sparse_status)
    call check(sparse_status%code == 0, "sampled load sizes the state space")
    allocate(constrained(size(load)), source=.false.)
    allocate(constrained_values(size(load)), source=0.0_dp)
    allocate(constrained_values_dot(size(load)), source=0.0_dp)
    constrained(1) = .true.
    constrained_values(1) = 0.17_dp
    constrained_values_dot(1) = -0.025_dp
    allocate(state(size(load)), state_dot(size(load)))

    call solve_triangle_nedelec_sampled_state( &
        mesh, order, quadrature_degree, curl_coefficient, mass_tensor, source, &
        constrained, constrained_values, state, status)
    call check(status == 0, "sampled triangle Nedelec state succeeds")
    call solve_triangle_nedelec_sampled_state_jvp( &
        mesh, order, quadrature_degree, curl_coefficient, mass_tensor, source, &
        source_gradients, constrained, constrained_values, vertices_dot, &
        curl_coefficient_dot, mass_tensor_dot, source_dot, &
        constrained_values_dot, state_dot, status)
    call check(status == 0, "sampled triangle Nedelec state JVP succeeds")
    allocate(state_plus(size(state)), state_minus(size(state)))
    call solve_triangle_nedelec_sampled_state( &
        plus_mesh, order, quadrature_degree, &
        curl_coefficient + step*curl_coefficient_dot, &
        mass_tensor + step*mass_tensor_dot, source_plus, constrained, &
        constrained_values + step*constrained_values_dot, state_plus, status)
    call solve_triangle_nedelec_sampled_state( &
        minus_mesh, order, quadrature_degree, &
        curl_coefficient - step*curl_coefficient_dot, &
        mass_tensor - step*mass_tensor_dot, source_minus, constrained, &
        constrained_values - step*constrained_values_dot, state_minus, status)
    call check(maxval(abs( &
        state_dot - (state_plus - state_minus)/(2.0_dp*step))) < 5.0e-8_dp, &
        "sampled triangle Nedelec state JVP matches full reassembly")

    allocate(state_bar(size(state)), source_bar(2, size(weights), 2))
    allocate(constrained_values_bar(size(state)))
    state_bar = [(0.03_dp*dof - 0.11_dp, dof = 1, size(state))]
    call solve_triangle_nedelec_sampled_state_vjp( &
        mesh, order, quadrature_degree, curl_coefficient, mass_tensor, source, &
        source_gradients, constrained, constrained_values, state, state_bar, &
        vertices_bar, curl_coefficient_bar, mass_tensor_bar, source_bar, &
        constrained_values_bar, status)
    call check(status == 0, "sampled triangle Nedelec state VJP succeeds")
    lhs = dot_product(state_bar, state_dot)
    rhs = sum(vertices_bar*vertices_dot) + &
        curl_coefficient_bar*curl_coefficient_dot + &
        sum(mass_tensor_bar*mass_tensor_dot) + sum(source_bar*source_dot) + &
        dot_product(constrained_values_bar, constrained_values_dot)
    call check(abs(lhs - rhs) < 8.0e-10_dp*max(1.0_dp, abs(lhs), abs(rhs)), &
        "sampled triangle Nedelec state products obey the adjoint identity")

    if (failures > 0) then
        write(error_unit, "(i0,a)") failures, " test(s) failed"
        stop 1
    end if
    write(*, "(a)") "PASS"

contains

    subroutine physical_sample( &
            coordinates, triangles, cell, point_index, point)
        real(dp), intent(in) :: coordinates(:, :)
        integer, intent(in) :: triangles(:, :), cell, point_index
        real(dp), intent(out) :: point(2)

        real(dp) :: jacobian(2, 2)

        jacobian(:, 1) = coordinates(:, triangles(2, cell)) - &
            coordinates(:, triangles(1, cell))
        jacobian(:, 2) = coordinates(:, triangles(3, cell)) - &
            coordinates(:, triangles(1, cell))
        point = coordinates(:, triangles(1, cell)) + &
            matmul(jacobian, [xi(point_index), eta(point_index)])
    end subroutine physical_sample

    subroutine check(condition, label)
        logical, intent(in) :: condition
        character(*), intent(in) :: label

        if (condition) return
        failures = failures + 1
        write(error_unit, "(a)") "FAIL: "//label
    end subroutine check

end program test_triangle_nedelec_sampled_state_ad
