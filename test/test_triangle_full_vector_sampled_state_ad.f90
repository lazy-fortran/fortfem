program test_triangle_full_vector_sampled_state_ad
    use, intrinsic :: iso_fortran_env, only: error_unit
    use fortfem_feec, only: &
        assemble_triangle_bdm_vector_load_samples, &
        solve_triangle_bdm_sampled_state, &
        solve_triangle_bdm_sampled_state_jvp, &
        solve_triangle_bdm_sampled_state_vjp, &
        solve_triangle_nedelec_second_sampled_state, &
        solve_triangle_nedelec_second_sampled_state_jvp, &
        solve_triangle_nedelec_second_sampled_state_vjp, &
        triangle_duffy_quadrature
    use fortfem_kinds, only: dp
    use fortfem_mesh_2d, only: mesh_2d_t
    use fortsparse, only: fortsparse_status_t
    implicit none

    integer, parameter :: degree = 1, quadrature_degree = 4
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
    real(dp) :: derivative_coefficient, derivative_coefficient_bar
    real(dp) :: derivative_coefficient_dot
    real(dp) :: mass_coefficient, mass_coefficient_bar, mass_coefficient_dot
    real(dp) :: lhs, point(2), point_dot(2), rhs
    integer :: cell, dof, failures, point_index, quadrature_status, status

    failures = 0
    vertices_dot = reshape([ &
        0.01_dp, -0.02_dp, -0.03_dp, 0.015_dp, &
        0.025_dp, -0.01_dp, 0.02_dp, 0.03_dp], [2, 4])
    derivative_coefficient = 1.4_dp
    mass_coefficient = 0.9_dp
    derivative_coefficient_dot = -0.11_dp
    mass_coefficient_dot = 0.06_dp
    call triangle_duffy_quadrature( &
        quadrature_degree, xi, eta, weights, quadrature_status)
    call check(quadrature_status == 0, "triangle quadrature succeeds")
    allocate(source(2, size(weights), 2), source_dot(2, size(weights), 2))
    allocate(source_minus(2, size(weights), 2))
    allocate(source_plus(2, size(weights), 2))
    allocate(source_gradients(2, 2, size(weights), 2))

    call check_family(.true., "BDM")
    call check_family(.false., "second-kind Nedelec")
    if (failures > 0) then
        write(error_unit, "(i0,a)") failures, " test(s) failed"
        stop 1
    end if
    write(*, "(a)") "PASS"

contains

    subroutine check_family(normal_family, label)
        logical, intent(in) :: normal_family
        character(*), intent(in) :: label

        call initialize_case()
        call initialize_state_arrays()
        if (normal_family) then
            call solve_triangle_bdm_sampled_state( &
                mesh, degree, quadrature_degree, derivative_coefficient, &
                mass_coefficient, source, constrained, constrained_values, &
                state, status)
            call solve_triangle_bdm_sampled_state_jvp( &
                mesh, degree, quadrature_degree, derivative_coefficient, &
                mass_coefficient, source, source_gradients, constrained, &
                constrained_values, vertices_dot, derivative_coefficient_dot, &
                mass_coefficient_dot, source_dot, constrained_values_dot, &
                state_dot, status)
            call solve_triangle_bdm_sampled_state( &
                plus_mesh, degree, quadrature_degree, &
                derivative_coefficient + step*derivative_coefficient_dot, &
                mass_coefficient + step*mass_coefficient_dot, source_plus, &
                constrained, constrained_values + &
                step*constrained_values_dot, state_plus, status)
            call solve_triangle_bdm_sampled_state( &
                minus_mesh, degree, quadrature_degree, &
                derivative_coefficient - step*derivative_coefficient_dot, &
                mass_coefficient - step*mass_coefficient_dot, source_minus, &
                constrained, constrained_values - &
                step*constrained_values_dot, state_minus, status)
        else
            call solve_triangle_nedelec_second_sampled_state( &
                mesh, degree, quadrature_degree, derivative_coefficient, &
                mass_coefficient, source, constrained, constrained_values, &
                state, status)
            call solve_triangle_nedelec_second_sampled_state_jvp( &
                mesh, degree, quadrature_degree, derivative_coefficient, &
                mass_coefficient, source, source_gradients, constrained, &
                constrained_values, vertices_dot, derivative_coefficient_dot, &
                mass_coefficient_dot, source_dot, constrained_values_dot, &
                state_dot, status)
            call solve_triangle_nedelec_second_sampled_state( &
                plus_mesh, degree, quadrature_degree, &
                derivative_coefficient + step*derivative_coefficient_dot, &
                mass_coefficient + step*mass_coefficient_dot, source_plus, &
                constrained, constrained_values + &
                step*constrained_values_dot, state_plus, status)
            call solve_triangle_nedelec_second_sampled_state( &
                minus_mesh, degree, quadrature_degree, &
                derivative_coefficient - step*derivative_coefficient_dot, &
                mass_coefficient - step*mass_coefficient_dot, source_minus, &
                constrained, constrained_values - &
                step*constrained_values_dot, state_minus, status)
        end if
        call check(status == 0, label//" sampled state products succeed")
        call check(maxval(abs( &
            state_dot - (state_plus - state_minus)/(2.0_dp*step))) < &
            6.0e-8_dp, label//" state JVP matches full reassembly")
        state_bar = [(0.03_dp*dof - 0.11_dp, dof = 1, size(state))]
        if (normal_family) then
            call solve_triangle_bdm_sampled_state_vjp( &
                mesh, degree, quadrature_degree, derivative_coefficient, &
                mass_coefficient, source, source_gradients, constrained, &
                constrained_values, state, state_bar, vertices_bar, &
                derivative_coefficient_bar, mass_coefficient_bar, source_bar, &
                constrained_values_bar, status)
        else
            call solve_triangle_nedelec_second_sampled_state_vjp( &
                mesh, degree, quadrature_degree, derivative_coefficient, &
                mass_coefficient, source, source_gradients, constrained, &
                constrained_values, state, state_bar, vertices_bar, &
                derivative_coefficient_bar, mass_coefficient_bar, source_bar, &
                constrained_values_bar, status)
        end if
        lhs = dot_product(state_bar, state_dot)
        rhs = sum(vertices_bar*vertices_dot) + &
            derivative_coefficient_bar*derivative_coefficient_dot + &
            mass_coefficient_bar*mass_coefficient_dot + &
            sum(source_bar*source_dot) + &
            dot_product(constrained_values_bar, constrained_values_dot)
        call check(status == 0, label//" sampled state VJP succeeds")
        call check( &
            abs(lhs - rhs) < 9.0e-10_dp*max(1.0_dp, abs(lhs), abs(rhs)), &
            label//" state products obey the adjoint identity")
    end subroutine check_family

    subroutine initialize_case()
        call mesh%create_rectangular(2, 2, 0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp)
        call plus_mesh%create_rectangular( &
            2, 2, 0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp)
        call minus_mesh%create_rectangular( &
            2, 2, 0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp)
        plus_mesh%vertices = plus_mesh%vertices + step*vertices_dot
        minus_mesh%vertices = minus_mesh%vertices - step*vertices_dot
        do cell = 1, 2
            do point_index = 1, size(weights)
                source_gradients(:, :, point_index, cell) = reshape([ &
                    0.2_dp*cell, -0.1_dp, 0.15_dp, 0.25_dp*cell], [2, 2])
                call physical_sample( &
                    mesh%vertices, cell, point_index, point)
                call physical_sample( &
                    vertices_dot, cell, point_index, point_dot)
                source(:, point_index, cell) = [0.4_dp, -0.3_dp*cell] + &
                    matmul(source_gradients(:, :, point_index, cell), point)
                source_dot(:, point_index, cell) = [ &
                    0.01_dp*cell, -0.012_dp*cell] + 0.0004_dp*point_index
                source_plus(:, point_index, cell) = &
                    source(:, point_index, cell) + step*( &
                    source_dot(:, point_index, cell) + matmul( &
                    source_gradients(:, :, point_index, cell), point_dot))
                source_minus(:, point_index, cell) = &
                    source(:, point_index, cell) - step*( &
                    source_dot(:, point_index, cell) + matmul( &
                    source_gradients(:, :, point_index, cell), point_dot))
            end do
        end do
    end subroutine initialize_case

    subroutine initialize_state_arrays()
        if (allocated(load)) deallocate(load)
        call assemble_triangle_bdm_vector_load_samples( &
            mesh, degree, quadrature_degree, source, load, sparse_status)
        call check(sparse_status%code == 0, "sampled load sizes state space")
        if (allocated(constrained)) deallocate(constrained)
        if (allocated(constrained_values)) deallocate(constrained_values)
        if (allocated(constrained_values_dot)) &
            deallocate(constrained_values_dot)
        if (allocated(constrained_values_bar)) &
            deallocate(constrained_values_bar)
        if (allocated(state)) deallocate(state)
        if (allocated(state_dot)) deallocate(state_dot)
        if (allocated(state_plus)) deallocate(state_plus)
        if (allocated(state_minus)) deallocate(state_minus)
        if (allocated(state_bar)) deallocate(state_bar)
        if (allocated(source_bar)) deallocate(source_bar)
        allocate(constrained(size(load)), source=.false.)
        allocate(constrained_values(size(load)), source=0.0_dp)
        allocate(constrained_values_dot(size(load)), source=0.0_dp)
        allocate(constrained_values_bar(size(load)))
        allocate(state(size(load)), state_dot(size(load)))
        allocate(state_plus(size(load)), state_minus(size(load)))
        allocate(state_bar(size(load)), source_bar(2, size(weights), 2))
        constrained(1) = .true.
        constrained_values(1) = 0.17_dp
        constrained_values_dot(1) = -0.025_dp
    end subroutine initialize_state_arrays

    subroutine physical_sample(coordinates, cell, point_index, point)
        real(dp), intent(in) :: coordinates(:, :)
        integer, intent(in) :: cell, point_index
        real(dp), intent(out) :: point(2)

        real(dp) :: jacobian(2, 2)

        jacobian(:, 1) = coordinates(:, mesh%triangles(2, cell)) - &
            coordinates(:, mesh%triangles(1, cell))
        jacobian(:, 2) = coordinates(:, mesh%triangles(3, cell)) - &
            coordinates(:, mesh%triangles(1, cell))
        point = coordinates(:, mesh%triangles(1, cell)) + &
            matmul(jacobian, [xi(point_index), eta(point_index)])
    end subroutine physical_sample

    subroutine check(condition, label)
        logical, intent(in) :: condition
        character(*), intent(in) :: label

        if (condition) return
        failures = failures + 1
        write(error_unit, "(a)") "FAIL: "//label
    end subroutine check

end program test_triangle_full_vector_sampled_state_ad
