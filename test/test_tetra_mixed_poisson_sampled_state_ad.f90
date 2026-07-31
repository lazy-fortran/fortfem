program test_tetra_mixed_poisson_sampled_state_ad
    use, intrinsic :: iso_fortran_env, only: error_unit
    use fortfem_api, only: &
        solve_tetra_mixed_poisson_sampled_state, &
        solve_tetra_mixed_poisson_sampled_state_jvp, &
        solve_tetra_mixed_poisson_sampled_state_vjp
    use fortfem_kinds, only: dp
    use fortfem_tetra_duffy_quadrature, only: tetra_duffy_quadrature
    implicit none

    integer, parameter :: degree = 1, quadrature_degree = 4
    integer, parameter :: tetrahedra(4, 2) = reshape([ &
        1, 2, 3, 4, 2, 3, 4, 5], [4, 2])
    real(dp), parameter :: step = 2.0e-7_dp
    real(dp) :: vertices(3, 5), vertices_dot(3, 5), vertices_bar(3, 5)
    real(dp), allocatable :: flux(:), flux_bar(:), flux_dot(:)
    real(dp), allocatable :: flux_minus(:), flux_plus(:)
    real(dp), allocatable :: pressure(:), pressure_bar(:), pressure_dot(:)
    real(dp), allocatable :: pressure_minus(:), pressure_plus(:)
    real(dp), allocatable :: source(:,:), source_bar(:,:), source_dot(:,:)
    real(dp), allocatable :: source_gradients(:,:,:)
    real(dp), allocatable :: source_minus(:,:), source_plus(:,:)
    real(dp), allocatable :: weights(:), x(:), y(:), z(:)
    real(dp) :: coefficient, coefficient_bar, coefficient_dot
    real(dp) :: lhs, point(3), point_dot(3), rhs
    integer :: cell, dof, failures, point_index, quadrature_status, status

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
    coefficient = 1.6_dp
    coefficient_dot = -0.13_dp
    call tetra_duffy_quadrature( &
        quadrature_degree, x, y, z, weights, quadrature_status)
    call check(quadrature_status == 0, "tetra quadrature succeeds")
    allocate(source(size(weights), 2), source_dot(size(weights), 2))
    allocate(source_minus(size(weights), 2), source_plus(size(weights), 2))
    allocate(source_gradients(3, size(weights), 2))
    source_gradients(:, :, 1) = spread([0.5_dp, -0.2_dp, 0.3_dp], 2, &
        size(weights))
    source_gradients(:, :, 2) = spread([-0.4_dp, 0.1_dp, 0.6_dp], 2, &
        size(weights))
    do cell = 1, 2
        do point_index = 1, size(weights)
            call physical_sample(vertices, cell, point_index, point)
            call physical_sample(vertices_dot, cell, point_index, point_dot)
            source(point_index, cell) = 0.8_dp - 0.25_dp*cell + &
                dot_product(source_gradients(:, point_index, cell), point)
            source_dot(point_index, cell) = &
                0.02_dp*cell - 0.001_dp*point_index
            source_plus(point_index, cell) = source(point_index, cell) + &
                step*(source_dot(point_index, cell) + dot_product( &
                source_gradients(:, point_index, cell), point_dot))
            source_minus(point_index, cell) = source(point_index, cell) - &
                step*(source_dot(point_index, cell) + dot_product( &
                source_gradients(:, point_index, cell), point_dot))
        end do
    end do

    call solve_tetra_mixed_poisson_sampled_state( &
        vertices, tetrahedra, degree, quadrature_degree, coefficient, source, &
        flux, pressure, status)
    call check(status == 0, "sampled mixed Poisson state succeeds")
    allocate(flux_bar(size(flux)), pressure_bar(size(pressure)))
    flux_bar = [(0.015_dp*dof - 0.08_dp, dof = 1, size(flux))]
    pressure_bar = [(0.03_dp*dof - 0.09_dp, dof = 1, size(pressure))]
    call solve_tetra_mixed_poisson_sampled_state_jvp( &
        vertices, tetrahedra, degree, quadrature_degree, coefficient, source, &
        source_gradients, vertices_dot, coefficient_dot, source_dot, &
        flux_dot, pressure_dot, status)
    call check(status == 0, "sampled mixed Poisson state JVP succeeds")
    call solve_tetra_mixed_poisson_sampled_state( &
        vertices + step*vertices_dot, tetrahedra, degree, quadrature_degree, &
        coefficient + step*coefficient_dot, source_plus, flux_plus, &
        pressure_plus, status)
    call solve_tetra_mixed_poisson_sampled_state( &
        vertices - step*vertices_dot, tetrahedra, degree, quadrature_degree, &
        coefficient - step*coefficient_dot, source_minus, flux_minus, &
        pressure_minus, status)
    call check(maxval(abs( &
        flux_dot - (flux_plus - flux_minus)/(2.0_dp*step))) < 4.0e-8_dp, &
        "sampled mixed flux JVP matches full perturbed reassembly")
    call check(maxval(abs( &
        pressure_dot - (pressure_plus - pressure_minus)/(2.0_dp*step))) < &
        4.0e-8_dp, &
        "sampled mixed pressure JVP matches full perturbed reassembly")

    allocate(source_bar(size(source, 1), size(source, 2)))
    call solve_tetra_mixed_poisson_sampled_state_vjp( &
        vertices, tetrahedra, degree, quadrature_degree, coefficient, source, &
        source_gradients, flux, pressure, flux_bar, pressure_bar, &
        vertices_bar, coefficient_bar, source_bar, status)
    call check(status == 0, "sampled mixed Poisson state VJP succeeds")
    lhs = dot_product(flux_bar, flux_dot) + &
        dot_product(pressure_bar, pressure_dot)
    rhs = sum(vertices_bar*vertices_dot) + coefficient_bar*coefficient_dot + &
        sum(source_bar*source_dot)
    call check(abs(lhs - rhs) < 5.0e-10_dp*max(1.0_dp, abs(lhs), abs(rhs)), &
        "sampled mixed Poisson state products obey the adjoint identity")

    if (failures > 0) then
        write(error_unit, "(i0,a)") failures, " test(s) failed"
        stop 1
    end if
    write(*, "(a)") "PASS"

contains

    subroutine physical_sample(coordinates, cell, point_index, point)
        real(dp), intent(in) :: coordinates(:, :)
        integer, intent(in) :: cell, point_index
        real(dp), intent(out) :: point(3)

        real(dp) :: jacobian(3, 3)
        integer :: local_node

        do local_node = 1, 3
            jacobian(:, local_node) = &
                coordinates(:, tetrahedra(local_node + 1, cell)) - &
                coordinates(:, tetrahedra(1, cell))
        end do
        point = coordinates(:, tetrahedra(1, cell)) + matmul( &
            jacobian, [x(point_index), y(point_index), z(point_index)])
    end subroutine physical_sample

    subroutine check(condition, label)
        logical, intent(in) :: condition
        character(*), intent(in) :: label

        if (condition) return
        failures = failures + 1
        write(error_unit, "(a)") "FAIL: "//label
    end subroutine check

end program test_tetra_mixed_poisson_sampled_state_ad
