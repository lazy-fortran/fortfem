program test_tetra_lagrange_load_samples_ad
    use, intrinsic :: iso_fortran_env, only: error_unit
    use fortfem_api, only: assemble_tetra_lagrange_scalar_load_samples, &
        assemble_tetra_lagrange_scalar_load_samples_jvp, &
        assemble_tetra_lagrange_scalar_load_samples_vjp, &
        tetra_duffy_quadrature
    use fortfem_kinds, only: dp
    use fortsparse, only: fortsparse_status_t
    implicit none

    integer, parameter :: degree = 2, quadrature_degree = 6
    integer, parameter :: tetrahedra(4, 2) = reshape([ &
        1, 2, 3, 4, 2, 3, 4, 5], [4, 2])
    real(dp), parameter :: step = 2.0e-7_dp
    real(dp) :: vertices(3, 5), vertices_dot(3, 5), vertices_bar(3, 5)
    real(dp), allocatable :: eta(:), load(:), load_bar(:), load_dot(:)
    real(dp), allocatable :: load_minus(:), load_plus(:), weights(:), xi(:)
    real(dp), allocatable :: source_gradients(:, :, :)
    real(dp), allocatable :: source_parameter_dot(:, :)
    real(dp), allocatable :: source_values(:, :), source_values_bar(:, :)
    real(dp), allocatable :: source_values_minus(:, :)
    real(dp), allocatable :: source_values_plus(:, :), zeta(:)
    type(fortsparse_status_t) :: sparse_status
    real(dp) :: lhs, parameter, parameter_dot, relative_error, rhs
    integer :: dof, failures, status

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
    parameter = 0.7_dp
    parameter_dot = -0.11_dp
    call tetra_duffy_quadrature( &
        quadrature_degree, xi, eta, zeta, weights, status)
    call check(status == 0, "quadrature fixture succeeds")
    allocate(source_values(size(weights), size(tetrahedra, 2)))
    allocate(source_values_plus, mold=source_values)
    allocate(source_values_minus, mold=source_values)
    allocate(source_parameter_dot(size(weights), size(tetrahedra, 2)))
    allocate(source_gradients(3, size(weights), size(tetrahedra, 2)))
    call sample_source(vertices, parameter, source_values, source_gradients)
    call sample_source( &
        vertices + step*vertices_dot, parameter + step*parameter_dot, &
        source_values_plus, source_gradients)
    call sample_source( &
        vertices - step*vertices_dot, parameter - step*parameter_dot, &
        source_values_minus, source_gradients)
    source_parameter_dot = parameter_dot
    source_gradients = spread(spread([1.0_dp, 2.0_dp, -0.5_dp], 2, &
        size(weights)), 3, size(tetrahedra, 2))

    call assemble_tetra_lagrange_scalar_load_samples( &
        vertices, tetrahedra, degree, quadrature_degree, source_values, load, &
        sparse_status)
    call check(sparse_status%code == 0, "sampled scalar load succeeds")
    allocate(load_bar(size(load)))
    do dof = 1, size(load)
        load_bar(dof) = 0.03_dp*dof - 0.002_dp*dof**2
    end do
    call assemble_tetra_lagrange_scalar_load_samples_jvp( &
        vertices, tetrahedra, degree, quadrature_degree, source_values, &
        source_gradients, vertices_dot, source_parameter_dot, load_dot, &
        sparse_status)
    call check(sparse_status%code == 0, "sampled scalar load JVP succeeds")
    call assemble_tetra_lagrange_scalar_load_samples( &
        vertices + step*vertices_dot, tetrahedra, degree, quadrature_degree, &
        source_values_plus, load_plus, sparse_status)
    call assemble_tetra_lagrange_scalar_load_samples( &
        vertices - step*vertices_dot, tetrahedra, degree, quadrature_degree, &
        source_values_minus, load_minus, sparse_status)
    relative_error = maxval(abs( &
        load_dot - (load_plus - load_minus)/(2.0_dp*step)))/ &
        max(1.0_dp, maxval(abs(load_dot)))
    call check(relative_error < 5.0e-8_dp, &
        "sampled load JVP matches independent Eulerian reassembly")

    allocate(source_values_bar, mold=source_values)
    call assemble_tetra_lagrange_scalar_load_samples_vjp( &
        vertices, tetrahedra, degree, quadrature_degree, source_values, &
        source_gradients, load_bar, vertices_bar, source_values_bar, &
        sparse_status)
    call check(sparse_status%code == 0, "sampled scalar load VJP succeeds")
    lhs = dot_product(load_bar, load_dot)
    rhs = sum(vertices_bar*vertices_dot) + &
        sum(source_values_bar*source_parameter_dot)
    call check(abs(lhs - rhs) < 2.0e-10_dp*max(1.0_dp, abs(lhs), abs(rhs)), &
        "sampled load products satisfy the Eulerian adjoint identity")

    if (failures > 0) then
        write(error_unit, "(i0,a)") failures, " test(s) failed"
        stop 1
    end if
    write(*, "(a)") "PASS"

contains

    subroutine sample_source( &
            mesh_vertices, source_parameter, samples, gradients)
        real(dp), intent(in) :: mesh_vertices(:, :), source_parameter
        real(dp), intent(out) :: samples(:, :), gradients(:, :, :)

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
                gradients(:, quadrature_point, cell) = &
                    [1.0_dp, 2.0_dp, -0.5_dp]
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

end program test_tetra_lagrange_load_samples_ad
