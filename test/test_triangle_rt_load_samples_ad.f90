program test_triangle_rt_load_samples_ad
    use, intrinsic :: iso_fortran_env, only: error_unit
    use fortfem_api, only: &
        assemble_triangle_rt_vector_load_samples, &
        assemble_triangle_rt_vector_load_samples_jvp, &
        assemble_triangle_rt_vector_load_samples_vjp, &
        triangle_duffy_quadrature
    use fortfem_kinds, only: dp
    use fortfem_mesh_2d, only: mesh_2d_t
    use fortsparse, only: fortsparse_status_t
    implicit none

    integer, parameter :: degree = 1, quadrature_degree = 5
    real(dp), parameter :: step = 2.0e-7_dp
    type(mesh_2d_t) :: mesh, minus_mesh, plus_mesh
    real(dp) :: vertices_dot(2, 4), vertices_bar(2, 4)
    real(dp), allocatable :: load(:), load_bar(:), load_dot(:)
    real(dp), allocatable :: load_minus(:), load_plus(:)
    real(dp), allocatable :: source(:,:,:), source_bar(:,:,:)
    real(dp), allocatable :: source_dot(:,:,:), source_gradients(:,:,:,:)
    real(dp), allocatable :: source_minus(:,:,:), source_plus(:,:,:)
    real(dp), allocatable :: eta(:), weights(:), xi(:)
    type(fortsparse_status_t) :: status
    real(dp) :: lhs, point(2), point_dot(2), rhs
    integer :: cell, dof, failures, point_index, quadrature_status

    failures = 0
    call mesh%create_rectangular(2, 2, 0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp)
    call plus_mesh%create_rectangular(2, 2, 0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp)
    call minus_mesh%create_rectangular(2, 2, 0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp)
    vertices_dot = reshape([ &
        0.01_dp, -0.02_dp, -0.03_dp, 0.015_dp, &
        0.025_dp, -0.01_dp, 0.02_dp, 0.03_dp], [2, 4])
    plus_mesh%vertices = plus_mesh%vertices + step*vertices_dot
    minus_mesh%vertices = minus_mesh%vertices - step*vertices_dot
    call triangle_duffy_quadrature( &
        quadrature_degree, xi, eta, weights, quadrature_status)
    call check(quadrature_status == 0, "triangle quadrature succeeds")
    allocate(source(2, size(weights), 2), source_dot(2, size(weights), 2))
    allocate(source_minus(2, size(weights), 2))
    allocate(source_plus(2, size(weights), 2))
    allocate(source_gradients(2, 2, size(weights), 2))
    source_gradients(:, :, 1, 1) = reshape([0.4_dp, -0.2_dp, &
        0.15_dp, 0.3_dp], [2, 2])
    source_gradients(:, :, 1, 2) = reshape([-0.25_dp, 0.35_dp, &
        0.2_dp, 0.1_dp], [2, 2])
    do cell = 1, 2
        do point_index = 1, size(weights)
            source_gradients(:, :, point_index, cell) = &
                source_gradients(:, :, 1, cell)
            call physical_sample( &
                mesh%vertices, mesh%triangles, cell, point_index, point)
            call physical_sample( &
                vertices_dot, mesh%triangles, cell, point_index, point_dot)
            source(:, point_index, cell) = [0.3_dp*cell, -0.2_dp*cell] + &
                matmul(source_gradients(:, :, point_index, cell), point)
            source_dot(:, point_index, cell) = [ &
                0.01_dp*cell, -0.015_dp*cell] + &
                0.0005_dp*point_index
            source_plus(:, point_index, cell) = source(:, point_index, cell) + &
                step*(source_dot(:, point_index, cell) + matmul( &
                source_gradients(:, :, point_index, cell), point_dot))
            source_minus(:, point_index, cell) = source(:, point_index, cell) - &
                step*(source_dot(:, point_index, cell) + matmul( &
                source_gradients(:, :, point_index, cell), point_dot))
        end do
    end do

    call assemble_triangle_rt_vector_load_samples( &
        mesh, degree, quadrature_degree, source, load, status)
    call check(status%code == 0, "sampled triangle RT load succeeds")
    call assemble_triangle_rt_vector_load_samples_jvp( &
        mesh, degree, quadrature_degree, source, source_gradients, &
        vertices_dot, source_dot, load_dot, status)
    call check(status%code == 0, "sampled triangle RT load JVP succeeds")
    call assemble_triangle_rt_vector_load_samples( &
        plus_mesh, degree, quadrature_degree, source_plus, load_plus, status)
    call assemble_triangle_rt_vector_load_samples( &
        minus_mesh, degree, quadrature_degree, source_minus, load_minus, status)
    call check(maxval(abs( &
        load_dot - (load_plus - load_minus)/(2.0_dp*step))) < 3.0e-9_dp, &
        "sampled triangle RT load JVP matches full Eulerian reassembly")

    allocate(load_bar(size(load)), source_bar(2, size(weights), 2))
    load_bar = [(0.04_dp*dof - 0.13_dp, dof = 1, size(load))]
    call assemble_triangle_rt_vector_load_samples_vjp( &
        mesh, degree, quadrature_degree, source, source_gradients, load_bar, &
        vertices_bar, source_bar, status)
    call check(status%code == 0, "sampled triangle RT load VJP succeeds")
    lhs = dot_product(load_bar, load_dot)
    rhs = sum(vertices_bar*vertices_dot) + sum(source_bar*source_dot)
    call check(abs(lhs - rhs) < 3.0e-11_dp*max(1.0_dp, abs(lhs), abs(rhs)), &
        "sampled triangle RT load products obey the adjoint identity")

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

        jacobian(:, 1) = &
            coordinates(:, triangles(2, cell)) - coordinates(:, triangles(1, cell))
        jacobian(:, 2) = &
            coordinates(:, triangles(3, cell)) - coordinates(:, triangles(1, cell))
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

end program test_triangle_rt_load_samples_ad
