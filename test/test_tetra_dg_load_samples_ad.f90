program test_tetra_dg_load_samples_ad
    use, intrinsic :: iso_fortran_env, only: error_unit
    use fortfem_feec, only: &
        assemble_tetra_dg_source_load_samples, &
        assemble_tetra_dg_source_load_samples_jvp, &
        assemble_tetra_dg_source_load_samples_vjp
    use fortfem_kinds, only: dp
    use fortfem_tetra_duffy_quadrature, only: tetra_duffy_quadrature
    use fortsparse, only: fortsparse_status_t
    implicit none

    integer, parameter :: degree = 2, quadrature_degree = 4
    integer, parameter :: tetrahedra(4, 2) = reshape([ &
        1, 2, 3, 4, 2, 3, 4, 5], [4, 2])
    real(dp), parameter :: step = 2.0e-7_dp
    real(dp) :: vertices(3, 5), vertices_dot(3, 5), vertices_bar(3, 5)
    real(dp), allocatable :: load(:), load_bar(:), load_dot(:)
    real(dp), allocatable :: load_minus(:), load_plus(:)
    real(dp), allocatable :: source(:,:), source_bar(:,:)
    real(dp), allocatable :: source_dot(:,:), source_gradients(:,:,:)
    real(dp), allocatable :: source_minus(:,:), source_plus(:,:)
    real(dp), allocatable :: weights(:), x(:), y(:), z(:)
    real(dp) :: lhs, point(3), point_dot(3), rhs
    type(fortsparse_status_t) :: status
    integer :: cell, failures, node, point_index, quadrature_status

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
    call tetra_duffy_quadrature( &
        quadrature_degree, x, y, z, weights, quadrature_status)
    call check(quadrature_status == 0, "tetra quadrature succeeds")
    allocate(source(size(weights), 2), source_dot(size(weights), 2))
    allocate(source_minus(size(weights), 2), source_plus(size(weights), 2))
    allocate(source_gradients(3, size(weights), 2))
    source_gradients(:, :, 1) = spread([0.7_dp, -0.4_dp, 0.2_dp], 2, &
        size(weights))
    source_gradients(:, :, 2) = spread([-0.3_dp, 0.6_dp, 0.5_dp], 2, &
        size(weights))
    do cell = 1, 2
        do point_index = 1, size(weights)
            call physical_sample(vertices, cell, point_index, point)
            source(point_index, cell) = 0.4_dp*cell + &
                dot_product(source_gradients(:, point_index, cell), point)
            source_dot(point_index, cell) = &
                0.03_dp*cell - 0.002_dp*point_index
            call physical_sample(vertices_dot, cell, point_index, point_dot)
            source_plus(point_index, cell) = source(point_index, cell) + &
                step*(source_dot(point_index, cell) + dot_product( &
                source_gradients(:, point_index, cell), point_dot))
            source_minus(point_index, cell) = source(point_index, cell) - &
                step*(source_dot(point_index, cell) + dot_product( &
                source_gradients(:, point_index, cell), point_dot))
        end do
    end do

    call assemble_tetra_dg_source_load_samples( &
        vertices, tetrahedra, degree, quadrature_degree, source, load, status)
    call check(status%code == 0, "sampled tetrahedral DG load succeeds")
    call assemble_tetra_dg_source_load_samples_jvp( &
        vertices, tetrahedra, degree, quadrature_degree, source, &
        source_gradients, vertices_dot, source_dot, load_dot, status)
    call check(status%code == 0, "sampled tetrahedral DG load JVP succeeds")
    call assemble_tetra_dg_source_load_samples( &
        vertices + step*vertices_dot, tetrahedra, degree, quadrature_degree, &
        source_plus, load_plus, status)
    call assemble_tetra_dg_source_load_samples( &
        vertices - step*vertices_dot, tetrahedra, degree, quadrature_degree, &
        source_minus, load_minus, status)
    call check(maxval(abs( &
        load_dot - (load_plus - load_minus)/(2.0_dp*step))) < 2.0e-9_dp, &
        "sampled DG load JVP matches independent full reassembly")

    allocate(load_bar(size(load)), source_bar(size(source, 1), 2))
    load_bar = [(0.04_dp*node - 0.13_dp, node = 1, size(load))]
    call assemble_tetra_dg_source_load_samples_vjp( &
        vertices, tetrahedra, degree, quadrature_degree, source, &
        source_gradients, load_bar, vertices_bar, source_bar, status)
    call check(status%code == 0, "sampled tetrahedral DG load VJP succeeds")
    lhs = dot_product(load_bar, load_dot)
    rhs = sum(vertices_bar*vertices_dot) + sum(source_bar*source_dot)
    call check(abs(lhs - rhs) < 2.0e-11_dp*max(1.0_dp, abs(lhs), abs(rhs)), &
        "sampled DG load products obey the adjoint identity")

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

end program test_tetra_dg_load_samples_ad
