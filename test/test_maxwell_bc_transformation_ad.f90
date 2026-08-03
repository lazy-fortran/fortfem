program test_maxwell_bc_transformation_ad
    use check, only: check_condition, check_summary
    use fortfem_feec, only: &
        barycentric_refine_surface_mesh, &
        build_maxwell_bc_transformation, &
        differentiate_maxwell_bc_transformation_jvp, &
        differentiate_maxwell_bc_transformation_vjp
    use fortfem_core, only: generate_torus_surface_mesh
    use fortfem_kinds, only: dp
    implicit none

    real(dp), parameter :: major_radius = 2.0_dp, minor_radius = 0.6_dp
    real(dp), parameter :: step = 2.0e-6_dp
    real(dp) :: lhs, rhs
    real(dp), allocatable :: vertices(:, :), refined_vertices(:, :)
    real(dp), allocatable :: vertices_dot(:, :)
    real(dp), allocatable :: refined_vertices_dot(:, :)
    real(dp), allocatable :: refined_vertices_bar(:, :)
    real(dp), allocatable :: refined_vertices_minus(:, :)
    real(dp), allocatable :: refined_vertices_plus(:, :)
    real(dp), allocatable :: transformation(:, :), transformation_dot(:, :)
    real(dp), allocatable :: transformation_bar(:, :)
    real(dp), allocatable :: transformation_minus(:, :)
    real(dp), allocatable :: transformation_plus(:, :)
    real(dp), allocatable :: parameters(:, :)
    integer, allocatable :: triangles(:, :), refined_triangles(:, :)
    integer, allocatable :: refined_triangles_minus(:, :)
    integer, allocatable :: refined_triangles_plus(:, :)
    integer :: status, status_minus, status_plus, vertex
    logical :: all_passed

    all_passed = .true.
    call generate_torus_surface_mesh( &
        major_radius, minor_radius, 4, 5, vertices, triangles, parameters)
    call build_maxwell_bc_transformation( &
        vertices, triangles, refined_vertices, refined_triangles, &
        transformation, status)
    allocate( &
        vertices_dot(size(vertices, 1), size(vertices, 2)), &
        refined_vertices_dot(size(refined_vertices, 1), size(refined_vertices, 2)))
    do vertex = 1, size(vertices_dot, 2)
        vertices_dot(:, vertex) = [ &
            0.003_dp*sin(real(vertex, dp)), &
            -0.002_dp*cos(real(2*vertex, dp)), &
            0.004_dp*sin(real(3*vertex, dp))]
    end do
    call barycentric_refine_surface_mesh( &
        vertices + step*vertices_dot, triangles, refined_vertices_plus, &
        refined_triangles_plus)
    call barycentric_refine_surface_mesh( &
        vertices - step*vertices_dot, triangles, refined_vertices_minus, &
        refined_triangles_minus)
    refined_vertices_dot = (refined_vertices_plus - refined_vertices_minus)/(2.0_dp*step)
    call differentiate_maxwell_bc_transformation_jvp( &
        vertices, triangles, refined_vertices, refined_triangles, &
        refined_vertices_dot, transformation_dot, status)
    call build_maxwell_bc_transformation( &
        vertices + step*vertices_dot, triangles, refined_vertices_plus, &
        refined_triangles_plus, transformation_plus, status_plus)
    call build_maxwell_bc_transformation( &
        vertices - step*vertices_dot, triangles, refined_vertices_minus, &
        refined_triangles_minus, transformation_minus, status_minus)
    call record_condition(status == 0, "BC transformation JVP succeeds")
    call record_condition(status_plus == 0 .and. status_minus == 0, &
        "perturbed BC transformations succeed")
    call record_condition(maxval(abs(transformation_dot - &
        (transformation_plus - transformation_minus)/(2.0_dp*step))) < 5.0e-8_dp, &
        "BC transformation JVP matches reassembly")

    allocate( &
        transformation_bar(size(transformation, 1), size(transformation, 2)), &
        refined_vertices_bar(size(refined_vertices, 1), size(refined_vertices, 2)))
    do vertex = 1, size(transformation_bar, 2)
        transformation_bar(:, vertex) = [ &
            0.007_dp*sin(real(vertex, dp)), &
            -0.005_dp*cos(real(2*vertex, dp)), &
            0.006_dp*sin(real(3*vertex, dp))]
    end do
    call differentiate_maxwell_bc_transformation_vjp( &
        vertices, triangles, refined_vertices, refined_triangles, &
        transformation_bar, refined_vertices_bar, status)
    lhs = sum(transformation_bar*transformation_dot)
    rhs = sum(refined_vertices_bar*refined_vertices_dot)
    call record_condition(status == 0, "BC transformation VJP succeeds")
    call record_condition(abs(lhs - rhs) < 5.0e-12_dp*max(1.0_dp, abs(lhs), abs(rhs)), &
        "BC transformation products obey the adjoint identity")

    call check_summary("Differentiable Buffa--Christiansen transformation")
    if (.not. all_passed) error stop 1

contains

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_maxwell_bc_transformation_ad
