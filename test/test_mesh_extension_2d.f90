program test_mesh_extension_2d
    use check, only: check_condition, check_summary
    use fortfem_kinds, only: dp
    use fortfem_mesh_2d, only: mesh_2d_t
    use fortfem_mesh_extension_2d, only: extend_triangle_mesh
    implicit none

    real(dp), parameter :: outer_vertices(2, 4) = reshape([ &
        1.0_dp, -1.0_dp, 4.0_dp, -1.0_dp, &
        4.0_dp, 2.0_dp, 1.0_dp, 2.0_dp], [2, 4])
    type(mesh_2d_t) :: core, extended
    real(dp), allocatable :: areas(:)
    integer :: edge, status
    logical :: all_passed

    all_passed = .true.
    core%n_vertices = 4
    core%n_triangles = 2
    core%has_triangles = .true.
    allocate(core%vertices(2, 4), core%triangles(3, 2))
    core%vertices(:, 1) = [2.0_dp, 0.0_dp]
    core%vertices(:, 2) = [3.0_dp, 0.0_dp]
    core%vertices(:, 3) = [3.0_dp, 1.0_dp]
    core%vertices(:, 4) = [2.0_dp, 1.0_dp]
    core%triangles(:, 1) = [1, 2, 3]
    core%triangles(:, 2) = [1, 3, 4]

    call extend_triangle_mesh(core, outer_vertices, extended, status)
    call record_condition(status == 0, &
        "Triangle mesh extension succeeds for nested squares")
    call record_condition(extended%n_vertices == 8 .and. &
        extended%n_triangles == 10 .and. extended%n_edges == 17, &
        "Extended square mesh satisfies planar topology")
    call record_condition(extended%n_boundary_edges == 4, &
        "Only the supplied outer loop remains on the boundary")
    call record_condition(all(extended%triangles(:, 1:2) == core%triangles), &
        "Mesh extension preserves the core triangle connectivity")

    areas = extended%compute_areas()
    call record_condition(abs(sum(areas) - 9.0_dp) < 1.0e-13_dp, &
        "Extended mesh covers the exact outer-square area")
    call record_condition(all_positive_triangles(extended), &
        "Extended mesh contains only counter-clockwise triangles")

    do edge = 1, 4
        call record_condition(.not. pair_is_boundary(extended, edge, &
            mod(edge, 4) + 1), &
            "Former core boundary edge becomes an interior edge")
    end do

    call check_summary("Triangle mesh extension")
    if (.not. all_passed) error stop 1

contains

    logical function all_positive_triangles(mesh) result(positive)
        type(mesh_2d_t), intent(in) :: mesh

        real(dp) :: determinant
        integer :: triangle

        positive = .true.
        do triangle = 1, mesh%n_triangles
            determinant = &
                (mesh%vertices(1, mesh%triangles(2, triangle)) - &
                mesh%vertices(1, mesh%triangles(1, triangle))) * &
                (mesh%vertices(2, mesh%triangles(3, triangle)) - &
                mesh%vertices(2, mesh%triangles(1, triangle))) - &
                (mesh%vertices(1, mesh%triangles(3, triangle)) - &
                mesh%vertices(1, mesh%triangles(1, triangle))) * &
                (mesh%vertices(2, mesh%triangles(2, triangle)) - &
                mesh%vertices(2, mesh%triangles(1, triangle)))
            positive = positive .and. determinant > 0.0_dp
        end do
    end function all_positive_triangles

    logical function pair_is_boundary(mesh, vertex_a, vertex_b) result(found)
        type(mesh_2d_t), intent(in) :: mesh
        integer, intent(in) :: vertex_a, vertex_b

        integer :: edge

        found = .false.
        do edge = 1, mesh%n_edges
            if (all(mesh%edges(:, edge) == [ &
                min(vertex_a, vertex_b), max(vertex_a, vertex_b)])) then
                found = mesh%is_boundary_edge(edge)
                return
            end if
        end do
    end function pair_is_boundary

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_mesh_extension_2d
