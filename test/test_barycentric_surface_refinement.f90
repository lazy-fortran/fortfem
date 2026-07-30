program test_barycentric_surface_refinement
    use check, only: check_condition, check_summary
    use fortfem_api, only: barycentric_refine_surface_mesh
    use fortfem_kinds, only: dp
    implicit none

    integer, allocatable :: refined_triangles(:, :)
    integer, allocatable :: edge_midpoints(:), face_centroids(:)
    integer, allocatable :: parent_faces(:), primal_edges(:, :), sectors(:)
    real(dp), allocatable :: refined_vertices(:, :)
    real(dp) :: normal(3), vertices(3, 4)
    integer :: boundary_triangles(3, 4), local_sector, triangle
    logical :: all_passed

    all_passed = .true.
    vertices(:, 1) = [0.0_dp, 0.0_dp, 0.0_dp]
    vertices(:, 2) = [1.0_dp, 0.0_dp, 0.0_dp]
    vertices(:, 3) = [0.0_dp, 1.0_dp, 0.0_dp]
    vertices(:, 4) = [0.0_dp, 0.0_dp, 1.0_dp]
    boundary_triangles(:, 1) = [1, 3, 2]
    boundary_triangles(:, 2) = [1, 2, 4]
    boundary_triangles(:, 3) = [1, 4, 3]
    boundary_triangles(:, 4) = [2, 3, 4]
    call barycentric_refine_surface_mesh( &
        vertices, boundary_triangles, refined_vertices, refined_triangles, &
        primal_edges, edge_midpoints, face_centroids, parent_faces, sectors)
    call record_condition(size(refined_vertices, 2) == 14 .and. &
        size(refined_triangles, 2) == 24, &
        "Barycentric tetrahedral surface has V+E+F vertices and 6F panels")
    do triangle = 1, size(refined_triangles, 2)
        normal = cross_product( &
            refined_vertices(:, refined_triangles(2, triangle)) - &
            refined_vertices(:, refined_triangles(1, triangle)), &
            refined_vertices(:, refined_triangles(3, triangle)) - &
            refined_vertices(:, refined_triangles(1, triangle)))
        call record_condition(dot_product(normal, sum( &
            refined_vertices(:, refined_triangles(:, triangle)), dim=2) - &
            3.0_dp*[0.25_dp, 0.25_dp, 0.25_dp]) > 0.0_dp, &
            "Barycentric surface refinement preserves outward orientation")
    end do
    call record_condition(minval(sqrt(sum(( &
        refined_vertices(:, 5:) - spread( &
        [0.5_dp, 0.0_dp, 0.0_dp], 2, 10))**2, dim=1))) < 2.0e-14_dp, &
        "Barycentric refinement contains the analytical edge midpoint")
    call record_condition(size(primal_edges, 2) == 6 .and. &
        all(edge_midpoints == [(triangle, triangle=5, 10)]) .and. &
        all(face_centroids == [(triangle, triangle=11, 14)]), &
        "Barycentric topology exposes primal-edge and dual-vertex numbering")
    call record_condition(all(parent_faces == [( &
        (triangle, local_sector=1, 6), triangle=1, 4)]) .and. &
        all(sectors == [(local_sector, local_sector=1, 6), &
        (local_sector, local_sector=1, 6), &
        (local_sector, local_sector=1, 6), &
        (local_sector, local_sector=1, 6)]), &
        "Barycentric topology records each parent face and local sector")
    call check_summary("Barycentric surface refinement")
    if (.not. all_passed) error stop 1

contains

    pure function cross_product(first, second) result(product)
        real(dp), intent(in) :: first(3), second(3)
        real(dp) :: product(3)

        product = [ &
            first(2)*second(3) - first(3)*second(2), &
            first(3)*second(1) - first(1)*second(3), &
            first(1)*second(2) - first(2)*second(1)]
    end function cross_product

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_barycentric_surface_refinement
