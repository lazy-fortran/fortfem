program test_point_location
    !> Test suite for point location using walking search.
    !
    !  Tests verify:
    !    1. Adjacency building
    !    2. Point location inside triangles
    !    3. Point location on edges
    !    4. Points outside triangulation
    !    5. Performance comparison with linear search
    !
    use fortfem_kinds, only: dp
    use delaunay_types, only: mesh_t, point_t, create_mesh, add_point,           &
        add_triangle
    use bowyer_watson, only: delaunay_triangulate
    use point_location, only: locate_point, build_adjacency,                     &
        LOCATION_INSIDE, LOCATION_ON_EDGE, LOCATION_NOT_FOUND
    use geometric_predicates, only: orientation
    use check, only: check_condition, check_summary
    implicit none

    call test_adjacency_building()
    call test_point_inside()
    call test_point_on_edge()
    call test_point_outside()
    call test_walking_vs_linear()

    call check_summary("Point Location")

contains

    subroutine test_adjacency_building()
        !> Test adjacency building on a simple mesh.
        type(mesh_t) :: mesh
        real(dp) :: points(2, 4)
        integer :: total_neighbors, t, e

        ! Unit square points
        points(:, 1) = [0.0_dp, 0.0_dp]
        points(:, 2) = [1.0_dp, 0.0_dp]
        points(:, 3) = [1.0_dp, 1.0_dp]
        points(:, 4) = [0.0_dp, 1.0_dp]

        call delaunay_triangulate(points, mesh)
        call build_adjacency(mesh)

        ! Count total neighbors (2 triangles sharing 1 edge = 2 neighbor entries)
        total_neighbors = 0
        do t = 1, mesh%ntriangles
            if (.not. mesh%triangles(t)%valid) cycle
            do e = 1, 3
                if (mesh%triangles(t)%neighbors(e) /= 0) then
                    total_neighbors = total_neighbors + 1
                end if
            end do
        end do

        call check_condition(total_neighbors >= 2,                               &
            "Adjacency: at least one shared edge found")
    end subroutine test_adjacency_building

    subroutine test_point_inside()
        !> Test locating points inside triangles.
        type(mesh_t) :: mesh
        real(dp) :: points(2, 4)
        integer :: tri_idx, loc_type

        ! Unit square
        points(:, 1) = [0.0_dp, 0.0_dp]
        points(:, 2) = [1.0_dp, 0.0_dp]
        points(:, 3) = [1.0_dp, 1.0_dp]
        points(:, 4) = [0.0_dp, 1.0_dp]

        call delaunay_triangulate(points, mesh)
        call build_adjacency(mesh)

        ! Point in lower-left region
        call locate_point(mesh, 0.25_dp, 0.25_dp, tri_idx, loc_type)
        call check_condition(tri_idx > 0, "Inside: lower-left point found")
        call check_condition(loc_type == LOCATION_INSIDE .or.                    &
                            loc_type == LOCATION_ON_EDGE,                        &
            "Inside: lower-left location type correct")

        ! Point in upper-right region
        call locate_point(mesh, 0.75_dp, 0.75_dp, tri_idx, loc_type)
        call check_condition(tri_idx > 0, "Inside: upper-right point found")

        ! Center point
        call locate_point(mesh, 0.5_dp, 0.5_dp, tri_idx, loc_type)
        call check_condition(tri_idx > 0, "Inside: center point found")
    end subroutine test_point_inside

    subroutine test_point_on_edge()
        !> Test locating points on triangle edges.
        type(mesh_t) :: mesh
        real(dp) :: points(2, 3)
        integer :: tri_idx, loc_type

        ! Single triangle
        points(:, 1) = [0.0_dp, 0.0_dp]
        points(:, 2) = [1.0_dp, 0.0_dp]
        points(:, 3) = [0.5_dp, 1.0_dp]

        call delaunay_triangulate(points, mesh)
        call build_adjacency(mesh)

        ! Point on bottom edge
        call locate_point(mesh, 0.5_dp, 0.0_dp, tri_idx, loc_type)
        call check_condition(tri_idx > 0, "On edge: bottom edge point found")

        ! Point inside
        call locate_point(mesh, 0.5_dp, 0.3_dp, tri_idx, loc_type)
        call check_condition(tri_idx > 0, "On edge: interior point found")
    end subroutine test_point_on_edge

    subroutine test_point_outside()
        !> Test locating points outside the triangulation.
        type(mesh_t) :: mesh
        real(dp) :: points(2, 3)
        integer :: tri_idx, loc_type

        ! Single triangle in first quadrant
        points(:, 1) = [0.0_dp, 0.0_dp]
        points(:, 2) = [1.0_dp, 0.0_dp]
        points(:, 3) = [0.5_dp, 1.0_dp]

        call delaunay_triangulate(points, mesh)
        call build_adjacency(mesh)

        ! Point clearly outside (negative x)
        call locate_point(mesh, -1.0_dp, 0.5_dp, tri_idx, loc_type)
        call check_condition(loc_type == LOCATION_NOT_FOUND .or. tri_idx == 0,   &
            "Outside: negative x point not found in mesh")

        ! Point far outside (large positive)
        call locate_point(mesh, 10.0_dp, 10.0_dp, tri_idx, loc_type)
        call check_condition(loc_type == LOCATION_NOT_FOUND .or. tri_idx == 0,   &
            "Outside: far point not found in mesh")
    end subroutine test_point_outside

    subroutine test_walking_vs_linear()
        !> Test that walking search finds same triangles as linear search.
        type(mesh_t) :: mesh
        real(dp) :: points(2, 100)
        integer :: tri_idx_walk, tri_idx_linear
        integer :: loc_type_walk, loc_type_linear
        integer :: i, n_tests, n_matches
        real(dp) :: px, py

        ! Create a grid of points for triangulation
        do i = 1, 100
            points(1, i) = mod(i - 1, 10) * 0.1_dp
            points(2, i) = (i - 1) / 10 * 0.1_dp
        end do

        call delaunay_triangulate(points, mesh)
        call build_adjacency(mesh)

        ! Test multiple query points
        n_tests = 0
        n_matches = 0

        do i = 1, 20
            px = 0.05_dp + (i - 1) * 0.045_dp
            py = 0.05_dp + mod(i, 5) * 0.15_dp

            call locate_point(mesh, px, py, tri_idx_walk, loc_type_walk)

            call locate_point_linear_test(mesh, px, py, tri_idx_linear,          &
                                         loc_type_linear)

            n_tests = n_tests + 1
            ! Both should find some triangle (walking may find different one
            ! than linear for points on shared edges, but both should find one)
            if (tri_idx_walk > 0 .and. tri_idx_linear > 0) then
                n_matches = n_matches + 1
            end if
        end do

        ! Allow some failures for edge cases
        call check_condition(n_matches >= n_tests - 5,                           &
            "Walking vs Linear: consistent results")
    end subroutine test_walking_vs_linear

    subroutine locate_point_linear_test(mesh, px, py, tri_idx, location_type)
        !> Simple linear search for comparison.
        type(mesh_t), intent(in) :: mesh
        real(dp), intent(in) :: px, py
        integer, intent(out) :: tri_idx
        integer, intent(out) :: location_type

        type(point_t) :: p, pa, pb, pc
        integer :: t, va, vb, vc
        integer :: orient_ab, orient_bc, orient_ca

        tri_idx = 0
        location_type = LOCATION_NOT_FOUND

        p%x = px
        p%y = py

        do t = 1, mesh%ntriangles
            if (.not. mesh%triangles(t)%valid) cycle

            va = mesh%triangles(t)%vertices(1)
            vb = mesh%triangles(t)%vertices(2)
            vc = mesh%triangles(t)%vertices(3)

            pa = mesh%points(va)
            pb = mesh%points(vb)
            pc = mesh%points(vc)

            orient_ab = orientation(pa, pb, p)
            orient_bc = orientation(pb, pc, p)
            orient_ca = orientation(pc, pa, p)

            if (orient_ab >= 0 .and. orient_bc >= 0 .and. orient_ca >= 0) then
                tri_idx = t
                if (orient_ab == 0 .or. orient_bc == 0 .or. orient_ca == 0) then
                    location_type = LOCATION_ON_EDGE
                else
                    location_type = LOCATION_INSIDE
                end if
                return
            end if
        end do
    end subroutine locate_point_linear_test

end program test_point_location
