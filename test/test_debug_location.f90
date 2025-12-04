program test_debug_location
    use fortfem_kinds, only: dp
    use delaunay_types, only: mesh_t
    use bowyer_watson, only: delaunay_triangulate
    use point_location, only: build_adjacency, locate_point,                     &
        LOCATION_INSIDE, LOCATION_ON_EDGE, LOCATION_NOT_FOUND
    implicit none

    type(mesh_t) :: mesh
    real(dp) :: points(2, 4)
    integer :: t, e, n_valid, tri_idx, loc_type, i, total_neighbors

    points(:, 1) = [0.0_dp, 0.0_dp]
    points(:, 2) = [1.0_dp, 0.0_dp]
    points(:, 3) = [1.0_dp, 1.0_dp]
    points(:, 4) = [0.0_dp, 1.0_dp]

    call delaunay_triangulate(points, mesh)

    n_valid = 0
    do t = 1, mesh%ntriangles
        if (mesh%triangles(t)%valid) n_valid = n_valid + 1
    end do
    print *, "Valid triangles after triangulate:", n_valid

    call build_adjacency(mesh)

    print *, ""
    print *, "After build_adjacency - valid triangles with neighbors:"
    total_neighbors = 0
    do t = 1, mesh%ntriangles
        if (mesh%triangles(t)%valid) then
            print *, "Triangle", t, "vertices:", mesh%triangles(t)%vertices,     &
                     "neighbors:", mesh%triangles(t)%neighbors
            do e = 1, 3
                if (mesh%triangles(t)%neighbors(e) /= 0) then
                    total_neighbors = total_neighbors + 1
                end if
            end do
        end if
    end do
    print *, "Total neighbor entries:", total_neighbors

    print *, ""
    print *, "Testing locate_point(0.25, 0.25):"
    call locate_point(mesh, 0.25_dp, 0.25_dp, tri_idx, loc_type)
    print *, "  tri_idx:", tri_idx, "loc_type:", loc_type

    print *, ""
    print *, "Testing locate_point(0.75, 0.75):"
    call locate_point(mesh, 0.75_dp, 0.75_dp, tri_idx, loc_type)
    print *, "  tri_idx:", tri_idx, "loc_type:", loc_type

    print *, ""
    print *, "Testing locate_point(0.5, 0.5):"
    call locate_point(mesh, 0.5_dp, 0.5_dp, tri_idx, loc_type)
    print *, "  tri_idx:", tri_idx, "loc_type:", loc_type
end program test_debug_location
