program test_debug_refinement
    !> Debug test for mesh refinement issues
    use fortfem_kinds, only: dp
    use delaunay_types, only: mesh_t, point_t, add_point, add_triangle,       &
        create_mesh, destroy_mesh
    use bowyer_watson, only: insert_point, find_cavity, delaunay_triangulate
    use geometric_predicates, only: circumcenter, triangle_area,              &
        triangle_angles, disable_robust_predicates, is_robust_mode,           &
        init_robust_predicates
    use constrained_delaunay, only: constrained_delaunay_triangulate
    use mesh_refinement, only: refine_delaunay
    use point_location, only: build_adjacency
    implicit none

    type(mesh_t) :: mesh
    real(dp) :: points(2, 5)
    integer :: segments(2, 4)
    integer, allocatable :: cdt_segments(:,:)
    type(point_t) :: pa, pb, pc, center
    integer :: t, v1, v2, v3, new_vertex, ncavity
    integer :: cavity_triangles(64)
    real(dp) :: area, angles(3), tri_min

    points(:,1) = [0.0_dp, 0.0_dp]
    points(:,2) = [2.0_dp, 0.0_dp]
    points(:,3) = [2.0_dp, 2.0_dp]
    points(:,4) = [0.0_dp, 2.0_dp]
    points(:,5) = [0.5_dp, 0.15_dp]

    segments(:,1) = [1, 2]
    segments(:,2) = [2, 3]
    segments(:,3) = [3, 4]
    segments(:,4) = [4, 1]

    write(*,*) "=== Debug Refinement Test ==="
    write(*,*) ""
    write(*,*) "--- Investigating incircle issue ---"

    call constrained_delaunay_triangulate(points, segments, mesh,             &
                                          final_segments=cdt_segments)
    write(*,'(A,I5)') " After CDT: npoints =", mesh%npoints
    write(*,'(A,I5)') " After CDT: ntriangles =", mesh%ntriangles
    write(*,'(A,L5)') " Robust mode after CDT:", is_robust_mode()

    ! Count valid triangles and show super-triangle status
    ncavity = 0
    do t = 1, mesh%ntriangles
        if (mesh%triangles(t)%valid) then
            ncavity = ncavity + 1
            write(*,'(A,I3,A,3I5)') "   Valid tri ", t, " vertices:",       &
                mesh%triangles(t)%vertices(1), mesh%triangles(t)%vertices(2),&
                mesh%triangles(t)%vertices(3)
        end if
    end do
    write(*,'(A,I5)') " Valid triangles after CDT:", ncavity
    write(*,'(A,3I5)') " Super vertices:", mesh%super_vertices(1),          &
        mesh%super_vertices(2), mesh%super_vertices(3)

    ! Show CDT segments
    write(*,'(A,I5)') " CDT segments count:", size(cdt_segments, 2)
    do t = 1, size(cdt_segments, 2)
        v1 = cdt_segments(1, t)
        v2 = cdt_segments(2, t)
        write(*,'(A,I3,A,2I5,A,4F8.3)') "   Seg ", t, " v:", v1, v2,        &
            " coords:", mesh%points(v1)%x, mesh%points(v1)%y,               &
            mesh%points(v2)%x, mesh%points(v2)%y
    end do
    write(*,*) ""

    ! Instead of calling refine_delaunay, do single step manually
    write(*,*) "--- Manual single refinement step ---"
    call disable_robust_predicates()
    call build_adjacency(mesh)

    ! Find first bad triangle after CDT
    do t = 1, mesh%ntriangles
        if (.not. mesh%triangles(t)%valid) cycle

        v1 = mesh%triangles(t)%vertices(1)
        v2 = mesh%triangles(t)%vertices(2)
        v3 = mesh%triangles(t)%vertices(3)

        if (v1 < 1 .or. v1 > mesh%npoints) cycle
        if (v2 < 1 .or. v2 > mesh%npoints) cycle
        if (v3 < 1 .or. v3 > mesh%npoints) cycle
        if (.not. mesh%points(v1)%valid) cycle
        if (.not. mesh%points(v2)%valid) cycle
        if (.not. mesh%points(v3)%valid) cycle

        pa = mesh%points(v1)
        pb = mesh%points(v2)
        pc = mesh%points(v3)

        area = triangle_area(pa, pb, pc)
        if (area <= 1.0e-14_dp) cycle

        angles = triangle_angles(pa, pb, pc)
        tri_min = min(angles(1), min(angles(2), angles(3)))

        write(*,'(A,I3,A,3I5,A,F8.2)') " Triangle ", t, " vertices:", v1,   &
            v2, v3, " min_angle:", tri_min

        if (tri_min >= 20.0_dp) cycle

        write(*,'(A,I3,A)') " Found bad triangle ", t, ":"
        write(*,'(A,2F10.4)') "   pa:", pa%x, pa%y
        write(*,'(A,2F10.4)') "   pb:", pb%x, pb%y
        write(*,'(A,2F10.4)') "   pc:", pc%x, pc%y

        center = circumcenter(pa, pb, pc)
        write(*,'(A,2F10.4)') "   circumcenter:", center%x, center%y

        ! Add the point
        new_vertex = add_point(mesh, center%x, center%y)
        write(*,'(A,I5)') " Added point as vertex:", new_vertex

        ! Try to find cavity
        call find_cavity(mesh, new_vertex, cavity_triangles, ncavity)
        write(*,'(A,I5)') " find_cavity returned ncavity:", ncavity

        if (ncavity == 0) then
            write(*,*) " ERROR: No cavity found for circumcenter!"
        else
            write(*,'(A)', advance='no') "   cavity triangles:"
            write(*,*) cavity_triangles(1:ncavity)
            call insert_point(mesh, new_vertex)
        end if

        exit
    end do

    ! Count valid triangles
    ncavity = 0
    do t = 1, mesh%ntriangles
        if (mesh%triangles(t)%valid) ncavity = ncavity + 1
    end do
    write(*,'(A,I5)') " Valid triangles after one step:", ncavity

    ! Try a second iteration
    write(*,*) ""
    write(*,*) "--- Second refinement step ---"
    call build_adjacency(mesh)

    do t = 1, mesh%ntriangles
        if (.not. mesh%triangles(t)%valid) cycle

        v1 = mesh%triangles(t)%vertices(1)
        v2 = mesh%triangles(t)%vertices(2)
        v3 = mesh%triangles(t)%vertices(3)

        if (v1 < 1 .or. v1 > mesh%npoints) cycle
        if (v2 < 1 .or. v2 > mesh%npoints) cycle
        if (v3 < 1 .or. v3 > mesh%npoints) cycle
        if (.not. mesh%points(v1)%valid) cycle
        if (.not. mesh%points(v2)%valid) cycle
        if (.not. mesh%points(v3)%valid) cycle

        pa = mesh%points(v1)
        pb = mesh%points(v2)
        pc = mesh%points(v3)

        area = triangle_area(pa, pb, pc)
        if (area <= 1.0e-14_dp) cycle

        angles = triangle_angles(pa, pb, pc)
        tri_min = min(angles(1), min(angles(2), angles(3)))

        if (tri_min >= 20.0_dp) cycle

        write(*,'(A,I3,A)') " Found bad triangle ", t, ":"
        write(*,'(A,3I5)') "   vertices:", v1, v2, v3
        write(*,'(A,2F10.4)') "   pa:", pa%x, pa%y
        write(*,'(A,2F10.4)') "   pb:", pb%x, pb%y
        write(*,'(A,2F10.4)') "   pc:", pc%x, pc%y

        center = circumcenter(pa, pb, pc)
        write(*,'(A,2F10.4)') "   circumcenter:", center%x, center%y

        new_vertex = add_point(mesh, center%x, center%y)
        write(*,'(A,I5)') " Added point as vertex:", new_vertex

        call find_cavity(mesh, new_vertex, cavity_triangles, ncavity)
        write(*,'(A,I5)') " find_cavity returned ncavity:", ncavity

        if (ncavity == 0) then
            write(*,*) " ERROR: No cavity found for circumcenter!"
            write(*,*) " This is the bug - point added but not triangulated"
        else
            write(*,'(A)', advance='no') "   cavity triangles:"
            write(*,*) cavity_triangles(1:ncavity)
            call insert_point(mesh, new_vertex)
        end if

        exit
    end do

    ncavity = 0
    do t = 1, mesh%ntriangles
        if (mesh%triangles(t)%valid) ncavity = ncavity + 1
    end do
    write(*,'(A,I5)') " Valid triangles after second step:", ncavity

    call destroy_mesh(mesh)

    write(*,*) ""
    write(*,*) "--- Testing with plain Delaunay (like debug above) ---"
    call delaunay_triangulate(points, mesh)
    write(*,'(A,I5)') " Initial triangulation: npoints =", mesh%npoints
    write(*,'(A,I5)') " Initial triangulation: ntriangles =", mesh%ntriangles
    write(*,'(A,L5)') " Robust mode after triangulate:", is_robust_mode()
    write(*,*) ""

    ! Disable robust predicates as refinement does
    call disable_robust_predicates()
    write(*,'(A,L5)') " Robust mode after disable:", is_robust_mode()

    ! Rebuild adjacency
    call build_adjacency(mesh)

    ! Find first bad triangle
    do t = 1, mesh%ntriangles
        if (.not. mesh%triangles(t)%valid) cycle

        v1 = mesh%triangles(t)%vertices(1)
        v2 = mesh%triangles(t)%vertices(2)
        v3 = mesh%triangles(t)%vertices(3)

        if (v1 < 1 .or. v1 > mesh%npoints) cycle
        if (v2 < 1 .or. v2 > mesh%npoints) cycle
        if (v3 < 1 .or. v3 > mesh%npoints) cycle
        if (.not. mesh%points(v1)%valid) cycle
        if (.not. mesh%points(v2)%valid) cycle
        if (.not. mesh%points(v3)%valid) cycle

        pa = mesh%points(v1)
        pb = mesh%points(v2)
        pc = mesh%points(v3)

        area = triangle_area(pa, pb, pc)
        if (area <= 1.0e-14_dp) cycle

        angles = triangle_angles(pa, pb, pc)
        tri_min = min(angles(1), min(angles(2), angles(3)))

        if (tri_min >= 20.0_dp) cycle

        write(*,'(A,I3,A)') " Found bad triangle ", t, ":"
        write(*,'(A,3I5)') "   vertices:", v1, v2, v3
        write(*,'(A,F10.4)') "   min angle:", tri_min
        write(*,'(A,2F10.4)') "   pa:", pa%x, pa%y
        write(*,'(A,2F10.4)') "   pb:", pb%x, pb%y
        write(*,'(A,2F10.4)') "   pc:", pc%x, pc%y

        center = circumcenter(pa, pb, pc)
        write(*,'(A,2F10.4)') "   circumcenter:", center%x, center%y
        write(*,*) ""

        ! Add the point
        new_vertex = add_point(mesh, center%x, center%y)
        write(*,'(A,I5)') " Added point as vertex:", new_vertex
        write(*,'(A,I5)') " mesh%npoints now:", mesh%npoints

        ! Try to find cavity
        call find_cavity(mesh, new_vertex, cavity_triangles, ncavity)
        write(*,'(A,I5)') " find_cavity returned ncavity:", ncavity

        if (ncavity == 0) then
            write(*,*) " ERROR: No cavity found for circumcenter!"
            write(*,*) " This is the bug - point added but not triangulated"
        else
            write(*,'(A)', advance='no') "   cavity triangles:"
            write(*,*) cavity_triangles(1:ncavity)

            ! Actually insert the point
            call insert_point(mesh, new_vertex)
            write(*,'(A,I5)') " After insert: ntriangles =", mesh%ntriangles
        end if

        exit
    end do

    call destroy_mesh(mesh)
    write(*,*) ""
    write(*,*) "=== Debug Test Complete ==="

end program test_debug_refinement
