program test_cdt_4pt_debug
    use fortfem_kinds
    use delaunay_types
    use constrained_delaunay, only: constrained_delaunay_triangulate
    implicit none

    type(mesh_t) :: mesh
    real(dp) :: points(2, 4)
    integer, allocatable :: segments(:,:)
    integer :: i, valid_count

    write(*,*) "=== CDT 4-Point Debug ==="

    points(:, 1) = [0.0_dp, 0.0_dp]
    points(:, 2) = [1.0_dp, 0.0_dp]
    points(:, 3) = [1.0_dp, 1.0_dp]
    points(:, 4) = [0.0_dp, 1.0_dp]

    allocate(segments(2, 0))

    write(*,*) "Calling constrained_delaunay_triangulate..."
    call constrained_delaunay_triangulate(points, segments, mesh)

    write(*,*) "After CDT:"
    write(*,*) "  Total points:", mesh%npoints
    write(*,*) "  Total triangles:", mesh%ntriangles

    valid_count = 0
    do i = 1, mesh%ntriangles
        if (mesh%triangles(i)%valid) then
            valid_count = valid_count + 1
            write(*,*) "  Valid triangle", i, "vertices:", mesh%triangles(i)%vertices
        end if
    end do

    write(*,*) "Valid triangles:", valid_count

    if (valid_count == 2) then
        write(*,*) "SUCCESS: Got expected 2 triangles"
    else
        write(*,*) "FAILURE: Expected 2 triangles"
        stop 1
    end if
end program test_cdt_4pt_debug
