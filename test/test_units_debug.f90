program test_units_debug
    use fortfem_kinds
    use delaunay_types
    use constrained_delaunay, only: constrained_delaunay_triangulate
    implicit none

    type(mesh_t) :: mesh
    real(dp) :: points(2, 4)
    integer, allocatable :: segments(:,:)
    integer :: i, count

    write(*,*) "=== Replicating test_triangulation_units Test 2 ==="

    points(:, 1) = [0.0_dp, 0.0_dp]
    points(:, 2) = [1.0_dp, 0.0_dp]
    points(:, 3) = [1.0_dp, 1.0_dp]
    points(:, 4) = [0.0_dp, 1.0_dp]

    allocate(segments(2, 0))
    call constrained_delaunay_triangulate(points, segments, mesh)

    write(*,*) "mesh%npoints =", mesh%npoints
    write(*,*) "mesh%ntriangles =", mesh%ntriangles

    count = 0
    do i = 1, mesh%ntriangles
        write(*,*) "Triangle", i, "valid=", mesh%triangles(i)%valid,            &
            "vertices=", mesh%triangles(i)%vertices
        if (mesh%triangles(i)%valid) count = count + 1
    end do

    write(*,*) "count_valid_triangles =", count

    if (count == 0) then
        write(*,*) "FAILURE: 0 valid triangles"
        stop 1
    else
        write(*,*) "SUCCESS:", count, "valid triangles"
    end if
end program test_units_debug
