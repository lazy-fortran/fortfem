program test_state_debug
    use fortfem_kinds
    use delaunay_types
    use constrained_delaunay, only: constrained_delaunay_triangulate
    implicit none

    type(mesh_t) :: mesh1, mesh2
    real(dp) :: points3(2, 3), points4(2, 4)
    integer, allocatable :: segments(:,:)
    integer :: valid1, valid2

    write(*,*) "=== Testing Global State Between Triangulations ==="

    points3(:, 1) = [0.0_dp, 0.0_dp]
    points3(:, 2) = [1.0_dp, 0.0_dp]
    points3(:, 3) = [0.5_dp, 1.0_dp]

    points4(:, 1) = [0.0_dp, 0.0_dp]
    points4(:, 2) = [1.0_dp, 0.0_dp]
    points4(:, 3) = [1.0_dp, 1.0_dp]
    points4(:, 4) = [0.0_dp, 1.0_dp]

    allocate(segments(2, 0))

    write(*,*) "First: 3-point triangulation..."
    call constrained_delaunay_triangulate(points3, segments, mesh1)
    valid1 = count_valid(mesh1)
    write(*,*) "  Valid triangles:", valid1

    write(*,*) "Second: 4-point triangulation..."
    call constrained_delaunay_triangulate(points4, segments, mesh2)
    valid2 = count_valid(mesh2)
    write(*,*) "  Valid triangles:", valid2

    if (valid1 == 1 .and. valid2 == 2) then
        write(*,*) "SUCCESS: Both triangulations correct"
    else
        write(*,*) "FAILURE: Expected 1 and 2 triangles"
        stop 1
    end if

contains

    integer function count_valid(mesh)
        type(mesh_t), intent(in) :: mesh
        integer :: i
        count_valid = 0
        do i = 1, mesh%ntriangles
            if (mesh%triangles(i)%valid) count_valid = count_valid + 1
        end do
    end function

end program test_state_debug
