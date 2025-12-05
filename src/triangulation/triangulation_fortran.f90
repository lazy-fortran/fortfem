module triangulation_fortran
    use fortfem_kinds, only: dp
    use delaunay_types
    use bowyer_watson
    use constrained_delaunay
    use mesh_refinement, only: refine_delaunay
    implicit none

    private
    public :: triangulation_result_t, triangulate_fortran
    public :: triangulate_with_hole_fortran
    public :: triangulate_with_quality_fortran, triangulate_triangle_lib
    public :: cleanup_triangulation

    type :: triangulation_result_t
        integer :: npoints
        integer :: ntriangles
        integer :: nsegments
        real(dp), allocatable :: points(:,:)
        integer, allocatable :: triangles(:,:)
        integer, allocatable :: segments(:,:)
        integer, allocatable :: neighbors(:,:)
    end type triangulation_result_t

contains

    subroutine triangulate_fortran(points, segments, result, status)
        real(dp), intent(in) :: points(:,:)
        integer, intent(in) :: segments(:,:)
        type(triangulation_result_t), intent(out) :: result
        integer, intent(out), optional :: status

        type(mesh_t) :: mesh

        if (present(status)) status = 0

        call constrained_delaunay_triangulate(points, segments, mesh)
        call mesh_to_result(mesh, segments, result)
        call destroy_mesh(mesh)
    end subroutine triangulate_fortran

    subroutine triangulate_with_hole_fortran(points, segments, hole_points,   &
                                             result, status)
        !> Triangulate with one or more holes.
        !
        !  hole_points can be either:
        !    - A 1D array of length 2: [x, y] for a single hole
        !    - A 2D array of shape (2, n_holes): each column is a hole point
        !
        real(dp), intent(in) :: points(:,:)
        integer, intent(in) :: segments(:,:)
        real(dp), intent(in) :: hole_points(:,:)
        type(triangulation_result_t), intent(out) :: result
        integer, intent(out), optional :: status

        type(mesh_t) :: mesh

        if (size(hole_points, 1) /= 2) then
            if (present(status)) status = 1
            return
        end if

        if (present(status)) status = 0

        call constrained_delaunay_triangulate(points, segments, mesh,         &
                                             hole_points)
        call mesh_to_result(mesh, segments, result)
        call destroy_mesh(mesh)
    end subroutine triangulate_with_hole_fortran

    subroutine triangulate_with_quality_fortran(points, segments, min_angle,  &
                                                result, status)
        real(dp), intent(in) :: points(:,:)
        integer, intent(in) :: segments(:,:)
        real(dp), intent(in) :: min_angle
        type(triangulation_result_t), intent(out) :: result
        integer, intent(out), optional :: status

        type(mesh_t) :: mesh
        integer, allocatable :: cdt_segments(:,:)
        integer, allocatable :: empty_segments(:,:)
        real(dp) :: min_angle_measured
        real(dp) :: max_area

        if (present(status)) status = 0

        ! Build constrained Delaunay triangulation and capture the final
        ! constraint segment set (including any Steiner splits).
        call constrained_delaunay_triangulate(points, segments, mesh,         &
                                             final_segments=cdt_segments)

        ! Use angle-based refinement only for now; disable area constraint
        ! by setting max_area to a very large value.
        max_area = huge(1.0_dp)

        if (allocated(cdt_segments)) then
            call refine_delaunay(mesh, cdt_segments, min_angle, max_area)
        else
            ! Unconstrained case: refine using an empty segment set so that
            ! only angle/area criteria drive refinement.
            allocate(empty_segments(2, 0))
            call refine_delaunay(mesh, empty_segments, min_angle, max_area)
            deallocate(empty_segments)
        end if

        call mesh_to_result(mesh, segments, result)
        call destroy_mesh(mesh)

        if (allocated(cdt_segments)) deallocate(cdt_segments)

        min_angle_measured = compute_min_triangle_angle(result)

        if (present(status)) then
            if (min_angle_measured >= min_angle) then
                status = 0
            else
                status = 1
            end if
        end if
    end subroutine triangulate_with_quality_fortran

    subroutine triangulate_triangle_lib(points, segments, result, status)
        real(dp), intent(in) :: points(:,:)
        integer, intent(in) :: segments(:,:)
        type(triangulation_result_t), intent(out) :: result
        integer, intent(out), optional :: status

        call triangulate_fortran(points, segments, result, status)
    end subroutine triangulate_triangle_lib

    subroutine allocate_result(result, npoints, ntriangles, nsegments)
        type(triangulation_result_t), intent(out) :: result
        integer, intent(in) :: npoints, ntriangles, nsegments

        result%npoints = npoints
        result%ntriangles = ntriangles
        result%nsegments = nsegments

        allocate(result%points(2, npoints))
        allocate(result%triangles(3, ntriangles))
        allocate(result%segments(2, nsegments))
        allocate(result%neighbors(3, ntriangles))
    end subroutine allocate_result

    subroutine mesh_to_result(mesh, input_segments, result)
        type(mesh_t), intent(in) :: mesh
        integer, intent(in) :: input_segments(:,:)
        type(triangulation_result_t), intent(out) :: result

        integer :: i, valid_points, valid_triangles

        valid_points = 0
        do i = 1, mesh%npoints
            if (mesh%points(i)%valid) valid_points = valid_points + 1
        end do

        valid_triangles = 0
        do i = 1, mesh%ntriangles
            if (mesh%triangles(i)%valid) then
                if (all(mesh%triangles(i)%vertices >= 1 .and.                 &
                        mesh%triangles(i)%vertices <= mesh%npoints)) then
                    if (all(mesh%points(mesh%triangles(i)%vertices)%valid)) then
                        valid_triangles = valid_triangles + 1
                    end if
                end if
            end if
        end do

        call allocate_result(result, valid_points, valid_triangles,           &
                            size(input_segments, 2))

        valid_points = 0
        do i = 1, mesh%npoints
            if (mesh%points(i)%valid) then
                valid_points = valid_points + 1
                result%points(1, valid_points) = mesh%points(i)%x
                result%points(2, valid_points) = mesh%points(i)%y
            end if
        end do

        valid_triangles = 0
        do i = 1, mesh%ntriangles
            if (mesh%triangles(i)%valid) then
                if (all(mesh%triangles(i)%vertices >= 1 .and.                 &
                        mesh%triangles(i)%vertices <= mesh%npoints)) then
                    if (all(mesh%points(mesh%triangles(i)%vertices)%valid)) then
                        valid_triangles = valid_triangles + 1
                        result%triangles(1, valid_triangles) =                &
                            remap_vertex_index(mesh, mesh%triangles(i)%vertices(1))
                        result%triangles(2, valid_triangles) =                &
                            remap_vertex_index(mesh, mesh%triangles(i)%vertices(2))
                        result%triangles(3, valid_triangles) =                &
                            remap_vertex_index(mesh, mesh%triangles(i)%vertices(3))
                    end if
                end if
            end if
        end do

        result%segments = input_segments
        result%npoints = valid_points
        result%ntriangles = valid_triangles
    end subroutine mesh_to_result

    integer function remap_vertex_index(mesh, original_idx)
        type(mesh_t), intent(in) :: mesh
        integer, intent(in) :: original_idx

        integer :: i, valid_count

        if (original_idx < 1 .or. original_idx > mesh%npoints) then
            remap_vertex_index = -1
            return
        end if

        if (.not. mesh%points(original_idx)%valid) then
            remap_vertex_index = -1
            return
        end if

        valid_count = 0
        do i = 1, original_idx
            if (mesh%points(i)%valid) then
                valid_count = valid_count + 1
            end if
        end do

        remap_vertex_index = valid_count
    end function remap_vertex_index

    real(dp) function compute_min_triangle_angle(result) result(min_angle_deg)
        type(triangulation_result_t), intent(in) :: result

        integer :: i
        real(dp) :: ax, ay, bx, by, cx, cy
        real(dp) :: angle1, angle2, angle3, tri_min
        real(dp), parameter :: pi = acos(-1.0_dp)

        min_angle_deg = 180.0_dp

        do i = 1, result%ntriangles
            ax = result%points(1, result%triangles(1, i))
            ay = result%points(2, result%triangles(1, i))
            bx = result%points(1, result%triangles(2, i))
            by = result%points(2, result%triangles(2, i))
            cx = result%points(1, result%triangles(3, i))
            cy = result%points(2, result%triangles(3, i))

            angle1 = interior_angle(ax, ay, bx, by, cx, cy, pi)
            angle2 = interior_angle(bx, by, cx, cy, ax, ay, pi)
            angle3 = interior_angle(cx, cy, ax, ay, bx, by, pi)

            tri_min = min(angle1, min(angle2, angle3))
            if (tri_min < min_angle_deg) min_angle_deg = tri_min
        end do
    end function compute_min_triangle_angle

    real(dp) function interior_angle(px, py, qx, qy, rx, ry, pi)              &
        result(angle_deg)
        real(dp), intent(in) :: px, py, qx, qy, rx, ry, pi
        real(dp) :: v1x, v1y, v2x, v2y, dotpr, norm1, norm2, cos_theta

        v1x = px - qx
        v1y = py - qy
        v2x = rx - qx
        v2y = ry - qy

        norm1 = sqrt(v1x*v1x + v1y*v1y)
        norm2 = sqrt(v2x*v2x + v2y*v2y)

        if (norm1 <= 0.0_dp .or. norm2 <= 0.0_dp) then
            angle_deg = 0.0_dp
            return
        end if

        dotpr = v1x*v2x + v1y*v2y
        cos_theta = max(-1.0_dp, min(1.0_dp, dotpr / (norm1*norm2)))
        angle_deg = acos(cos_theta) * 180.0_dp / pi
    end function interior_angle

    subroutine cleanup_triangulation(result)
        type(triangulation_result_t), intent(inout) :: result

        if (allocated(result%points)) deallocate(result%points)
        if (allocated(result%triangles)) deallocate(result%triangles)
        if (allocated(result%segments)) deallocate(result%segments)
        if (allocated(result%neighbors)) deallocate(result%neighbors)

        result%npoints = 0
        result%ntriangles = 0
        result%nsegments = 0
    end subroutine cleanup_triangulation

end module triangulation_fortran
