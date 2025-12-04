module polygon_triangulation
    ! Utilities for robust triangulation of simple polygon boundaries.
    use fortfem_kinds, only: dp
    use delaunay_types, only: point_t
    use geometric_predicates, only: orientation, ORIENTATION_CCW
    implicit none

    private
    public :: is_simple_boundary_chain
    public :: triangulate_simple_polygon

contains

    logical function is_simple_boundary_chain(boundary_points, segments) &
        result(is_simple)
        !> Detect boundaries created by boundary_t: segments form the chain
        !  (1,2), (2,3), ..., (n-1,n), (n,1) with no additional edges.
        real(dp), intent(in) :: boundary_points(:,:)
        integer, intent(in) :: segments(:,:)

        integer :: n_points, n_segments, i

        n_points = size(boundary_points, 2)
        n_segments = size(segments, 2)

        if (n_points <= 2) then
            is_simple = .false.
            return
        end if

        if (n_segments /= n_points) then
            is_simple = .false.
            return
        end if

        do i = 1, n_points - 1
            if (segments(1, i) /= i) then
                is_simple = .false.
                return
            end if
            if (segments(2, i) /= i + 1) then
                is_simple = .false.
                return
            end if
        end do

        if (segments(1, n_points) /= n_points) then
            is_simple = .false.
            return
        end if
        if (segments(2, n_points) /= 1) then
            is_simple = .false.
            return
        end if

        is_simple = .true.
    end function is_simple_boundary_chain

    subroutine triangulate_simple_polygon(boundary_points, mesh_points, &
                                          mesh_triangles, n_points,    &
                                          n_triangles)
        !> Ear-clipping triangulation for a single simple polygon whose
        !  vertices are given in boundary order. This produces a valid,
        !  non-overlapping triangulation of the polygon interior.
        real(dp), intent(in) :: boundary_points(:,:)
        real(dp), allocatable, intent(out) :: mesh_points(:,:)
        integer, allocatable, intent(out) :: mesh_triangles(:,:)
        integer, intent(out) :: n_points, n_triangles

        integer :: n_vertices
        integer, allocatable :: polygon(:)
        integer :: i, prev_idx, next_idx
        integer :: t
        logical :: ear_found
        real(dp) :: signed_area
        logical :: polygon_ccw
        type(point_t) :: pa, pb, pc

        n_vertices = size(boundary_points, 2)
        n_points = n_vertices

        if (n_vertices < 3) then
            n_triangles = 0
            allocate(mesh_points(2, n_points))
            allocate(mesh_triangles(3, 0))
            mesh_points = boundary_points(:, 1:n_points)
            return
        end if

        allocate(polygon(n_vertices))
        do i = 1, n_vertices
            polygon(i) = i
        end do

        signed_area = compute_polygon_signed_area(boundary_points, polygon, &
                                                  n_vertices)
        polygon_ccw = signed_area > 0.0_dp

        n_triangles = n_vertices - 2
        allocate(mesh_points(2, n_points))
        allocate(mesh_triangles(3, n_triangles))
        mesh_points = boundary_points(:, 1:n_points)

        t = 0
        do while (n_vertices > 3)
            ear_found = .false.

            do i = 1, n_vertices
                prev_idx = modulo_index(i - 1, n_vertices)
                next_idx = modulo_index(i + 1, n_vertices)

                pa%x = boundary_points(1, polygon(prev_idx))
                pa%y = boundary_points(2, polygon(prev_idx))
                pb%x = boundary_points(1, polygon(i))
                pb%y = boundary_points(2, polygon(i))
                pc%x = boundary_points(1, polygon(next_idx))
                pc%y = boundary_points(2, polygon(next_idx))

                if (.not. is_convex_vertex(pa, pb, pc, polygon_ccw)) cycle
                if (has_point_in_triangle(boundary_points, polygon,      &
                                          n_vertices, prev_idx, i,      &
                                          next_idx, polygon_ccw)) cycle

                t = t + 1
                mesh_triangles(1, t) = polygon(prev_idx)
                mesh_triangles(2, t) = polygon(i)
                mesh_triangles(3, t) = polygon(next_idx)

                call remove_vertex_from_polygon(polygon, n_vertices, i)
                ear_found = .true.
                exit
            end do

            if (.not. ear_found) then
                exit
            end if
        end do

        if (n_vertices == 3) then
            t = t + 1
            mesh_triangles(1, t) = polygon(1)
            mesh_triangles(2, t) = polygon(2)
            mesh_triangles(3, t) = polygon(3)
        end if

        if (t < n_triangles) then
            n_triangles = t
            if (n_triangles > 0) then
                mesh_triangles = mesh_triangles(:, 1:n_triangles)
            else
                deallocate(mesh_triangles)
                allocate(mesh_triangles(3, 0))
            end if
        end if

        deallocate(polygon)
    end subroutine triangulate_simple_polygon

    real(dp) function compute_polygon_signed_area(points, polygon,       &
                                                  n_vertices) result(area)
        !> Compute signed area (2 * area) of polygon defined by vertex
        !  indices.
        real(dp), intent(in) :: points(:,:)
        integer, intent(in) :: polygon(:)
        integer, intent(in) :: n_vertices

        integer :: i, j
        real(dp) :: x1, y1, x2, y2

        area = 0.0_dp
        do i = 1, n_vertices
            j = modulo_index(i + 1, n_vertices)
            x1 = points(1, polygon(i))
            y1 = points(2, polygon(i))
            x2 = points(1, polygon(j))
            y2 = points(2, polygon(j))
            area = area + (x1 * y2 - x2 * y1)
        end do
    end function compute_polygon_signed_area

    logical function is_convex_vertex(pa, pb, pc, polygon_ccw)           &
        result(is_convex)
        !> Check if vertex B is convex with respect to polygon orientation.
        type(point_t), intent(in) :: pa, pb, pc
        logical, intent(in) :: polygon_ccw

        integer :: orient_val

        orient_val = orientation(pa, pb, pc)

        if (polygon_ccw) then
            is_convex = (orient_val == ORIENTATION_CCW)
        else
            is_convex = (orient_val /= ORIENTATION_CCW)
        end if
    end function is_convex_vertex

    logical function has_point_in_triangle(points, polygon, n_vertices,  &
                                           i_prev, i_curr, i_next,       &
                                           polygon_ccw) result(has_point)
        !> Check if any polygon vertex other than the ear vertices lies
        !  inside the candidate ear triangle.
        real(dp), intent(in) :: points(:,:)
        integer, intent(in) :: polygon(:)
        integer, intent(in) :: n_vertices
        integer, intent(in) :: i_prev, i_curr, i_next
        logical, intent(in) :: polygon_ccw

        integer :: k
        integer :: v_prev, v_curr, v_next, v_test
        type(point_t) :: pa, pb, pc, pd
        integer :: o1, o2, o3

        v_prev = polygon(i_prev)
        v_curr = polygon(i_curr)
        v_next = polygon(i_next)

        pa%x = points(1, v_prev)
        pa%y = points(2, v_prev)
        pb%x = points(1, v_curr)
        pb%y = points(2, v_curr)
        pc%x = points(1, v_next)
        pc%y = points(2, v_next)

        has_point = .false.
        do k = 1, n_vertices
            if (k == i_prev .or. k == i_curr .or. k == i_next) cycle

            v_test = polygon(k)
            pd%x = points(1, v_test)
            pd%y = points(2, v_test)

            o1 = orientation(pa, pb, pd)
            o2 = orientation(pb, pc, pd)
            o3 = orientation(pc, pa, pd)

            if (polygon_ccw) then
                if (o1 == ORIENTATION_CCW .and. o2 == ORIENTATION_CCW .and. &
                    o3 == ORIENTATION_CCW) then
                    has_point = .true.
                    return
                end if
            else
                if (o1 /= ORIENTATION_CCW .and. o2 /= ORIENTATION_CCW .and. &
                    o3 /= ORIENTATION_CCW) then
                    has_point = .true.
                    return
                end if
            end if
        end do
    end function has_point_in_triangle

    subroutine remove_vertex_from_polygon(polygon, n_vertices, idx)
        !> Remove vertex at position idx from polygon list.
        integer, intent(inout) :: polygon(:)
        integer, intent(inout) :: n_vertices
        integer, intent(in) :: idx

        integer :: i

        if (idx < 1 .or. idx > n_vertices) return

        do i = idx, n_vertices - 1
            polygon(i) = polygon(i + 1)
        end do
        n_vertices = n_vertices - 1
    end subroutine remove_vertex_from_polygon

    integer function modulo_index(i, n)
        !> 1-based cyclic index modulo n, robust for any integer i.
        integer, intent(in) :: i, n
        integer :: k

        if (n <= 0) then
            modulo_index = 1
            return
        end if

        k = i
        do while (k < 1)
            k = k + n
        end do
        do while (k > n)
            k = k - n
        end do
        modulo_index = k
    end function modulo_index

end module polygon_triangulation
