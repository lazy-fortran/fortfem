module fortfem_api_mesh
    use fortfem_kinds, only: dp
    use fortfem_mesh_2d, only: mesh_2d_t
    use fortfem_boundary, only: boundary_t
    use triangle_io, only: read_triangle_mesh
    use fortfem_api_types, only: mesh_t, function_t
    use fortfem_api_forms, only: init_measures
    implicit none

    private

    public :: mesh_t
    public :: unit_square_mesh
    public :: rectangle_mesh
    public :: unit_disk_mesh
    public :: circle_boundary
    public :: rectangle_boundary
    public :: line_segment
    public :: arc_segment
    public :: l_shape_boundary
    public :: mesh_from_boundary
    public :: mesh_from_arrays
    public :: mesh_from_triangle_files
    public :: mesh_from_domain
    public :: structured_quad_mesh
    public :: refine_uniform
    public :: refine_adaptive
    public :: refine_adaptive_markers
    public :: refine_adaptive_solution
    public :: compute_gradient_indicators
    public :: find_triangle_edges

    interface refine_adaptive
        module procedure refine_adaptive_markers
        module procedure refine_adaptive_solution
    end interface refine_adaptive

contains

    subroutine compute_gradient_indicators(mesh, solution, indicators)
        type(mesh_t), intent(in) :: mesh
        type(function_t), intent(in) :: solution
        real(dp), intent(out) :: indicators(:)

        integer :: e

        if (size(indicators) /= mesh%data%n_triangles) then
            error stop "compute_gradient_indicators: size mismatch"
        end if

        do e = 1, mesh%data%n_triangles
            call compute_element_gradient_indicator(mesh, solution, e,      &
                indicators(e))
        end do
    end subroutine compute_gradient_indicators

    pure subroutine compute_element_gradient_indicator(mesh, solution,      &
        element_index, indicator)
        type(mesh_t), intent(in) :: mesh
        type(function_t), intent(in) :: solution
        integer, intent(in) :: element_index
        real(dp), intent(out) :: indicator

        integer :: v1, v2, v3
        real(dp) :: x1, y1, x2, y2, x3, y3
        real(dp) :: a(2,2), det_a
        real(dp) :: b(3), c(3)
        real(dp) :: u1, u2, u3
        real(dp) :: gradx, grady

        v1 = mesh%data%triangles(1, element_index)
        v2 = mesh%data%triangles(2, element_index)
        v3 = mesh%data%triangles(3, element_index)

        x1 = mesh%data%vertices(1, v1)
        y1 = mesh%data%vertices(2, v1)
        x2 = mesh%data%vertices(1, v2)
        y2 = mesh%data%vertices(2, v2)
        x3 = mesh%data%vertices(1, v3)
        y3 = mesh%data%vertices(2, v3)

        a(1,1) = x2 - x1
        a(1,2) = x3 - x1
        a(2,1) = y2 - y1
        a(2,2) = y3 - y1

        det_a = a(1,1) * a(2,2) - a(1,2) * a(2,1)

        if (abs(det_a) < 1.0e-14_dp) then
            indicator = 0.0_dp
            return
        end if

        b(1) = (-a(2,2) + a(2,1)) / det_a
        c(1) = ( a(1,2) - a(1,1)) / det_a
        b(2) = a(2,2) / det_a
        c(2) = -a(1,2) / det_a
        b(3) = -a(2,1) / det_a
        c(3) = a(1,1) / det_a

        u1 = solution%values(v1)
        u2 = solution%values(v2)
        u3 = solution%values(v3)

        gradx = b(1) * u1 + b(2) * u2 + b(3) * u3
        grady = c(1) * u1 + c(2) * u2 + c(3) * u3

        indicator = sqrt(gradx * gradx + grady * grady)
    end subroutine compute_element_gradient_indicator

    function unit_square_mesh(n) result(mesh)
        integer, intent(in) :: n
        type(mesh_t) :: mesh

        call init_measures()

        call mesh%data%create_rectangular(nx=n, ny=n,                      &
                                          x_min=0.0_dp, x_max=1.0_dp,      &
                                          y_min=0.0_dp, y_max=1.0_dp)
        call mesh%data%build_connectivity()
        call mesh%data%find_boundary()
    end function unit_square_mesh

    function rectangle_mesh(nx, ny, domain) result(mesh)
        integer, intent(in) :: nx, ny
        real(dp), intent(in) :: domain(4)
        type(mesh_t) :: mesh

        call init_measures()
        call mesh%data%create_rectangular(nx=nx, ny=ny,                    &
                                          x_min=domain(1), x_max=domain(2),&
                                          y_min=domain(3), y_max=domain(4))
        call mesh%data%build_connectivity()
        call mesh%data%find_boundary()
    end function rectangle_mesh

    function unit_disk_mesh(resolution) result(mesh)
        real(dp), intent(in), optional :: resolution
        type(mesh_t) :: mesh
        real(dp) :: h

        h = 0.1_dp
        if (present(resolution)) h = resolution

        call init_measures()
        call mesh%data%create_unit_disk(h)
        call mesh%data%build_connectivity()
        call mesh%data%find_boundary()
    end function unit_disk_mesh

    function circle_boundary(center, radius, n) result(boundary)
        real(dp), intent(in) :: center(2), radius
        integer, intent(in) :: n
        type(boundary_t) :: boundary
        integer :: i
        real(dp) :: theta

        boundary%n_points = n
        allocate(boundary%points(2, n))
        allocate(boundary%labels(n-1))

        do i = 1, n
            theta = 2.0_dp * acos(-1.0_dp) * real(i-1, dp) / real(n, dp)
            boundary%points(1, i) = center(1) + radius * cos(theta)
            boundary%points(2, i) = center(2) + radius * sin(theta)
        end do

        boundary%labels = 1
        boundary%is_closed = .true.
    end function circle_boundary

    function rectangle_boundary(domain, n) result(boundary)
        real(dp), intent(in) :: domain(4)
        integer, intent(in) :: n
        type(boundary_t) :: boundary
        integer :: i, idx
        real(dp) :: t

        boundary%n_points = 4 * n
        allocate(boundary%points(2, 4*n))
        allocate(boundary%labels(4*n-1))

        idx = 0

        do i = 1, n
            idx = idx + 1
            t = real(i-1, dp) / real(n-1, dp)
            boundary%points(1, idx) = domain(1) + t * (domain(2) - domain(1))
            boundary%points(2, idx) = domain(3)
        end do

        do i = 1, n
            idx = idx + 1
            t = real(i-1, dp) / real(n-1, dp)
            boundary%points(1, idx) = domain(2)
            boundary%points(2, idx) = domain(3) + t * (domain(4) - domain(3))
        end do

        do i = 1, n
            idx = idx + 1
            t = real(i-1, dp) / real(n-1, dp)
            boundary%points(1, idx) = domain(2) - t * (domain(2) - domain(1))
            boundary%points(2, idx) = domain(4)
        end do

        do i = 1, n
            idx = idx + 1
            t = real(i-1, dp) / real(n-1, dp)
            boundary%points(1, idx) = domain(1)
            boundary%points(2, idx) = domain(4) - t * (domain(4) - domain(3))
        end do

        call set_rectangle_boundary_labels(boundary, n)

        boundary%is_closed = .true.
    end function rectangle_boundary

    pure subroutine set_rectangle_boundary_labels(boundary, n)
        type(boundary_t), intent(inout) :: boundary
        integer, intent(in) :: n
        integer :: i

        do i = 1, n - 1
            boundary%labels(i) = 1
        end do
        do i = n, 2*n - 2
            boundary%labels(i) = 2
        end do
        do i = 2*n - 1, 3*n - 3
            boundary%labels(i) = 3
        end do
        do i = 3*n - 2, 4*n - 1
            boundary%labels(i) = 4
        end do
    end subroutine set_rectangle_boundary_labels

    function line_segment(p1, p2, n) result(boundary)
        real(dp), intent(in) :: p1(2), p2(2)
        integer, intent(in) :: n
        type(boundary_t) :: boundary
        integer :: i
        real(dp) :: t

        boundary%n_points = n
        allocate(boundary%points(2, n))
        allocate(boundary%labels(n-1))

        do i = 1, n
            t = real(i-1, dp) / real(n-1, dp)
            boundary%points(:, i) = p1 + t * (p2 - p1)
        end do

        boundary%labels = 1
        boundary%is_closed = .false.
    end function line_segment

    function arc_segment(p1, p2, center, n) result(boundary)
        real(dp), intent(in) :: p1(2), p2(2), center(2)
        integer, intent(in) :: n
        type(boundary_t) :: boundary
        integer :: i
        real(dp) :: radius1, radius2, radius
        real(dp) :: theta1, theta2, dtheta, angle
        real(dp), parameter :: pi = acos(-1.0_dp)

        boundary%n_points = n
        allocate(boundary%points(2, n))
        allocate(boundary%labels(n-1))

        radius1 = sqrt((p1(1) - center(1))**2 + (p1(2) - center(2))**2)
        radius2 = sqrt((p2(1) - center(1))**2 + (p2(2) - center(2))**2)
        radius = 0.5_dp * (radius1 + radius2)

        theta1 = atan2(p1(2) - center(2), p1(1) - center(1))
        theta2 = atan2(p2(2) - center(2), p2(1) - center(1))
        dtheta = theta2 - theta1
        if (dtheta <= 0.0_dp) dtheta = dtheta + 2.0_dp * pi

        do i = 1, n
            angle = theta1 + dtheta * real(i-1, dp) /                        &
                real(max(n-1, 1), dp)
            boundary%points(1, i) = center(1) + radius * cos(angle)
            boundary%points(2, i) = center(2) + radius * sin(angle)
        end do

        boundary%labels = 1
        boundary%is_closed = .false.
    end function arc_segment

    function l_shape_boundary(size, n) result(boundary)
        real(dp), intent(in) :: size
        integer, intent(in) :: n
        type(boundary_t) :: boundary
        integer :: idx, n_per_segment
        real(dp) :: s

        s = size

        n_per_segment = max(n - 1, 1)
        boundary%n_points = 6 * n_per_segment
        allocate(boundary%points(2, boundary%n_points))
        allocate(boundary%labels(boundary%n_points))

        idx = 0

        call add_l_shape_segment(boundary%points, idx, n_per_segment,       &
                                 0.0_dp, 0.0_dp, s, 0.0_dp)
        call add_l_shape_segment(boundary%points, idx, n_per_segment,       &
                                 s, 0.0_dp, 0.0_dp, s)
        call add_l_shape_segment(boundary%points, idx, n_per_segment,       &
                                 s, s, s, 0.0_dp)
        call add_l_shape_segment(boundary%points, idx, n_per_segment,       &
                                 2.0_dp * s, s, 0.0_dp, s)
        call add_l_shape_segment(boundary%points, idx, n_per_segment,       &
                                 2.0_dp * s, 2.0_dp * s, -2.0_dp * s, 0.0_dp)
        call add_l_shape_segment(boundary%points, idx, n_per_segment,       &
                                 0.0_dp, 2.0_dp * s, 0.0_dp, -2.0_dp * s)

        boundary%labels = 1
        boundary%is_closed = .true.
    end function l_shape_boundary

    pure subroutine add_l_shape_segment(points, idx, n_per_segment,        &
        start_x, start_y, dx, dy)
        real(dp), intent(inout) :: points(:,:)
        integer, intent(inout) :: idx
        integer, intent(in) :: n_per_segment
        real(dp), intent(in) :: start_x, start_y, dx, dy
        integer :: i
        real(dp) :: t

        do i = 0, n_per_segment - 1
            idx = idx + 1
            t = real(i, dp) / real(n_per_segment, dp)
            points(1, idx) = start_x + dx * t
            points(2, idx) = start_y + dy * t
        end do
    end subroutine add_l_shape_segment

    function mesh_from_boundary(boundary, resolution) result(mesh)
        type(boundary_t), intent(in) :: boundary
        real(dp), intent(in), optional :: resolution
        type(mesh_t) :: mesh
        real(dp) :: h

        h = 0.1_dp
        if (present(resolution)) h = resolution

        call init_measures()
        call mesh%data%create_from_boundary(boundary, h)
        call mesh%data%build_connectivity()
        call mesh%data%find_boundary()
    end function mesh_from_boundary

    function mesh_from_arrays(vertices, triangles) result(mesh)
        real(dp), intent(in) :: vertices(:,:)
        integer, intent(in) :: triangles(:,:)
        type(mesh_t) :: mesh

        call init_measures()

        call set_triangle_mesh_metadata(mesh, size(vertices, 2),            &
                                        size(triangles, 2))

        allocate(mesh%data%vertices(2, mesh%data%n_vertices))
        allocate(mesh%data%triangles(3, mesh%data%n_triangles))

        mesh%data%vertices = vertices
        mesh%data%triangles = triangles

        call mesh%data%build_connectivity()
        call mesh%data%find_boundary()
    end function mesh_from_arrays

    function mesh_from_triangle_files(basename) result(mesh)
        character(len=*), intent(in) :: basename
        type(mesh_t) :: mesh

        real(dp), allocatable :: vertices(:,:)
        integer, allocatable :: triangles(:,:)
        integer :: n_vertices, n_triangles, stat

        call init_measures()

        call read_triangle_mesh(basename, vertices, triangles,              &
                                n_vertices, n_triangles, stat)

        if (stat /= 0) then
            mesh%data%n_vertices = 0
            mesh%data%n_triangles = 0
            return
        end if

        mesh%data%n_vertices = n_vertices
        mesh%data%n_triangles = n_triangles
        mesh%data%n_quads = 0
        mesh%data%has_triangles = .true.
        mesh%data%has_quads = .false.
        mesh%data%has_mixed_elements = .false.

        call move_alloc(vertices, mesh%data%vertices)
        call move_alloc(triangles, mesh%data%triangles)

        call mesh%data%build_connectivity()
        call mesh%data%find_boundary()
    end function mesh_from_triangle_files

    pure subroutine set_triangle_mesh_metadata(mesh, n_vertices,           &
        n_triangles)
        type(mesh_t), intent(inout) :: mesh
        integer, intent(in) :: n_vertices, n_triangles

        mesh%data%n_vertices = n_vertices
        mesh%data%n_triangles = n_triangles
        mesh%data%n_quads = 0
        mesh%data%has_triangles = .true.
        mesh%data%has_quads = .false.
        mesh%data%has_mixed_elements = .false.
    end subroutine set_triangle_mesh_metadata

    function mesh_from_domain(vertices, segments, hole_points, min_angle)   &
        result(mesh)
        use triangulation_fortran, only: triangulation_result_t,            &
            triangulate_with_hole_fortran,                                  &
            triangulate_with_quality_fortran, cleanup_triangulation
        real(dp), intent(in) :: vertices(:,:)
        integer, intent(in) :: segments(:,:)
        real(dp), intent(in), optional :: hole_points(:,:)
        real(dp), intent(in), optional :: min_angle
        type(mesh_t) :: mesh

        type(triangulation_result_t) :: result
        real(dp) :: angle
        integer :: stat

        call init_measures()

        angle = 20.0_dp
        if (present(min_angle)) angle = min_angle

        if (present(hole_points)) then
            call triangulate_with_hole_fortran(vertices, segments,          &
                                               hole_points, result, stat)
        else
            call triangulate_with_quality_fortran(vertices, segments,       &
                                                  angle, result, stat)
        end if

        if (result%ntriangles == 0) then
            mesh%data%n_vertices = 0
            mesh%data%n_triangles = 0
            return
        end if

        call set_triangle_mesh_metadata(mesh, result%npoints,               &
                                        result%ntriangles)

        allocate(mesh%data%vertices(2, mesh%data%n_vertices))
        allocate(mesh%data%triangles(3, mesh%data%n_triangles))

        mesh%data%vertices = result%points
        mesh%data%triangles = result%triangles

        call cleanup_triangulation(result)

        call mesh%data%build_connectivity()
        call mesh%data%find_boundary()
    end function mesh_from_domain

    function structured_quad_mesh(nx, ny, x0, x1, y0, y1) result(mesh)
        integer, intent(in) :: nx, ny
        real(dp), intent(in) :: x0, x1, y0, y1
        type(mesh_t) :: mesh

        call init_measures()
        call mesh%data%create_structured_quads(nx, ny, x0, x1, y0, y1)
        call mesh%data%build_connectivity()
        call mesh%data%find_boundary()
    end function structured_quad_mesh

    function refine_uniform(mesh, levels) result(refined_mesh)
        type(mesh_t), intent(in) :: mesh
        integer, intent(in), optional :: levels
        type(mesh_t) :: refined_mesh
        type(mesh_t) :: current_mesh, next_mesh
        integer :: l, nlevels

        nlevels = 1
        if (present(levels)) then
            if (levels > 1) nlevels = levels
        end if

        current_mesh = mesh
        do l = 1, nlevels
            call current_mesh%data%refine_uniform(next_mesh%data)
            if (l < nlevels) then
                call current_mesh%destroy()
                current_mesh = next_mesh
            else
                refined_mesh = next_mesh
            end if
        end do
    end function refine_uniform

    function refine_adaptive_markers(mesh, refine_markers) result(refined_mesh)
        type(mesh_t), intent(in) :: mesh
        logical, intent(in) :: refine_markers(:)
        type(mesh_t) :: refined_mesh

        call mesh%data%refine_adaptive(refine_markers, refined_mesh%data)
    end function refine_adaptive_markers

    function refine_adaptive_solution(mesh, solution, tolerance)            &
        result(refined_mesh)
        type(mesh_t), intent(in) :: mesh
        type(function_t), intent(in) :: solution
        real(dp), intent(in) :: tolerance
        type(mesh_t) :: refined_mesh
        real(dp), allocatable :: indicators(:)
        logical, allocatable :: refine_markers(:)
        real(dp) :: max_eta, threshold
        integer :: n

        n = mesh%data%n_triangles
        if (n <= 0) then
            refined_mesh = mesh
            return
        end if

        allocate(indicators(n))
        allocate(refine_markers(n))

        call compute_gradient_indicators(mesh, solution, indicators)

        max_eta = maxval(indicators)
        if (max_eta <= 0.0_dp) then
            refine_markers = .false.
        else
            threshold = tolerance * max_eta
            refine_markers = indicators >= threshold
        end if

        call mesh%data%refine_adaptive(refine_markers, refined_mesh%data)

        deallocate(indicators, refine_markers)
    end function refine_adaptive_solution

    subroutine find_triangle_edges(mesh, triangle_idx, edge1, edge2, edge3)
        type(mesh_2d_t), intent(in) :: mesh
        integer, intent(in) :: triangle_idx
        integer, intent(out) :: edge1, edge2, edge3

        integer :: v1, v2, v3, i
        integer :: e1_v1, e1_v2, e2_v1, e2_v2, e3_v1, e3_v2

        v1 = mesh%triangles(1, triangle_idx)
        v2 = mesh%triangles(2, triangle_idx)
        v3 = mesh%triangles(3, triangle_idx)

        e1_v1 = min(v1, v2)
        e1_v2 = max(v1, v2)
        e2_v1 = min(v2, v3)
        e2_v2 = max(v2, v3)
        e3_v1 = min(v3, v1)
        e3_v2 = max(v3, v1)

        edge1 = 0
        edge2 = 0
        edge3 = 0

        do i = 1, mesh%n_edges
            if (mesh%edges(1, i) == e1_v1 .and.                              &
                mesh%edges(2, i) == e1_v2) then
                edge1 = i
            else if (mesh%edges(1, i) == e2_v1 .and.                         &
                     mesh%edges(2, i) == e2_v2) then
                edge2 = i
            else if (mesh%edges(1, i) == e3_v1 .and.                         &
                     mesh%edges(2, i) == e3_v2) then
                edge3 = i
            end if
        end do

        if (edge1 == 0 .or. edge2 == 0 .or. edge3 == 0) then
            error stop "find_triangle_edges: edges not found for triangle"
        end if
    end subroutine find_triangle_edges

end module fortfem_api_mesh
