module fortfem_adaptive_surface_bem
    use fortfem_kinds, only: dp
    use fortfem_laplace_galerkin_3d, only: &
        assemble_laplace_single_layer_p0_3d
    implicit none
    private

    public :: estimate_laplace_p0_two_level_residual_3d
    public :: mark_bem_dorfler
    public :: refine_surface_mesh_marked

contains

    subroutine refine_surface_mesh_marked( &
            vertices, triangles, marked, refined_vertices, refined_triangles, &
            parent, status)
        real(dp), intent(in) :: vertices(:, :)
        integer, intent(in) :: triangles(:, :)
        logical, intent(in) :: marked(:)
        real(dp), allocatable, intent(out) :: refined_vertices(:, :)
        integer, allocatable, intent(out) :: refined_triangles(:, :)
        integer, allocatable, intent(out) :: parent(:)
        integer, intent(out) :: status

        integer, allocatable :: edge_midpoint(:), edge_vertices(:, :)
        integer :: a, b, c, child, edge, edge_count, local, midpoint(3)
        integer :: new_vertex_count, split_count, triangle, unsplit

        status = 1
        if (size(vertices, 1) /= 3 .or. size(triangles, 1) /= 3) return
        if (size(marked) /= size(triangles, 2)) return
        if (any(triangles < 1) .or. any(triangles > size(vertices, 2))) return
        allocate(edge_vertices(2, 3*size(triangles, 2)))
        allocate(edge_midpoint(3*size(triangles, 2)), source=0)
        edge_count = 0
        do triangle = 1, size(triangles, 2)
            do local = 1, 3
                a = triangles(local, triangle)
                b = triangles(modulo(local, 3) + 1, triangle)
                call register_edge( &
                    min(a, b), max(a, b), edge_vertices, edge_count, edge)
                if (marked(triangle)) edge_midpoint(edge) = -1
            end do
        end do

        new_vertex_count = size(vertices, 2) + count(edge_midpoint(:edge_count) < 0)
        allocate(refined_vertices(3, new_vertex_count))
        refined_vertices(:, :size(vertices, 2)) = vertices
        new_vertex_count = size(vertices, 2)
        do edge = 1, edge_count
            if (edge_midpoint(edge) >= 0) cycle
            new_vertex_count = new_vertex_count + 1
            edge_midpoint(edge) = new_vertex_count
            refined_vertices(:, new_vertex_count) = 0.5_dp*( &
                vertices(:, edge_vertices(1, edge)) + &
                vertices(:, edge_vertices(2, edge)))
        end do

        child = 0
        do triangle = 1, size(triangles, 2)
            call triangle_midpoints( &
                triangles(:, triangle), edge_vertices(:, :edge_count), &
                edge_midpoint(:edge_count), midpoint)
            child = child + 1 + count(midpoint > 0)
        end do
        allocate(refined_triangles(3, child), parent(child))
        child = 0
        do triangle = 1, size(triangles, 2)
            call triangle_midpoints( &
                triangles(:, triangle), edge_vertices(:, :edge_count), &
                edge_midpoint(:edge_count), midpoint)
            split_count = count(midpoint > 0)
            select case (split_count)
            case (0)
                call add_child(triangles(:, triangle))
            case (1)
                local = findloc(midpoint > 0, .true., dim=1)
                a = triangles(local, triangle)
                b = triangles(modulo(local, 3) + 1, triangle)
                c = triangles(modulo(local + 1, 3) + 1, triangle)
                call add_child([a, midpoint(local), c])
                call add_child([midpoint(local), b, c])
            case (2)
                unsplit = findloc(midpoint == 0, .true., dim=1)
                a = triangles(unsplit, triangle)
                b = triangles(modulo(unsplit, 3) + 1, triangle)
                c = triangles(modulo(unsplit + 1, 3) + 1, triangle)
                call add_child([ &
                    a, b, midpoint(modulo(unsplit + 1, 3) + 1)])
                call add_child([ &
                    b, midpoint(modulo(unsplit, 3) + 1), &
                    midpoint(modulo(unsplit + 1, 3) + 1)])
                call add_child([ &
                    midpoint(modulo(unsplit, 3) + 1), c, &
                    midpoint(modulo(unsplit + 1, 3) + 1)])
            case (3)
                a = triangles(1, triangle)
                b = triangles(2, triangle)
                c = triangles(3, triangle)
                call add_child([a, midpoint(1), midpoint(3)])
                call add_child([midpoint(1), b, midpoint(2)])
                call add_child([midpoint(3), midpoint(2), c])
                call add_child([midpoint(1), midpoint(2), midpoint(3)])
            end select
        end do
        status = 0

    contains

        subroutine add_child(child_vertices)
            integer, intent(in) :: child_vertices(3)

            child = child + 1
            refined_triangles(:, child) = child_vertices
            parent(child) = triangle
        end subroutine add_child

    end subroutine refine_surface_mesh_marked

    subroutine estimate_laplace_p0_two_level_residual_3d( &
            vertices, triangles, density, boundary_value, quadrature_degree, &
            indicators, status)
        real(dp), intent(in) :: vertices(:, :), density(:), boundary_value
        integer, intent(in) :: triangles(:, :), quadrature_degree
        real(dp), allocatable, intent(out) :: indicators(:)
        integer, intent(out) :: status

        integer, allocatable :: parent(:), refined_triangles(:, :)
        real(dp), allocatable :: matrix(:, :), prolonged(:), residual(:)
        real(dp), allocatable :: refined_vertices(:, :), right_hand_side(:)
        logical, allocatable :: marked(:)
        real(dp) :: area
        integer :: child

        status = 1
        if (size(density) /= size(triangles, 2)) return
        if (quadrature_degree < 1) return
        allocate(marked(size(triangles, 2)), source=.true.)
        call refine_surface_mesh_marked( &
            vertices, triangles, marked, refined_vertices, refined_triangles, &
            parent, status)
        if (status /= 0) return
        call assemble_laplace_single_layer_p0_3d( &
            refined_vertices, refined_triangles, quadrature_degree, matrix, &
            status)
        if (status /= 0) return
        allocate(prolonged(size(parent)), residual(size(parent)))
        allocate(right_hand_side(size(parent)))
        do child = 1, size(parent)
            prolonged(child) = density(parent(child))
            area = triangle_area( &
                refined_vertices(:, refined_triangles(:, child)))
            right_hand_side(child) = boundary_value*area
        end do
        residual = right_hand_side - matmul(matrix, prolonged)
        allocate(indicators(size(triangles, 2)), source=0.0_dp)
        do child = 1, size(parent)
            area = triangle_area( &
                refined_vertices(:, refined_triangles(:, child)))
            indicators(parent(child)) = indicators(parent(child)) + &
                residual(child)**2/max(area, tiny(1.0_dp))
        end do
        indicators = sqrt(indicators)
        status = 0
    end subroutine estimate_laplace_p0_two_level_residual_3d

    subroutine mark_bem_dorfler(indicators, theta, marked, status)
        real(dp), intent(in) :: indicators(:), theta
        logical, allocatable, intent(out) :: marked(:)
        integer, intent(out) :: status

        real(dp) :: accumulated, target, value
        integer :: panel

        status = 1
        if (size(indicators) < 1 .or. any(indicators < 0.0_dp)) return
        if (theta <= 0.0_dp .or. theta > 1.0_dp) return
        allocate(marked(size(indicators)), source=.false.)
        target = theta*sum(indicators**2)
        if (target <= tiny(1.0_dp)) then
            status = 0
            return
        end if
        accumulated = 0.0_dp
        do while (accumulated < target)
            value = -1.0_dp
            panel = 0
            call find_largest_unmarked(indicators, marked, panel, value)
            if (panel == 0) return
            marked(panel) = .true.
            accumulated = accumulated + value**2
        end do
        status = 0
    end subroutine mark_bem_dorfler

    pure subroutine find_largest_unmarked( &
            indicators, marked, panel, value)
        real(dp), intent(in) :: indicators(:)
        logical, intent(in) :: marked(:)
        integer, intent(out) :: panel
        real(dp), intent(out) :: value
        integer :: candidate

        panel = 0
        value = -1.0_dp
        do candidate = 1, size(indicators)
            if (marked(candidate)) cycle
            if (indicators(candidate) <= value) cycle
            panel = candidate
            value = indicators(candidate)
        end do
    end subroutine find_largest_unmarked

    subroutine register_edge( &
            first, second, edges, edge_count, edge)
        integer, intent(in) :: first, second
        integer, intent(inout) :: edges(:, :), edge_count
        integer, intent(out) :: edge

        do edge = 1, edge_count
            if (edges(1, edge) == first .and. edges(2, edge) == second) return
        end do
        edge_count = edge_count + 1
        edge = edge_count
        edges(:, edge) = [first, second]
    end subroutine register_edge

    subroutine triangle_midpoints( &
            triangle, edges, edge_midpoint, midpoint)
        integer, intent(in) :: triangle(3), edges(:, :), edge_midpoint(:)
        integer, intent(out) :: midpoint(3)
        integer :: edge, first, local, second

        midpoint = 0
        do local = 1, 3
            first = min(triangle(local), triangle(modulo(local, 3) + 1))
            second = max(triangle(local), triangle(modulo(local, 3) + 1))
            do edge = 1, size(edges, 2)
                if (edges(1, edge) /= first .or. &
                    edges(2, edge) /= second) cycle
                midpoint(local) = edge_midpoint(edge)
                exit
            end do
        end do
    end subroutine triangle_midpoints

    pure function triangle_area(vertices) result(area)
        real(dp), intent(in) :: vertices(3, 3)
        real(dp) :: area

        area = 0.5_dp*norm2(cross_product( &
            vertices(:, 2) - vertices(:, 1), &
            vertices(:, 3) - vertices(:, 1)))
    end function triangle_area

    pure function cross_product(first, second) result(product)
        real(dp), intent(in) :: first(3), second(3)
        real(dp) :: product(3)

        product = [ &
            first(2)*second(3) - first(3)*second(2), &
            first(3)*second(1) - first(1)*second(3), &
            first(1)*second(2) - first(2)*second(1)]
    end function cross_product

end module fortfem_adaptive_surface_bem
