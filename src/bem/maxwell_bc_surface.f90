module fortfem_maxwell_bc_surface
    !! Buffa--Christiansen basis transformation on closed triangular surfaces.
    !!
    !! The columns express one BC function per primal RWG edge in the
    !! element-local RWG basis of the barycentric refinement. The construction
    !! follows the coefficient map of Buffa--Christiansen and Bempp-cl.
    use fortfem_barycentric_surface_refinement, only: &
        barycentric_refine_surface_mesh
    use fortfem_kinds, only: dp
    use fortfem_maxwell_localized_rwg_surface, only: &
        evaluate_maxwell_localized_rwg_basis
    use fortfem_maxwell_rwg_surface, only: build_maxwell_rwg_surface_space
    use fortfem_maxwell_rwg_surface, only: evaluate_maxwell_rwg_basis
    use fortfem_triangle_duffy_quadrature, only: triangle_duffy_quadrature
    implicit none
    private

    public :: build_maxwell_bc_transformation
    public :: assemble_maxwell_rwg_rbc_pairing

contains

    subroutine assemble_maxwell_rwg_rbc_pairing( &
            vertices, triangles, quadrature_degree, matrix, status)
        real(dp), intent(in) :: vertices(:, :)
        integer, intent(in) :: triangles(:, :), quadrature_degree
        real(dp), allocatable, intent(out) :: matrix(:, :)
        integer, intent(out) :: status

        integer, allocatable :: edge_panels(:, :), edge_vertices(:, :)
        integer, allocatable :: refined_triangles(:, :)
        real(dp), allocatable :: eta(:), refined_vertices(:, :)
        real(dp), allocatable :: transformation(:, :), weights(:), xi(:)
        real(dp), allocatable :: bc_values(:, :), rwg_values(:, :)
        real(dp) :: divergence, jacobian, local_value(3), normal(3), point(3)
        integer :: basis, local_edge, node, parent, refined_panel, row
        integer :: test_basis

        call build_maxwell_bc_transformation( &
            vertices, triangles, refined_vertices, refined_triangles, &
            transformation, status)
        if (status /= 0) return
        call build_maxwell_rwg_surface_space( &
            vertices, triangles, edge_vertices, edge_panels, status)
        if (status /= 0) return
        call triangle_duffy_quadrature( &
            quadrature_degree, xi, eta, weights, status)
        if (status /= 0) return
        allocate(matrix(size(edge_vertices, 2), size(edge_vertices, 2)))
        allocate( &
            bc_values(3, size(edge_vertices, 2)), &
            rwg_values(3, size(edge_vertices, 2)))
        matrix = 0.0_dp
        do refined_panel = 1, size(refined_triangles, 2)
            parent = (refined_panel - 1)/6 + 1
            normal = triangle_normal( &
                refined_vertices(:, refined_triangles(:, refined_panel)))
            jacobian = 2.0_dp*triangle_area( &
                refined_vertices(:, refined_triangles(:, refined_panel)))
            do node = 1, size(weights)
                point = triangle_point( &
                    refined_vertices(:, refined_triangles(:, refined_panel)), &
                    xi(node), eta(node))
                bc_values = 0.0_dp
                do local_edge = 1, 3
                    call evaluate_maxwell_localized_rwg_basis( &
                        refined_vertices, refined_triangles, refined_panel, &
                        local_edge, point, local_value, divergence, status)
                    if (status /= 0) return
                    row = 3*(refined_panel - 1) + local_edge
                    do test_basis = 1, size(edge_vertices, 2)
                        bc_values(:, test_basis) = &
                            bc_values(:, test_basis) + &
                            transformation(row, test_basis)*local_value
                    end do
                end do
                rwg_values = 0.0_dp
                do basis = 1, size(edge_vertices, 2)
                    if (.not. any(edge_panels(:, basis) == parent)) cycle
                    call evaluate_maxwell_rwg_basis( &
                        vertices, triangles, edge_vertices, edge_panels, &
                        basis, parent, point, rwg_values(:, basis), divergence, &
                        status)
                    if (status /= 0) return
                end do
                do test_basis = 1, size(edge_vertices, 2)
                    do basis = 1, size(edge_vertices, 2)
                        matrix(test_basis, basis) = matrix(test_basis, basis) + &
                            jacobian*weights(node)*dot_product( &
                            cross_product(normal, bc_values(:, test_basis)), &
                            rwg_values(:, basis))
                    end do
                end do
            end do
        end do
        status = 0
    end subroutine assemble_maxwell_rwg_rbc_pairing

    subroutine build_maxwell_bc_transformation( &
            vertices, triangles, refined_vertices, refined_triangles, &
            transformation, status, sphere_radius)
        real(dp), intent(in) :: vertices(:, :)
        integer, intent(in) :: triangles(:, :)
        real(dp), allocatable, intent(out) :: refined_vertices(:, :)
        integer, allocatable, intent(out) :: refined_triangles(:, :)
        real(dp), allocatable, intent(out) :: transformation(:, :)
        integer, intent(out) :: status
        real(dp), optional, intent(in) :: sphere_radius

        integer, allocatable :: edge_panels(:, :), primal_edges(:, :)
        integer, allocatable :: refined_element_edges(:, :)
        integer, allocatable :: refined_edges(:, :)
        integer :: basis, lower, upper, vertex1, vertex2
        integer :: local_vertex1, local_vertex2
        integer :: upper_minus, upper_plus, lower_minus, lower_plus
        integer :: refined_vertex

        status = 1
        call build_maxwell_rwg_surface_space( &
            vertices, triangles, primal_edges, edge_panels, status)
        if (status /= 0) return
        if (size(primal_edges, 2) /= 3*size(triangles, 2)/2) then
            status = 2
            return
        end if
        call barycentric_refine_surface_mesh( &
            vertices, triangles, refined_vertices, refined_triangles)
        if (present(sphere_radius)) then
            if (sphere_radius <= 0.0_dp) return
            do refined_vertex = 1, size(refined_vertices, 2)
                if (norm2(refined_vertices(:, refined_vertex)) <= &
                    tiny(1.0_dp)) return
                refined_vertices(:, refined_vertex) = sphere_radius* &
                    refined_vertices(:, refined_vertex)/ &
                    norm2(refined_vertices(:, refined_vertex))
            end do
        end if
        call enumerate_local_edges( &
            refined_triangles, refined_edges, refined_element_edges)
        allocate(transformation(3*size(refined_triangles, 2), &
            size(primal_edges, 2)))
        transformation = 0.0_dp
        do basis = 1, size(primal_edges, 2)
            call upper_lower_panels( &
                triangles, primal_edges(:, basis), edge_panels(:, basis), &
                upper, lower)
            vertex1 = primal_edges(1, basis)
            vertex2 = primal_edges(2, basis)
            local_vertex1 = find_local_vertex(triangles(:, upper), vertex1)
            if (triangles(modulo(local_vertex1 - 2, 3) + 1, upper) == &
                vertex2) then
                vertex1 = primal_edges(2, basis)
                vertex2 = primal_edges(1, basis)
            end if
            local_vertex1 = find_local_vertex(triangles(:, upper), vertex1)
            local_vertex2 = find_local_vertex(triangles(:, lower), vertex2)
            upper_minus = 6*(upper - 1) + 2*(local_vertex1 - 1) + 1
            upper_plus = upper_minus + 1
            lower_minus = 6*(lower - 1) + 2*(local_vertex2 - 1) + 1
            lower_plus = lower_minus + 1
            call add_vertex_ring( &
                vertex1, upper_minus, -1.0_dp, basis, refined_vertices, &
                refined_triangles, refined_element_edges, &
                transformation, status)
            if (status /= 0) return
            call add_vertex_ring( &
                vertex2, lower_minus, 1.0_dp, basis, refined_vertices, &
                refined_triangles, refined_element_edges, &
                transformation, status)
            if (status /= 0) return
            call add_reference_edge( &
                upper_minus, upper_plus, lower_minus, lower_plus, basis, &
                refined_vertices, refined_triangles, refined_element_edges, &
                transformation)
        end do
        status = 0
    end subroutine build_maxwell_bc_transformation

    subroutine upper_lower_panels( &
            triangles, edge, panels, upper, lower)
        integer, intent(in) :: triangles(:, :), edge(2), panels(2)
        integer, intent(out) :: upper, lower

        integer :: local, next, orientation

        orientation = 0
        do local = 1, 3
            next = modulo(local, 3) + 1
            if (triangles(local, panels(1)) == edge(1) .and. &
                triangles(next, panels(1)) == edge(2)) orientation = 1
            if (triangles(local, panels(1)) == edge(2) .and. &
                triangles(next, panels(1)) == edge(1)) orientation = -1
        end do
        if (orientation > 0) then
            lower = panels(1)
            upper = panels(2)
        else
            lower = panels(2)
            upper = panels(1)
        end if
    end subroutine upper_lower_panels

    subroutine add_vertex_ring( &
            vertex, start_element, initial_sign, basis, vertices, triangles, &
            element_edges, transformation, status)
        integer, intent(in) :: vertex, start_element, basis
        real(dp), intent(in) :: initial_sign, vertices(:, :)
        integer, intent(in) :: triangles(:, :), element_edges(:, :)
        real(dp), intent(inout) :: transformation(:, :)
        integer, intent(out) :: status

        integer, allocatable :: ring_elements(:), ring_first(:), ring_second(:)
        integer, allocatable :: entry_elements(:), entry_edges(:)
        integer :: count_index, entry, item, local_edge, nc, row
        real(dp) :: edge_length, sign

        call ordered_vertex_ring( &
            vertex, start_element, triangles, element_edges, ring_elements, &
            ring_first, ring_second, status)
        if (status /= 0) return
        allocate( &
            entry_elements(2*size(ring_elements) - 2), &
            entry_edges(2*size(ring_elements) - 2))
        entry = 0
        do item = 1, size(ring_elements)
            if (item == 1) then
                entry = entry + 1
                entry_elements(entry) = ring_elements(item)
                entry_edges(entry) = ring_second(item)
            else if (item == size(ring_elements)) then
                entry = entry + 1
                entry_elements(entry) = ring_elements(item)
                entry_edges(entry) = ring_first(item)
            else
                entry = entry + 1
                entry_elements(entry) = ring_elements(item)
                entry_edges(entry) = ring_first(item)
                entry = entry + 1
                entry_elements(entry) = ring_elements(item)
                entry_edges(entry) = ring_second(item)
            end if
        end do
        nc = size(ring_elements)/2
        count_index = 0
        sign = initial_sign
        do entry = 1, size(entry_elements)
            if (modulo(entry, 2) == 1) count_index = count_index + 1
            item = entry_elements(entry)
            local_edge = entry_edges(entry)
            call add_ring_entry()
        end do
        status = 0

    contains

        subroutine add_ring_entry()
            edge_length = local_edge_length( &
                vertices, triangles(:, item), local_edge)
            row = 3*(item - 1) + local_edge
            transformation(row, basis) = transformation(row, basis) + &
                sign*real(nc - count_index, dp)/ &
                (2.0_dp*real(nc, dp)*edge_length)
            sign = -sign
        end subroutine add_ring_entry

    end subroutine add_vertex_ring

    subroutine ordered_vertex_ring( &
            vertex, start_element, triangles, element_edges, ring_elements, &
            ring_first, ring_second, status)
        integer, intent(in) :: vertex, start_element
        integer, intent(in) :: triangles(:, :), element_edges(:, :)
        integer, allocatable, intent(out) :: ring_elements(:)
        integer, allocatable, intent(out) :: ring_first(:), ring_second(:)
        integer, intent(out) :: status

        integer, allocatable :: candidates(:), first_edges(:), second_edges(:)
        logical, allocatable :: used(:)
        integer :: candidate, current, local, next, ring_count

        ring_count = count(triangles == vertex)
        allocate( &
            candidates(ring_count), first_edges(ring_count), &
            second_edges(ring_count), used(ring_count), &
            ring_elements(ring_count), ring_first(ring_count), &
            ring_second(ring_count))
        ring_count = 0
        do candidate = 1, size(triangles, 2)
            local = find_local_vertex(triangles(:, candidate), vertex)
            if (local == 0) cycle
            ring_count = ring_count + 1
            candidates(ring_count) = candidate
            call incident_local_edges(local, first_edges(ring_count), &
                second_edges(ring_count))
        end do
        used = .false.
        current = find_index(candidates, start_element)
        status = 1
        if (current == 0) return
        do local = 1, ring_count
            ring_elements(local) = candidates(current)
            ring_first(local) = first_edges(current)
            ring_second(local) = second_edges(current)
            used(current) = .true.
            if (local == ring_count) exit
            next = 0
            do candidate = 1, ring_count
                if (used(candidate)) cycle
                if (element_edges(first_edges(candidate), candidates(candidate)) &
                    == element_edges(second_edges(current), &
                    candidates(current))) then
                    next = candidate
                    exit
                end if
            end do
            if (next == 0) return
            current = next
        end do
        status = 0
    end subroutine ordered_vertex_ring

    pure subroutine incident_local_edges(local_vertex, first, second)
        integer, intent(in) :: local_vertex
        integer, intent(out) :: first, second

        select case (local_vertex)
        case (1)
            first = 1
            second = 2
        case (2)
            first = 3
            second = 1
        case default
            first = 2
            second = 3
        end select
    end subroutine incident_local_edges

    subroutine add_reference_edge( &
            upper_minus, upper_plus, lower_minus, lower_plus, basis, vertices, &
            triangles, element_edges, transformation)
        integer, intent(in) :: upper_minus, upper_plus, lower_minus, lower_plus
        integer, intent(in) :: basis, triangles(:, :), element_edges(:, :)
        real(dp), intent(in) :: vertices(:, :)
        real(dp), intent(inout) :: transformation(:, :)

        real(dp) :: lower_length, upper_length

        associate(unused => element_edges)
        end associate
        upper_length = local_edge_length( &
            vertices, triangles(:, upper_minus), 3)
        lower_length = local_edge_length( &
            vertices, triangles(:, lower_minus), 3)
        transformation(3*(upper_minus - 1) + 3, basis) = &
            1.0_dp/(2.0_dp*upper_length)
        transformation(3*(upper_plus - 1) + 3, basis) = &
            -1.0_dp/(2.0_dp*upper_length)
        transformation(3*(lower_minus - 1) + 3, basis) = &
            -1.0_dp/(2.0_dp*lower_length)
        transformation(3*(lower_plus - 1) + 3, basis) = &
            1.0_dp/(2.0_dp*lower_length)
    end subroutine add_reference_edge

    subroutine enumerate_local_edges(triangles, edges, element_edges)
        integer, intent(in) :: triangles(:, :)
        integer, allocatable, intent(out) :: edges(:, :), element_edges(:, :)

        integer :: candidate(2), edge_count, element, local_edge, found
        integer, allocatable :: storage(:, :)

        allocate(storage(2, 3*size(triangles, 2)))
        allocate(element_edges(3, size(triangles, 2)))
        edge_count = 0
        do element = 1, size(triangles, 2)
            do local_edge = 1, 3
                candidate = local_edge_vertex_ids( &
                    triangles(:, element), local_edge)
                candidate = [minval(candidate), maxval(candidate)]
                found = find_edge(storage, edge_count, candidate)
                if (found == 0) then
                    edge_count = edge_count + 1
                    storage(:, edge_count) = candidate
                    found = edge_count
                end if
                element_edges(local_edge, element) = found
            end do
        end do
        allocate(edges(2, edge_count))
        edges = storage(:, :edge_count)
    end subroutine enumerate_local_edges

    pure function local_edge_vertex_ids(element, local_edge) result(ids)
        integer, intent(in) :: element(3), local_edge
        integer :: ids(2)

        select case (local_edge)
        case (1)
            ids = element([1, 2])
        case (2)
            ids = element([3, 1])
        case default
            ids = element([2, 3])
        end select
    end function local_edge_vertex_ids

    pure real(dp) function local_edge_length( &
            vertices, element, local_edge) result(length)
        real(dp), intent(in) :: vertices(:, :)
        integer, intent(in) :: element(3), local_edge
        integer :: ids(2)

        ids = local_edge_vertex_ids(element, local_edge)
        length = norm2(vertices(:, ids(2)) - vertices(:, ids(1)))
    end function local_edge_length

    pure integer function find_local_vertex(element, vertex) result(local)
        integer, intent(in) :: element(3), vertex

        do local = 1, 3
            if (element(local) == vertex) return
        end do
        local = 0
    end function find_local_vertex

    pure integer function find_edge(edges, edge_count, target) result(index)
        integer, intent(in) :: edges(:, :), edge_count, target(2)

        do index = 1, edge_count
            if (all(edges(:, index) == target)) return
        end do
        index = 0
    end function find_edge

    pure integer function find_index(values, target) result(index)
        integer, intent(in) :: values(:), target

        do index = 1, size(values)
            if (values(index) == target) return
        end do
        index = 0
    end function find_index

    pure function triangle_normal(points) result(normal)
        real(dp), intent(in) :: points(3, 3)
        real(dp) :: normal(3)

        normal = cross_product( &
            points(:, 2) - points(:, 1), points(:, 3) - points(:, 1))
        normal = normal/norm2(normal)
    end function triangle_normal

    pure function triangle_area(points) result(area)
        real(dp), intent(in) :: points(3, 3)
        real(dp) :: area

        area = 0.5_dp*norm2(cross_product( &
            points(:, 2) - points(:, 1), points(:, 3) - points(:, 1)))
    end function triangle_area

    pure function triangle_point(points, xi, eta) result(point)
        real(dp), intent(in) :: points(3, 3), xi, eta
        real(dp) :: point(3)

        point = points(:, 1) + xi*(points(:, 2) - points(:, 1)) + &
            eta*(points(:, 3) - points(:, 1))
    end function triangle_point

    pure function cross_product(first, second) result(product)
        real(dp), intent(in) :: first(3), second(3)
        real(dp) :: product(3)

        product = [ &
            first(2)*second(3) - first(3)*second(2), &
            first(3)*second(1) - first(1)*second(3), &
            first(1)*second(2) - first(2)*second(1)]
    end function cross_product

end module fortfem_maxwell_bc_surface
