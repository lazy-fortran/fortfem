module fortfem_maxwell_rwg_surface
    !! Lowest-order Rao--Wilton--Glisson surface-current space on an
    !! orientable triangular surface. One basis function is created for each
    !! edge shared by exactly two panels.
    use fortfem_kinds, only: dp
    use fortfem_tetra_edge_dof_map, only: build_tetra_edge_dof_map
    use fortfem_triangle_duffy_quadrature, only: triangle_duffy_quadrature
    implicit none
    private

    public :: build_maxwell_rwg_surface_space
    public :: assemble_maxwell_rwg_mass_matrix
    public :: evaluate_maxwell_rwg_basis
    public :: map_maxwell_rwg_to_tetra_nedelec_edges

contains

    subroutine assemble_maxwell_rwg_mass_matrix( &
            vertices, triangles, quadrature_degree, matrix, status)
        real(dp), intent(in) :: vertices(:, :)
        integer, intent(in) :: triangles(:, :), quadrature_degree
        real(dp), allocatable, intent(out) :: matrix(:, :)
        integer, intent(out) :: status

        integer, allocatable :: edge_triangles(:, :), edge_vertices(:, :)
        real(dp), allocatable :: eta(:), weights(:), xi(:)
        real(dp), allocatable :: values(:, :)
        real(dp) :: basis_divergence, jacobian, point(3)
        integer :: basis, node, panel, test_basis

        status = 1
        call build_maxwell_rwg_surface_space( &
            vertices, triangles, edge_vertices, edge_triangles, status)
        if (status /= 0 .or. size(edge_vertices, 2) == 0) return
        call triangle_duffy_quadrature( &
            quadrature_degree, xi, eta, weights, status)
        if (status /= 0) return
        allocate( &
            matrix(size(edge_vertices, 2), size(edge_vertices, 2)), &
            values(3, size(edge_vertices, 2)))
        matrix = 0.0_dp
        do panel = 1, size(triangles, 2)
            jacobian = 2.0_dp*triangle_area( &
                vertices(:, triangles(:, panel)))
            do node = 1, size(weights)
                point = triangle_point( &
                    vertices(:, triangles(:, panel)), xi(node), eta(node))
                values = 0.0_dp
                do basis = 1, size(edge_vertices, 2)
                    if (.not. any(edge_triangles(:, basis) == panel)) cycle
                    call evaluate_maxwell_rwg_basis( &
                        vertices, triangles, edge_vertices, edge_triangles, &
                        basis, panel, point, values(:, basis), basis_divergence, &
                        status)
                    if (status /= 0) return
                end do
                do basis = 1, size(edge_vertices, 2)
                    if (.not. any(edge_triangles(:, basis) == panel)) cycle
                    do test_basis = 1, size(edge_vertices, 2)
                        if (.not. any( &
                            edge_triangles(:, test_basis) == panel)) cycle
                        matrix(test_basis, basis) = matrix(test_basis, basis) + &
                            jacobian*weights(node)*dot_product( &
                            values(:, test_basis), values(:, basis))
                    end do
                end do
            end do
        end do
        status = 0
    end subroutine assemble_maxwell_rwg_mass_matrix

    subroutine map_maxwell_rwg_to_tetra_nedelec_edges( &
            vertices, tetrahedra, rwg_edges, nedelec_dofs, scales, status)
        real(dp), intent(in) :: vertices(:, :)
        integer, intent(in) :: tetrahedra(:, :), rwg_edges(:, :)
        integer, allocatable, intent(out) :: nedelec_dofs(:)
        real(dp), allocatable, intent(out) :: scales(:)
        integer, intent(out) :: status

        integer, allocatable :: edge_dofs(:, :), edge_orientations(:, :)
        integer, allocatable :: volume_edges(:, :)
        integer :: edge, volume_edge

        status = 1
        if (size(vertices, 1) /= 3 .or. size(rwg_edges, 1) /= 2) return
        if (any(rwg_edges < 1) .or. any(rwg_edges > size(vertices, 2))) return
        call build_tetra_edge_dof_map( &
            tetrahedra, volume_edges, edge_dofs, edge_orientations, status)
        if (status /= 0) return
        allocate(nedelec_dofs(size(rwg_edges, 2)), scales(size(rwg_edges, 2)))
        do edge = 1, size(rwg_edges, 2)
            nedelec_dofs(edge) = 0
            do volume_edge = 1, size(volume_edges, 2)
                if (all(volume_edges(:, volume_edge) == rwg_edges(:, edge))) then
                    nedelec_dofs(edge) = volume_edge
                    exit
                end if
            end do
            if (nedelec_dofs(edge) == 0) then
                status = 2
                return
            end if
            scales(edge) = 1.0_dp/norm2( &
                vertices(:, rwg_edges(2, edge)) - &
                vertices(:, rwg_edges(1, edge)))
        end do
        status = 0
    end subroutine map_maxwell_rwg_to_tetra_nedelec_edges

    subroutine build_maxwell_rwg_surface_space( &
            vertices, triangles, edge_vertices, edge_triangles, status)
        real(dp), intent(in) :: vertices(:, :)
        integer, intent(in) :: triangles(:, :)
        integer, allocatable, intent(out) :: edge_vertices(:, :)
        integer, allocatable, intent(out) :: edge_triangles(:, :)
        integer, intent(out) :: status

        integer, allocatable :: candidate_edges(:, :), candidate_triangles(:, :)
        integer, allocatable :: occurrence_count(:)
        integer :: candidate_count, edge, first, local, second, triangle
        integer :: vertex_count

        status = 1
        if (size(vertices, 1) /= 3 .or. size(triangles, 1) /= 3) return
        vertex_count = size(vertices, 2)
        if (vertex_count < 3 .or. size(triangles, 2) < 1) return
        if (any(triangles < 1) .or. any(triangles > vertex_count)) return
        allocate( &
            candidate_edges(2, 3*size(triangles, 2)), &
            candidate_triangles(2, 3*size(triangles, 2)), &
            occurrence_count(3*size(triangles, 2)))
        candidate_count = 0
        occurrence_count = 0
        candidate_triangles = 0
        do triangle = 1, size(triangles, 2)
            if (triangle_area(vertices(:, triangles(:, triangle))) <= &
                tiny(1.0_dp)) return
            do local = 1, 3
                first = triangles(local, triangle)
                second = triangles(modulo(local, 3) + 1, triangle)
                if (first == second) return
                call register_edge( &
                    min(first, second), max(first, second), triangle, &
                    candidate_edges, candidate_triangles, occurrence_count, &
                    candidate_count, status)
                if (status /= 0) return
            end do
        end do
        allocate( &
            edge_vertices(2, count(occurrence_count(:candidate_count) == 2)), &
            edge_triangles(2, count(occurrence_count(:candidate_count) == 2)))
        edge = 0
        do first = 1, candidate_count
            if (occurrence_count(first) /= 2) cycle
            edge = edge + 1
            edge_vertices(:, edge) = candidate_edges(:, first)
            edge_triangles(:, edge) = candidate_triangles(:, first)
        end do
        status = 0
    end subroutine build_maxwell_rwg_surface_space

    subroutine register_edge( &
            first, second, triangle, edges, adjacent, occurrences, count, status)
        integer, intent(in) :: first, second, triangle
        integer, intent(inout) :: edges(:, :), adjacent(:, :), occurrences(:)
        integer, intent(inout) :: count
        integer, intent(out) :: status

        integer :: edge

        status = 0
        do edge = 1, count
            if (edges(1, edge) /= first .or. edges(2, edge) /= second) cycle
            occurrences(edge) = occurrences(edge) + 1
            if (occurrences(edge) > 2) then
                status = 2
                return
            end if
            adjacent(occurrences(edge), edge) = triangle
            return
        end do
        count = count + 1
        edges(:, count) = [first, second]
        occurrences(count) = 1
        adjacent(1, count) = triangle
    end subroutine register_edge

    subroutine evaluate_maxwell_rwg_basis( &
            vertices, triangles, edge_vertices, edge_triangles, basis, &
            triangle, point, value, surface_divergence, status)
        real(dp), intent(in) :: vertices(:, :), point(3)
        integer, intent(in) :: triangles(:, :), edge_vertices(:, :)
        integer, intent(in) :: edge_triangles(:, :), basis, triangle
        real(dp), intent(out) :: value(3), surface_divergence
        integer, intent(out) :: status

        real(dp) :: area, edge_length, panel(3, 3)
        integer :: local, next, opposite, orientation

        value = 0.0_dp
        surface_divergence = 0.0_dp
        status = 1
        if (basis < 1 .or. basis > size(edge_vertices, 2)) return
        if (triangle < 1 .or. triangle > size(triangles, 2)) return
        if (.not. any(edge_triangles(:, basis) == triangle)) return
        panel = vertices(:, triangles(:, triangle))
        if (.not. point_in_triangle(panel, point)) return
        orientation = 0
        opposite = 0
        do local = 1, 3
            next = modulo(local, 3) + 1
            if (triangles(local, triangle) == edge_vertices(1, basis) .and. &
                triangles(next, triangle) == edge_vertices(2, basis)) then
                orientation = 1
                opposite = modulo(next, 3) + 1
                exit
            end if
            if (triangles(local, triangle) == edge_vertices(2, basis) .and. &
                triangles(next, triangle) == edge_vertices(1, basis)) then
                orientation = -1
                opposite = modulo(next, 3) + 1
                exit
            end if
        end do
        if (orientation == 0) return
        area = triangle_area(panel)
        edge_length = norm2( &
            vertices(:, edge_vertices(2, basis)) - &
            vertices(:, edge_vertices(1, basis)))
        value = real(orientation, dp)*edge_length/(2.0_dp*area)* &
            (point - panel(:, opposite))
        surface_divergence = real(orientation, dp)*edge_length/area
        status = 0
    end subroutine evaluate_maxwell_rwg_basis

    pure logical function point_in_triangle(vertices, point) result(inside)
        real(dp), intent(in) :: vertices(3, 3), point(3)

        real(dp) :: displacement(3), dot00, dot01, dot02, dot11, dot12
        real(dp) :: denominator, normal(3), u, v

        normal = cross_product( &
            vertices(:, 2) - vertices(:, 1), &
            vertices(:, 3) - vertices(:, 1))
        displacement = point - vertices(:, 1)
        if (abs(dot_product(displacement, normal)) > &
            128.0_dp*epsilon(1.0_dp)*max(1.0_dp, norm2(normal))) then
            inside = .false.
            return
        end if
        dot00 = dot_product( &
            vertices(:, 2) - vertices(:, 1), &
            vertices(:, 2) - vertices(:, 1))
        dot01 = dot_product( &
            vertices(:, 2) - vertices(:, 1), &
            vertices(:, 3) - vertices(:, 1))
        dot02 = dot_product(vertices(:, 2) - vertices(:, 1), displacement)
        dot11 = dot_product( &
            vertices(:, 3) - vertices(:, 1), &
            vertices(:, 3) - vertices(:, 1))
        dot12 = dot_product(vertices(:, 3) - vertices(:, 1), displacement)
        denominator = dot00*dot11 - dot01*dot01
        if (denominator <= tiny(1.0_dp)) then
            inside = .false.
            return
        end if
        u = (dot11*dot02 - dot01*dot12)/denominator
        v = (dot00*dot12 - dot01*dot02)/denominator
        inside = u >= -128.0_dp*epsilon(1.0_dp) .and. &
            v >= -128.0_dp*epsilon(1.0_dp) .and. &
            u + v <= 1.0_dp + 128.0_dp*epsilon(1.0_dp)
    end function point_in_triangle

    pure function triangle_area(vertices) result(area)
        real(dp), intent(in) :: vertices(3, 3)
        real(dp) :: area

        area = 0.5_dp*norm2(cross_product( &
            vertices(:, 2) - vertices(:, 1), &
            vertices(:, 3) - vertices(:, 1)))
    end function triangle_area

    pure function triangle_point(vertices, xi, eta) result(point)
        real(dp), intent(in) :: vertices(3, 3), xi, eta
        real(dp) :: point(3)

        point = vertices(:, 1) + xi*(vertices(:, 2) - vertices(:, 1)) + &
            eta*(vertices(:, 3) - vertices(:, 1))
    end function triangle_point

    pure function cross_product(first, second) result(product)
        real(dp), intent(in) :: first(3), second(3)
        real(dp) :: product(3)

        product = [ &
            first(2)*second(3) - first(3)*second(2), &
            first(3)*second(1) - first(1)*second(3), &
            first(1)*second(2) - first(2)*second(1)]
    end function cross_product

end module fortfem_maxwell_rwg_surface
