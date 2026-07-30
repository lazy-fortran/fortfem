module fortfem_maxwell_surface_rt
    !! Arbitrary-order divergence-conforming currents on parametric panels.
    use fortfem_kinds, only: dp
    use fortfem_triangle_duffy_quadrature, only: triangle_duffy_quadrature
    use fortfem_triangle_rt_arbitrary_order, only: &
        evaluate_triangle_raviart_thomas, initialize_triangle_raviart_thomas, &
        triangle_rt_basis_t, triangle_rt_dof_count
    implicit none
    private

    public :: assemble_maxwell_surface_rt_mass_matrix
    public :: build_maxwell_surface_rt_dof_map
    public :: evaluate_maxwell_surface_rt_basis
    public :: evaluate_maxwell_surface_rt_global_basis

contains

    subroutine assemble_maxwell_surface_rt_mass_matrix( &
            vertices, triangles, degree, quadrature_degree, matrix, status)
        real(dp), intent(in) :: vertices(:, :)
        integer, intent(in) :: triangles(:, :), degree, quadrature_degree
        real(dp), allocatable, intent(out) :: matrix(:, :)
        integer, intent(out) :: status

        type(triangle_rt_basis_t) :: basis
        integer, allocatable :: edge_vertices(:, :), global_dofs(:, :)
        integer, allocatable :: transforms(:, :)
        real(dp), allocatable :: divergences(:), eta(:), values(:, :)
        real(dp), allocatable :: weights(:), xi(:)
        real(dp) :: jacobian, tangent_eta(3), tangent_xi(3)
        integer :: column, global_column, global_row, local_dof_count
        integer :: node, row, triangle

        status = 1
        if (size(vertices, 1) /= 3 .or. size(triangles, 1) /= 3) return
        if (degree < 0 .or. quadrature_degree < 0) return
        if (any(triangles < 1) .or. any(triangles > size(vertices, 2))) return
        call initialize_triangle_raviart_thomas(degree, basis, status)
        if (status /= 0) return
        call build_maxwell_surface_rt_dof_map( &
            triangles, degree, edge_vertices, global_dofs, transforms, status)
        if (status /= 0) return
        call triangle_duffy_quadrature( &
            quadrature_degree, xi, eta, weights, status)
        if (status /= 0) return
        local_dof_count = triangle_rt_dof_count(basis)
        allocate( &
            matrix(maxval(global_dofs), maxval(global_dofs)), &
            values(3, local_dof_count), divergences(local_dof_count))
        matrix = 0.0_dp
        do triangle = 1, size(triangles, 2)
            tangent_xi = &
                vertices(:, triangles(2, triangle)) - &
                vertices(:, triangles(1, triangle))
            tangent_eta = &
                vertices(:, triangles(3, triangle)) - &
                vertices(:, triangles(1, triangle))
            jacobian = norm2(vector_cross_product(tangent_xi, tangent_eta))
            if (jacobian <= tiny(1.0_dp)) return
            do node = 1, size(weights)
                call evaluate_maxwell_surface_rt_basis( &
                    basis, xi(node), eta(node), tangent_xi, tangent_eta, &
                    jacobian, values, divergences, status)
                if (status /= 0) return
                do column = 1, local_dof_count
                    global_column = global_dofs(column, triangle)
                    do row = 1, local_dof_count
                        global_row = global_dofs(row, triangle)
                        matrix(global_row, global_column) = &
                            matrix(global_row, global_column) + &
                            jacobian*weights(node)* &
                            real( &
                            transforms(row, triangle)* &
                            transforms(column, triangle), dp)* &
                            dot_product(values(:, row), values(:, column))
                    end do
                end do
            end do
        end do
        status = 0
    end subroutine assemble_maxwell_surface_rt_mass_matrix

    subroutine build_maxwell_surface_rt_dof_map( &
            triangles, degree, edge_vertices, global_dofs, transforms, status)
        integer, intent(in) :: triangles(:, :), degree
        integer, allocatable, intent(out) :: edge_vertices(:, :)
        integer, allocatable, intent(out) :: global_dofs(:, :), transforms(:, :)
        integer, intent(out) :: status

        integer, allocatable :: candidates(:, :), occurrences(:)
        integer :: cell_dof, cell_dof_count, edge, edge_count
        integer :: edge_index, edge_moment_count, first, local_dof
        integer :: local_dof_count, moment, second, triangle

        status = 1
        if (degree < 0 .or. size(triangles, 1) /= 3) return
        if (size(triangles, 2) < 1 .or. any(triangles < 1)) return
        edge_moment_count = degree + 1
        cell_dof_count = degree*(degree + 1)
        local_dof_count = edge_moment_count*(degree + 3)
        allocate( &
            candidates(2, 3*size(triangles, 2)), &
            occurrences(3*size(triangles, 2)))
        edge_count = 0
        occurrences = 0
        do triangle = 1, size(triangles, 2)
            do edge = 1, 3
                first = triangles(edge, triangle)
                second = triangles(modulo(edge, 3) + 1, triangle)
                if (first == second) return
                call register_surface_edge( &
                    min(first, second), max(first, second), candidates, &
                    occurrences, edge_count, edge_index)
                if (occurrences(edge_index) > 2) return
            end do
        end do
        allocate( &
            edge_vertices(2, edge_count), &
            global_dofs(local_dof_count, size(triangles, 2)), &
            transforms(local_dof_count, size(triangles, 2)))
        edge_vertices = candidates(:, :edge_count)
        global_dofs = 0
        transforms = 1
        do triangle = 1, size(triangles, 2)
            do edge = 1, 3
                first = triangles(edge, triangle)
                second = triangles(modulo(edge, 3) + 1, triangle)
                edge_index = find_surface_edge( &
                    min(first, second), max(first, second), edge_vertices)
                if (edge_index == 0) return
                do moment = 1, edge_moment_count
                    local_dof = (edge - 1)*edge_moment_count + moment
                    global_dofs(local_dof, triangle) = &
                        (edge_index - 1)*edge_moment_count + moment
                    if (first > second .and. mod(moment, 2) == 1) &
                        transforms(local_dof, triangle) = -1
                end do
            end do
            do cell_dof = 1, cell_dof_count
                local_dof = 3*edge_moment_count + cell_dof
                global_dofs(local_dof, triangle) = &
                    edge_count*edge_moment_count + &
                    (triangle - 1)*cell_dof_count + cell_dof
            end do
        end do
        status = 0
    end subroutine build_maxwell_surface_rt_dof_map

    subroutine register_surface_edge( &
            first, second, edges, occurrences, edge_count, edge_index)
        integer, intent(in) :: first, second
        integer, intent(inout) :: edges(:, :), occurrences(:), edge_count
        integer, intent(out) :: edge_index

        edge_index = find_surface_edge( &
            first, second, edges(:, :edge_count))
        if (edge_index == 0) then
            edge_count = edge_count + 1
            edge_index = edge_count
            edges(:, edge_index) = [first, second]
        end if
        occurrences(edge_index) = occurrences(edge_index) + 1
    end subroutine register_surface_edge

    pure integer function find_surface_edge( &
            first, second, edges) result(edge_index)
        integer, intent(in) :: first, second, edges(:, :)

        integer :: edge

        edge_index = 0
        do edge = 1, size(edges, 2)
            if (edges(1, edge) /= first .or. edges(2, edge) /= second) cycle
            edge_index = edge
            return
        end do
    end function find_surface_edge

    pure function vector_cross_product(first, second) result(product)
        real(dp), intent(in) :: first(3), second(3)
        real(dp) :: product(3)

        product = [ &
            first(2)*second(3) - first(3)*second(2), &
            first(3)*second(1) - first(1)*second(3), &
            first(1)*second(2) - first(2)*second(1)]
    end function vector_cross_product

    subroutine evaluate_maxwell_surface_rt_basis( &
            basis, xi, eta, tangent_xi, tangent_eta, surface_jacobian, &
            values, surface_divergences, status)
        type(triangle_rt_basis_t), intent(in) :: basis
        real(dp), intent(in) :: xi, eta, tangent_xi(3), tangent_eta(3)
        real(dp), intent(in) :: surface_jacobian
        real(dp), intent(out) :: values(:, :), surface_divergences(:)
        integer, intent(out) :: status

        real(dp), allocatable :: reference_divergences(:)
        real(dp), allocatable :: reference_values(:, :)
        integer :: dof, dof_count

        status = 1
        values = 0.0_dp
        surface_divergences = 0.0_dp
        dof_count = triangle_rt_dof_count(basis)
        if (dof_count < 1 .or. surface_jacobian <= 0.0_dp) return
        if (size(values, 1) /= 3 .or. size(values, 2) /= dof_count) return
        if (size(surface_divergences) /= dof_count) return
        allocate( &
            reference_values(2, dof_count), &
            reference_divergences(dof_count))
        call evaluate_triangle_raviart_thomas( &
            basis, xi, eta, reference_values, reference_divergences, status)
        if (status /= 0) return
        do dof = 1, dof_count
            values(:, dof) = ( &
                tangent_xi*reference_values(1, dof) + &
                tangent_eta*reference_values(2, dof))/surface_jacobian
        end do
        surface_divergences = reference_divergences/surface_jacobian
        status = 0
    end subroutine evaluate_maxwell_surface_rt_basis

    subroutine evaluate_maxwell_surface_rt_global_basis( &
            vertices, triangles, degree, panel, global_basis, point, value, &
            surface_divergence, status)
        real(dp), intent(in) :: vertices(:, :), point(3)
        integer, intent(in) :: triangles(:, :), degree, panel, global_basis
        real(dp), intent(out) :: value(3), surface_divergence
        integer, intent(out) :: status

        type(triangle_rt_basis_t) :: basis
        integer, allocatable :: edge_vertices(:, :), global_dofs(:, :)
        integer, allocatable :: transforms(:, :)
        real(dp), allocatable :: divergences(:), values(:, :)
        real(dp) :: determinant, displacement(3), eta, gram11, gram12
        real(dp) :: gram22, jacobian, scale, tangent_eta(3), tangent_xi(3), xi
        integer :: local_dof, local_dof_count

        status = 1
        value = 0.0_dp
        surface_divergence = 0.0_dp
        if (size(vertices, 1) /= 3 .or. size(triangles, 1) /= 3) return
        if (degree < 0 .or. panel < 1 .or. panel > size(triangles, 2)) return
        if (any(triangles < 1) .or. any(triangles > size(vertices, 2))) return
        call initialize_triangle_raviart_thomas(degree, basis, status)
        if (status /= 0) return
        call build_maxwell_surface_rt_dof_map( &
            triangles, degree, edge_vertices, global_dofs, transforms, status)
        if (status /= 0) return
        if (global_basis < 1 .or. global_basis > maxval(global_dofs)) return
        local_dof = 0
        do local_dof_count = 1, size(global_dofs, 1)
            if (global_dofs(local_dof_count, panel) /= global_basis) cycle
            local_dof = local_dof_count
            exit
        end do
        if (local_dof == 0) then
            status = 2
            return
        end if
        tangent_xi = &
            vertices(:, triangles(2, panel)) - &
            vertices(:, triangles(1, panel))
        tangent_eta = &
            vertices(:, triangles(3, panel)) - &
            vertices(:, triangles(1, panel))
        displacement = point - vertices(:, triangles(1, panel))
        gram11 = dot_product(tangent_xi, tangent_xi)
        gram12 = dot_product(tangent_xi, tangent_eta)
        gram22 = dot_product(tangent_eta, tangent_eta)
        determinant = gram11*gram22 - gram12*gram12
        scale = max(1.0_dp, gram11*gram22)
        if (determinant <= 64.0_dp*epsilon(1.0_dp)*scale) return
        xi = ( &
            gram22*dot_product(tangent_xi, displacement) - &
            gram12*dot_product(tangent_eta, displacement))/determinant
        eta = ( &
            gram11*dot_product(tangent_eta, displacement) - &
            gram12*dot_product(tangent_xi, displacement))/determinant
        if (xi < -128.0_dp*epsilon(1.0_dp) .or. &
            eta < -128.0_dp*epsilon(1.0_dp) .or. &
            xi + eta > 1.0_dp + 256.0_dp*epsilon(1.0_dp)) return
        jacobian = sqrt(determinant)
        local_dof_count = triangle_rt_dof_count(basis)
        allocate( &
            values(3, local_dof_count), divergences(local_dof_count))
        call evaluate_maxwell_surface_rt_basis( &
            basis, xi, eta, tangent_xi, tangent_eta, jacobian, values, &
            divergences, status)
        if (status /= 0) return
        value = real(transforms(local_dof, panel), dp)*values(:, local_dof)
        surface_divergence = &
            real(transforms(local_dof, panel), dp)*divergences(local_dof)
        status = 0
    end subroutine evaluate_maxwell_surface_rt_global_basis

end module fortfem_maxwell_surface_rt
