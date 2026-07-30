module fortfem_maxwell_surface_rt_efie_3d
    !! Galerkin EFIE for arbitrary-order surface Raviart--Thomas currents.
    use fortfem_kinds, only: dp
    use fortfem_helmholtz_galerkin_3d, only: &
        assemble_helmholtz_single_layer_p0_adaptive_3d
    use fortfem_maxwell_surface_rt, only: &
        build_maxwell_surface_rt_dof_map, evaluate_maxwell_surface_rt_basis
    use fortfem_triangle_duffy_quadrature, only: triangle_duffy_quadrature
    use fortfem_triangle_rt_arbitrary_order, only: &
        initialize_triangle_raviart_thomas, triangle_rt_basis_t, &
        triangle_rt_dof_count
    use fortnum_quadrature, only: gauss_legendre_ab
    implicit none
    private

    public :: assemble_maxwell_surface_rt_efie_3d

contains

    subroutine assemble_maxwell_surface_rt_efie_3d( &
            vertices, triangles, degree, wave_number, impedance, &
            quadrature_degree, tolerance, max_depth, matrix, status)
        real(dp), intent(in) :: vertices(:, :), wave_number, impedance
        integer, intent(in) :: triangles(:, :), degree, quadrature_degree
        real(dp), intent(in) :: tolerance
        integer, intent(in) :: max_depth
        complex(dp), allocatable, intent(out) :: matrix(:, :)
        integer, intent(out) :: status

        type(triangle_rt_basis_t) :: basis
        integer, allocatable :: edges(:, :), global_dofs(:, :), transforms(:, :)
        real(dp), allocatable :: eta(:), line_nodes(:), line_weights(:)
        real(dp), allocatable :: weights(:), xi(:)
        complex(dp), allocatable :: scalar_block(:, :), vector_block(:, :)
        complex(dp), allocatable :: scalar_potential(:, :)
        complex(dp), allocatable :: vector_potential(:, :)
        complex(dp), allocatable :: panel_operator(:, :)
        real(dp), allocatable :: panel_divergences(:, :)
        real(dp) :: divergence(3), values(3, 3)
        real(dp) :: jacobian, tangent_eta(3), tangent_xi(3)
        integer :: first, i, j, node_count, second

        status = 1
        if (size(vertices, 1) /= 3 .or. size(triangles, 1) /= 3) return
        if (degree < 0 .or. wave_number <= 0.0_dp .or. impedance <= 0.0_dp) return
        if (quadrature_degree < 0 .or. tolerance <= 0.0_dp) return
        if (max_depth < 0 .or. any(triangles < 1) .or. &
            any(triangles > size(vertices, 2))) return
        call initialize_triangle_raviart_thomas(degree, basis, status)
        if (status /= 0) return
        call build_maxwell_surface_rt_dof_map( &
            triangles, degree, edges, global_dofs, transforms, status)
        if (status /= 0) return
        call triangle_duffy_quadrature( &
            quadrature_degree, xi, eta, weights, status)
        if (status /= 0) return
        ! The radial singularity cancels, but the remaining angular factor
        ! 1/|d(t)| is not polynomial and needs more nodes than smooth terms.
        node_count = max(2, (quadrature_degree + 3)/2)
        allocate(line_nodes(node_count), line_weights(node_count))
        call gauss_legendre_ab( &
            node_count, 0.0_dp, 1.0_dp, line_nodes, line_weights)
        allocate( &
            vector_potential(maxval(global_dofs), maxval(global_dofs)), &
            scalar_potential(maxval(global_dofs), maxval(global_dofs)))
        vector_potential = cmplx(0.0_dp, 0.0_dp, dp)
        scalar_potential = cmplx(0.0_dp, 0.0_dp, dp)
        do first = 1, size(triangles, 2)
            do second = 1, size(triangles, 2)
                call integrate_panel_pair( &
                    basis, vertices(:, triangles(:, first)), &
                    vertices(:, triangles(:, second)), first == second, &
                    panels_touch(triangles(:, first), triangles(:, second)), &
                    wave_number, xi, eta, weights, line_nodes, line_weights, &
                    tolerance, max_depth, vector_block, scalar_block, status)
                if (status /= 0) return
                do j = 1, size(global_dofs, 1)
                    do i = 1, size(global_dofs, 1)
                        vector_potential( &
                            global_dofs(i, first), global_dofs(j, second)) = &
                            vector_potential( &
                            global_dofs(i, first), global_dofs(j, second)) + &
                            real(transforms(i, first)*transforms(j, second), dp)* &
                            vector_block(i, j)
                        scalar_potential( &
                            global_dofs(i, first), global_dofs(j, second)) = &
                            scalar_potential( &
                            global_dofs(i, first), global_dofs(j, second)) + &
                            real(transforms(i, first)*transforms(j, second), dp)* &
                            scalar_block(i, j)
                    end do
                end do
            end do
        end do
        if (degree == 0) then
            call assemble_helmholtz_single_layer_p0_adaptive_3d( &
                vertices, triangles, wave_number, quadrature_degree, &
                tolerance, max_depth, panel_operator, status)
            if (status /= 0) return
            scalar_potential = cmplx(0.0_dp, 0.0_dp, dp)
            allocate(panel_divergences(3, size(triangles, 2)))
            do first = 1, size(triangles, 2)
                tangent_xi = vertices(:, triangles(2, first)) - &
                    vertices(:, triangles(1, first))
                tangent_eta = vertices(:, triangles(3, first)) - &
                    vertices(:, triangles(1, first))
                jacobian = norm2(cross_product(tangent_xi, tangent_eta))
                call evaluate_maxwell_surface_rt_basis( &
                    basis, 1.0_dp/3.0_dp, 1.0_dp/3.0_dp, tangent_xi, &
                    tangent_eta, jacobian, values, divergence, status)
                if (status /= 0) return
                panel_divergences(:, first) = &
                    real(transforms(:, first), dp)*divergence
            end do
            do first = 1, size(triangles, 2)
                do second = 1, size(triangles, 2)
                    do j = 1, 3
                        do i = 1, 3
                            scalar_potential( &
                                global_dofs(i, first), &
                                global_dofs(j, second)) = &
                                scalar_potential( &
                                global_dofs(i, first), &
                                global_dofs(j, second)) + &
                                panel_divergences(i, first)* &
                                panel_operator(first, second)* &
                                panel_divergences(j, second)
                        end do
                    end do
                end do
            end do
        end if
        vector_potential = 0.5_dp*( &
            vector_potential + transpose(vector_potential))
        scalar_potential = 0.5_dp*( &
            scalar_potential + transpose(scalar_potential))
        matrix = cmplx(0.0_dp, wave_number*impedance, dp)*vector_potential - &
            cmplx(0.0_dp, impedance/wave_number, dp)*scalar_potential
        status = 0
    end subroutine assemble_maxwell_surface_rt_efie_3d

    subroutine integrate_panel_pair( &
            basis, first_parent, second_parent, coincident, touching, &
            wave_number, xi, eta, weights, line_nodes, line_weights, &
            tolerance, max_depth, vector_block, scalar_block, status)
        type(triangle_rt_basis_t), intent(in) :: basis
        real(dp), intent(in) :: first_parent(3, 3), second_parent(3, 3)
        logical, intent(in) :: coincident, touching
        real(dp), intent(in) :: wave_number, xi(:), eta(:), weights(:)
        real(dp), intent(in) :: line_nodes(:), line_weights(:), tolerance
        integer, intent(in) :: max_depth
        complex(dp), allocatable, intent(out) :: vector_block(:, :)
        complex(dp), allocatable, intent(out) :: scalar_block(:, :)
        integer, intent(out) :: status

        integer :: count

        count = triangle_rt_dof_count(basis)
        allocate(vector_block(count, count), scalar_block(count, count))
        if (coincident) then
            call integrate_coincident_pair( &
                basis, first_parent, wave_number, xi, eta, weights, &
                line_nodes, line_weights, vector_block, scalar_block, status)
        else if (touching) then
            call integrate_adaptive_pair( &
                basis, first_parent, second_parent, first_parent, second_parent, &
                wave_number, xi, eta, weights, tolerance, 0, max_depth, &
                vector_block, scalar_block, status)
        else
            call integrate_regular_pair( &
                basis, first_parent, second_parent, first_parent, second_parent, &
                wave_number, xi, eta, weights, vector_block, scalar_block, &
                status)
        end if
    end subroutine integrate_panel_pair

    recursive subroutine integrate_adaptive_pair( &
            basis, first_parent, second_parent, first_panel, second_panel, &
            wave_number, xi, eta, weights, tolerance, depth, max_depth, &
            vector_block, scalar_block, status)
        type(triangle_rt_basis_t), intent(in) :: basis
        real(dp), intent(in) :: first_parent(3, 3), second_parent(3, 3)
        real(dp), intent(in) :: first_panel(3, 3), second_panel(3, 3)
        real(dp), intent(in) :: wave_number, xi(:), eta(:), weights(:)
        real(dp), intent(in) :: tolerance
        integer, intent(in) :: depth, max_depth
        complex(dp), intent(out) :: vector_block(:, :), scalar_block(:, :)
        integer, intent(out) :: status

        real(dp) :: first_children(3, 3, 4), second_children(3, 3, 4)
        complex(dp) :: coarse_v(size(vector_block, 1), size(vector_block, 2))
        complex(dp) :: coarse_s(size(scalar_block, 1), size(scalar_block, 2))
        complex(dp) :: part_v(size(vector_block, 1), size(vector_block, 2))
        complex(dp) :: part_s(size(scalar_block, 1), size(scalar_block, 2))
        complex(dp) :: refined_v(size(vector_block, 1), size(vector_block, 2))
        complex(dp) :: refined_s(size(scalar_block, 1), size(scalar_block, 2))
        integer :: first_child, second_child

        call integrate_regular_pair( &
            basis, first_parent, second_parent, first_panel, second_panel, &
            wave_number, xi, eta, weights, coarse_v, coarse_s, status)
        if (status /= 0) return
        call subdivide_triangle(first_panel, first_children)
        call subdivide_triangle(second_panel, second_children)
        refined_v = cmplx(0.0_dp, 0.0_dp, dp)
        refined_s = cmplx(0.0_dp, 0.0_dp, dp)
        do first_child = 1, 4
            do second_child = 1, 4
                call integrate_regular_pair( &
                    basis, first_parent, second_parent, &
                    first_children(:, :, first_child), &
                    second_children(:, :, second_child), wave_number, xi, eta, &
                    weights, part_v, part_s, status)
                if (status /= 0) return
                refined_v = refined_v + part_v
                refined_s = refined_s + part_s
            end do
        end do
        if (depth + 1 >= max_depth .or. converged_blocks( &
            coarse_v, refined_v, coarse_s, refined_s, tolerance)) then
            vector_block = refined_v
            scalar_block = refined_s
            status = 0
            return
        end if
        vector_block = cmplx(0.0_dp, 0.0_dp, dp)
        scalar_block = cmplx(0.0_dp, 0.0_dp, dp)
        do first_child = 1, 4
            do second_child = 1, 4
                if (geometric_panels_touch( &
                    first_children(:, :, first_child), &
                    second_children(:, :, second_child))) then
                    call integrate_adaptive_pair( &
                        basis, first_parent, second_parent, &
                        first_children(:, :, first_child), &
                        second_children(:, :, second_child), wave_number, xi, &
                        eta, weights, tolerance, depth + 1, max_depth, &
                        part_v, part_s, status)
                else
                    call integrate_regular_pair( &
                        basis, first_parent, second_parent, &
                        first_children(:, :, first_child), &
                        second_children(:, :, second_child), wave_number, xi, &
                        eta, weights, part_v, part_s, status)
                end if
                if (status /= 0) return
                vector_block = vector_block + part_v
                scalar_block = scalar_block + part_s
            end do
        end do
        status = 0
    end subroutine integrate_adaptive_pair

    subroutine integrate_regular_pair( &
            basis, first_parent, second_parent, first_panel, second_panel, &
            wave_number, xi, eta, weights, vector_block, scalar_block, status)
        type(triangle_rt_basis_t), intent(in) :: basis
        real(dp), intent(in) :: first_parent(3, 3), second_parent(3, 3)
        real(dp), intent(in) :: first_panel(3, 3), second_panel(3, 3)
        real(dp), intent(in) :: wave_number, xi(:), eta(:), weights(:)
        complex(dp), intent(out) :: vector_block(:, :), scalar_block(:, :)
        integer, intent(out) :: status

        real(dp) :: first_div(size(vector_block, 1))
        real(dp) :: first_values(3, size(vector_block, 1))
        real(dp) :: second_div(size(vector_block, 2))
        real(dp) :: second_values(3, size(vector_block, 2))
        real(dp) :: first_jacobian, first_point(3), radius
        real(dp) :: second_jacobian, second_point(3)
        complex(dp) :: factor
        integer :: first_node, i, j, second_node

        vector_block = cmplx(0.0_dp, 0.0_dp, dp)
        scalar_block = cmplx(0.0_dp, 0.0_dp, dp)
        first_jacobian = surface_jacobian(first_panel)
        second_jacobian = surface_jacobian(second_panel)
        do first_node = 1, size(weights)
            first_point = triangle_point( &
                first_panel, xi(first_node), eta(first_node))
            call evaluate_on_parent( &
                basis, first_parent, first_point, first_values, first_div, status)
            if (status /= 0) return
            do second_node = 1, size(weights)
                second_point = triangle_point( &
                    second_panel, xi(second_node), eta(second_node))
                call evaluate_on_parent( &
                    basis, second_parent, second_point, second_values, &
                    second_div, status)
                if (status /= 0) return
                radius = norm2(first_point - second_point)
                if (radius <= tiny(1.0_dp)) return
                factor = first_jacobian*second_jacobian*weights(first_node)* &
                    weights(second_node)*helmholtz_green(wave_number, radius)
                do j = 1, size(vector_block, 2)
                    do i = 1, size(vector_block, 1)
                        vector_block(i, j) = vector_block(i, j) + factor* &
                            dot_product(first_values(:, i), second_values(:, j))
                        scalar_block(i, j) = scalar_block(i, j) + factor* &
                            first_div(i)*second_div(j)
                    end do
                end do
            end do
        end do
        status = 0
    end subroutine integrate_regular_pair

    subroutine integrate_coincident_pair( &
            basis, panel, wave_number, xi, eta, weights, line_nodes, &
            line_weights, vector_block, scalar_block, status)
        type(triangle_rt_basis_t), intent(in) :: basis
        real(dp), intent(in) :: panel(3, 3), wave_number
        real(dp), intent(in) :: xi(:), eta(:), weights(:)
        real(dp), intent(in) :: line_nodes(:), line_weights(:)
        complex(dp), intent(out) :: vector_block(:, :), scalar_block(:, :)
        integer, intent(out) :: status

        real(dp) :: first_div(size(vector_block, 1))
        real(dp) :: first_values(3, size(vector_block, 1))
        real(dp) :: second_div(size(vector_block, 2))
        real(dp) :: second_values(3, size(vector_block, 2))
        real(dp) :: direction(3), first_point(3), jacobian, radius, rho
        real(dp) :: second_point(3), t, wedge_first(3), wedge_jacobian
        real(dp) :: wedge_second(3)
        complex(dp) :: factor
        integer :: first_node, i, j, radial_node, tangential_node, wedge

        vector_block = cmplx(0.0_dp, 0.0_dp, dp)
        scalar_block = cmplx(0.0_dp, 0.0_dp, dp)
        jacobian = surface_jacobian(panel)
        do first_node = 1, size(weights)
            first_point = triangle_point( &
                panel, xi(first_node), eta(first_node))
            call evaluate_on_parent( &
                basis, panel, first_point, first_values, first_div, status)
            if (status /= 0) return
            do wedge = 1, 3
                wedge_first = panel(:, wedge) - first_point
                wedge_second = panel(:, modulo(wedge, 3) + 1) - first_point
                wedge_jacobian = norm2(cross_product(wedge_first, wedge_second))
                do radial_node = 1, size(line_nodes)
                    rho = line_nodes(radial_node)
                    do tangential_node = 1, size(line_nodes)
                        t = line_nodes(tangential_node)
                        direction = &
                            (1.0_dp - t)*wedge_first + t*wedge_second
                        radius = rho*norm2(direction)
                        second_point = first_point + rho*direction
                        call evaluate_on_parent( &
                            basis, panel, second_point, second_values, &
                            second_div, status)
                        if (status /= 0) return
                        factor = jacobian*weights(first_node)* &
                            line_weights(radial_node)* &
                            line_weights(tangential_node)*rho*wedge_jacobian* &
                            helmholtz_green(wave_number, radius)
                        do j = 1, size(vector_block, 2)
                            do i = 1, size(vector_block, 1)
                                vector_block(i, j) = vector_block(i, j) + &
                                    factor*dot_product( &
                                    first_values(:, i), second_values(:, j))
                                scalar_block(i, j) = scalar_block(i, j) + &
                                    factor*first_div(i)*second_div(j)
                            end do
                        end do
                    end do
                end do
            end do
        end do
        status = 0
    end subroutine integrate_coincident_pair

    subroutine evaluate_on_parent( &
            basis, parent, point, values, divergences, status)
        type(triangle_rt_basis_t), intent(in) :: basis
        real(dp), intent(in) :: parent(3, 3), point(3)
        real(dp), intent(out) :: values(:, :), divergences(:)
        integer, intent(out) :: status

        real(dp) :: determinant, displacement(3), eta, gram11, gram12, gram22
        real(dp) :: jacobian, tangent_eta(3), tangent_xi(3), xi

        tangent_xi = parent(:, 2) - parent(:, 1)
        tangent_eta = parent(:, 3) - parent(:, 1)
        displacement = point - parent(:, 1)
        gram11 = dot_product(tangent_xi, tangent_xi)
        gram12 = dot_product(tangent_xi, tangent_eta)
        gram22 = dot_product(tangent_eta, tangent_eta)
        determinant = gram11*gram22 - gram12*gram12
        status = 1
        if (determinant <= tiny(1.0_dp)) return
        xi = (gram22*dot_product(tangent_xi, displacement) - &
            gram12*dot_product(tangent_eta, displacement))/determinant
        eta = (gram11*dot_product(tangent_eta, displacement) - &
            gram12*dot_product(tangent_xi, displacement))/determinant
        jacobian = sqrt(determinant)
        call evaluate_maxwell_surface_rt_basis( &
            basis, xi, eta, tangent_xi, tangent_eta, jacobian, values, &
            divergences, status)
    end subroutine evaluate_on_parent

    pure logical function converged_blocks( &
            coarse_v, refined_v, coarse_s, refined_s, tolerance) result(done)
        complex(dp), intent(in) :: coarse_v(:, :), refined_v(:, :)
        complex(dp), intent(in) :: coarse_s(:, :), refined_s(:, :)
        real(dp), intent(in) :: tolerance

        real(dp) :: error, scale

        error = max(maxval(abs(refined_v - coarse_v)), &
            maxval(abs(refined_s - coarse_s)))
        scale = max(tiny(1.0_dp), maxval(abs(refined_v)), &
            maxval(abs(refined_s)))
        done = error <= tolerance*scale
    end function converged_blocks

    pure function helmholtz_green(wave_number, radius) result(value)
        real(dp), intent(in) :: wave_number, radius
        complex(dp) :: value

        value = exp(cmplx(0.0_dp, wave_number*radius, dp))/ &
            (4.0_dp*acos(-1.0_dp)*radius)
    end function helmholtz_green

    pure function triangle_point(vertices, xi, eta) result(point)
        real(dp), intent(in) :: vertices(3, 3), xi, eta
        real(dp) :: point(3)

        point = vertices(:, 1) + xi*(vertices(:, 2) - vertices(:, 1)) + &
            eta*(vertices(:, 3) - vertices(:, 1))
    end function triangle_point

    pure function surface_jacobian(vertices) result(jacobian)
        real(dp), intent(in) :: vertices(3, 3)
        real(dp) :: jacobian

        jacobian = norm2(cross_product( &
            vertices(:, 2) - vertices(:, 1), &
            vertices(:, 3) - vertices(:, 1)))
    end function surface_jacobian

    pure logical function panels_touch(first, second) result(touch)
        integer, intent(in) :: first(3), second(3)
        integer :: vertex

        touch = .false.
        do vertex = 1, 3
            if (any(second == first(vertex))) touch = .true.
        end do
    end function panels_touch

    pure logical function geometric_panels_touch(first, second) result(touch)
        real(dp), intent(in) :: first(3, 3), second(3, 3)
        real(dp) :: scale
        integer :: i, j

        scale = max(1.0_dp, maxval(abs(first)), maxval(abs(second)))
        touch = .false.
        do i = 1, 3
            do j = 1, 3
                if (norm2(first(:, i) - second(:, j)) <= &
                    128.0_dp*epsilon(1.0_dp)*scale) touch = .true.
            end do
        end do
    end function geometric_panels_touch

    pure subroutine subdivide_triangle(vertices, children)
        real(dp), intent(in) :: vertices(3, 3)
        real(dp), intent(out) :: children(3, 3, 4)
        real(dp) :: midpoint_12(3), midpoint_23(3), midpoint_31(3)

        midpoint_12 = 0.5_dp*(vertices(:, 1) + vertices(:, 2))
        midpoint_23 = 0.5_dp*(vertices(:, 2) + vertices(:, 3))
        midpoint_31 = 0.5_dp*(vertices(:, 3) + vertices(:, 1))
        children(:, :, 1) = reshape( &
            [vertices(:, 1), midpoint_12, midpoint_31], [3, 3])
        children(:, :, 2) = reshape( &
            [midpoint_12, vertices(:, 2), midpoint_23], [3, 3])
        children(:, :, 3) = reshape( &
            [midpoint_31, midpoint_23, vertices(:, 3)], [3, 3])
        children(:, :, 4) = reshape( &
            [midpoint_12, midpoint_23, midpoint_31], [3, 3])
    end subroutine subdivide_triangle

    pure function cross_product(first, second) result(product)
        real(dp), intent(in) :: first(3), second(3)
        real(dp) :: product(3)

        product = [ &
            first(2)*second(3) - first(3)*second(2), &
            first(3)*second(1) - first(1)*second(3), &
            first(1)*second(2) - first(2)*second(1)]
    end function cross_product

end module fortfem_maxwell_surface_rt_efie_3d
