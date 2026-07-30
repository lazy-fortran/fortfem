module fortfem_maxwell_efie_rwg_3d
    !! Electric-field integral equation in a lowest-order RWG space.
    !!
    !! The scalar-potential block is a Galerkin product D^T V D because RWG
    !! surface divergences are panelwise constant. The vector-potential block
    !! is a true affine-RWG product integral; a radial Duffy map cancels the
    !! coincident-panel singularity analytically. Keeping the two blocks
    !! separate exposes low-frequency scaling and compatible preconditioning.
    use fortfem_helmholtz_galerkin_3d, only: &
        assemble_helmholtz_single_layer_p0_adaptive_3d
    use fortfem_kinds, only: dp
    use fortfem_maxwell_rwg_surface, only: &
        build_maxwell_rwg_surface_space, evaluate_maxwell_rwg_basis
    use fortfem_triangle_duffy_quadrature, only: triangle_duffy_quadrature
    use fortnum_linalg, only: dense_solve
    use fortnum_quadrature, only: gauss_legendre_ab
    implicit none
    private

    public :: assemble_maxwell_efie_rwg_3d
    public :: assemble_maxwell_plane_wave_rhs_rwg_3d
    public :: assemble_maxwell_rwg_potential_operators_3d
    public :: evaluate_maxwell_efie_field_rwg_3d
    public :: evaluate_maxwell_efie_far_field_rwg_3d
    public :: solve_maxwell_pec_efie_rwg_3d

contains

    subroutine evaluate_maxwell_efie_far_field_rwg_3d( &
            vertices, triangles, coefficients, direction, wave_number, &
            impedance, quadrature_degree, far_field, status)
        real(dp), intent(in) :: vertices(:, :), direction(3)
        complex(dp), intent(in) :: coefficients(:)
        real(dp), intent(in) :: wave_number, impedance
        integer, intent(in) :: triangles(:, :), quadrature_degree
        complex(dp), intent(out) :: far_field(3)
        integer, intent(out) :: status

        integer, allocatable :: edge_triangles(:, :), edge_vertices(:, :)
        real(dp), allocatable :: eta(:), weights(:), xi(:)
        real(dp) :: basis_divergence, basis_value(3), jacobian, point(3)
        complex(dp) :: phase, surface_current(3), transverse_current(3)
        integer :: basis, node, panel

        far_field = cmplx(0.0_dp, 0.0_dp, dp)
        status = 1
        if (wave_number <= 0.0_dp .or. impedance <= 0.0_dp) return
        if (abs(norm2(direction) - 1.0_dp) > 128.0_dp*epsilon(1.0_dp)) return
        call build_maxwell_rwg_surface_space( &
            vertices, triangles, edge_vertices, edge_triangles, status)
        if (status /= 0 .or. size(coefficients) /= size(edge_vertices, 2)) return
        call triangle_duffy_quadrature( &
            quadrature_degree, xi, eta, weights, status)
        if (status /= 0) return
        do panel = 1, size(triangles, 2)
            jacobian = 2.0_dp*triangle_area( &
                vertices(:, triangles(:, panel)))
            do node = 1, size(weights)
                point = triangle_point( &
                    vertices(:, triangles(:, panel)), xi(node), eta(node))
                surface_current = 0.0_dp
                do basis = 1, size(edge_vertices, 2)
                    if (.not. any(edge_triangles(:, basis) == panel)) cycle
                    call evaluate_maxwell_rwg_basis( &
                        vertices, triangles, edge_vertices, edge_triangles, &
                        basis, panel, point, basis_value, basis_divergence, status)
                    if (status /= 0) return
                    surface_current = surface_current + &
                        coefficients(basis)*basis_value
                end do
                transverse_current = surface_current - &
                    direction*sum(direction*surface_current)
                phase = exp(cmplx( &
                    0.0_dp, -wave_number*dot_product(direction, point), dp))
                far_field = far_field + jacobian*weights(node)*phase* &
                    transverse_current
            end do
        end do
        far_field = cmplx(0.0_dp, wave_number*impedance/ &
            (4.0_dp*acos(-1.0_dp)), dp)*far_field
        status = 0
    end subroutine evaluate_maxwell_efie_far_field_rwg_3d

    subroutine solve_maxwell_pec_efie_rwg_3d( &
            vertices, triangles, direction, polarization, wave_number, &
            impedance, quadrature_degree, tolerance, max_depth, density, status)
        real(dp), intent(in) :: vertices(:, :), direction(3), wave_number
        real(dp), intent(in) :: impedance, tolerance
        complex(dp), intent(in) :: polarization(3)
        integer, intent(in) :: triangles(:, :), quadrature_degree, max_depth
        complex(dp), allocatable, intent(out) :: density(:)
        integer, intent(out) :: status

        complex(dp), allocatable :: matrix(:, :), right_hand_side(:)
        integer :: info

        status = 1
        call assemble_maxwell_efie_rwg_3d( &
            vertices, triangles, wave_number, impedance, quadrature_degree, &
            tolerance, max_depth, matrix, status)
        if (status /= 0) return
        call assemble_maxwell_plane_wave_rhs_rwg_3d( &
            vertices, triangles, direction, polarization, wave_number, &
            quadrature_degree, right_hand_side, status)
        if (status /= 0) return
        allocate(density(size(right_hand_side)))
        call dense_solve(matrix, right_hand_side, density, info)
        if (info /= 0) then
            status = 2
            return
        end if
        status = 0
    end subroutine solve_maxwell_pec_efie_rwg_3d

    subroutine assemble_maxwell_plane_wave_rhs_rwg_3d( &
            vertices, triangles, direction, polarization, wave_number, &
            quadrature_degree, right_hand_side, status)
        real(dp), intent(in) :: vertices(:, :), direction(3), wave_number
        complex(dp), intent(in) :: polarization(3)
        integer, intent(in) :: triangles(:, :), quadrature_degree
        complex(dp), allocatable, intent(out) :: right_hand_side(:)
        integer, intent(out) :: status

        integer, allocatable :: edge_triangles(:, :), edge_vertices(:, :)
        real(dp), allocatable :: eta(:), weights(:), xi(:)
        real(dp) :: basis_divergence, basis_value(3), jacobian, point(3)
        complex(dp) :: incident_field(3)
        integer :: basis, node, panel

        status = 1
        if (wave_number < 0.0_dp) return
        if (abs(norm2(direction) - 1.0_dp) > 128.0_dp*epsilon(1.0_dp)) return
        if (abs(sum(polarization*direction)) > &
            128.0_dp*epsilon(1.0_dp)*max(1.0_dp, &
            sqrt(sum(abs(polarization)**2)))) return
        if (sqrt(sum(abs(polarization)**2)) <= tiny(1.0_dp)) return
        call build_maxwell_rwg_surface_space( &
            vertices, triangles, edge_vertices, edge_triangles, status)
        if (status /= 0 .or. size(edge_vertices, 2) == 0) return
        call triangle_duffy_quadrature( &
            quadrature_degree, xi, eta, weights, status)
        if (status /= 0) return
        allocate(right_hand_side(size(edge_vertices, 2)))
        right_hand_side = cmplx(0.0_dp, 0.0_dp, dp)
        do basis = 1, size(edge_vertices, 2)
            do panel = 1, 2
                jacobian = 2.0_dp*triangle_area(vertices(:, triangles(:, &
                    edge_triangles(panel, basis))))
                do node = 1, size(weights)
                    point = triangle_point( &
                        vertices(:, triangles(:, edge_triangles(panel, basis))), &
                        xi(node), eta(node))
                    call evaluate_maxwell_rwg_basis( &
                        vertices, triangles, edge_vertices, edge_triangles, &
                        basis, edge_triangles(panel, basis), point, basis_value, &
                        basis_divergence, status)
                    if (status /= 0) return
                    incident_field = polarization*exp(cmplx( &
                        0.0_dp, wave_number*dot_product(direction, point), dp))
                    right_hand_side(basis) = right_hand_side(basis) - &
                        jacobian*weights(node)* &
                        sum(cmplx(basis_value, 0.0_dp, dp)*incident_field)
                end do
            end do
        end do
        status = 0
    end subroutine assemble_maxwell_plane_wave_rhs_rwg_3d

    subroutine evaluate_maxwell_efie_field_rwg_3d( &
            vertices, triangles, coefficients, target, wave_number, impedance, &
            quadrature_degree, field, status)
        real(dp), intent(in) :: vertices(:, :), target(3)
        complex(dp), intent(in) :: coefficients(:)
        real(dp), intent(in) :: wave_number, impedance
        integer, intent(in) :: triangles(:, :), quadrature_degree
        complex(dp), intent(out) :: field(3)
        integer, intent(out) :: status

        integer, allocatable :: edge_triangles(:, :), edge_vertices(:, :)
        real(dp), allocatable :: eta(:), weights(:), xi(:)
        real(dp) :: basis_divergence, basis_value(3), displacement(3)
        real(dp) :: jacobian, point(3), radius
        complex(dp) :: gradient_green(3), green, surface_current(3)
        complex(dp) :: surface_divergence
        integer :: basis, node, panel

        field = cmplx(0.0_dp, 0.0_dp, dp)
        status = 1
        if (wave_number <= 0.0_dp .or. impedance <= 0.0_dp) return
        call build_maxwell_rwg_surface_space( &
            vertices, triangles, edge_vertices, edge_triangles, status)
        if (status /= 0 .or. size(coefficients) /= size(edge_vertices, 2)) return
        call triangle_duffy_quadrature( &
            quadrature_degree, xi, eta, weights, status)
        if (status /= 0) return
        do panel = 1, size(triangles, 2)
            jacobian = 2.0_dp*triangle_area( &
                vertices(:, triangles(:, panel)))
            do node = 1, size(weights)
                point = triangle_point( &
                    vertices(:, triangles(:, panel)), xi(node), eta(node))
                surface_current = 0.0_dp
                surface_divergence = 0.0_dp
                do basis = 1, size(edge_vertices, 2)
                    if (.not. any(edge_triangles(:, basis) == panel)) cycle
                    call evaluate_maxwell_rwg_basis( &
                        vertices, triangles, edge_vertices, edge_triangles, &
                        basis, panel, point, basis_value, basis_divergence, status)
                    if (status /= 0) return
                    surface_current = surface_current + &
                        coefficients(basis)*basis_value
                    surface_divergence = surface_divergence + &
                        coefficients(basis)*basis_divergence
                end do
                displacement = target - point
                radius = norm2(displacement)
                if (radius <= 0.0_dp) return
                green = helmholtz_green(wave_number, radius)
                gradient_green = green* &
                    cmplx(-1.0_dp, wave_number*radius, dp)* &
                    displacement/radius**2
                field = field + jacobian*weights(node)*( &
                    cmplx(0.0_dp, wave_number*impedance, dp)*green* &
                    surface_current + &
                    cmplx(0.0_dp, impedance/wave_number, dp)* &
                    gradient_green*surface_divergence)
            end do
        end do
        status = 0
    end subroutine evaluate_maxwell_efie_field_rwg_3d

    subroutine assemble_maxwell_efie_rwg_3d( &
            vertices, triangles, wave_number, impedance, quadrature_degree, &
            tolerance, max_depth, matrix, status)
        real(dp), intent(in) :: vertices(:, :), wave_number, impedance
        integer, intent(in) :: triangles(:, :), quadrature_degree, max_depth
        real(dp), intent(in) :: tolerance
        complex(dp), allocatable, intent(out) :: matrix(:, :)
        integer, intent(out) :: status

        complex(dp), allocatable :: scalar_potential(:, :)
        complex(dp), allocatable :: vector_potential(:, :)

        status = 1
        if (wave_number <= 0.0_dp .or. impedance <= 0.0_dp) return
        call assemble_maxwell_rwg_potential_operators_3d( &
            vertices, triangles, wave_number, quadrature_degree, tolerance, &
            max_depth, vector_potential, scalar_potential, status)
        if (status /= 0) return
        matrix = cmplx(0.0_dp, wave_number*impedance, dp)*vector_potential - &
            cmplx(0.0_dp, impedance/wave_number, dp)*scalar_potential
        status = 0
    end subroutine assemble_maxwell_efie_rwg_3d

    subroutine assemble_maxwell_rwg_potential_operators_3d( &
            vertices, triangles, wave_number, quadrature_degree, tolerance, &
            max_depth, vector_potential, scalar_potential, status)
        real(dp), intent(in) :: vertices(:, :), wave_number, tolerance
        integer, intent(in) :: triangles(:, :), quadrature_degree, max_depth
        complex(dp), allocatable, intent(out) :: vector_potential(:, :)
        complex(dp), allocatable, intent(out) :: scalar_potential(:, :)
        integer, intent(out) :: status

        complex(dp), allocatable :: panel_operator(:, :)
        integer, allocatable :: edge_triangles(:, :), edge_vertices(:, :)
        real(dp), allocatable :: eta(:), line_nodes(:), line_weights(:)
        real(dp), allocatable :: weights(:), xi(:)
        real(dp) :: first_divergence, first_value(3)
        real(dp) :: second_divergence, second_value(3)
        complex(dp) :: vector_entry
        integer :: first, first_panel, first_slot, second, second_panel
        integer :: node_count, second_slot

        status = 1
        if (wave_number < 0.0_dp) return
        call build_maxwell_rwg_surface_space( &
            vertices, triangles, edge_vertices, edge_triangles, status)
        if (status /= 0 .or. size(edge_vertices, 2) == 0) return
        call assemble_helmholtz_single_layer_p0_adaptive_3d( &
            vertices, triangles, wave_number, quadrature_degree, tolerance, &
            max_depth, panel_operator, status)
        if (status /= 0) return
        call triangle_duffy_quadrature( &
            quadrature_degree, xi, eta, weights, status)
        if (status /= 0) return
        node_count = max(2, (quadrature_degree + 3)/2)
        allocate(line_nodes(node_count), line_weights(node_count))
        call gauss_legendre_ab( &
            node_count, 0.0_dp, 1.0_dp, line_nodes, line_weights)
        allocate( &
            vector_potential(size(edge_vertices, 2), size(edge_vertices, 2)), &
            scalar_potential(size(edge_vertices, 2), size(edge_vertices, 2)))
        vector_potential = cmplx(0.0_dp, 0.0_dp, dp)
        scalar_potential = cmplx(0.0_dp, 0.0_dp, dp)
        do first = 1, size(edge_vertices, 2)
            do second = 1, first
                do first_slot = 1, 2
                    first_panel = edge_triangles(first_slot, first)
                    call evaluate_maxwell_rwg_basis( &
                        vertices, triangles, edge_vertices, edge_triangles, &
                        first, first_panel, panel_centroid( &
                        vertices(:, triangles(:, first_panel))), first_value, &
                        first_divergence, status)
                    if (status /= 0) return
                    do second_slot = 1, 2
                        second_panel = edge_triangles(second_slot, second)
                        call evaluate_maxwell_rwg_basis( &
                            vertices, triangles, edge_vertices, edge_triangles, &
                            second, second_panel, panel_centroid( &
                            vertices(:, triangles(:, second_panel))), &
                            second_value, &
                            second_divergence, status)
                        if (status /= 0) return
                        call integrate_rwg_vector_panel_pair( &
                            vertices, triangles, edge_vertices, edge_triangles, &
                            first, first_panel, second, second_panel, wave_number, &
                            xi, eta, weights, line_nodes, line_weights, &
                            tolerance, max_depth, vector_entry, status)
                        if (status /= 0) return
                        vector_potential(first, second) = &
                            vector_potential(first, second) + vector_entry
                        scalar_potential(first, second) = &
                            scalar_potential(first, second) + &
                            first_divergence*second_divergence* &
                            panel_operator(first_panel, second_panel)
                    end do
                end do
                vector_potential(second, first) = vector_potential(first, second)
                scalar_potential(second, first) = scalar_potential(first, second)
            end do
        end do
        status = 0
    end subroutine assemble_maxwell_rwg_potential_operators_3d

    subroutine integrate_rwg_vector_panel_pair( &
            vertices, triangles, edge_vertices, edge_triangles, first_basis, &
            first_panel, second_basis, second_panel, wave_number, xi, eta, &
            weights, line_nodes, line_weights, tolerance, max_depth, value, &
            status)
        real(dp), intent(in) :: vertices(:, :), wave_number
        integer, intent(in) :: triangles(:, :), edge_vertices(:, :)
        integer, intent(in) :: edge_triangles(:, :)
        integer, intent(in) :: first_basis, first_panel, second_basis
        integer, intent(in) :: second_panel
        real(dp), intent(in) :: xi(:), eta(:), weights(:)
        real(dp), intent(in) :: line_nodes(:), line_weights(:)
        real(dp), intent(in) :: tolerance
        integer, intent(in) :: max_depth
        complex(dp), intent(out) :: value
        integer, intent(out) :: status

        real(dp) :: first_vertices(3, 3), second_vertices(3, 3)

        first_vertices = vertices(:, triangles(:, first_panel))
        second_vertices = vertices(:, triangles(:, second_panel))
        if (first_panel == second_panel) then
            call integrate_coincident_rwg_pair( &
                vertices, triangles, edge_vertices, edge_triangles, &
                first_basis, first_panel, second_basis, first_vertices, &
                wave_number, xi, eta, weights, line_nodes, line_weights, &
                value, status)
        else if (panels_touch( &
                triangles(:, first_panel), triangles(:, second_panel))) then
            call integrate_adaptive_rwg_pair( &
                vertices, triangles, edge_vertices, edge_triangles, &
                first_basis, first_panel, second_basis, second_panel, &
                first_vertices, second_vertices, wave_number, xi, eta, &
                weights, tolerance, 0, max_depth, value, status)
        else
            call integrate_regular_rwg_pair( &
                vertices, triangles, edge_vertices, edge_triangles, &
                first_basis, first_panel, second_basis, second_panel, &
                first_vertices, second_vertices, wave_number, xi, eta, &
                weights, value, status)
        end if
    end subroutine integrate_rwg_vector_panel_pair

    recursive subroutine integrate_adaptive_rwg_pair( &
            vertices, triangles, edge_vertices, edge_triangles, first_basis, &
            first_panel, second_basis, second_panel, first_vertices, &
            second_vertices, wave_number, xi, eta, weights, tolerance, depth, &
            max_depth, value, status)
        real(dp), intent(in) :: vertices(:, :), first_vertices(3, 3)
        real(dp), intent(in) :: second_vertices(3, 3), wave_number, tolerance
        integer, intent(in) :: triangles(:, :), edge_vertices(:, :)
        integer, intent(in) :: edge_triangles(:, :)
        integer, intent(in) :: first_basis, first_panel, second_basis
        integer, intent(in) :: second_panel, depth, max_depth
        real(dp), intent(in) :: xi(:), eta(:), weights(:)
        complex(dp), intent(out) :: value
        integer, intent(out) :: status

        real(dp) :: first_children(3, 3, 4), second_children(3, 3, 4)
        complex(dp) :: coarse, refined, contribution
        integer :: first_child, second_child

        call integrate_regular_rwg_pair( &
            vertices, triangles, edge_vertices, edge_triangles, first_basis, &
            first_panel, second_basis, second_panel, first_vertices, &
            second_vertices, wave_number, xi, eta, weights, coarse, status)
        if (status /= 0) return
        call subdivide_triangle(first_vertices, first_children)
        call subdivide_triangle(second_vertices, second_children)
        refined = cmplx(0.0_dp, 0.0_dp, dp)
        do first_child = 1, 4
            do second_child = 1, 4
                call integrate_regular_rwg_pair( &
                    vertices, triangles, edge_vertices, edge_triangles, &
                    first_basis, first_panel, second_basis, second_panel, &
                    first_children(:, :, first_child), &
                    second_children(:, :, second_child), wave_number, xi, eta, &
                    weights, contribution, status)
                if (status /= 0) return
                refined = refined + contribution
            end do
        end do
        if (depth + 1 >= max_depth .or. abs(refined - coarse) <= &
            tolerance*max(tiny(1.0_dp), abs(refined))) then
            value = refined
            status = 0
            return
        end if
        value = cmplx(0.0_dp, 0.0_dp, dp)
        do first_child = 1, 4
            do second_child = 1, 4
                if (geometric_panels_touch( &
                    first_children(:, :, first_child), &
                    second_children(:, :, second_child))) then
                    call integrate_adaptive_rwg_pair( &
                        vertices, triangles, edge_vertices, edge_triangles, &
                        first_basis, first_panel, second_basis, second_panel, &
                        first_children(:, :, first_child), &
                        second_children(:, :, second_child), wave_number, xi, &
                        eta, weights, tolerance, depth + 1, max_depth, &
                        contribution, status)
                else
                    call integrate_regular_rwg_pair( &
                        vertices, triangles, edge_vertices, edge_triangles, &
                        first_basis, first_panel, second_basis, second_panel, &
                        first_children(:, :, first_child), &
                        second_children(:, :, second_child), wave_number, xi, &
                        eta, weights, contribution, status)
                end if
                if (status /= 0) return
                value = value + contribution
            end do
        end do
        status = 0
    end subroutine integrate_adaptive_rwg_pair

    subroutine integrate_regular_rwg_pair( &
            vertices, triangles, edge_vertices, edge_triangles, first_basis, &
            first_panel, second_basis, second_panel, first_vertices, &
            second_vertices, wave_number, xi, eta, weights, value, status)
        real(dp), intent(in) :: vertices(:, :), first_vertices(3, 3)
        real(dp), intent(in) :: second_vertices(3, 3), wave_number
        integer, intent(in) :: triangles(:, :), edge_vertices(:, :)
        integer, intent(in) :: edge_triangles(:, :)
        integer, intent(in) :: first_basis, first_panel, second_basis
        integer, intent(in) :: second_panel
        real(dp), intent(in) :: xi(:), eta(:), weights(:)
        complex(dp), intent(out) :: value
        integer, intent(out) :: status

        real(dp) :: first_divergence, first_point(3), first_value(3)
        real(dp) :: first_jacobian, radius, second_divergence
        real(dp) :: second_jacobian, second_point(3), second_value(3)
        integer :: first_node, second_node

        value = cmplx(0.0_dp, 0.0_dp, dp)
        first_jacobian = 2.0_dp*triangle_area(first_vertices)
        second_jacobian = 2.0_dp*triangle_area(second_vertices)
        do first_node = 1, size(weights)
            first_point = triangle_point( &
                first_vertices, xi(first_node), eta(first_node))
            call evaluate_maxwell_rwg_basis( &
                vertices, triangles, edge_vertices, edge_triangles, &
                first_basis, first_panel, first_point, first_value, &
                first_divergence, status)
            if (status /= 0) return
            do second_node = 1, size(weights)
                second_point = triangle_point( &
                    second_vertices, xi(second_node), eta(second_node))
                call evaluate_maxwell_rwg_basis( &
                    vertices, triangles, edge_vertices, edge_triangles, &
                    second_basis, second_panel, second_point, second_value, &
                    second_divergence, status)
                if (status /= 0) return
                radius = norm2(first_point - second_point)
                value = value + first_jacobian*second_jacobian* &
                    weights(first_node)*weights(second_node)* &
                    helmholtz_green(wave_number, radius)* &
                    dot_product(first_value, second_value)
            end do
        end do
        status = 0
    end subroutine integrate_regular_rwg_pair

    subroutine integrate_coincident_rwg_pair( &
            vertices, triangles, edge_vertices, edge_triangles, first_basis, &
            panel, second_basis, panel_vertices, wave_number, xi, eta, &
            weights, line_nodes, line_weights, value, status)
        real(dp), intent(in) :: vertices(:, :), panel_vertices(3, 3)
        real(dp), intent(in) :: wave_number
        integer, intent(in) :: triangles(:, :), edge_vertices(:, :)
        integer, intent(in) :: edge_triangles(:, :)
        integer, intent(in) :: first_basis, panel, second_basis
        real(dp), intent(in) :: xi(:), eta(:), weights(:)
        real(dp), intent(in) :: line_nodes(:), line_weights(:)
        complex(dp), intent(out) :: value
        integer, intent(out) :: status

        real(dp) :: first_divergence, first_point(3), first_value(3)
        real(dp) :: direction(3), jacobian, radius, rho, second_divergence
        real(dp) :: second_point(3), second_value(3), t
        real(dp) :: wedge_first(3), wedge_jacobian, wedge_second(3)
        integer :: first_node, radial_node, tangential_node, wedge

        value = cmplx(0.0_dp, 0.0_dp, dp)
        jacobian = 2.0_dp*triangle_area(panel_vertices)
        do first_node = 1, size(weights)
            first_point = triangle_point( &
                panel_vertices, xi(first_node), eta(first_node))
            call evaluate_maxwell_rwg_basis( &
                vertices, triangles, edge_vertices, edge_triangles, &
                first_basis, panel, first_point, first_value, &
                first_divergence, status)
            if (status /= 0) return
            do wedge = 1, 3
                wedge_first = panel_vertices(:, wedge) - first_point
                wedge_second = panel_vertices(:, modulo(wedge, 3) + 1) - &
                    first_point
                wedge_jacobian = norm2( &
                    cross_product(wedge_first, wedge_second))
                do radial_node = 1, size(line_nodes)
                    rho = line_nodes(radial_node)
                    do tangential_node = 1, size(line_nodes)
                        t = line_nodes(tangential_node)
                        direction = (1.0_dp - t)*wedge_first + t*wedge_second
                        radius = rho*norm2(direction)
                        second_point = first_point + rho*direction
                        call evaluate_maxwell_rwg_basis( &
                            vertices, triangles, edge_vertices, edge_triangles, &
                            second_basis, panel, second_point, second_value, &
                            second_divergence, status)
                        if (status /= 0) return
                        value = value + jacobian*weights(first_node)* &
                            line_weights(radial_node)* &
                            line_weights(tangential_node)*rho*wedge_jacobian* &
                            helmholtz_green(wave_number, radius)* &
                            dot_product(first_value, second_value)
                    end do
                end do
            end do
        end do
        status = 0
    end subroutine integrate_coincident_rwg_pair

    pure function helmholtz_green(wave_number, radius) result(value)
        real(dp), intent(in) :: wave_number, radius
        complex(dp) :: value

        value = exp(cmplx(0.0_dp, wave_number*radius, dp))/ &
            (4.0_dp*acos(-1.0_dp)*radius)
    end function helmholtz_green

    pure function panel_centroid(vertices) result(point)
        real(dp), intent(in) :: vertices(3, 3)
        real(dp) :: point(3)

        point = sum(vertices, dim=2)/3.0_dp
    end function panel_centroid

    pure function triangle_point(vertices, xi, eta) result(point)
        real(dp), intent(in) :: vertices(3, 3), xi, eta
        real(dp) :: point(3)

        point = vertices(:, 1) + xi*(vertices(:, 2) - vertices(:, 1)) + &
            eta*(vertices(:, 3) - vertices(:, 1))
    end function triangle_point

    pure function triangle_area(vertices) result(area)
        real(dp), intent(in) :: vertices(3, 3)
        real(dp) :: area

        area = 0.5_dp*norm2(cross_product( &
            vertices(:, 2) - vertices(:, 1), &
            vertices(:, 3) - vertices(:, 1)))
    end function triangle_area

    pure logical function panels_touch(first, second) result(touch)
        integer, intent(in) :: first(3), second(3)

        integer :: first_vertex

        touch = .false.
        do first_vertex = 1, 3
            if (any(second == first(first_vertex))) then
                touch = .true.
                return
            end if
        end do
    end function panels_touch

    pure logical function geometric_panels_touch(first, second) result(touch)
        real(dp), intent(in) :: first(3, 3), second(3, 3)

        real(dp) :: scale
        integer :: first_vertex, second_vertex

        scale = max(1.0_dp, maxval(abs(first)), maxval(abs(second)))
        touch = .false.
        do first_vertex = 1, 3
            do second_vertex = 1, 3
                if (norm2(first(:, first_vertex) - &
                    second(:, second_vertex)) <= &
                    128.0_dp*epsilon(1.0_dp)*scale) then
                    touch = .true.
                    return
                end if
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

end module fortfem_maxwell_efie_rwg_3d
