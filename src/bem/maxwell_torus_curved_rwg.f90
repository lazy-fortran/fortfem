module fortfem_maxwell_torus_curved_rwg
    !! Surface-Piola RWG basis on exact parametric torus panels.
    use fortfem_kinds, only: dp
    use fortfem_maxwell_bc_surface, only: build_maxwell_bc_transformation
    use fortfem_maxwell_rwg_surface, only: build_maxwell_rwg_surface_space
    use fortfem_torus_curved_panel, only: evaluate_torus_curved_panel
    use fortfem_triangle_duffy_quadrature, only: triangle_duffy_quadrature
    use fortnum_quadrature, only: gauss_legendre_ab
    implicit none
    private

    public :: assemble_maxwell_torus_curved_rwg_mass_matrix
    public :: assemble_maxwell_torus_curved_rwg_rbc_pairing
    public :: assemble_maxwell_torus_curved_plane_wave_rhs_rwg_3d
    public :: assemble_maxwell_torus_curved_efie_rwg_3d
    public :: assemble_maxwell_torus_curved_potential_operators_rwg_3d
    public :: evaluate_maxwell_torus_curved_far_field_rwg_3d
    public :: evaluate_maxwell_torus_curved_localized_rwg_basis
    public :: evaluate_maxwell_torus_curved_rwg_basis
    public :: integrate_maxwell_torus_curved_adjacent_rwg_pair_3d
    public :: integrate_maxwell_torus_curved_coincident_rwg_pair_3d

contains

    subroutine assemble_maxwell_torus_curved_rwg_rbc_pairing( &
            vertices, triangles, parameters, major_radius, minor_radius, &
            quadrature_degree, matrix, status)
        real(dp), intent(in) :: vertices(:, :), parameters(:, :)
        real(dp), intent(in) :: major_radius, minor_radius
        integer, intent(in) :: triangles(:, :), quadrature_degree
        real(dp), allocatable, intent(out) :: matrix(:, :)
        integer, intent(out) :: status

        integer, allocatable :: edge_triangles(:, :), edge_vertices(:, :)
        integer, allocatable :: refined_triangles(:, :)
        real(dp), allocatable :: eta(:), refined_parameters(:, :)
        real(dp), allocatable :: refined_vertices(:, :)
        real(dp), allocatable :: transformation(:, :), weights(:), xi(:)
        real(dp), allocatable :: bc_values(:, :), rwg_values(:, :)
        real(dp) :: coarse_eta, coarse_jacobian, coarse_xi, divergence
        real(dp) :: local_value(3)
        real(dp) :: normal(3), point(3), refined_jacobian, rotated_bc(3)
        integer :: basis, local_edge, node, parent, refined_panel, row
        integer :: test_basis

        status = 1
        if (allocated(matrix)) deallocate(matrix)
        call build_maxwell_bc_transformation( &
            vertices, triangles, refined_vertices, refined_triangles, &
            transformation, status, torus_parameters=parameters, &
            torus_major_radius=major_radius, torus_minor_radius=minor_radius, &
            refined_torus_parameters=refined_parameters)
        if (status /= 0) return
        call build_maxwell_rwg_surface_space( &
            vertices, triangles, edge_vertices, edge_triangles, status)
        if (status /= 0) return
        call triangle_duffy_quadrature( &
            quadrature_degree, xi, eta, weights, status)
        if (status /= 0) return
        allocate( &
            matrix(size(edge_vertices, 2), size(edge_vertices, 2)), &
            bc_values(3, size(edge_vertices, 2)), &
            rwg_values(3, size(edge_vertices, 2)))
        matrix = 0.0_dp
        do refined_panel = 1, size(refined_triangles, 2)
            parent = (refined_panel - 1)/6 + 1
            do node = 1, size(weights)
                bc_values = 0.0_dp
                do local_edge = 1, 3
                    call evaluate_maxwell_torus_curved_localized_rwg_basis( &
                        refined_vertices, refined_triangles, &
                        refined_parameters, refined_panel, local_edge, &
                        major_radius, minor_radius, xi(node), eta(node), point, &
                        local_value, divergence, refined_jacobian, status)
                    if (status /= 0) return
                    row = 3*(refined_panel - 1) + local_edge
                    do test_basis = 1, size(edge_vertices, 2)
                        bc_values(:, test_basis) = &
                            bc_values(:, test_basis) + &
                            transformation(row, test_basis)*local_value
                    end do
                end do
                call map_refined_torus_point_to_parent( &
                    refined_parameters(:, refined_triangles(:, refined_panel)), &
                    parameters(:, triangles(:, parent)), xi(node), eta(node), &
                    coarse_xi, coarse_eta, status)
                if (status /= 0) return
                rwg_values = 0.0_dp
                do basis = 1, size(edge_vertices, 2)
                    if (.not. any(edge_triangles(:, basis) == parent)) cycle
                    call evaluate_maxwell_torus_curved_rwg_basis( &
                        vertices, triangles, parameters, edge_vertices, &
                        edge_triangles, basis, parent, major_radius, &
                        minor_radius, coarse_xi, coarse_eta, point, &
                        rwg_values(:, basis), divergence, coarse_jacobian, status)
                    if (status /= 0) return
                end do
                normal = torus_unit_normal(point, major_radius)
                do test_basis = 1, size(edge_vertices, 2)
                    rotated_bc = real_cross_product( &
                        normal, bc_values(:, test_basis))
                    do basis = 1, size(edge_vertices, 2)
                        matrix(test_basis, basis) = &
                            matrix(test_basis, basis) + &
                            refined_jacobian*weights(node)*dot_product( &
                            rotated_bc, rwg_values(:, basis))
                    end do
                end do
            end do
        end do
        status = 0
    end subroutine assemble_maxwell_torus_curved_rwg_rbc_pairing

    pure subroutine evaluate_maxwell_torus_curved_localized_rwg_basis( &
            vertices, triangles, parameters, panel, local_edge, major_radius, &
            minor_radius, xi, eta, point, value, surface_divergence, &
            surface_jacobian, status)
        real(dp), intent(in) :: vertices(:, :), parameters(:, :)
        real(dp), intent(in) :: major_radius, minor_radius, xi, eta
        integer, intent(in) :: triangles(:, :), panel, local_edge
        real(dp), intent(out) :: point(3), value(3), surface_divergence
        real(dp), intent(out) :: surface_jacobian
        integer, intent(out) :: status

        real(dp) :: edge_length, opposite_coordinates(2)
        real(dp) :: panel_parameters(2, 3), tangent_eta(3), tangent_xi(3)
        integer :: edge_local_vertices(2), local, opposite

        point = 0.0_dp
        value = 0.0_dp
        surface_divergence = 0.0_dp
        surface_jacobian = 0.0_dp
        status = 1
        if (panel < 1 .or. panel > size(triangles, 2)) return
        if (local_edge < 1 .or. local_edge > 3) return
        select case (local_edge)
        case (1)
            edge_local_vertices = [1, 2]
            opposite = 3
        case (2)
            edge_local_vertices = [3, 1]
            opposite = 2
        case (3)
            edge_local_vertices = [2, 3]
            opposite = 1
        end select
        do local = 1, 3
            panel_parameters(:, local) = &
                parameters(:, triangles(local, panel))
        end do
        call evaluate_torus_curved_panel( &
            panel_parameters, major_radius, minor_radius, xi, eta, point, &
            tangent_xi, tangent_eta, surface_jacobian, status)
        if (status /= 0) return
        select case (opposite)
        case (1)
            opposite_coordinates = [0.0_dp, 0.0_dp]
        case (2)
            opposite_coordinates = [1.0_dp, 0.0_dp]
        case (3)
            opposite_coordinates = [0.0_dp, 1.0_dp]
        end select
        edge_length = norm2( &
            vertices(:, triangles(edge_local_vertices(2), panel)) - &
            vertices(:, triangles(edge_local_vertices(1), panel)))
        value = edge_length/surface_jacobian*( &
            (xi - opposite_coordinates(1))*tangent_xi + &
            (eta - opposite_coordinates(2))*tangent_eta)
        surface_divergence = 2.0_dp*edge_length/surface_jacobian
        status = 0
    end subroutine evaluate_maxwell_torus_curved_localized_rwg_basis

    subroutine assemble_maxwell_torus_curved_efie_rwg_3d( &
            vertices, triangles, parameters, major_radius, minor_radius, &
            wave_number, impedance, quadrature_degree, tolerance, max_depth, &
            matrix, status)
        real(dp), intent(in) :: vertices(:, :), parameters(:, :)
        real(dp), intent(in) :: major_radius, minor_radius, wave_number
        real(dp), intent(in) :: impedance, tolerance
        integer, intent(in) :: triangles(:, :), quadrature_degree, max_depth
        complex(dp), allocatable, intent(out) :: matrix(:, :)
        integer, intent(out) :: status

        complex(dp), allocatable :: scalar_potential(:, :)
        complex(dp), allocatable :: vector_potential(:, :)

        status = 1
        if (allocated(matrix)) deallocate(matrix)
        if (wave_number <= 0.0_dp .or. impedance <= 0.0_dp) return
        call assemble_maxwell_torus_curved_potential_operators_rwg_3d( &
            vertices, triangles, parameters, major_radius, minor_radius, &
            wave_number, quadrature_degree, tolerance, max_depth, &
            vector_potential, scalar_potential, status)
        if (status /= 0) return
        matrix = cmplx(0.0_dp, wave_number*impedance, dp)*vector_potential - &
            cmplx(0.0_dp, impedance/wave_number, dp)*scalar_potential
        status = 0
    end subroutine assemble_maxwell_torus_curved_efie_rwg_3d

    subroutine assemble_maxwell_torus_curved_potential_operators_rwg_3d( &
            vertices, triangles, parameters, major_radius, minor_radius, &
            wave_number, quadrature_degree, tolerance, max_depth, &
            vector_potential, scalar_potential, status, decaying_kernel)
        real(dp), intent(in) :: vertices(:, :), parameters(:, :)
        real(dp), intent(in) :: major_radius, minor_radius, wave_number
        real(dp), intent(in) :: tolerance
        integer, intent(in) :: triangles(:, :), quadrature_degree, max_depth
        complex(dp), allocatable, intent(out) :: vector_potential(:, :)
        complex(dp), allocatable, intent(out) :: scalar_potential(:, :)
        integer, intent(out) :: status
        logical, optional, intent(in) :: decaying_kernel

        integer, allocatable :: edge_triangles(:, :), edge_vertices(:, :)
        real(dp), allocatable :: reference_divergence(:, :)
        complex(dp) :: contribution, scalar_green
        real(dp) :: divergence, jacobian, point(3), rwg_value(3)
        integer :: basis, first, first_panel, first_slot, second, second_panel
        integer :: second_slot
        logical :: use_decaying_kernel

        status = 1
        if (allocated(vector_potential)) deallocate(vector_potential)
        if (allocated(scalar_potential)) deallocate(scalar_potential)
        use_decaying_kernel = .false.
        if (present(decaying_kernel)) use_decaying_kernel = decaying_kernel
        call build_maxwell_rwg_surface_space( &
            vertices, triangles, edge_vertices, edge_triangles, status)
        if (status /= 0) return
        if (size(edge_vertices, 2) == 0) then
            status = 1
            return
        end if
        allocate( &
            vector_potential(size(edge_vertices, 2), size(edge_vertices, 2)), &
            scalar_potential(size(edge_vertices, 2), size(edge_vertices, 2)), &
            reference_divergence(2, size(edge_vertices, 2)))
        vector_potential = cmplx(0.0_dp, 0.0_dp, dp)
        scalar_potential = cmplx(0.0_dp, 0.0_dp, dp)
        do basis = 1, size(edge_vertices, 2)
            do first_slot = 1, 2
                call evaluate_maxwell_torus_curved_rwg_basis( &
                    vertices, triangles, parameters, edge_vertices, &
                    edge_triangles, basis, edge_triangles(first_slot, basis), &
                    major_radius, minor_radius, 1.0_dp/3.0_dp, 1.0_dp/3.0_dp, &
                    point, rwg_value, divergence, jacobian, status)
                if (status /= 0) return
                reference_divergence(first_slot, basis) = divergence*jacobian
            end do
        end do
        do first = 1, size(edge_vertices, 2)
            do second = 1, first
                do first_slot = 1, 2
                    first_panel = edge_triangles(first_slot, first)
                    do second_slot = 1, 2
                        second_panel = edge_triangles(second_slot, second)
                        if (first_panel == second_panel) then
                            call &
                                integrate_maxwell_torus_curved_coincident_rwg_pair_3d( &
                                vertices, triangles, parameters, edge_vertices, &
                                edge_triangles, first, first_panel, second, &
                                major_radius, minor_radius, wave_number, &
                                quadrature_degree, contribution, status, &
                                scalar_green, use_decaying_kernel)
                        else
                            call &
                                integrate_maxwell_torus_curved_adjacent_rwg_pair_3d( &
                                vertices, triangles, parameters, edge_vertices, &
                                edge_triangles, first, first_panel, second, &
                                second_panel, major_radius, minor_radius, &
                                wave_number, quadrature_degree, tolerance, &
                                max_depth, contribution, status, scalar_green, &
                                use_decaying_kernel)
                        end if
                        if (status /= 0) return
                        vector_potential(first, second) = &
                            vector_potential(first, second) + contribution
                        scalar_potential(first, second) = &
                            scalar_potential(first, second) + &
                            reference_divergence(first_slot, first)* &
                            reference_divergence(second_slot, second)* &
                            scalar_green
                    end do
                end do
                vector_potential(second, first) = vector_potential(first, second)
                scalar_potential(second, first) = scalar_potential(first, second)
            end do
        end do
        status = 0
    end subroutine assemble_maxwell_torus_curved_potential_operators_rwg_3d

    subroutine integrate_maxwell_torus_curved_adjacent_rwg_pair_3d( &
            vertices, triangles, parameters, edge_vertices, edge_triangles, &
            first_basis, first_panel, second_basis, second_panel, major_radius, &
            minor_radius, wave_number, quadrature_degree, tolerance, max_depth, &
            value, status, scalar_value, decaying_kernel)
        real(dp), intent(in) :: vertices(:, :), parameters(:, :)
        real(dp), intent(in) :: major_radius, minor_radius, wave_number
        real(dp), intent(in) :: tolerance
        integer, intent(in) :: triangles(:, :), edge_vertices(:, :)
        integer, intent(in) :: edge_triangles(:, :)
        integer, intent(in) :: first_basis, first_panel, second_basis
        integer, intent(in) :: second_panel, quadrature_degree, max_depth
        complex(dp), intent(out) :: value
        integer, intent(out) :: status
        complex(dp), optional, intent(out) :: scalar_value
        logical, optional, intent(in) :: decaying_kernel

        real(dp), allocatable :: eta(:), weights(:), xi(:)
        real(dp) :: reference_triangle(2, 3)
        complex(dp) :: scalar_integral
        logical :: use_decaying_kernel

        value = cmplx(0.0_dp, 0.0_dp, dp)
        scalar_integral = cmplx(0.0_dp, 0.0_dp, dp)
        if (present(scalar_value)) scalar_value = scalar_integral
        use_decaying_kernel = .false.
        if (present(decaying_kernel)) use_decaying_kernel = decaying_kernel
        status = 1
        if (major_radius <= minor_radius .or. minor_radius <= 0.0_dp) return
        if (wave_number < 0.0_dp .or. tolerance <= 0.0_dp) return
        if (max_depth < 1 .or. first_panel == second_panel) return
        if (.not. any(edge_triangles(:, first_basis) == first_panel)) return
        if (.not. any(edge_triangles(:, second_basis) == second_panel)) return
        call triangle_duffy_quadrature( &
            quadrature_degree, xi, eta, weights, status)
        if (status /= 0) return
        reference_triangle = reshape( &
            [0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, 1.0_dp], [2, 3])
        call integrate_adaptive_torus_curved_rwg_pair( &
            vertices, triangles, parameters, edge_vertices, edge_triangles, &
            first_basis, first_panel, second_basis, second_panel, major_radius, &
            minor_radius, wave_number, reference_triangle, reference_triangle, &
            xi, eta, weights, tolerance, present(scalar_value), &
            use_decaying_kernel, 0, max_depth, value, scalar_integral, status)
        if (present(scalar_value)) scalar_value = scalar_integral
    end subroutine integrate_maxwell_torus_curved_adjacent_rwg_pair_3d

    recursive subroutine integrate_adaptive_torus_curved_rwg_pair( &
            vertices, triangles, parameters, edge_vertices, edge_triangles, &
            first_basis, first_panel, second_basis, second_panel, major_radius, &
            minor_radius, wave_number, first_reference, second_reference, xi, &
            eta, weights, tolerance, need_scalar, decaying_kernel, depth, &
            max_depth, value, scalar_value, status)
        real(dp), intent(in) :: vertices(:, :), parameters(:, :)
        real(dp), intent(in) :: major_radius, minor_radius, wave_number
        integer, intent(in) :: triangles(:, :), edge_vertices(:, :)
        integer, intent(in) :: edge_triangles(:, :)
        integer, intent(in) :: first_basis, first_panel, second_basis
        integer, intent(in) :: second_panel, depth, max_depth
        real(dp), intent(in) :: first_reference(2, 3)
        real(dp), intent(in) :: second_reference(2, 3)
        real(dp), intent(in) :: xi(:), eta(:), weights(:), tolerance
        logical, intent(in) :: need_scalar, decaying_kernel
        complex(dp), intent(out) :: value, scalar_value
        integer, intent(out) :: status

        real(dp) :: first_children(2, 3, 4), second_children(2, 3, 4)
        complex(dp) :: coarse, coarse_scalar, contribution
        complex(dp) :: contribution_scalar, refined, refined_scalar
        integer :: first_child, second_child

        call integrate_regular_torus_curved_rwg_pair( &
            vertices, triangles, parameters, edge_vertices, edge_triangles, &
            first_basis, first_panel, second_basis, second_panel, major_radius, &
            minor_radius, wave_number, first_reference, second_reference, xi, &
            eta, weights, coarse, coarse_scalar, decaying_kernel, status)
        if (status /= 0) return
        call subdivide_reference_triangle(first_reference, first_children)
        call subdivide_reference_triangle(second_reference, second_children)
        refined = cmplx(0.0_dp, 0.0_dp, dp)
        refined_scalar = cmplx(0.0_dp, 0.0_dp, dp)
        do first_child = 1, 4
            do second_child = 1, 4
                call integrate_regular_torus_curved_rwg_pair( &
                    vertices, triangles, parameters, edge_vertices, &
                    edge_triangles, first_basis, first_panel, second_basis, &
                    second_panel, major_radius, minor_radius, wave_number, &
                    first_children(:, :, first_child), &
                    second_children(:, :, second_child), xi, eta, weights, &
                    contribution, contribution_scalar, decaying_kernel, status)
                if (status /= 0) return
                refined = refined + contribution
                refined_scalar = refined_scalar + contribution_scalar
            end do
        end do
        if (depth + 1 >= max_depth .or. (abs(refined - coarse) <= &
            tolerance*max(tiny(1.0_dp), abs(refined)) .and. &
            (.not. need_scalar .or. abs(refined_scalar - coarse_scalar) <= &
            tolerance*max(tiny(1.0_dp), abs(refined_scalar))))) then
            value = refined
            scalar_value = refined_scalar
            status = 0
            return
        end if
        value = cmplx(0.0_dp, 0.0_dp, dp)
        scalar_value = cmplx(0.0_dp, 0.0_dp, dp)
        do first_child = 1, 4
            do second_child = 1, 4
                if (torus_reference_children_touch( &
                    parameters, triangles, first_panel, second_panel, &
                    major_radius, minor_radius, &
                    first_children(:, :, first_child), &
                    second_children(:, :, second_child))) then
                    call integrate_adaptive_torus_curved_rwg_pair( &
                        vertices, triangles, parameters, edge_vertices, &
                        edge_triangles, first_basis, first_panel, second_basis, &
                        second_panel, major_radius, minor_radius, wave_number, &
                        first_children(:, :, first_child), &
                        second_children(:, :, second_child), xi, eta, weights, &
                        tolerance, need_scalar, decaying_kernel, depth + 1, &
                        max_depth, contribution, contribution_scalar, status)
                else
                    call integrate_regular_torus_curved_rwg_pair( &
                        vertices, triangles, parameters, edge_vertices, &
                        edge_triangles, first_basis, first_panel, second_basis, &
                        second_panel, major_radius, minor_radius, wave_number, &
                        first_children(:, :, first_child), &
                        second_children(:, :, second_child), xi, eta, weights, &
                        contribution, contribution_scalar, decaying_kernel, &
                        status)
                end if
                if (status /= 0) return
                value = value + contribution
                scalar_value = scalar_value + contribution_scalar
            end do
        end do
        status = 0
    end subroutine integrate_adaptive_torus_curved_rwg_pair

    subroutine integrate_regular_torus_curved_rwg_pair( &
            vertices, triangles, parameters, edge_vertices, edge_triangles, &
            first_basis, first_panel, second_basis, second_panel, major_radius, &
            minor_radius, wave_number, first_reference, second_reference, xi, &
            eta, weights, value, scalar_value, decaying_kernel, status)
        real(dp), intent(in) :: vertices(:, :), parameters(:, :)
        real(dp), intent(in) :: major_radius, minor_radius, wave_number
        integer, intent(in) :: triangles(:, :), edge_vertices(:, :)
        integer, intent(in) :: edge_triangles(:, :)
        integer, intent(in) :: first_basis, first_panel, second_basis
        integer, intent(in) :: second_panel
        real(dp), intent(in) :: first_reference(2, 3)
        real(dp), intent(in) :: second_reference(2, 3)
        real(dp), intent(in) :: xi(:), eta(:), weights(:)
        complex(dp), intent(out) :: value, scalar_value
        logical, intent(in) :: decaying_kernel
        integer, intent(out) :: status

        real(dp) :: first_divergence, first_jacobian, first_point(3)
        real(dp) :: first_reference_jacobian, first_value(3), first_xi_eta(2)
        real(dp) :: physical_distance, second_divergence, second_jacobian
        real(dp) :: second_point(3), second_reference_jacobian
        real(dp) :: second_value(3), second_xi_eta(2)
        complex(dp) :: green
        integer :: first_node, second_node

        value = cmplx(0.0_dp, 0.0_dp, dp)
        scalar_value = cmplx(0.0_dp, 0.0_dp, dp)
        first_reference_jacobian = reference_triangle_jacobian(first_reference)
        second_reference_jacobian = &
            reference_triangle_jacobian(second_reference)
        do first_node = 1, size(weights)
            first_xi_eta = reference_point( &
                first_reference, xi(first_node), eta(first_node))
            call evaluate_maxwell_torus_curved_rwg_basis( &
                vertices, triangles, parameters, edge_vertices, &
                edge_triangles, first_basis, first_panel, major_radius, &
                minor_radius, first_xi_eta(1), first_xi_eta(2), first_point, &
                first_value, first_divergence, first_jacobian, status)
            if (status /= 0) return
            do second_node = 1, size(weights)
                second_xi_eta = reference_point( &
                    second_reference, xi(second_node), eta(second_node))
                call evaluate_maxwell_torus_curved_rwg_basis( &
                    vertices, triangles, parameters, edge_vertices, &
                    edge_triangles, second_basis, second_panel, major_radius, &
                    minor_radius, second_xi_eta(1), second_xi_eta(2), &
                    second_point, second_value, second_divergence, &
                    second_jacobian, status)
                if (status /= 0) return
                physical_distance = norm2(first_point - second_point)
                green = torus_boundary_green( &
                    wave_number, physical_distance, decaying_kernel)
                value = value + first_reference_jacobian* &
                    second_reference_jacobian*weights(first_node)* &
                    weights(second_node)*first_jacobian*second_jacobian*green* &
                    dot_product(first_value, second_value)
                scalar_value = scalar_value + first_reference_jacobian* &
                    second_reference_jacobian*weights(first_node)* &
                    weights(second_node)*green
            end do
        end do
        status = 0
    end subroutine integrate_regular_torus_curved_rwg_pair

    subroutine integrate_maxwell_torus_curved_coincident_rwg_pair_3d( &
            vertices, triangles, parameters, edge_vertices, edge_triangles, &
            first_basis, panel, second_basis, major_radius, minor_radius, &
            wave_number, quadrature_degree, value, status, scalar_value, &
            decaying_kernel)
        real(dp), intent(in) :: vertices(:, :), parameters(:, :)
        real(dp), intent(in) :: major_radius, minor_radius, wave_number
        integer, intent(in) :: triangles(:, :), edge_vertices(:, :)
        integer, intent(in) :: edge_triangles(:, :)
        integer, intent(in) :: first_basis, panel, second_basis
        integer, intent(in) :: quadrature_degree
        complex(dp), intent(out) :: value
        integer, intent(out) :: status
        complex(dp), optional, intent(out) :: scalar_value
        logical, optional, intent(in) :: decaying_kernel

        real(dp), allocatable :: eta(:), line_nodes(:), line_weights(:)
        real(dp), allocatable :: weights(:), xi(:)
        real(dp) :: direction(2), first_divergence, first_jacobian
        real(dp) :: first_point(3), first_value(3), reference_vertices(2, 3)
        real(dp) :: physical_distance, rho, second_divergence, second_eta
        real(dp) :: second_jacobian, second_point(3), second_value(3)
        real(dp) :: second_xi, t, wedge_first(2), wedge_jacobian
        real(dp) :: wedge_second(2)
        complex(dp) :: green, scalar_integral
        logical :: use_decaying_kernel
        integer :: first_node, line_count, radial_node, tangential_node, wedge

        value = cmplx(0.0_dp, 0.0_dp, dp)
        scalar_integral = cmplx(0.0_dp, 0.0_dp, dp)
        if (present(scalar_value)) scalar_value = scalar_integral
        use_decaying_kernel = .false.
        if (present(decaying_kernel)) use_decaying_kernel = decaying_kernel
        status = 1
        if (major_radius <= minor_radius .or. minor_radius <= 0.0_dp) return
        if (wave_number < 0.0_dp) return
        if (.not. any(edge_triangles(:, first_basis) == panel)) return
        if (.not. any(edge_triangles(:, second_basis) == panel)) return
        call triangle_duffy_quadrature( &
            quadrature_degree, xi, eta, weights, status)
        if (status /= 0) return
        line_count = max(2, quadrature_degree)
        allocate(line_nodes(line_count), line_weights(line_count))
        call gauss_legendre_ab( &
            line_count, 0.0_dp, 1.0_dp, line_nodes, line_weights)
        reference_vertices = reshape( &
            [0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, 1.0_dp], [2, 3])
        do first_node = 1, size(weights)
            call evaluate_maxwell_torus_curved_rwg_basis( &
                vertices, triangles, parameters, edge_vertices, &
                edge_triangles, first_basis, panel, major_radius, &
                minor_radius, xi(first_node), eta(first_node), first_point, &
                first_value, first_divergence, first_jacobian, status)
            if (status /= 0) return
            do wedge = 1, 3
                wedge_first = &
                    reference_vertices(:, wedge) - &
                    [xi(first_node), eta(first_node)]
                wedge_second = &
                    reference_vertices(:, modulo(wedge, 3) + 1) - &
                    [xi(first_node), eta(first_node)]
                wedge_jacobian = abs( &
                    wedge_first(1)*wedge_second(2) - &
                    wedge_first(2)*wedge_second(1))
                do radial_node = 1, line_count
                    rho = line_nodes(radial_node)
                    do tangential_node = 1, line_count
                        t = line_nodes(tangential_node)
                        direction = (1.0_dp - t)*wedge_first + t*wedge_second
                        second_xi = xi(first_node) + rho*direction(1)
                        second_eta = eta(first_node) + rho*direction(2)
                        call evaluate_maxwell_torus_curved_rwg_basis( &
                            vertices, triangles, parameters, edge_vertices, &
                            edge_triangles, second_basis, panel, major_radius, &
                            minor_radius, second_xi, second_eta, second_point, &
                            second_value, second_divergence, second_jacobian, &
                            status)
                        if (status /= 0) return
                        physical_distance = norm2(first_point - second_point)
                        green = torus_boundary_green( &
                            wave_number, physical_distance, &
                            use_decaying_kernel)
                        value = value + weights(first_node)*first_jacobian* &
                            line_weights(radial_node)* &
                            line_weights(tangential_node)*rho*wedge_jacobian* &
                            second_jacobian*green* &
                            dot_product(first_value, second_value)
                        scalar_integral = scalar_integral + &
                            weights(first_node)*line_weights(radial_node)* &
                            line_weights(tangential_node)*rho*wedge_jacobian* &
                            green
                    end do
                end do
            end do
        end do
        if (present(scalar_value)) scalar_value = scalar_integral
        status = 0
    end subroutine integrate_maxwell_torus_curved_coincident_rwg_pair_3d

    subroutine evaluate_maxwell_torus_curved_far_field_rwg_3d( &
            vertices, triangles, parameters, major_radius, minor_radius, &
            coefficients, direction, wave_number, impedance, &
            quadrature_degree, far_field, status)
        real(dp), intent(in) :: vertices(:, :), parameters(:, :)
        real(dp), intent(in) :: major_radius, minor_radius, direction(3)
        complex(dp), intent(in) :: coefficients(:)
        real(dp), intent(in) :: wave_number, impedance
        integer, intent(in) :: triangles(:, :), quadrature_degree
        complex(dp), intent(out) :: far_field(3)
        integer, intent(out) :: status

        integer, allocatable :: edge_triangles(:, :), edge_vertices(:, :)
        real(dp), allocatable :: eta(:), weights(:), xi(:)
        real(dp) :: basis_value(3), divergence, jacobian, point(3)
        complex(dp) :: phase, surface_current(3), transverse_current(3)
        integer :: basis, node, panel

        far_field = cmplx(0.0_dp, 0.0_dp, dp)
        status = 1
        if (major_radius <= minor_radius .or. minor_radius <= 0.0_dp) return
        if (wave_number <= 0.0_dp .or. impedance <= 0.0_dp) return
        if (abs(norm2(direction) - 1.0_dp) > &
            128.0_dp*epsilon(1.0_dp)) return
        call build_maxwell_rwg_surface_space( &
            vertices, triangles, edge_vertices, edge_triangles, status)
        if (status /= 0) return
        if (size(coefficients) /= size(edge_vertices, 2)) then
            status = 1
            return
        end if
        call triangle_duffy_quadrature( &
            quadrature_degree, xi, eta, weights, status)
        if (status /= 0) return
        do panel = 1, size(triangles, 2)
            do node = 1, size(weights)
                surface_current = cmplx(0.0_dp, 0.0_dp, dp)
                do basis = 1, size(edge_vertices, 2)
                    if (.not. any(edge_triangles(:, basis) == panel)) cycle
                    call evaluate_maxwell_torus_curved_rwg_basis( &
                        vertices, triangles, parameters, edge_vertices, &
                        edge_triangles, basis, panel, major_radius, &
                        minor_radius, xi(node), eta(node), point, basis_value, &
                        divergence, jacobian, status)
                    if (status /= 0) return
                    surface_current = surface_current + &
                        coefficients(basis)*basis_value
                end do
                transverse_current = surface_current - &
                    direction*sum(direction*surface_current)
                phase = exp(cmplx( &
                    0.0_dp, -wave_number*dot_product(direction, point), dp))
                far_field = far_field + weights(node)*jacobian*phase* &
                    transverse_current
            end do
        end do
        far_field = cmplx(0.0_dp, wave_number*impedance/ &
            (4.0_dp*acos(-1.0_dp)), dp)*far_field
        status = 0
    end subroutine evaluate_maxwell_torus_curved_far_field_rwg_3d

    subroutine assemble_maxwell_torus_curved_plane_wave_rhs_rwg_3d( &
            vertices, triangles, parameters, major_radius, minor_radius, &
            direction, polarization, wave_number, quadrature_degree, &
            right_hand_side, status)
        real(dp), intent(in) :: vertices(:, :), parameters(:, :)
        real(dp), intent(in) :: major_radius, minor_radius, direction(3)
        complex(dp), intent(in) :: polarization(3)
        real(dp), intent(in) :: wave_number
        integer, intent(in) :: triangles(:, :), quadrature_degree
        complex(dp), allocatable, intent(out) :: right_hand_side(:)
        integer, intent(out) :: status

        integer, allocatable :: edge_triangles(:, :), edge_vertices(:, :)
        real(dp), allocatable :: eta(:), weights(:), xi(:)
        real(dp) :: basis_value(3), divergence, jacobian, point(3)
        complex(dp) :: incident_field(3)
        integer :: basis, node, panel

        status = 1
        if (allocated(right_hand_side)) deallocate(right_hand_side)
        if (major_radius <= minor_radius .or. minor_radius <= 0.0_dp) return
        if (wave_number < 0.0_dp) return
        if (abs(norm2(direction) - 1.0_dp) > &
            128.0_dp*epsilon(1.0_dp)) return
        if (abs(sum(polarization*direction)) > &
            128.0_dp*epsilon(1.0_dp)*max(1.0_dp, &
            sqrt(sum(abs(polarization)**2)))) return
        if (sqrt(sum(abs(polarization)**2)) <= tiny(1.0_dp)) return
        call build_maxwell_rwg_surface_space( &
            vertices, triangles, edge_vertices, edge_triangles, status)
        if (status /= 0) return
        if (size(edge_vertices, 2) == 0) then
            status = 1
            return
        end if
        call triangle_duffy_quadrature( &
            quadrature_degree, xi, eta, weights, status)
        if (status /= 0) return
        allocate(right_hand_side(size(edge_vertices, 2)))
        right_hand_side = cmplx(0.0_dp, 0.0_dp, dp)
        do basis = 1, size(edge_vertices, 2)
            do panel = 1, 2
                do node = 1, size(weights)
                    call evaluate_maxwell_torus_curved_rwg_basis( &
                        vertices, triangles, parameters, edge_vertices, &
                        edge_triangles, basis, edge_triangles(panel, basis), &
                        major_radius, minor_radius, xi(node), eta(node), point, &
                        basis_value, divergence, jacobian, status)
                    if (status /= 0) return
                    incident_field = polarization*exp(cmplx( &
                        0.0_dp, wave_number*dot_product(direction, point), dp))
                    right_hand_side(basis) = right_hand_side(basis) - &
                        weights(node)*jacobian* &
                        sum(cmplx(basis_value, 0.0_dp, dp)*incident_field)
                end do
            end do
        end do
        status = 0
    end subroutine assemble_maxwell_torus_curved_plane_wave_rhs_rwg_3d

    subroutine assemble_maxwell_torus_curved_rwg_mass_matrix( &
            vertices, triangles, parameters, major_radius, minor_radius, &
            quadrature_degree, matrix, status)
        real(dp), intent(in) :: vertices(:, :), parameters(:, :)
        real(dp), intent(in) :: major_radius, minor_radius
        integer, intent(in) :: triangles(:, :), quadrature_degree
        real(dp), allocatable, intent(out) :: matrix(:, :)
        integer, intent(out) :: status

        integer, allocatable :: edge_triangles(:, :), edge_vertices(:, :)
        real(dp), allocatable :: eta(:), values(:, :), weights(:), xi(:)
        real(dp) :: divergence, jacobian, point(3)
        integer :: basis, node, panel, test_basis

        status = 1
        if (allocated(matrix)) deallocate(matrix)
        if (size(vertices, 1) /= 3) return
        if (size(parameters, 1) /= 2) return
        if (size(parameters, 2) /= size(vertices, 2)) return
        if (major_radius <= minor_radius .or. minor_radius <= 0.0_dp) return
        call build_maxwell_rwg_surface_space( &
            vertices, triangles, edge_vertices, edge_triangles, status)
        if (status /= 0) return
        if (size(edge_vertices, 2) == 0) then
            status = 1
            return
        end if
        call triangle_duffy_quadrature( &
            quadrature_degree, xi, eta, weights, status)
        if (status /= 0) return
        allocate( &
            matrix(size(edge_vertices, 2), size(edge_vertices, 2)), &
            values(3, size(edge_vertices, 2)))
        matrix = 0.0_dp
        do panel = 1, size(triangles, 2)
            do node = 1, size(weights)
                values = 0.0_dp
                do basis = 1, size(edge_vertices, 2)
                    if (.not. any(edge_triangles(:, basis) == panel)) cycle
                    call evaluate_maxwell_torus_curved_rwg_basis( &
                        vertices, triangles, parameters, edge_vertices, &
                        edge_triangles, basis, panel, major_radius, &
                        minor_radius, xi(node), eta(node), point, &
                        values(:, basis), divergence, jacobian, status)
                    if (status /= 0) return
                end do
                do basis = 1, size(edge_vertices, 2)
                    if (.not. any(edge_triangles(:, basis) == panel)) cycle
                    do test_basis = 1, size(edge_vertices, 2)
                        if (.not. any( &
                            edge_triangles(:, test_basis) == panel)) cycle
                        matrix(test_basis, basis) = &
                            matrix(test_basis, basis) + &
                            weights(node)*jacobian*dot_product( &
                            values(:, test_basis), values(:, basis))
                    end do
                end do
            end do
        end do
        status = 0
    end subroutine assemble_maxwell_torus_curved_rwg_mass_matrix

    pure subroutine evaluate_maxwell_torus_curved_rwg_basis( &
            vertices, triangles, parameters, edge_vertices, edge_triangles, &
            basis, panel, major_radius, minor_radius, xi, eta, point, value, &
            surface_divergence, surface_jacobian, status)
        real(dp), intent(in) :: vertices(:, :), parameters(:, :)
        real(dp), intent(in) :: major_radius, minor_radius, xi, eta
        integer, intent(in) :: triangles(:, :), edge_vertices(:, :)
        integer, intent(in) :: edge_triangles(:, :), basis, panel
        real(dp), intent(out) :: point(3), value(3), surface_divergence
        real(dp), intent(out) :: surface_jacobian
        integer, intent(out) :: status

        real(dp) :: edge_length, opposite_coordinates(2)
        real(dp) :: panel_parameters(2, 3), tangent_eta(3), tangent_xi(3)
        integer :: local, next, opposite, orientation

        point = 0.0_dp
        value = 0.0_dp
        surface_divergence = 0.0_dp
        surface_jacobian = 0.0_dp
        status = 1
        if (size(vertices, 1) /= 3 .or. size(parameters, 1) /= 2) return
        if (size(parameters, 2) /= size(vertices, 2)) return
        if (basis < 1 .or. basis > size(edge_vertices, 2)) return
        if (panel < 1 .or. panel > size(triangles, 2)) return
        if (.not. any(edge_triangles(:, basis) == panel)) return
        orientation = 0
        opposite = 0
        do local = 1, 3
            next = modulo(local, 3) + 1
            if (triangles(local, panel) == edge_vertices(1, basis) .and. &
                triangles(next, panel) == edge_vertices(2, basis)) then
                orientation = 1
                opposite = modulo(next, 3) + 1
                exit
            end if
            if (triangles(local, panel) == edge_vertices(2, basis) .and. &
                triangles(next, panel) == edge_vertices(1, basis)) then
                orientation = -1
                opposite = modulo(next, 3) + 1
                exit
            end if
        end do
        if (orientation == 0) return
        do local = 1, 3
            panel_parameters(:, local) = &
                parameters(:, triangles(local, panel))
        end do
        call evaluate_torus_curved_panel( &
            panel_parameters, major_radius, minor_radius, xi, eta, point, &
            tangent_xi, tangent_eta, surface_jacobian, status)
        if (status /= 0) return
        select case (opposite)
        case (1)
            opposite_coordinates = [0.0_dp, 0.0_dp]
        case (2)
            opposite_coordinates = [1.0_dp, 0.0_dp]
        case (3)
            opposite_coordinates = [0.0_dp, 1.0_dp]
        end select
        edge_length = norm2( &
            vertices(:, edge_vertices(2, basis)) - &
            vertices(:, edge_vertices(1, basis)))
        value = real(orientation, dp)*edge_length/surface_jacobian*( &
            (xi - opposite_coordinates(1))*tangent_xi + &
            (eta - opposite_coordinates(2))*tangent_eta)
        surface_divergence = &
            2.0_dp*real(orientation, dp)*edge_length/surface_jacobian
        status = 0
    end subroutine evaluate_maxwell_torus_curved_rwg_basis

    pure logical function torus_reference_children_touch( &
            parameters, triangles, first_panel, second_panel, major_radius, &
            minor_radius, first_reference, second_reference) result(touch)
        real(dp), intent(in) :: parameters(:, :)
        integer, intent(in) :: triangles(:, :), first_panel, second_panel
        real(dp), intent(in) :: major_radius, minor_radius
        real(dp), intent(in) :: first_reference(2, 3)
        real(dp), intent(in) :: second_reference(2, 3)

        real(dp) :: first_panel_parameters(2, 3), first_point(3), jacobian
        real(dp) :: scale, second_panel_parameters(2, 3), second_point(3)
        real(dp) :: tangent_eta(3), tangent_xi(3)
        integer :: first_vertex, second_vertex, status

        scale = max(1.0_dp, major_radius + minor_radius)
        touch = .false.
        do first_vertex = 1, 3
            first_panel_parameters(:, first_vertex) = &
                parameters(:, triangles(first_vertex, first_panel))
            second_panel_parameters(:, first_vertex) = &
                parameters(:, triangles(first_vertex, second_panel))
        end do
        do first_vertex = 1, 3
            call evaluate_torus_curved_panel( &
                first_panel_parameters, major_radius, minor_radius, &
                first_reference(1, first_vertex), &
                first_reference(2, first_vertex), first_point, tangent_xi, &
                tangent_eta, jacobian, status)
            if (status /= 0) return
            do second_vertex = 1, 3
                call evaluate_torus_curved_panel( &
                    second_panel_parameters, major_radius, minor_radius, &
                    second_reference(1, second_vertex), &
                    second_reference(2, second_vertex), second_point, &
                    tangent_xi, tangent_eta, jacobian, status)
                if (status /= 0) return
                if (norm2(first_point - second_point) <= &
                    512.0_dp*epsilon(1.0_dp)*scale) then
                    touch = .true.
                    return
                end if
            end do
        end do
    end function torus_reference_children_touch

    pure subroutine map_refined_torus_point_to_parent( &
            refined_panel_parameters, parent_panel_parameters, xi, eta, &
            parent_xi, parent_eta, status)
        real(dp), intent(in) :: refined_panel_parameters(2, 3)
        real(dp), intent(in) :: parent_panel_parameters(2, 3), xi, eta
        real(dp), intent(out) :: parent_xi, parent_eta
        integer, intent(out) :: status

        real(dp) :: determinant, parent_parameters(2, 3)
        real(dp) :: refined_parameters(2, 3), right_hand_side(2), target(2)

        parent_xi = 0.0_dp
        parent_eta = 0.0_dp
        status = 1
        parent_parameters = parent_panel_parameters
        refined_parameters = refined_panel_parameters
        call unwrap_torus_parameter( &
            parent_parameters(:, 1), parent_parameters(:, 2))
        call unwrap_torus_parameter( &
            parent_parameters(:, 1), parent_parameters(:, 3))
        call unwrap_torus_parameter( &
            refined_parameters(:, 1), refined_parameters(:, 2))
        call unwrap_torus_parameter( &
            refined_parameters(:, 1), refined_parameters(:, 3))
        target = refined_parameters(:, 1) + &
            xi*(refined_parameters(:, 2) - refined_parameters(:, 1)) + &
            eta*(refined_parameters(:, 3) - refined_parameters(:, 1))
        call unwrap_torus_parameter(parent_parameters(:, 1), target)
        right_hand_side = target - parent_parameters(:, 1)
        determinant = &
            (parent_parameters(1, 2) - parent_parameters(1, 1))* &
            (parent_parameters(2, 3) - parent_parameters(2, 1)) - &
            (parent_parameters(1, 3) - parent_parameters(1, 1))* &
            (parent_parameters(2, 2) - parent_parameters(2, 1))
        if (abs(determinant) <= tiny(1.0_dp)) return
        parent_xi = ( &
            right_hand_side(1)* &
            (parent_parameters(2, 3) - parent_parameters(2, 1)) - &
            right_hand_side(2)* &
            (parent_parameters(1, 3) - parent_parameters(1, 1)))/determinant
        parent_eta = ( &
            (parent_parameters(1, 2) - parent_parameters(1, 1))* &
            right_hand_side(2) - &
            (parent_parameters(2, 2) - parent_parameters(2, 1))* &
            right_hand_side(1))/determinant
        if (parent_xi < -256.0_dp*epsilon(1.0_dp)) return
        if (parent_eta < -256.0_dp*epsilon(1.0_dp)) return
        if (parent_xi + parent_eta > 1.0_dp + &
            256.0_dp*epsilon(1.0_dp)) return
        parent_xi = max(0.0_dp, parent_xi)
        parent_eta = max(0.0_dp, parent_eta)
        status = 0
    end subroutine map_refined_torus_point_to_parent

    pure subroutine unwrap_torus_parameter(reference, value)
        real(dp), intent(in) :: reference(2)
        real(dp), intent(inout) :: value(2)

        integer :: coordinate

        do coordinate = 1, 2
            do while (value(coordinate) - reference(coordinate) > &
                    acos(-1.0_dp))
                value(coordinate) = value(coordinate) - 2.0_dp*acos(-1.0_dp)
            end do
            do while (value(coordinate) - reference(coordinate) < &
                    -acos(-1.0_dp))
                value(coordinate) = value(coordinate) + 2.0_dp*acos(-1.0_dp)
            end do
        end do
    end subroutine unwrap_torus_parameter

    pure function torus_unit_normal(point, major_radius) result(normal)
        real(dp), intent(in) :: point(3), major_radius
        real(dp) :: normal(3)
        real(dp) :: cylindrical_radius

        cylindrical_radius = sqrt(point(1)**2 + point(2)**2)
        normal = [ &
            (cylindrical_radius - major_radius)*point(1)/cylindrical_radius, &
            (cylindrical_radius - major_radius)*point(2)/cylindrical_radius, &
            point(3)]
        normal = normal/norm2(normal)
    end function torus_unit_normal

    pure function real_cross_product(first, second) result(product)
        real(dp), intent(in) :: first(3), second(3)
        real(dp) :: product(3)

        product = [ &
            first(2)*second(3) - first(3)*second(2), &
            first(3)*second(1) - first(1)*second(3), &
            first(1)*second(2) - first(2)*second(1)]
    end function real_cross_product

    pure function torus_boundary_green( &
            wave_number, distance, decaying_kernel) result(value)
        real(dp), intent(in) :: wave_number, distance
        logical, intent(in) :: decaying_kernel
        complex(dp) :: value

        if (decaying_kernel) then
            value = cmplx(exp(-wave_number*distance)/ &
                (4.0_dp*acos(-1.0_dp)*distance), 0.0_dp, dp)
        else
            value = exp(cmplx(0.0_dp, wave_number*distance, dp))/ &
                (4.0_dp*acos(-1.0_dp)*distance)
        end if
    end function torus_boundary_green

    pure function reference_point(vertices, xi, eta) result(point)
        real(dp), intent(in) :: vertices(2, 3), xi, eta
        real(dp) :: point(2)

        point = vertices(:, 1) + xi*(vertices(:, 2) - vertices(:, 1)) + &
            eta*(vertices(:, 3) - vertices(:, 1))
    end function reference_point

    pure function reference_triangle_jacobian(vertices) result(jacobian)
        real(dp), intent(in) :: vertices(2, 3)
        real(dp) :: jacobian

        jacobian = abs( &
            (vertices(1, 2) - vertices(1, 1))* &
            (vertices(2, 3) - vertices(2, 1)) - &
            (vertices(2, 2) - vertices(2, 1))* &
            (vertices(1, 3) - vertices(1, 1)))
    end function reference_triangle_jacobian

    pure subroutine subdivide_reference_triangle(vertices, children)
        real(dp), intent(in) :: vertices(2, 3)
        real(dp), intent(out) :: children(2, 3, 4)

        real(dp) :: midpoint_12(2), midpoint_23(2), midpoint_31(2)

        midpoint_12 = 0.5_dp*(vertices(:, 1) + vertices(:, 2))
        midpoint_23 = 0.5_dp*(vertices(:, 2) + vertices(:, 3))
        midpoint_31 = 0.5_dp*(vertices(:, 3) + vertices(:, 1))
        children(:, 1, 1) = vertices(:, 1)
        children(:, 2, 1) = midpoint_12
        children(:, 3, 1) = midpoint_31
        children(:, 1, 2) = midpoint_12
        children(:, 2, 2) = vertices(:, 2)
        children(:, 3, 2) = midpoint_23
        children(:, 1, 3) = midpoint_31
        children(:, 2, 3) = midpoint_23
        children(:, 3, 3) = vertices(:, 3)
        children(:, 1, 4) = midpoint_12
        children(:, 2, 4) = midpoint_23
        children(:, 3, 4) = midpoint_31
    end subroutine subdivide_reference_triangle

end module fortfem_maxwell_torus_curved_rwg
