module fortfem_maxwell_sphere_curved_rwg
    !! Surface-Piola image of the affine RWG basis on a radial sphere panel.
    use fortfem_kinds, only: dp
    use fortfem_maxwell_rwg_surface, only: build_maxwell_rwg_surface_space
    use fortfem_sphere_curved_panel, only: evaluate_sphere_curved_panel
    use fortfem_triangle_duffy_quadrature, only: triangle_duffy_quadrature
    use fortnum_quadrature, only: gauss_legendre_ab
    implicit none
    private

    public :: evaluate_maxwell_sphere_curved_rwg_basis
    public :: assemble_maxwell_sphere_curved_rwg_mass_matrix
    public :: assemble_maxwell_sphere_curved_plane_wave_rhs_rwg_3d
    public :: evaluate_maxwell_sphere_curved_far_field_rwg_3d
    public :: integrate_maxwell_sphere_curved_coincident_rwg_pair_3d
    public :: integrate_maxwell_sphere_curved_adjacent_rwg_pair_3d
    public :: assemble_maxwell_sphere_curved_vector_potential_rwg_3d
    public :: assemble_maxwell_sphere_curved_potential_operators_rwg_3d
    public :: assemble_maxwell_sphere_curved_efie_rwg_3d
    public :: solve_maxwell_pec_sphere_curved_efie_rwg_3d
    public :: evaluate_maxwell_sphere_curved_magnetic_field_rwg_3d

    interface
        subroutine zgesv(n, nrhs, a, lda, ipiv, b, ldb, info)
            import :: dp
            integer, intent(in) :: n, nrhs, lda, ldb
            complex(dp), intent(inout) :: a(lda, *)
            integer, intent(out) :: ipiv(*)
            complex(dp), intent(inout) :: b(ldb, *)
            integer, intent(out) :: info
        end subroutine zgesv
    end interface

contains

    subroutine evaluate_maxwell_sphere_curved_magnetic_field_rwg_3d( &
            vertices, triangles, radius, coefficients, observation, &
            wave_number, quadrature_degree, magnetic_field, status)
        real(dp), intent(in) :: vertices(:, :), radius, observation(3)
        complex(dp), intent(in) :: coefficients(:)
        real(dp), intent(in) :: wave_number
        integer, intent(in) :: triangles(:, :), quadrature_degree
        complex(dp), intent(out) :: magnetic_field(3)
        integer, intent(out) :: status

        integer, allocatable :: edge_triangles(:, :), edge_vertices(:, :)
        real(dp), allocatable :: eta(:), weights(:), xi(:)
        real(dp) :: basis_value(3), displacement(3), divergence, jacobian
        real(dp) :: point(3), distance
        complex(dp) :: gradient_green(3), green, surface_current(3)
        integer :: basis, node, panel

        magnetic_field = cmplx(0.0_dp, 0.0_dp, dp)
        status = 1
        if (radius <= 0.0_dp .or. wave_number < 0.0_dp) return
        call build_maxwell_rwg_surface_space( &
            vertices, triangles, edge_vertices, edge_triangles, status)
        if (status /= 0 .or. size(coefficients) /= size(edge_vertices, 2)) return
        call triangle_duffy_quadrature( &
            quadrature_degree, xi, eta, weights, status)
        if (status /= 0) return
        do panel = 1, size(triangles, 2)
            do node = 1, size(weights)
                surface_current = cmplx(0.0_dp, 0.0_dp, dp)
                do basis = 1, size(edge_vertices, 2)
                    if (.not. any(edge_triangles(:, basis) == panel)) cycle
                    call evaluate_maxwell_sphere_curved_rwg_basis( &
                        vertices, triangles, edge_vertices, edge_triangles, &
                        basis, panel, radius, xi(node), eta(node), point, &
                        basis_value, divergence, jacobian, status)
                    if (status /= 0) return
                    surface_current = surface_current + &
                        coefficients(basis)*basis_value
                end do
                displacement = observation - point
                distance = norm2(displacement)
                if (distance <= 128.0_dp*epsilon(1.0_dp)*radius) return
                green = exp(cmplx(0.0_dp, wave_number*distance, dp))/ &
                    (4.0_dp*acos(-1.0_dp)*distance)
                gradient_green = green* &
                    (cmplx(0.0_dp, wave_number, dp) - 1.0_dp/distance)* &
                    displacement/distance
                magnetic_field = magnetic_field + weights(node)*jacobian* &
                    complex_cross_product(gradient_green, surface_current)
            end do
        end do
        status = 0
    end subroutine evaluate_maxwell_sphere_curved_magnetic_field_rwg_3d

    subroutine solve_maxwell_pec_sphere_curved_efie_rwg_3d( &
            vertices, triangles, radius, direction, polarization, wave_number, &
            impedance, quadrature_degree, tolerance, max_depth, density, status)
        real(dp), intent(in) :: vertices(:, :), radius, direction(3)
        complex(dp), intent(in) :: polarization(3)
        real(dp), intent(in) :: wave_number, impedance, tolerance
        integer, intent(in) :: triangles(:, :), quadrature_degree, max_depth
        complex(dp), allocatable, intent(out) :: density(:)
        integer, intent(out) :: status

        complex(dp), allocatable :: matrix(:, :), right_hand_side(:)
        complex(dp), allocatable :: right_hand_side_matrix(:, :)
        integer, allocatable :: pivots(:)
        integer :: info

        status = 1
        call assemble_maxwell_sphere_curved_efie_rwg_3d( &
            vertices, triangles, radius, wave_number, impedance, &
            quadrature_degree, tolerance, max_depth, matrix, status)
        if (status /= 0) return
        call assemble_maxwell_sphere_curved_plane_wave_rhs_rwg_3d( &
            vertices, triangles, radius, direction, polarization, wave_number, &
            quadrature_degree, right_hand_side, status)
        if (status /= 0) return
        allocate( &
            density(size(right_hand_side)), &
            right_hand_side_matrix(size(right_hand_side), 1), &
            pivots(size(right_hand_side)))
        right_hand_side_matrix(:, 1) = right_hand_side
        call zgesv( &
            size(right_hand_side), 1, matrix, size(right_hand_side), pivots, &
            right_hand_side_matrix, size(right_hand_side), info)
        if (info /= 0) then
            status = 2
            return
        end if
        density = right_hand_side_matrix(:, 1)
        status = 0
    end subroutine solve_maxwell_pec_sphere_curved_efie_rwg_3d

    subroutine assemble_maxwell_sphere_curved_efie_rwg_3d( &
            vertices, triangles, radius, wave_number, impedance, &
            quadrature_degree, tolerance, max_depth, matrix, status)
        real(dp), intent(in) :: vertices(:, :), radius, wave_number, impedance
        integer, intent(in) :: triangles(:, :), quadrature_degree, max_depth
        real(dp), intent(in) :: tolerance
        complex(dp), allocatable, intent(out) :: matrix(:, :)
        integer, intent(out) :: status

        complex(dp), allocatable :: scalar_potential(:, :)
        complex(dp), allocatable :: vector_potential(:, :)

        status = 1
        if (radius <= 0.0_dp .or. wave_number <= 0.0_dp .or. &
            impedance <= 0.0_dp) return
        call assemble_maxwell_sphere_curved_potential_operators_rwg_3d( &
            vertices, triangles, radius, wave_number, quadrature_degree, &
            tolerance, max_depth, vector_potential, scalar_potential, status)
        if (status /= 0) return
        matrix = cmplx(0.0_dp, wave_number*impedance, dp)*vector_potential - &
            cmplx(0.0_dp, impedance/wave_number, dp)*scalar_potential
        status = 0
    end subroutine assemble_maxwell_sphere_curved_efie_rwg_3d

    subroutine assemble_maxwell_sphere_curved_potential_operators_rwg_3d( &
            vertices, triangles, radius, wave_number, quadrature_degree, &
            tolerance, max_depth, vector_potential, scalar_potential, status)
        real(dp), intent(in) :: vertices(:, :), radius, wave_number, tolerance
        integer, intent(in) :: triangles(:, :), quadrature_degree, max_depth
        complex(dp), allocatable, intent(out) :: vector_potential(:, :)
        complex(dp), allocatable, intent(out) :: scalar_potential(:, :)
        integer, intent(out) :: status

        integer, allocatable :: edge_triangles(:, :), edge_vertices(:, :)
        real(dp), allocatable :: reference_divergence(:, :)
        complex(dp) :: contribution, scalar_green
        real(dp) :: divergence, jacobian, point(3), rwg_value(3)
        integer :: basis, first, first_panel, first_slot, second, second_panel
        integer :: second_slot

        status = 1
        call build_maxwell_rwg_surface_space( &
            vertices, triangles, edge_vertices, edge_triangles, status)
        if (status /= 0 .or. size(edge_vertices, 2) == 0) return
        allocate( &
            vector_potential(size(edge_vertices, 2), size(edge_vertices, 2)), &
            scalar_potential(size(edge_vertices, 2), size(edge_vertices, 2)), &
            reference_divergence(2, size(edge_vertices, 2)))
        vector_potential = cmplx(0.0_dp, 0.0_dp, dp)
        scalar_potential = cmplx(0.0_dp, 0.0_dp, dp)
        do basis = 1, size(edge_vertices, 2)
            do first_slot = 1, 2
                call evaluate_maxwell_sphere_curved_rwg_basis( &
                    vertices, triangles, edge_vertices, edge_triangles, basis, &
                    edge_triangles(first_slot, basis), radius, 1.0_dp/3.0_dp, &
                    1.0_dp/3.0_dp, point, rwg_value, divergence, jacobian, &
                    status)
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
                                integrate_maxwell_sphere_curved_coincident_rwg_pair_3d( &
                                vertices, triangles, edge_vertices, &
                                edge_triangles, first, first_panel, second, &
                                radius, wave_number, quadrature_degree, &
                                contribution, status, scalar_green)
                        else
                            call &
                                integrate_maxwell_sphere_curved_adjacent_rwg_pair_3d( &
                                vertices, triangles, edge_vertices, &
                                edge_triangles, first, first_panel, second, &
                                second_panel, radius, wave_number, &
                                quadrature_degree, tolerance, max_depth, &
                                contribution, status, scalar_green)
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
    end subroutine assemble_maxwell_sphere_curved_potential_operators_rwg_3d

    subroutine assemble_maxwell_sphere_curved_vector_potential_rwg_3d( &
            vertices, triangles, radius, wave_number, quadrature_degree, &
            tolerance, max_depth, matrix, status)
        real(dp), intent(in) :: vertices(:, :), radius, wave_number, tolerance
        integer, intent(in) :: triangles(:, :), quadrature_degree, max_depth
        complex(dp), allocatable, intent(out) :: matrix(:, :)
        integer, intent(out) :: status

        integer, allocatable :: edge_triangles(:, :), edge_vertices(:, :)
        complex(dp) :: contribution
        integer :: first, first_panel, first_slot, second, second_panel
        integer :: second_slot

        status = 1
        call build_maxwell_rwg_surface_space( &
            vertices, triangles, edge_vertices, edge_triangles, status)
        if (status /= 0 .or. size(edge_vertices, 2) == 0) return
        allocate(matrix(size(edge_vertices, 2), size(edge_vertices, 2)))
        matrix = cmplx(0.0_dp, 0.0_dp, dp)
        do first = 1, size(edge_vertices, 2)
            do second = 1, first
                do first_slot = 1, 2
                    first_panel = edge_triangles(first_slot, first)
                    do second_slot = 1, 2
                        second_panel = edge_triangles(second_slot, second)
                        if (first_panel == second_panel) then
                            call &
                                integrate_maxwell_sphere_curved_coincident_rwg_pair_3d( &
                                vertices, triangles, edge_vertices, &
                                edge_triangles, first, first_panel, second, &
                                radius, wave_number, quadrature_degree, &
                                contribution, status)
                        else
                            call &
                                integrate_maxwell_sphere_curved_adjacent_rwg_pair_3d( &
                                vertices, triangles, edge_vertices, &
                                edge_triangles, first, first_panel, second, &
                                second_panel, radius, wave_number, &
                                quadrature_degree, tolerance, max_depth, &
                                contribution, status)
                        end if
                        if (status /= 0) return
                        matrix(first, second) = &
                            matrix(first, second) + contribution
                    end do
                end do
                matrix(second, first) = matrix(first, second)
            end do
        end do
        status = 0
    end subroutine assemble_maxwell_sphere_curved_vector_potential_rwg_3d

    subroutine integrate_maxwell_sphere_curved_adjacent_rwg_pair_3d( &
            vertices, triangles, edge_vertices, edge_triangles, first_basis, &
            first_panel, second_basis, second_panel, radius, wave_number, &
            quadrature_degree, tolerance, max_depth, value, status, scalar_value)
        real(dp), intent(in) :: vertices(:, :), radius, wave_number, tolerance
        integer, intent(in) :: triangles(:, :), edge_vertices(:, :)
        integer, intent(in) :: edge_triangles(:, :)
        integer, intent(in) :: first_basis, first_panel, second_basis
        integer, intent(in) :: second_panel, quadrature_degree, max_depth
        complex(dp), intent(out) :: value
        integer, intent(out) :: status
        complex(dp), optional, intent(out) :: scalar_value

        real(dp), allocatable :: eta(:), weights(:), xi(:)
        real(dp) :: reference_triangle(2, 3)
        complex(dp) :: scalar_integral

        value = cmplx(0.0_dp, 0.0_dp, dp)
        if (present(scalar_value)) scalar_value = cmplx(0.0_dp, 0.0_dp, dp)
        status = 1
        if (radius <= 0.0_dp .or. wave_number < 0.0_dp .or. &
            tolerance <= 0.0_dp .or. max_depth < 1) return
        if (first_panel == second_panel) return
        if (.not. any(edge_triangles(:, first_basis) == first_panel) .or. &
            .not. any(edge_triangles(:, second_basis) == second_panel)) return
        call triangle_duffy_quadrature( &
            quadrature_degree, xi, eta, weights, status)
        if (status /= 0) return
        reference_triangle = reshape( &
            [0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, 1.0_dp], [2, 3])
        call integrate_adaptive_curved_rwg_pair( &
            vertices, triangles, edge_vertices, edge_triangles, first_basis, &
            first_panel, second_basis, second_panel, radius, wave_number, &
            reference_triangle, reference_triangle, xi, eta, weights, &
            tolerance, present(scalar_value), 0, max_depth, value, &
            scalar_integral, status)
        if (present(scalar_value)) scalar_value = scalar_integral
    end subroutine integrate_maxwell_sphere_curved_adjacent_rwg_pair_3d

    recursive subroutine integrate_adaptive_curved_rwg_pair( &
            vertices, triangles, edge_vertices, edge_triangles, first_basis, &
            first_panel, second_basis, second_panel, radius, wave_number, &
            first_reference, second_reference, xi, eta, weights, tolerance, &
            need_scalar, depth, max_depth, value, scalar_value, status)
        real(dp), intent(in) :: vertices(:, :), radius, wave_number
        integer, intent(in) :: triangles(:, :), edge_vertices(:, :)
        integer, intent(in) :: edge_triangles(:, :)
        integer, intent(in) :: first_basis, first_panel, second_basis
        integer, intent(in) :: second_panel, depth, max_depth
        real(dp), intent(in) :: first_reference(2, 3)
        real(dp), intent(in) :: second_reference(2, 3)
        real(dp), intent(in) :: xi(:), eta(:), weights(:), tolerance
        logical, intent(in) :: need_scalar
        complex(dp), intent(out) :: value, scalar_value
        integer, intent(out) :: status

        real(dp) :: first_children(2, 3, 4), second_children(2, 3, 4)
        complex(dp) :: coarse, coarse_scalar, contribution
        complex(dp) :: contribution_scalar, refined, refined_scalar
        integer :: first_child, second_child

        call integrate_regular_curved_rwg_pair( &
            vertices, triangles, edge_vertices, edge_triangles, first_basis, &
            first_panel, second_basis, second_panel, radius, wave_number, &
            first_reference, second_reference, xi, eta, weights, coarse, &
            coarse_scalar, status)
        if (status /= 0) return
        call subdivide_reference_triangle(first_reference, first_children)
        call subdivide_reference_triangle(second_reference, second_children)
        refined = cmplx(0.0_dp, 0.0_dp, dp)
        refined_scalar = cmplx(0.0_dp, 0.0_dp, dp)
        do first_child = 1, 4
            do second_child = 1, 4
                call integrate_regular_curved_rwg_pair( &
                    vertices, triangles, edge_vertices, edge_triangles, &
                    first_basis, first_panel, second_basis, second_panel, &
                    radius, wave_number, first_children(:, :, first_child), &
                    second_children(:, :, second_child), xi, eta, weights, &
                    contribution, contribution_scalar, status)
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
                if (curved_reference_children_touch( &
                    vertices, triangles, first_panel, second_panel, radius, &
                    first_children(:, :, first_child), &
                    second_children(:, :, second_child))) then
                    call integrate_adaptive_curved_rwg_pair( &
                        vertices, triangles, edge_vertices, edge_triangles, &
                        first_basis, first_panel, second_basis, second_panel, &
                        radius, wave_number, first_children(:, :, first_child), &
                        second_children(:, :, second_child), xi, eta, weights, &
                        tolerance, need_scalar, depth + 1, max_depth, &
                        contribution, &
                        contribution_scalar, status)
                else
                    call integrate_regular_curved_rwg_pair( &
                        vertices, triangles, edge_vertices, edge_triangles, &
                        first_basis, first_panel, second_basis, second_panel, &
                        radius, wave_number, first_children(:, :, first_child), &
                        second_children(:, :, second_child), xi, eta, weights, &
                        contribution, contribution_scalar, status)
                end if
                if (status /= 0) return
                value = value + contribution
                scalar_value = scalar_value + contribution_scalar
            end do
        end do
        status = 0
    end subroutine integrate_adaptive_curved_rwg_pair

    subroutine integrate_regular_curved_rwg_pair( &
            vertices, triangles, edge_vertices, edge_triangles, first_basis, &
            first_panel, second_basis, second_panel, radius, wave_number, &
            first_reference, second_reference, xi, eta, weights, value, &
            scalar_value, status)
        real(dp), intent(in) :: vertices(:, :), radius, wave_number
        integer, intent(in) :: triangles(:, :), edge_vertices(:, :)
        integer, intent(in) :: edge_triangles(:, :)
        integer, intent(in) :: first_basis, first_panel, second_basis
        integer, intent(in) :: second_panel
        real(dp), intent(in) :: first_reference(2, 3)
        real(dp), intent(in) :: second_reference(2, 3)
        real(dp), intent(in) :: xi(:), eta(:), weights(:)
        complex(dp), intent(out) :: value, scalar_value
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
            call evaluate_maxwell_sphere_curved_rwg_basis( &
                vertices, triangles, edge_vertices, edge_triangles, &
                first_basis, first_panel, radius, first_xi_eta(1), &
                first_xi_eta(2), first_point, first_value, first_divergence, &
                first_jacobian, status)
            if (status /= 0) return
            do second_node = 1, size(weights)
                second_xi_eta = reference_point( &
                    second_reference, xi(second_node), eta(second_node))
                call evaluate_maxwell_sphere_curved_rwg_basis( &
                    vertices, triangles, edge_vertices, edge_triangles, &
                    second_basis, second_panel, radius, second_xi_eta(1), &
                    second_xi_eta(2), second_point, second_value, &
                    second_divergence, second_jacobian, status)
                if (status /= 0) return
                physical_distance = norm2(first_point - second_point)
                green = helmholtz_green(wave_number, physical_distance)
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
    end subroutine integrate_regular_curved_rwg_pair

    subroutine integrate_maxwell_sphere_curved_coincident_rwg_pair_3d( &
            vertices, triangles, edge_vertices, edge_triangles, first_basis, &
            panel, second_basis, radius, wave_number, quadrature_degree, &
            value, status, scalar_value)
        real(dp), intent(in) :: vertices(:, :), radius, wave_number
        integer, intent(in) :: triangles(:, :), edge_vertices(:, :)
        integer, intent(in) :: edge_triangles(:, :)
        integer, intent(in) :: first_basis, panel, second_basis
        integer, intent(in) :: quadrature_degree
        complex(dp), intent(out) :: value
        integer, intent(out) :: status
        complex(dp), optional, intent(out) :: scalar_value

        real(dp), allocatable :: eta(:), line_nodes(:), line_weights(:)
        real(dp), allocatable :: weights(:), xi(:)
        real(dp) :: direction(2), first_divergence, first_jacobian
        real(dp) :: first_point(3), first_value(3), reference_vertices(2, 3)
        real(dp) :: physical_distance, rho, second_divergence, second_eta
        real(dp) :: second_jacobian, second_point(3), second_value(3)
        real(dp) :: second_xi, t, wedge_first(2), wedge_jacobian
        real(dp) :: wedge_second(2)
        complex(dp) :: green, scalar_integral
        integer :: first_node, line_count, radial_node, tangential_node, wedge

        value = cmplx(0.0_dp, 0.0_dp, dp)
        scalar_integral = cmplx(0.0_dp, 0.0_dp, dp)
        if (present(scalar_value)) scalar_value = scalar_integral
        status = 1
        if (radius <= 0.0_dp .or. wave_number < 0.0_dp) return
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
            call evaluate_maxwell_sphere_curved_rwg_basis( &
                vertices, triangles, edge_vertices, edge_triangles, &
                first_basis, panel, radius, xi(first_node), eta(first_node), &
                first_point, first_value, first_divergence, first_jacobian, &
                status)
            if (status /= 0) return
            do wedge = 1, 3
                wedge_first(1) = reference_vertices(1, wedge) - xi(first_node)
                wedge_first(2) = reference_vertices(2, wedge) - eta(first_node)
                wedge_second(1) = &
                    reference_vertices(1, modulo(wedge, 3) + 1) - xi(first_node)
                wedge_second(2) = &
                    reference_vertices(2, modulo(wedge, 3) + 1) - eta(first_node)
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
                        call evaluate_maxwell_sphere_curved_rwg_basis( &
                            vertices, triangles, edge_vertices, edge_triangles, &
                            second_basis, panel, radius, second_xi, second_eta, &
                            second_point, second_value, second_divergence, &
                            second_jacobian, status)
                        if (status /= 0) return
                        physical_distance = norm2(first_point - second_point)
                        green = helmholtz_green( &
                            wave_number, physical_distance)
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
    end subroutine integrate_maxwell_sphere_curved_coincident_rwg_pair_3d

    subroutine evaluate_maxwell_sphere_curved_far_field_rwg_3d( &
            vertices, triangles, radius, coefficients, direction, wave_number, &
            impedance, quadrature_degree, far_field, status)
        real(dp), intent(in) :: vertices(:, :), radius, direction(3)
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
        if (radius <= 0.0_dp .or. wave_number <= 0.0_dp .or. &
            impedance <= 0.0_dp) return
        if (abs(norm2(direction) - 1.0_dp) > 128.0_dp*epsilon(1.0_dp)) return
        call build_maxwell_rwg_surface_space( &
            vertices, triangles, edge_vertices, edge_triangles, status)
        if (status /= 0 .or. size(coefficients) /= size(edge_vertices, 2)) return
        call triangle_duffy_quadrature( &
            quadrature_degree, xi, eta, weights, status)
        if (status /= 0) return
        do panel = 1, size(triangles, 2)
            do node = 1, size(weights)
                surface_current = cmplx(0.0_dp, 0.0_dp, dp)
                do basis = 1, size(edge_vertices, 2)
                    if (.not. any(edge_triangles(:, basis) == panel)) cycle
                    call evaluate_maxwell_sphere_curved_rwg_basis( &
                        vertices, triangles, edge_vertices, edge_triangles, &
                        basis, panel, radius, xi(node), eta(node), point, &
                        basis_value, divergence, jacobian, status)
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
    end subroutine evaluate_maxwell_sphere_curved_far_field_rwg_3d

    subroutine assemble_maxwell_sphere_curved_plane_wave_rhs_rwg_3d( &
            vertices, triangles, radius, direction, polarization, wave_number, &
            quadrature_degree, right_hand_side, status)
        real(dp), intent(in) :: vertices(:, :), radius, direction(3)
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
        if (radius <= 0.0_dp .or. wave_number < 0.0_dp) return
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
                do node = 1, size(weights)
                    call evaluate_maxwell_sphere_curved_rwg_basis( &
                        vertices, triangles, edge_vertices, edge_triangles, &
                        basis, edge_triangles(panel, basis), radius, xi(node), &
                        eta(node), point, basis_value, divergence, jacobian, &
                        status)
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
    end subroutine assemble_maxwell_sphere_curved_plane_wave_rhs_rwg_3d

    subroutine assemble_maxwell_sphere_curved_rwg_mass_matrix( &
            vertices, triangles, radius, quadrature_degree, matrix, status)
        real(dp), intent(in) :: vertices(:, :), radius
        integer, intent(in) :: triangles(:, :), quadrature_degree
        real(dp), allocatable, intent(out) :: matrix(:, :)
        integer, intent(out) :: status

        integer, allocatable :: edge_triangles(:, :), edge_vertices(:, :)
        real(dp), allocatable :: eta(:), weights(:), xi(:)
        real(dp), allocatable :: values(:, :)
        real(dp) :: divergence, jacobian, point(3)
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
            do node = 1, size(weights)
                values = 0.0_dp
                do basis = 1, size(edge_vertices, 2)
                    if (.not. any(edge_triangles(:, basis) == panel)) cycle
                    call evaluate_maxwell_sphere_curved_rwg_basis( &
                        vertices, triangles, edge_vertices, edge_triangles, &
                        basis, panel, radius, xi(node), eta(node), point, &
                        values(:, basis), divergence, jacobian, status)
                    if (status /= 0) return
                end do
                do basis = 1, size(edge_vertices, 2)
                    if (.not. any(edge_triangles(:, basis) == panel)) cycle
                    do test_basis = 1, size(edge_vertices, 2)
                        if (.not. any( &
                            edge_triangles(:, test_basis) == panel)) cycle
                        matrix(test_basis, basis) = matrix(test_basis, basis) + &
                            weights(node)*jacobian*dot_product( &
                            values(:, test_basis), values(:, basis))
                    end do
                end do
            end do
        end do
        status = 0
    end subroutine assemble_maxwell_sphere_curved_rwg_mass_matrix

    pure subroutine evaluate_maxwell_sphere_curved_rwg_basis( &
            vertices, triangles, edge_vertices, edge_triangles, basis, panel, &
            radius, xi, eta, point, value, surface_divergence, &
            surface_jacobian, status)
        real(dp), intent(in) :: vertices(:, :), radius, xi, eta
        integer, intent(in) :: triangles(:, :), edge_vertices(:, :)
        integer, intent(in) :: edge_triangles(:, :), basis, panel
        real(dp), intent(out) :: point(3), value(3), surface_divergence
        real(dp), intent(out) :: surface_jacobian
        integer, intent(out) :: status

        real(dp) :: edge_length, opposite_coordinates(2), panel_vertices(3, 3)
        real(dp) :: tangent_eta(3), tangent_xi(3)
        integer :: local, next, opposite, orientation

        point = 0.0_dp
        value = 0.0_dp
        surface_divergence = 0.0_dp
        surface_jacobian = 0.0_dp
        status = 1
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
            panel_vertices(:, local) = vertices(:, triangles(local, panel))
        end do
        call evaluate_sphere_curved_panel( &
            panel_vertices, radius, xi, eta, point, &
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
    end subroutine evaluate_maxwell_sphere_curved_rwg_basis

    pure function helmholtz_green(wave_number, radius) result(value)
        real(dp), intent(in) :: wave_number, radius
        complex(dp) :: value

        value = exp(cmplx(0.0_dp, wave_number*radius, dp))/ &
            (4.0_dp*acos(-1.0_dp)*radius)
    end function helmholtz_green

    pure function complex_cross_product(first, second) result(product)
        complex(dp), intent(in) :: first(3), second(3)
        complex(dp) :: product(3)

        product = [ &
            first(2)*second(3) - first(3)*second(2), &
            first(3)*second(1) - first(1)*second(3), &
            first(1)*second(2) - first(2)*second(1)]
    end function complex_cross_product

    pure logical function curved_reference_children_touch( &
            vertices, triangles, first_panel, second_panel, radius, &
            first_reference, second_reference) result(touch)
        real(dp), intent(in) :: vertices(:, :), radius
        integer, intent(in) :: triangles(:, :), first_panel, second_panel
        real(dp), intent(in) :: first_reference(2, 3)
        real(dp), intent(in) :: second_reference(2, 3)

        real(dp) :: first_panel_vertices(3, 3), first_point(3), jacobian
        real(dp) :: scale, second_panel_vertices(3, 3), second_point(3)
        real(dp) :: tangent_eta(3), tangent_xi(3)
        integer :: first_vertex, second_vertex, status

        scale = max(1.0_dp, radius)
        touch = .false.
        do first_vertex = 1, 3
            first_panel_vertices(:, first_vertex) = &
                vertices(:, triangles(first_vertex, first_panel))
            second_panel_vertices(:, first_vertex) = &
                vertices(:, triangles(first_vertex, second_panel))
        end do
        do first_vertex = 1, 3
            call evaluate_sphere_curved_panel( &
                first_panel_vertices, radius, &
                first_reference(1, first_vertex), &
                first_reference(2, first_vertex), first_point, tangent_xi, &
                tangent_eta, jacobian, status)
            if (status /= 0) return
            do second_vertex = 1, 3
                call evaluate_sphere_curved_panel( &
                    second_panel_vertices, radius, &
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
    end function curved_reference_children_touch

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

end module fortfem_maxwell_sphere_curved_rwg
