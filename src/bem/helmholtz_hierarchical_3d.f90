module fortfem_helmholtz_hierarchical_3d
    !! Pointer-free Barnes--Hut Helmholtz single- and combined-layer operators.
    !!
    !! Analytical Galerkin self panels are retained. Far clusters use complex
    !! area-weighted monopoles with the outgoing Helmholtz kernel; near leaves
    !! use panel-centroid interactions. Storage is O(N).
    use fortfem_helmholtz_galerkin_3d, only: &
        assemble_helmholtz_single_layer_p0_3d
    use fortfem_kinds, only: dp
    use fortfem_panel_cluster_tree_3d, only: build_panel_cluster_tree_3d
    use fortnum_krylov, only: complex_gmres_operator, KRYLOV_OK
    implicit none

    private

    public :: apply_helmholtz_single_layer_p0_hierarchical_3d
    public :: apply_helmholtz_cfie_p0_hierarchical_3d
    public :: solve_helmholtz_dirichlet_p0_hierarchical_3d
    public :: solve_helmholtz_cfie_p0_hierarchical_3d

contains

    subroutine apply_helmholtz_single_layer_p0_hierarchical_3d( &
            vertices, triangles, density, wave_number, opening_angle, &
            leaf_size, result, status, interaction_count)
        real(dp), intent(in) :: vertices(:, :), wave_number, opening_angle
        integer, intent(in) :: triangles(:, :), leaf_size
        complex(dp), intent(in) :: density(:)
        complex(dp), allocatable, intent(out) :: result(:)
        integer, intent(out) :: status, interaction_count

        integer, allocatable :: first_child(:), inverse_position(:)
        integer, allocatable :: last(:), first(:), permutation(:)
        integer, allocatable :: second_child(:)
        real(dp), allocatable :: areas(:), centers(:, :)
        real(dp), allocatable :: node_centers(:, :), radii(:)
        complex(dp), allocatable :: self_diagonal(:)
        integer :: node_count

        status = 1
        interaction_count = 0
        if (size(vertices, 1) /= 3 .or. size(vertices, 2) < 3) return
        if (size(triangles, 1) /= 3 .or. size(triangles, 2) < 1) return
        if (size(density) /= size(triangles, 2)) return
        if (any(triangles < 1) .or. &
            any(triangles > size(vertices, 2))) return
        if (wave_number < 0.0_dp) return
        if (opening_angle <= 0.0_dp .or. opening_angle >= 1.0_dp) return
        if (leaf_size < 1) return

        call build_panel_cluster_tree_3d( &
            vertices, triangles, leaf_size, areas, centers, permutation, &
            inverse_position, first, last, first_child, second_child, &
            node_centers, radii, node_count, status)
        if (status /= 0) return
        call precompute_self_diagonal( &
            vertices, triangles, wave_number, self_diagonal, status)
        if (status /= 0) return
        allocate(result(size(density)))
        call apply_prebuilt( &
            density, wave_number, opening_angle, areas, centers, permutation, &
            inverse_position, first, last, first_child, second_child, &
            node_centers, radii, node_count, self_diagonal, result, &
            interaction_count)
        status = 0
    end subroutine apply_helmholtz_single_layer_p0_hierarchical_3d

    subroutine apply_helmholtz_cfie_p0_hierarchical_3d( &
            vertices, triangles, density, wave_number, coupling, &
            opening_angle, leaf_size, result, status, interaction_count)
        real(dp), intent(in) :: vertices(:, :), wave_number, coupling
        real(dp), intent(in) :: opening_angle
        integer, intent(in) :: triangles(:, :), leaf_size
        complex(dp), intent(in) :: density(:)
        complex(dp), allocatable, intent(out) :: result(:)
        integer, intent(out) :: status, interaction_count

        integer, allocatable :: first_child(:), inverse_position(:)
        integer, allocatable :: last(:), first(:), permutation(:)
        integer, allocatable :: second_child(:)
        real(dp), allocatable :: areas(:), centers(:, :), normals(:, :)
        real(dp), allocatable :: node_centers(:, :), radii(:)
        complex(dp), allocatable :: self_diagonal(:)
        integer :: node_count

        status = 1
        interaction_count = 0
        if (.not. valid_geometry_and_parameters( &
            vertices, triangles, wave_number, opening_angle, leaf_size)) return
        if (size(density) /= size(triangles, 2) .or. coupling <= 0.0_dp) return
        call build_panel_cluster_tree_3d( &
            vertices, triangles, leaf_size, areas, centers, permutation, &
            inverse_position, first, last, first_child, second_child, &
            node_centers, radii, node_count, status)
        if (status /= 0) return
        call precompute_self_diagonal( &
            vertices, triangles, wave_number, self_diagonal, status)
        if (status /= 0) return
        call compute_panel_normals(vertices, triangles, normals, status)
        if (status /= 0) return
        allocate(result(size(density)))
        call apply_cfie_prebuilt( &
            density, wave_number, coupling, opening_angle, areas, centers, &
            normals, permutation, inverse_position, first, last, first_child, &
            second_child, node_centers, radii, node_count, self_diagonal, &
            result, interaction_count)
        status = 0
    end subroutine apply_helmholtz_cfie_p0_hierarchical_3d

    subroutine solve_helmholtz_dirichlet_p0_hierarchical_3d( &
            vertices, triangles, boundary_value, wave_number, opening_angle, &
            leaf_size, tolerance, max_iterations, restart, density, status, &
            iterations, residual_norm, interaction_count)
        real(dp), intent(in) :: vertices(:, :), wave_number, opening_angle
        real(dp), intent(in) :: tolerance
        integer, intent(in) :: triangles(:, :), leaf_size
        integer, intent(in) :: max_iterations, restart
        complex(dp), intent(in) :: boundary_value
        complex(dp), allocatable, intent(out) :: density(:)
        integer, intent(out) :: status, iterations, interaction_count
        real(dp), intent(out) :: residual_norm

        integer, allocatable :: first_child(:), inverse_position(:)
        integer, allocatable :: last(:), first(:), permutation(:)
        integer, allocatable :: second_child(:)
        real(dp), allocatable :: areas(:), centers(:, :)
        real(dp), allocatable :: node_centers(:, :), radii(:)
        complex(dp), allocatable :: right_hand_side(:), self_diagonal(:)
        integer :: krylov_status, node_count

        status = 1
        iterations = 0
        interaction_count = 0
        residual_norm = huge(1.0_dp)
        if (size(vertices, 1) /= 3 .or. size(vertices, 2) < 3) return
        if (size(triangles, 1) /= 3 .or. size(triangles, 2) < 1) return
        if (any(triangles < 1) .or. &
            any(triangles > size(vertices, 2))) return
        if (wave_number < 0.0_dp) return
        if (opening_angle <= 0.0_dp .or. opening_angle >= 1.0_dp) return
        if (leaf_size < 1 .or. tolerance <= 0.0_dp) return
        if (max_iterations < 1 .or. restart < 1) return
        if (restart > max_iterations) return

        call build_panel_cluster_tree_3d( &
            vertices, triangles, leaf_size, areas, centers, permutation, &
            inverse_position, first, last, first_child, second_child, &
            node_centers, radii, node_count, status)
        if (status /= 0) return
        call precompute_self_diagonal( &
            vertices, triangles, wave_number, self_diagonal, status)
        if (status /= 0) return
        allocate(density(size(areas)), right_hand_side(size(areas)))
        density = cmplx(0.0_dp, 0.0_dp, dp)
        right_hand_side = boundary_value*areas
        call complex_gmres_operator( &
            matvec, right_hand_side, density, tolerance, max_iterations, &
            restart, krylov_status, iterations, residual_norm)
        if (krylov_status /= KRYLOV_OK) then
            status = 2
            return
        end if
        status = 0

    contains

        subroutine matvec(input, output)
            complex(dp), intent(in) :: input(:)
            complex(dp), intent(out) :: output(:)

            integer :: matvec_interactions

            call apply_prebuilt( &
                input, wave_number, opening_angle, areas, centers, &
                permutation, inverse_position, first, last, first_child, &
                second_child, node_centers, radii, node_count, self_diagonal, &
                output, matvec_interactions)
            interaction_count = interaction_count + matvec_interactions
        end subroutine matvec

    end subroutine solve_helmholtz_dirichlet_p0_hierarchical_3d

    subroutine solve_helmholtz_cfie_p0_hierarchical_3d( &
            vertices, triangles, boundary_value, wave_number, coupling, &
            opening_angle, leaf_size, tolerance, max_iterations, restart, &
            density, status, iterations, residual_norm, interaction_count)
        real(dp), intent(in) :: vertices(:, :), wave_number, coupling
        real(dp), intent(in) :: opening_angle, tolerance
        integer, intent(in) :: triangles(:, :), leaf_size
        integer, intent(in) :: max_iterations, restart
        complex(dp), intent(in) :: boundary_value
        complex(dp), allocatable, intent(out) :: density(:)
        integer, intent(out) :: status, iterations, interaction_count
        real(dp), intent(out) :: residual_norm

        integer, allocatable :: first_child(:), inverse_position(:)
        integer, allocatable :: last(:), first(:), permutation(:)
        integer, allocatable :: second_child(:)
        real(dp), allocatable :: areas(:), centers(:, :), normals(:, :)
        real(dp), allocatable :: node_centers(:, :), radii(:)
        complex(dp), allocatable :: right_hand_side(:), self_diagonal(:)
        integer :: krylov_status, node_count

        status = 1
        iterations = 0
        interaction_count = 0
        residual_norm = huge(1.0_dp)
        if (.not. valid_geometry_and_parameters( &
            vertices, triangles, wave_number, opening_angle, leaf_size)) return
        if (coupling <= 0.0_dp .or. tolerance <= 0.0_dp) return
        if (max_iterations < 1 .or. restart < 1) return
        if (restart > max_iterations) return
        call build_panel_cluster_tree_3d( &
            vertices, triangles, leaf_size, areas, centers, permutation, &
            inverse_position, first, last, first_child, second_child, &
            node_centers, radii, node_count, status)
        if (status /= 0) return
        call precompute_self_diagonal( &
            vertices, triangles, wave_number, self_diagonal, status)
        if (status /= 0) return
        call compute_panel_normals(vertices, triangles, normals, status)
        if (status /= 0) return
        allocate(density(size(areas)), right_hand_side(size(areas)))
        density = cmplx(0.0_dp, 0.0_dp, dp)
        right_hand_side = boundary_value*areas
        call complex_gmres_operator( &
            matvec, right_hand_side, density, tolerance, max_iterations, &
            restart, krylov_status, iterations, residual_norm)
        if (krylov_status /= KRYLOV_OK) then
            status = 2
            return
        end if
        status = 0

    contains

        subroutine matvec(input, output)
            complex(dp), intent(in) :: input(:)
            complex(dp), intent(out) :: output(:)

            integer :: matvec_interactions

            call apply_cfie_prebuilt( &
                input, wave_number, coupling, opening_angle, areas, centers, &
                normals, permutation, inverse_position, first, last, &
                first_child, second_child, node_centers, radii, node_count, &
                self_diagonal, output, matvec_interactions)
            interaction_count = interaction_count + matvec_interactions
        end subroutine matvec

    end subroutine solve_helmholtz_cfie_p0_hierarchical_3d

    subroutine precompute_self_diagonal( &
            vertices, triangles, wave_number, self_diagonal, status)
        real(dp), intent(in) :: vertices(:, :), wave_number
        integer, intent(in) :: triangles(:, :)
        complex(dp), allocatable, intent(out) :: self_diagonal(:)
        integer, intent(out) :: status

        complex(dp), allocatable :: self_matrix(:, :)
        real(dp) :: local_vertices(3, 3)
        integer :: local_triangle(3, 1), panel

        status = 1
        allocate(self_diagonal(size(triangles, 2)))
        local_triangle(:, 1) = [1, 2, 3]
        do panel = 1, size(triangles, 2)
            local_vertices = vertices(:, triangles(:, panel))
            call assemble_helmholtz_single_layer_p0_3d( &
                local_vertices, local_triangle, wave_number, 8, &
                self_matrix, status)
            if (status /= 0) return
            self_diagonal(panel) = self_matrix(1, 1)
        end do
        status = 0
    end subroutine precompute_self_diagonal

    subroutine apply_prebuilt( &
            density, wave_number, opening_angle, areas, centers, permutation, &
            inverse_position, first, last, first_child, second_child, &
            node_centers, radii, node_count, self_diagonal, result, &
            interaction_count)
        complex(dp), intent(in) :: density(:), self_diagonal(:)
        real(dp), intent(in) :: wave_number, opening_angle
        real(dp), intent(in) :: areas(:), centers(:, :)
        integer, intent(in) :: permutation(:), inverse_position(:)
        integer, intent(in) :: first(:), last(:), first_child(:)
        integer, intent(in) :: second_child(:), node_count
        real(dp), intent(in) :: node_centers(:, :), radii(:)
        complex(dp), intent(out) :: result(:)
        integer, intent(out) :: interaction_count

        complex(dp), allocatable :: charges(:)
        integer :: node, panel

        allocate(charges(node_count))
        do node = node_count, 1, -1
            if (first_child(node) == 0) then
                charges(node) = sum( &
                    density(permutation(first(node):last(node)))* &
                    areas(permutation(first(node):last(node))))
            else
                charges(node) = &
                    charges(first_child(node)) + charges(second_child(node))
            end if
        end do
        interaction_count = 0
        result = self_diagonal*density
        do panel = 1, size(density)
            interaction_count = interaction_count + 1
            call apply_node( &
                1, panel, inverse_position(panel), wave_number, &
                opening_angle, centers, areas, density, permutation, first, &
                last, first_child, second_child, node_centers, radii, &
                charges, result(panel), interaction_count)
        end do
    end subroutine apply_prebuilt

    subroutine apply_cfie_prebuilt( &
            density, wave_number, coupling, opening_angle, areas, centers, &
            normals, permutation, inverse_position, first, last, first_child, &
            second_child, node_centers, radii, node_count, self_diagonal, &
            result, interaction_count)
        complex(dp), intent(in) :: density(:), self_diagonal(:)
        real(dp), intent(in) :: wave_number, coupling, opening_angle
        real(dp), intent(in) :: areas(:), centers(:, :), normals(:, :)
        integer, intent(in) :: permutation(:), inverse_position(:)
        integer, intent(in) :: first(:), last(:), first_child(:)
        integer, intent(in) :: second_child(:), node_count
        real(dp), intent(in) :: node_centers(:, :), radii(:)
        complex(dp), intent(out) :: result(:)
        integer, intent(out) :: interaction_count

        complex(dp), allocatable :: charges(:), dipoles(:, :)
        integer :: node, panel, position

        allocate(charges(node_count), dipoles(3, node_count))
        do node = node_count, 1, -1
            if (first_child(node) == 0) then
                charges(node) = cmplx(0.0_dp, 0.0_dp, dp)
                dipoles(:, node) = cmplx(0.0_dp, 0.0_dp, dp)
                do position = first(node), last(node)
                    panel = permutation(position)
                    charges(node) = charges(node) + &
                        density(panel)*areas(panel)
                    dipoles(:, node) = dipoles(:, node) + &
                        density(panel)*areas(panel)*normals(:, panel)
                end do
            else
                charges(node) = &
                    charges(first_child(node)) + charges(second_child(node))
                dipoles(:, node) = dipoles(:, first_child(node)) + &
                    dipoles(:, second_child(node))
            end if
        end do
        interaction_count = 0
        result = ( &
            0.5_dp*areas - cmplx(0.0_dp, coupling, dp)*self_diagonal)*density
        do panel = 1, size(density)
            interaction_count = interaction_count + 1
            call apply_cfie_node( &
                1, panel, inverse_position(panel), wave_number, coupling, &
                opening_angle, centers, areas, normals, density, permutation, &
                first, last, first_child, second_child, node_centers, radii, &
                charges, dipoles, result(panel), interaction_count)
        end do
    end subroutine apply_cfie_prebuilt

    recursive subroutine apply_cfie_node( &
            node, target, target_position, wave_number, coupling, &
            opening_angle, centers, areas, normals, density, permutation, &
            first, last, first_child, second_child, node_centers, radii, &
            charges, dipoles, value, interaction_count)
        integer, intent(in) :: node, target, target_position
        real(dp), intent(in) :: wave_number, coupling, opening_angle
        real(dp), intent(in) :: centers(:, :), areas(:), normals(:, :)
        complex(dp), intent(in) :: density(:), charges(:), dipoles(:, :)
        integer, intent(in) :: permutation(:), first(:), last(:)
        integer, intent(in) :: first_child(:), second_child(:)
        real(dp), intent(in) :: node_centers(:, :), radii(:)
        complex(dp), intent(inout) :: value
        integer, intent(inout) :: interaction_count

        real(dp) :: displacement(3), distance
        complex(dp) :: green, normal_factor
        integer :: panel, position
        logical :: contains_target

        contains_target = target_position >= first(node) .and. &
            target_position <= last(node)
        displacement = centers(:, target) - node_centers(:, node)
        distance = norm2(displacement)
        if (.not. contains_target .and. distance > 0.0_dp .and. &
            radii(node)/distance < opening_angle .and. &
            wave_number*radii(node) < opening_angle) then
            green = outgoing_kernel(wave_number, distance)
            normal_factor = cmplx(1.0_dp, -wave_number*distance, dp)* &
                dot_product(displacement, dipoles(:, node))/distance**2
            value = value + areas(target)*green*( &
                normal_factor - cmplx(0.0_dp, coupling, dp)*charges(node))
            interaction_count = interaction_count + 1
            return
        end if
        if (first_child(node) == 0) then
            do position = first(node), last(node)
                panel = permutation(position)
                if (panel == target) cycle
                displacement = centers(:, target) - centers(:, panel)
                distance = norm2(displacement)
                green = outgoing_kernel(wave_number, distance)
                normal_factor = cmplx(1.0_dp, -wave_number*distance, dp)* &
                    dot_product(displacement, normals(:, panel))/distance**2
                value = value + areas(target)*areas(panel)*density(panel)* &
                    green*(normal_factor - cmplx(0.0_dp, coupling, dp))
                interaction_count = interaction_count + 1
            end do
            return
        end if
        call apply_cfie_node( &
            first_child(node), target, target_position, wave_number, coupling, &
            opening_angle, centers, areas, normals, density, permutation, &
            first, last, first_child, second_child, node_centers, radii, &
            charges, dipoles, value, interaction_count)
        call apply_cfie_node( &
            second_child(node), target, target_position, wave_number, coupling, &
            opening_angle, centers, areas, normals, density, permutation, &
            first, last, first_child, second_child, node_centers, radii, &
            charges, dipoles, value, interaction_count)
    end subroutine apply_cfie_node

    subroutine compute_panel_normals(vertices, triangles, normals, status)
        real(dp), intent(in) :: vertices(:, :)
        integer, intent(in) :: triangles(:, :)
        real(dp), allocatable, intent(out) :: normals(:, :)
        integer, intent(out) :: status

        real(dp) :: first(3), length, second(3)
        integer :: panel

        status = 1
        allocate(normals(3, size(triangles, 2)))
        do panel = 1, size(triangles, 2)
            first = vertices(:, triangles(2, panel)) - &
                vertices(:, triangles(1, panel))
            second = vertices(:, triangles(3, panel)) - &
                vertices(:, triangles(1, panel))
            normals(:, panel) = [ &
                first(2)*second(3) - first(3)*second(2), &
                first(3)*second(1) - first(1)*second(3), &
                first(1)*second(2) - first(2)*second(1)]
            length = norm2(normals(:, panel))
            if (length <= 0.0_dp) return
            normals(:, panel) = normals(:, panel)/length
        end do
        status = 0
    end subroutine compute_panel_normals

    pure logical function valid_geometry_and_parameters( &
            vertices, triangles, wave_number, opening_angle, leaf_size)
        real(dp), intent(in) :: vertices(:, :), wave_number, opening_angle
        integer, intent(in) :: triangles(:, :), leaf_size

        valid_geometry_and_parameters = .false.
        if (size(vertices, 1) /= 3 .or. size(vertices, 2) < 3) return
        if (size(triangles, 1) /= 3 .or. size(triangles, 2) < 1) return
        if (any(triangles < 1)) return
        if (any(triangles > size(vertices, 2))) return
        if (wave_number < 0.0_dp) return
        if (opening_angle <= 0.0_dp .or. opening_angle >= 1.0_dp) return
        if (leaf_size < 1) return
        valid_geometry_and_parameters = .true.
    end function valid_geometry_and_parameters

    recursive subroutine apply_node( &
            node, target, target_position, wave_number, opening_angle, &
            centers, areas, density, permutation, first, last, first_child, &
            second_child, node_centers, radii, charges, value, &
            interaction_count)
        integer, intent(in) :: node, target, target_position
        real(dp), intent(in) :: wave_number, opening_angle
        real(dp), intent(in) :: centers(:, :), areas(:)
        complex(dp), intent(in) :: density(:), charges(:)
        real(dp), intent(in) :: node_centers(:, :), radii(:)
        integer, intent(in) :: permutation(:), first(:), last(:)
        integer, intent(in) :: first_child(:), second_child(:)
        complex(dp), intent(inout) :: value
        integer, intent(inout) :: interaction_count

        real(dp) :: distance
        integer :: panel, position
        logical :: contains_target

        contains_target = target_position >= first(node) .and. &
            target_position <= last(node)
        distance = norm2(centers(:, target) - node_centers(:, node))
        if (.not. contains_target .and. distance > 0.0_dp .and. &
            radii(node)/distance < opening_angle .and. &
            wave_number*radii(node) < opening_angle) then
            value = value + areas(target)*charges(node)* &
                outgoing_kernel(wave_number, distance)
            interaction_count = interaction_count + 1
            return
        end if
        if (first_child(node) == 0) then
            do position = first(node), last(node)
                panel = permutation(position)
                if (panel == target) cycle
                distance = norm2(centers(:, target) - centers(:, panel))
                value = value + areas(target)*areas(panel)*density(panel)* &
                    outgoing_kernel(wave_number, distance)
                interaction_count = interaction_count + 1
            end do
            return
        end if
        call apply_node( &
            first_child(node), target, target_position, wave_number, &
            opening_angle, centers, areas, density, permutation, first, last, &
            first_child, second_child, node_centers, radii, charges, value, &
            interaction_count)
        call apply_node( &
            second_child(node), target, target_position, wave_number, &
            opening_angle, centers, areas, density, permutation, first, last, &
            first_child, second_child, node_centers, radii, charges, value, &
            interaction_count)
    end subroutine apply_node

    pure function outgoing_kernel(wave_number, distance) result(value)
        real(dp), intent(in) :: wave_number, distance
        complex(dp) :: value

        value = exp(cmplx(0.0_dp, wave_number*distance, dp))/ &
            (4.0_dp*acos(-1.0_dp)*distance)
    end function outgoing_kernel

end module fortfem_helmholtz_hierarchical_3d
