module fortfem_maxwell_sphere_curved_rwg
    !! Surface-Piola image of the affine RWG basis on a radial sphere panel.
    use fortfem_kinds, only: dp
    use fortfem_maxwell_rwg_surface, only: build_maxwell_rwg_surface_space
    use fortfem_sphere_curved_panel, only: evaluate_sphere_curved_panel
    use fortfem_triangle_duffy_quadrature, only: triangle_duffy_quadrature
    implicit none
    private

    public :: evaluate_maxwell_sphere_curved_rwg_basis
    public :: assemble_maxwell_sphere_curved_rwg_mass_matrix
    public :: assemble_maxwell_sphere_curved_plane_wave_rhs_rwg_3d
    public :: evaluate_maxwell_sphere_curved_far_field_rwg_3d

contains

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

        real(dp) :: edge_length, opposite_coordinates(2)
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
        call evaluate_sphere_curved_panel( &
            vertices(:, triangles(:, panel)), radius, xi, eta, point, &
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

end module fortfem_maxwell_sphere_curved_rwg
