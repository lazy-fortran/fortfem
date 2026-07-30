program test_maxwell_sphere_curved_rwg
    use check, only: check_condition, check_summary
    use fortfem_api, only: &
        build_maxwell_rwg_surface_space, &
        evaluate_maxwell_sphere_curved_localized_rwg_basis, &
        evaluate_maxwell_sphere_curved_rwg_basis, &
        generate_sphere_surface_mesh
    use fortfem_kinds, only: dp
    use fortfem_triangle_duffy_quadrature, only: triangle_duffy_quadrature
    implicit none

    integer, allocatable :: edge_triangles(:, :), edge_vertices(:, :)
    integer, allocatable :: triangles(:, :)
    real(dp), allocatable :: eta(:), vertices(:, :), weights(:), xi(:)
    real(dp) :: divergence, divergence_integral, edge_length, jacobian
    real(dp) :: local_integral
    real(dp) :: point(3), radius, value(3)
    integer :: basis, local, local_edge, node, orientation, panel, status
    logical :: all_passed, geometry_ok

    all_passed = .true.
    geometry_ok = .true.
    radius = 1.7_dp
    call generate_sphere_surface_mesh(radius, 0, vertices, triangles)
    call build_maxwell_rwg_surface_space( &
        vertices, triangles, edge_vertices, edge_triangles, status)
    call record_condition(status == 0, "closed sphere has an RWG trace space")
    basis = 1
    panel = edge_triangles(1, basis)
    orientation = 0
    do local = 1, 3
        if (triangles(local, panel) == edge_vertices(1, basis) .and. &
            triangles(modulo(local, 3) + 1, panel) == &
            edge_vertices(2, basis)) orientation = 1
        if (triangles(local, panel) == edge_vertices(2, basis) .and. &
            triangles(modulo(local, 3) + 1, panel) == &
            edge_vertices(1, basis)) orientation = -1
    end do
    edge_length = norm2( &
        vertices(:, edge_vertices(2, basis)) - &
        vertices(:, edge_vertices(1, basis)))
    call triangle_duffy_quadrature(16, xi, eta, weights, status)
    divergence_integral = 0.0_dp
    do node = 1, size(weights)
        call evaluate_maxwell_sphere_curved_rwg_basis( &
            vertices, triangles, edge_vertices, edge_triangles, basis, panel, &
            radius, xi(node), eta(node), point, value, divergence, jacobian, &
            status)
        geometry_ok = geometry_ok .and. status == 0 .and. &
            abs(norm2(point) - radius) < 8.0e-15_dp .and. &
            abs(dot_product(point, value)) < 3.0e-14_dp
        divergence_integral = divergence_integral + &
            weights(node)*jacobian*divergence
    end do
    call record_condition(geometry_ok, &
        "curved RWG values are tangent to the exact sphere")
    call record_condition(abs(divergence_integral - &
        real(orientation, dp)*edge_length) < 2.0e-13_dp, &
        "surface Piola map preserves the analytical RWG edge flux")
    do local_edge = 1, 3
        local_integral = 0.0_dp
        do node = 1, size(weights)
            call evaluate_maxwell_sphere_curved_localized_rwg_basis( &
                vertices, triangles, panel, local_edge, radius, xi(node), &
                eta(node), point, value, divergence, jacobian, status)
            geometry_ok = geometry_ok .and. status == 0 .and. &
                abs(dot_product(point, value)) < 3.0e-14_dp
            local_integral = local_integral + &
                weights(node)*jacobian*divergence
        end do
        select case (local_edge)
        case (1)
            edge_length = norm2(vertices(:, triangles(2, panel)) - &
                vertices(:, triangles(1, panel)))
        case (2)
            edge_length = norm2(vertices(:, triangles(1, panel)) - &
                vertices(:, triangles(3, panel)))
        case (3)
            edge_length = norm2(vertices(:, triangles(3, panel)) - &
                vertices(:, triangles(2, panel)))
        end select
        geometry_ok = geometry_ok .and. &
            abs(local_integral - edge_length) < 2.0e-13_dp
    end do
    call record_condition(geometry_ok, &
        "all curved localized RWG functions preserve positive edge flux")

    call check_summary("Curved-sphere RWG surface Piola map")
    if (.not. all_passed) error stop 1

contains

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_maxwell_sphere_curved_rwg
