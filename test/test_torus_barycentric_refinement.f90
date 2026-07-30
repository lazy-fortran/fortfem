program test_torus_barycentric_refinement
    use check, only: check_condition, check_summary
    use fortfem_api, only: &
        barycentric_refine_torus_surface_mesh, evaluate_torus_curved_panel, &
        generate_torus_surface_mesh
    use fortfem_kinds, only: dp
    use fortfem_triangle_duffy_quadrature, only: triangle_duffy_quadrature
    implicit none

    real(dp), parameter :: major_radius = 2.0_dp, minor_radius = 0.6_dp
    integer, allocatable :: refined_triangles(:, :), triangles(:, :)
    real(dp), allocatable :: eta(:), parameters(:, :), refined_parameters(:, :)
    real(dp), allocatable :: refined_vertices(:, :), vertices(:, :)
    real(dp), allocatable :: weights(:), xi(:)
    real(dp) :: area, area_error, jacobian, level_error, point(3)
    real(dp) :: tangent_eta(3), tangent_xi(3)
    integer :: expected_edges, node, panel, status
    logical :: all_passed

    all_passed = .true.
    call generate_torus_surface_mesh( &
        major_radius, minor_radius, 5, 6, vertices, triangles, parameters)
    call barycentric_refine_torus_surface_mesh( &
        vertices, triangles, parameters, major_radius, minor_radius, &
        refined_vertices, refined_triangles, refined_parameters, status)
    expected_edges = 3*size(triangles, 2)/2
    call record_condition(status == 0 .and. &
        size(refined_vertices, 2) == size(vertices, 2) + expected_edges + &
        size(triangles, 2) .and. &
        size(refined_triangles, 2) == 6*size(triangles, 2), &
        "toroidal barycentric refinement has the exact closed-mesh topology")
    level_error = 0.0_dp
    do node = 1, size(refined_vertices, 2)
        level_error = max(level_error, abs(torus_level(refined_vertices(:, node))))
    end do
    call record_condition(level_error < 3.0e-14_dp, &
        "all refined vertices lie on the exact torus")

    call triangle_duffy_quadrature(12, xi, eta, weights, status)
    area = 0.0_dp
    do panel = 1, size(refined_triangles, 2)
        do node = 1, size(weights)
            call evaluate_torus_curved_panel( &
                refined_parameters(:, refined_triangles(:, panel)), &
                major_radius, minor_radius, xi(node), eta(node), point, &
                tangent_xi, tangent_eta, jacobian, status)
            area = area + weights(node)*jacobian
        end do
    end do
    area_error = abs(area - 4.0_dp*acos(-1.0_dp)**2* &
        major_radius*minor_radius)
    call record_condition(status == 0 .and. area_error < 3.0e-12_dp, &
        "refined exact panels recover the analytical torus area")

    call check_summary("Exact-torus barycentric refinement")
    if (.not. all_passed) error stop 1

contains

    pure function torus_level(position) result(value)
        real(dp), intent(in) :: position(3)
        real(dp) :: value

        value = (sqrt(position(1)**2 + position(2)**2) - major_radius)**2 + &
            position(3)**2 - minor_radius**2
    end function torus_level

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_torus_barycentric_refinement
