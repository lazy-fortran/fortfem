program test_sphere_curved_panel
    use check, only: check_condition, check_summary
    use fortfem_api, only: &
        evaluate_sphere_curved_panel, generate_sphere_surface_mesh
    use fortfem_kinds, only: dp
    use fortfem_triangle_duffy_quadrature, only: triangle_duffy_quadrature
    implicit none

    integer, allocatable :: triangles(:, :)
    real(dp), allocatable :: eta(:), vertices(:, :), weights(:), xi(:)
    real(dp) :: areas(3), errors(3), jacobian, point(3), radius
    real(dp) :: tangent_eta(3), tangent_xi(3)
    integer :: degree, degrees(3), node, panel, status
    logical :: all_passed, geometry_ok

    all_passed = .true.
    geometry_ok = .true.
    radius = 2.3_dp
    degrees = [2, 8, 16]
    call generate_sphere_surface_mesh(radius, 0, vertices, triangles)
    do degree = 1, 3
        call triangle_duffy_quadrature( &
            degrees(degree), xi, eta, weights, status)
        areas(degree) = 0.0_dp
        do panel = 1, size(triangles, 2)
            do node = 1, size(weights)
                call evaluate_sphere_curved_panel( &
                    vertices(:, triangles(:, panel)), radius, xi(node), &
                    eta(node), point, tangent_xi, tangent_eta, jacobian, status)
                areas(degree) = areas(degree) + weights(node)*jacobian
                geometry_ok = geometry_ok .and. status == 0 .and. &
                    abs(norm2(point) - radius) < 8.0e-15_dp .and. &
                    abs(dot_product(point, tangent_xi)) < 2.0e-14_dp .and. &
                    abs(dot_product(point, tangent_eta)) < 2.0e-14_dp
            end do
        end do
    end do
    errors = abs(areas - 4.0_dp*acos(-1.0_dp)*radius**2)
    call record_condition(geometry_ok, &
        "radial sphere panel has exact radius and tangent plane")
    if (errors(2) >= 0.02_dp*errors(1) .or. &
        errors(3) >= 0.1_dp*errors(2)) write (*, *) "sphere area errors", errors
    call record_condition(errors(2) < 0.02_dp*errors(1) .and. &
        errors(3) < 0.1_dp*errors(2), &
        "curved-panel sphere area converges rapidly under quadrature refinement")
    call record_condition(errors(3)/(4.0_dp*acos(-1.0_dp)*radius**2) < &
        1.0e-6_dp, &
        "curved octahedral panels recover the analytical sphere area")

    call check_summary("Exact curved sphere panel geometry")
    if (.not. all_passed) error stop 1

contains

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_sphere_curved_panel
