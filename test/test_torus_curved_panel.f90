program test_torus_curved_panel
    use check, only: check_condition, check_summary
    use fortfem_core, only: &
        evaluate_torus_curved_panel, generate_torus_surface_mesh
    use fortfem_kinds, only: dp
    use fortfem_feec, only: triangle_duffy_quadrature
    implicit none

    real(dp), parameter :: major_radius = 2.0_dp, minor_radius = 0.6_dp
    integer, allocatable :: triangles(:, :)
    real(dp), allocatable :: eta(:), parameters(:, :), vertices(:, :)
    real(dp), allocatable :: weights(:), xi(:)
    real(dp) :: areas(3), errors(3), jacobian, normal(3), point(3)
    real(dp) :: tangent_eta(3), tangent_xi(3)
    integer :: degree, degrees(3), node, panel, status
    logical :: all_passed, geometry_ok

    all_passed = .true.
    geometry_ok = .true.
    degrees = [2, 8, 16]
    call generate_torus_surface_mesh( &
        major_radius, minor_radius, 5, 6, vertices, triangles, parameters)
    do degree = 1, 3
        call triangle_duffy_quadrature( &
            degrees(degree), xi, eta, weights, status)
        areas(degree) = 0.0_dp
        do panel = 1, size(triangles, 2)
            do node = 1, size(weights)
                call evaluate_torus_curved_panel( &
                    parameters(:, triangles(:, panel)), major_radius, &
                    minor_radius, xi(node), eta(node), point, tangent_xi, &
                    tangent_eta, jacobian, status)
                normal = torus_normal(point)
                areas(degree) = areas(degree) + weights(node)*jacobian
                geometry_ok = geometry_ok .and. status == 0 .and. &
                    abs(torus_level_set(point)) < 2.0e-14_dp .and. &
                    abs(dot_product(normal, tangent_xi)) < 3.0e-14_dp .and. &
                    abs(dot_product(normal, tangent_eta)) < 3.0e-14_dp
            end do
        end do
    end do
    errors = abs(areas - 4.0_dp*acos(-1.0_dp)**2* &
        major_radius*minor_radius)
    call record_condition(geometry_ok, &
        "curved torus panels lie on the exact surface with tangent derivatives")
    call record_condition(maxval(errors) < 2.0e-12_dp, &
        "curved panels recover the analytical torus area")
    call check_summary("Exact curved torus panel geometry")
    if (.not. all_passed) error stop 1

contains

    pure function torus_level_set(position) result(value)
        real(dp), intent(in) :: position(3)
        real(dp) :: value

        value = (sqrt(position(1)**2 + position(2)**2) - major_radius)**2 + &
            position(3)**2 - minor_radius**2
    end function torus_level_set

    pure function torus_normal(position) result(normal)
        real(dp), intent(in) :: position(3)
        real(dp) :: normal(3)
        real(dp) :: cylindrical_radius

        cylindrical_radius = sqrt(position(1)**2 + position(2)**2)
        normal = [ &
            (cylindrical_radius - major_radius)*position(1)/ &
            cylindrical_radius, &
            (cylindrical_radius - major_radius)*position(2)/ &
            cylindrical_radius, position(3)]/minor_radius
    end function torus_normal

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_torus_curved_panel
