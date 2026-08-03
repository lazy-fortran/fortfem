program test_maxwell_magnetic_rwg_3d
    use check, only: check_condition, check_summary
    use fortfem_boundary, only: evaluate_maxwell_magnetic_field_rwg_3d
    use fortfem_feec, only: &
        build_maxwell_rwg_surface_space, &
        evaluate_maxwell_rwg_basis
    use fortfem_kinds, only: dp
    use fortfem_triangle_duffy_quadrature, only: triangle_duffy_quadrature
    implicit none

    complex(dp) :: coefficients(6), finite_difference_curl(3)
    complex(dp) :: magnetic_field(3), minus_potential(3), plus_potential(3)
    real(dp) :: observation(3), step, vertices(3, 4)
    integer :: axis, status, triangles(3, 4)
    logical :: all_passed

    all_passed = .true.
    vertices(:, 1) = [0.0_dp, 0.0_dp, 0.0_dp]
    vertices(:, 2) = [1.0_dp, 0.0_dp, 0.0_dp]
    vertices(:, 3) = [0.0_dp, 1.0_dp, 0.0_dp]
    vertices(:, 4) = [0.0_dp, 0.0_dp, 1.0_dp]
    triangles(:, 1) = [1, 3, 2]
    triangles(:, 2) = [1, 2, 4]
    triangles(:, 3) = [1, 4, 3]
    triangles(:, 4) = [2, 3, 4]
    coefficients = [ &
        cmplx(0.3_dp, -0.1_dp, dp), cmplx(-0.2_dp, 0.4_dp, dp), &
        cmplx(0.7_dp, 0.2_dp, dp), cmplx(-0.5_dp, -0.3_dp, dp), &
        cmplx(0.1_dp, 0.6_dp, dp), cmplx(-0.4_dp, 0.2_dp, dp)]
    observation = [1.7_dp, -0.8_dp, 1.3_dp]

    call evaluate_maxwell_magnetic_field_rwg_3d( &
        vertices, triangles, coefficients, observation, 1.4_dp, 6, &
        magnetic_field, status)
    step = 2.0e-5_dp
    finite_difference_curl = cmplx(0.0_dp, 0.0_dp, dp)
    do axis = 1, 3
        observation(axis) = observation(axis) + step
        call vector_potential(observation, plus_potential, status)
        observation(axis) = observation(axis) - 2.0_dp*step
        call vector_potential(observation, minus_potential, status)
        observation(axis) = observation(axis) + step
        select case (axis)
        case (1)
            finite_difference_curl(2) = finite_difference_curl(2) - &
                (plus_potential(3) - minus_potential(3))/(2.0_dp*step)
            finite_difference_curl(3) = finite_difference_curl(3) + &
                (plus_potential(2) - minus_potential(2))/(2.0_dp*step)
        case (2)
            finite_difference_curl(1) = finite_difference_curl(1) + &
                (plus_potential(3) - minus_potential(3))/(2.0_dp*step)
            finite_difference_curl(3) = finite_difference_curl(3) - &
                (plus_potential(1) - minus_potential(1))/(2.0_dp*step)
        case (3)
            finite_difference_curl(1) = finite_difference_curl(1) - &
                (plus_potential(2) - minus_potential(2))/(2.0_dp*step)
            finite_difference_curl(2) = finite_difference_curl(2) + &
                (plus_potential(1) - minus_potential(1))/(2.0_dp*step)
        end select
    end do
    call record_condition(status == 0 .and. &
        sqrt(sum(abs(magnetic_field - finite_difference_curl)**2)) < 2.0e-10_dp, &
        "RWG magnetic representation equals an independent curl of potential")

    call check_summary("Three-dimensional Maxwell magnetic representation")
    if (.not. all_passed) error stop 1

contains

    subroutine vector_potential(point, potential, local_status)
        real(dp), intent(in) :: point(3)
        complex(dp), intent(out) :: potential(3)
        integer, intent(out) :: local_status

        integer, allocatable :: edge_panels(:, :), edge_vertices(:, :)
        real(dp), allocatable :: eta(:), weights(:), xi(:)
        real(dp) :: basis_divergence, basis_value(3), jacobian, radius
        real(dp) :: source(3)
        complex(dp) :: current(3), green
        integer :: basis, node, panel

        potential = cmplx(0.0_dp, 0.0_dp, dp)
        call build_maxwell_rwg_surface_space( &
            vertices, triangles, edge_vertices, edge_panels, local_status)
        if (local_status /= 0) return
        call triangle_duffy_quadrature(6, xi, eta, weights, local_status)
        if (local_status /= 0) return
        do panel = 1, 4
            jacobian = 2.0_dp*triangle_area( &
                vertices(:, triangles(:, panel)))
            do node = 1, size(weights)
                source = triangle_point( &
                    vertices(:, triangles(:, panel)), xi(node), eta(node))
                current = cmplx(0.0_dp, 0.0_dp, dp)
                do basis = 1, 6
                    if (.not. any(edge_panels(:, basis) == panel)) cycle
                    call evaluate_maxwell_rwg_basis( &
                        vertices, triangles, edge_vertices, edge_panels, &
                        basis, panel, source, basis_value, basis_divergence, &
                        local_status)
                    if (local_status /= 0) return
                    current = current + coefficients(basis)*basis_value
                end do
                radius = norm2(point - source)
                green = exp(cmplx(0.0_dp, 1.4_dp*radius, dp))/ &
                    (4.0_dp*acos(-1.0_dp)*radius)
                potential = potential + jacobian*weights(node)*green*current
            end do
        end do
    end subroutine vector_potential

    pure function triangle_point(points, xi, eta) result(point)
        real(dp), intent(in) :: points(3, 3), xi, eta
        real(dp) :: point(3)

        point = points(:, 1) + xi*(points(:, 2) - points(:, 1)) + &
            eta*(points(:, 3) - points(:, 1))
    end function triangle_point

    pure function triangle_area(points) result(area)
        real(dp), intent(in) :: points(3, 3)
        real(dp) :: area

        area = 0.5_dp*norm2(cross_product( &
            points(:, 2) - points(:, 1), points(:, 3) - points(:, 1)))
    end function triangle_area

    pure function cross_product(first, second) result(product)
        real(dp), intent(in) :: first(3), second(3)
        real(dp) :: product(3)

        product = [ &
            first(2)*second(3) - first(3)*second(2), &
            first(3)*second(1) - first(1)*second(3), &
            first(1)*second(2) - first(2)*second(1)]
    end function cross_product

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_maxwell_magnetic_rwg_3d
