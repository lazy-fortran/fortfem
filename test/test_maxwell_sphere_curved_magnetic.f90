program test_maxwell_sphere_curved_magnetic
    use check, only: check_condition, check_summary
    use fortfem_api, only: &
        build_maxwell_rwg_surface_space, &
        evaluate_maxwell_sphere_curved_magnetic_field_rwg_3d, &
        evaluate_maxwell_sphere_curved_rwg_basis, generate_sphere_surface_mesh
    use fortfem_kinds, only: dp
    use fortfem_triangle_duffy_quadrature, only: triangle_duffy_quadrature
    implicit none

    complex(dp), allocatable :: coefficients(:)
    complex(dp) :: curl_oracle(3), field(3), scaled_field(3)
    integer, allocatable :: edge_triangles(:, :), edge_vertices(:, :)
    integer, allocatable :: triangles(:, :)
    real(dp), allocatable :: scaled_vertices(:, :), vertices(:, :)
    real(dp) :: error, observation(3)
    integer :: basis, status
    logical :: all_passed

    all_passed = .true.
    observation = [2.2_dp, 0.3_dp, -0.1_dp]
    call generate_sphere_surface_mesh(1.0_dp, 0, vertices, triangles)
    call build_maxwell_rwg_surface_space( &
        vertices, triangles, edge_vertices, edge_triangles, status)
    allocate(coefficients(size(edge_vertices, 2)))
    do basis = 1, size(coefficients)
        coefficients(basis) = cmplx( &
            cos(real(basis, dp)), sin(real(2*basis, dp)), dp)
    end do
    call evaluate_maxwell_sphere_curved_magnetic_field_rwg_3d( &
        vertices, triangles, 1.0_dp, coefficients, observation, 0.7_dp, 12, &
        field, status)
    call finite_difference_curl( &
        vertices, triangles, edge_vertices, edge_triangles, coefficients, &
        observation, 1.0_dp, 0.7_dp, curl_oracle, status)
    error = sqrt(sum(abs(field - curl_oracle)**2))/ &
        sqrt(sum(abs(curl_oracle)**2))
    call record_condition(status == 0 .and. error < 2.0e-8_dp, &
        "curved magnetic field matches finite-difference curl oracle")

    call generate_sphere_surface_mesh(2.0_dp, 0, scaled_vertices, triangles)
    call evaluate_maxwell_sphere_curved_magnetic_field_rwg_3d( &
        scaled_vertices, triangles, 2.0_dp, coefficients, 2.0_dp*observation, &
        0.35_dp, 12, scaled_field, status)
    error = maxval(abs(scaled_field - field))/maxval(abs(field))
    call record_condition(status == 0 .and. error < 8.0e-13_dp, &
        "curved magnetic potential obeys electromagnetic similarity")

    call check_summary("Curved-sphere Maxwell magnetic potential")
    if (.not. all_passed) error stop 1

contains

    subroutine finite_difference_curl( &
            vertices, triangles, edge_vertices, edge_triangles, coefficients, &
            observation, radius, wave_number, curl_value, status)
        real(dp), intent(in) :: vertices(:, :), observation(3)
        real(dp), intent(in) :: radius, wave_number
        integer, intent(in) :: triangles(:, :), edge_vertices(:, :)
        integer, intent(in) :: edge_triangles(:, :)
        complex(dp), intent(in) :: coefficients(:)
        complex(dp), intent(out) :: curl_value(3)
        integer, intent(out) :: status

        complex(dp) :: derivative(3, 3), minus_value(3), plus_value(3)
        real(dp) :: displacement(3), step
        integer :: coordinate

        step = 2.0e-5_dp
        derivative = cmplx(0.0_dp, 0.0_dp, dp)
        do coordinate = 1, 3
            displacement = 0.0_dp
            displacement(coordinate) = step
            call vector_potential( &
                vertices, triangles, edge_vertices, edge_triangles, &
                coefficients, observation + displacement, radius, wave_number, &
                plus_value, status)
            if (status /= 0) return
            call vector_potential( &
                vertices, triangles, edge_vertices, edge_triangles, &
                coefficients, observation - displacement, radius, wave_number, &
                minus_value, status)
            if (status /= 0) return
            derivative(:, coordinate) = (plus_value - minus_value)/(2.0_dp*step)
        end do
        curl_value = [ &
            derivative(3, 2) - derivative(2, 3), &
            derivative(1, 3) - derivative(3, 1), &
            derivative(2, 1) - derivative(1, 2)]
        status = 0
    end subroutine finite_difference_curl

    subroutine vector_potential( &
            vertices, triangles, edge_vertices, edge_triangles, coefficients, &
            observation, radius, wave_number, value, status)
        real(dp), intent(in) :: vertices(:, :), observation(3)
        real(dp), intent(in) :: radius, wave_number
        integer, intent(in) :: triangles(:, :), edge_vertices(:, :)
        integer, intent(in) :: edge_triangles(:, :)
        complex(dp), intent(in) :: coefficients(:)
        complex(dp), intent(out) :: value(3)
        integer, intent(out) :: status

        real(dp), allocatable :: eta(:), weights(:), xi(:)
        real(dp) :: basis_value(3), divergence, jacobian, point(3), distance
        complex(dp) :: current(3), green
        integer :: basis, node, panel

        value = cmplx(0.0_dp, 0.0_dp, dp)
        call triangle_duffy_quadrature(12, xi, eta, weights, status)
        if (status /= 0) return
        do panel = 1, size(triangles, 2)
            do node = 1, size(weights)
                current = cmplx(0.0_dp, 0.0_dp, dp)
                do basis = 1, size(coefficients)
                    if (.not. any(edge_triangles(:, basis) == panel)) cycle
                    call evaluate_maxwell_sphere_curved_rwg_basis( &
                        vertices, triangles, edge_vertices, edge_triangles, &
                        basis, panel, radius, xi(node), eta(node), point, &
                        basis_value, divergence, jacobian, status)
                    if (status /= 0) return
                    current = current + coefficients(basis)*basis_value
                end do
                distance = norm2(observation - point)
                green = exp(cmplx(0.0_dp, wave_number*distance, dp))/ &
                    (4.0_dp*acos(-1.0_dp)*distance)
                value = value + weights(node)*jacobian*green*current
            end do
        end do
        status = 0
    end subroutine vector_potential

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_maxwell_sphere_curved_magnetic
