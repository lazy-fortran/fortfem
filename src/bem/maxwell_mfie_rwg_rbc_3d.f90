module fortfem_maxwell_mfie_rwg_rbc_3d
    !! RBC-tested Maxwell magnetic-field boundary operator.
    use fortfem_kinds, only: dp
    use fortfem_maxwell_bc_surface, only: &
        assemble_maxwell_rwg_rbc_pairing, build_maxwell_bc_transformation
    use fortfem_maxwell_localized_rwg_surface, only: &
        evaluate_maxwell_localized_rwg_basis
    use fortfem_maxwell_rwg_surface, only: &
        build_maxwell_rwg_surface_space, evaluate_maxwell_rwg_basis
    use fortfem_triangle_duffy_quadrature, only: triangle_duffy_quadrature
    implicit none
    private

    public :: assemble_maxwell_mfie_rwg_rbc_3d

contains

    subroutine assemble_maxwell_mfie_rwg_rbc_3d( &
            vertices, triangles, wave_number, quadrature_degree, tolerance, &
            max_depth, principal_value, matrix, status)
        real(dp), intent(in) :: vertices(:, :), wave_number, tolerance
        integer, intent(in) :: triangles(:, :), quadrature_degree, max_depth
        complex(dp), allocatable, intent(out) :: principal_value(:, :)
        complex(dp), allocatable, intent(out) :: matrix(:, :)
        integer, intent(out) :: status

        integer, allocatable :: edge_panels(:, :), edge_vertices(:, :)
        integer, allocatable :: refined_triangles(:, :)
        real(dp), allocatable :: eta(:), gram(:, :), refined_vertices(:, :)
        real(dp), allocatable :: transformation(:, :), weights(:), xi(:)
        complex(dp) :: contribution
        integer :: source_panel, test_basis, test_panel, trial_basis

        status = 1
        if (wave_number < 0.0_dp .or. tolerance <= 0.0_dp) return
        if (max_depth < 1) return
        call build_maxwell_bc_transformation( &
            vertices, triangles, refined_vertices, refined_triangles, &
            transformation, status)
        if (status /= 0) return
        call build_maxwell_rwg_surface_space( &
            vertices, triangles, edge_vertices, edge_panels, status)
        if (status /= 0) return
        call triangle_duffy_quadrature( &
            quadrature_degree, xi, eta, weights, status)
        if (status /= 0) return
        allocate( &
            principal_value(size(edge_vertices, 2), size(edge_vertices, 2)), &
            matrix(size(edge_vertices, 2), size(edge_vertices, 2)))
        principal_value = cmplx(0.0_dp, 0.0_dp, dp)
        do test_basis = 1, size(edge_vertices, 2)
            do trial_basis = 1, size(edge_vertices, 2)
                do test_panel = 1, size(refined_triangles, 2)
                    if (all(abs(transformation( &
                        3*test_panel - 2:3*test_panel, test_basis)) <= &
                        tiny(1.0_dp))) cycle
                    do source_panel = 1, 2
                        if ((test_panel - 1)/6 + 1 == &
                            edge_panels(source_panel, trial_basis)) cycle
                        call integrate_panel_pair( &
                            vertices, triangles, refined_vertices, &
                            refined_triangles, transformation, edge_vertices, &
                            edge_panels, test_basis, test_panel, trial_basis, &
                            edge_panels(source_panel, trial_basis), wave_number, &
                            xi, eta, weights, tolerance, max_depth, contribution, &
                            status)
                        if (status /= 0) return
                        principal_value(test_basis, trial_basis) = &
                            principal_value(test_basis, trial_basis) + &
                            contribution
                    end do
                end do
            end do
        end do
        call assemble_maxwell_rwg_rbc_pairing( &
            vertices, triangles, quadrature_degree, gram, status)
        if (status /= 0) return
        matrix = 0.5_dp*cmplx(gram, 0.0_dp, dp) - principal_value
        status = 0
    end subroutine assemble_maxwell_mfie_rwg_rbc_3d

    subroutine integrate_panel_pair( &
            vertices, triangles, refined_vertices, refined_triangles, &
            transformation, edge_vertices, edge_panels, test_basis, &
            test_panel, trial_basis, source_panel, wave_number, xi, eta, &
            weights, tolerance, max_depth, value, status)
        real(dp), intent(in) :: vertices(:, :), refined_vertices(:, :)
        real(dp), intent(in) :: transformation(:, :), wave_number
        real(dp), intent(in) :: xi(:), eta(:), weights(:), tolerance
        integer, intent(in) :: triangles(:, :), refined_triangles(:, :)
        integer, intent(in) :: edge_vertices(:, :), edge_panels(:, :)
        integer, intent(in) :: test_basis, test_panel, trial_basis, source_panel
        integer, intent(in) :: max_depth
        complex(dp), intent(out) :: value
        integer, intent(out) :: status

        real(dp) :: source_vertices(3, 3), test_vertices(3, 3)

        test_vertices = refined_vertices(:, refined_triangles(:, test_panel))
        source_vertices = vertices(:, triangles(:, source_panel))
        if (geometric_panels_touch(test_vertices, source_vertices)) then
            call integrate_adaptive_pair( &
                vertices, triangles, refined_vertices, refined_triangles, &
                transformation, edge_vertices, edge_panels, test_basis, &
                test_panel, trial_basis, source_panel, test_vertices, &
                source_vertices, wave_number, xi, eta, weights, tolerance, 0, &
                max_depth, value, status)
        else
            call integrate_regular_pair( &
                vertices, triangles, refined_vertices, refined_triangles, &
                transformation, edge_vertices, edge_panels, test_basis, &
                test_panel, trial_basis, source_panel, test_vertices, &
                source_vertices, wave_number, xi, eta, weights, value, status)
        end if
    end subroutine integrate_panel_pair

    recursive subroutine integrate_adaptive_pair( &
            vertices, triangles, refined_vertices, refined_triangles, &
            transformation, edge_vertices, edge_panels, test_basis, &
            test_panel, trial_basis, source_panel, test_vertices, &
            source_vertices, wave_number, xi, eta, weights, tolerance, depth, &
            max_depth, value, status)
        real(dp), intent(in) :: vertices(:, :), refined_vertices(:, :)
        real(dp), intent(in) :: transformation(:, :), test_vertices(3, 3)
        real(dp), intent(in) :: source_vertices(3, 3), wave_number
        real(dp), intent(in) :: xi(:), eta(:), weights(:), tolerance
        integer, intent(in) :: triangles(:, :), refined_triangles(:, :)
        integer, intent(in) :: edge_vertices(:, :), edge_panels(:, :)
        integer, intent(in) :: test_basis, test_panel, trial_basis, source_panel
        integer, intent(in) :: depth, max_depth
        complex(dp), intent(out) :: value
        integer, intent(out) :: status

        real(dp) :: source_children(3, 3, 4), test_children(3, 3, 4)
        complex(dp) :: coarse, contribution, refined
        integer :: source_child, test_child

        call integrate_regular_pair( &
            vertices, triangles, refined_vertices, refined_triangles, &
            transformation, edge_vertices, edge_panels, test_basis, test_panel, &
            trial_basis, source_panel, test_vertices, source_vertices, &
            wave_number, xi, eta, weights, coarse, status)
        if (status /= 0) return
        call subdivide_triangle(test_vertices, test_children)
        call subdivide_triangle(source_vertices, source_children)
        refined = cmplx(0.0_dp, 0.0_dp, dp)
        do test_child = 1, 4
            do source_child = 1, 4
                call integrate_regular_pair( &
                    vertices, triangles, refined_vertices, refined_triangles, &
                    transformation, edge_vertices, edge_panels, test_basis, &
                    test_panel, trial_basis, source_panel, &
                    test_children(:, :, test_child), &
                    source_children(:, :, source_child), wave_number, xi, eta, &
                    weights, contribution, status)
                if (status /= 0) return
                refined = refined + contribution
            end do
        end do
        if (depth + 1 >= max_depth .or. abs(refined - coarse) <= &
            tolerance*max(tiny(1.0_dp), abs(refined))) then
            value = refined
            status = 0
            return
        end if
        value = cmplx(0.0_dp, 0.0_dp, dp)
        do test_child = 1, 4
            do source_child = 1, 4
                if (geometric_panels_touch( &
                    test_children(:, :, test_child), &
                    source_children(:, :, source_child))) then
                    call integrate_adaptive_pair( &
                        vertices, triangles, refined_vertices, &
                        refined_triangles, transformation, edge_vertices, &
                        edge_panels, test_basis, test_panel, trial_basis, &
                        source_panel, test_children(:, :, test_child), &
                        source_children(:, :, source_child), wave_number, xi, &
                        eta, weights, tolerance, depth + 1, max_depth, &
                        contribution, status)
                else
                    call integrate_regular_pair( &
                        vertices, triangles, refined_vertices, &
                        refined_triangles, transformation, edge_vertices, &
                        edge_panels, test_basis, test_panel, trial_basis, &
                        source_panel, test_children(:, :, test_child), &
                        source_children(:, :, source_child), wave_number, xi, &
                        eta, weights, contribution, status)
                end if
                if (status /= 0) return
                value = value + contribution
            end do
        end do
        status = 0
    end subroutine integrate_adaptive_pair

    subroutine integrate_regular_pair( &
            vertices, triangles, refined_vertices, refined_triangles, &
            transformation, edge_vertices, edge_panels, test_basis, &
            test_panel, trial_basis, source_panel, test_vertices, &
            source_vertices, wave_number, xi, eta, weights, value, status)
        real(dp), intent(in) :: vertices(:, :), refined_vertices(:, :)
        real(dp), intent(in) :: transformation(:, :), test_vertices(3, 3)
        real(dp), intent(in) :: source_vertices(3, 3), wave_number
        real(dp), intent(in) :: xi(:), eta(:), weights(:)
        integer, intent(in) :: triangles(:, :), refined_triangles(:, :)
        integer, intent(in) :: edge_vertices(:, :), edge_panels(:, :)
        integer, intent(in) :: test_basis, test_panel, trial_basis, source_panel
        complex(dp), intent(out) :: value
        integer, intent(out) :: status

        real(dp) :: bc_value(3), divergence, displacement(3), local_value(3)
        real(dp) :: source_jacobian, source_point(3), source_value(3)
        real(dp) :: test_jacobian, test_point(3), radius
        complex(dp) :: gradient_green(3), green, magnetic_field(3)
        integer :: local_edge, source_node, test_node

        associate(unused => refined_triangles)
        end associate
        value = cmplx(0.0_dp, 0.0_dp, dp)
        test_jacobian = 2.0_dp*triangle_area(test_vertices)
        source_jacobian = 2.0_dp*triangle_area(source_vertices)
        do test_node = 1, size(weights)
            test_point = triangle_point( &
                test_vertices, xi(test_node), eta(test_node))
            bc_value = 0.0_dp
            do local_edge = 1, 3
                call evaluate_maxwell_localized_rwg_basis( &
                    refined_vertices, refined_triangles, test_panel, &
                    local_edge, test_point, local_value, divergence, status)
                if (status /= 0) return
                bc_value = bc_value + transformation( &
                    3*(test_panel - 1) + local_edge, test_basis)*local_value
            end do
            do source_node = 1, size(weights)
                source_point = triangle_point( &
                    source_vertices, xi(source_node), eta(source_node))
                call evaluate_maxwell_rwg_basis( &
                    vertices, triangles, edge_vertices, edge_panels, &
                    trial_basis, source_panel, source_point, source_value, &
                    divergence, status)
                if (status /= 0) return
                displacement = test_point - source_point
                radius = norm2(displacement)
                if (radius <= tiny(1.0_dp)) cycle
                green = exp(cmplx(0.0_dp, wave_number*radius, dp))/ &
                    (4.0_dp*acos(-1.0_dp)*radius)
                gradient_green = green* &
                    (cmplx(0.0_dp, wave_number, dp) - 1.0_dp/radius)* &
                    displacement/radius
                magnetic_field = cross_product(gradient_green, source_value)
                value = value - test_jacobian*source_jacobian* &
                    weights(test_node)*weights(source_node)* &
                    sum(bc_value*magnetic_field)
            end do
        end do
        status = 0
    end subroutine integrate_regular_pair

    pure function triangle_point(vertices, xi, eta) result(point)
        real(dp), intent(in) :: vertices(3, 3), xi, eta
        real(dp) :: point(3)
        point = vertices(:, 1) + xi*(vertices(:, 2) - vertices(:, 1)) + &
            eta*(vertices(:, 3) - vertices(:, 1))
    end function triangle_point

    pure function triangle_area(vertices) result(area)
        real(dp), intent(in) :: vertices(3, 3)
        real(dp) :: area
        area = 0.5_dp*norm2(cross_product_real( &
            vertices(:, 2) - vertices(:, 1), &
            vertices(:, 3) - vertices(:, 1)))
    end function triangle_area

    pure logical function geometric_panels_touch(first, second) result(touch)
        real(dp), intent(in) :: first(3, 3), second(3, 3)
        real(dp) :: scale
        integer :: first_vertex, second_vertex

        scale = max(1.0_dp, maxval(abs(first)), maxval(abs(second)))
        touch = .false.
        do first_vertex = 1, 3
            do second_vertex = 1, 3
                if (norm2(first(:, first_vertex) - second(:, second_vertex)) <= &
                    128.0_dp*epsilon(1.0_dp)*scale) then
                    touch = .true.
                    return
                end if
            end do
        end do
    end function geometric_panels_touch

    pure subroutine subdivide_triangle(vertices, children)
        real(dp), intent(in) :: vertices(3, 3)
        real(dp), intent(out) :: children(3, 3, 4)
        real(dp) :: midpoint_12(3), midpoint_23(3), midpoint_31(3)

        midpoint_12 = 0.5_dp*(vertices(:, 1) + vertices(:, 2))
        midpoint_23 = 0.5_dp*(vertices(:, 2) + vertices(:, 3))
        midpoint_31 = 0.5_dp*(vertices(:, 3) + vertices(:, 1))
        children(:, :, 1) = reshape( &
            [vertices(:, 1), midpoint_12, midpoint_31], [3, 3])
        children(:, :, 2) = reshape( &
            [midpoint_12, vertices(:, 2), midpoint_23], [3, 3])
        children(:, :, 3) = reshape( &
            [midpoint_31, midpoint_23, vertices(:, 3)], [3, 3])
        children(:, :, 4) = reshape( &
            [midpoint_12, midpoint_23, midpoint_31], [3, 3])
    end subroutine subdivide_triangle

    pure function cross_product(first, second) result(product)
        complex(dp), intent(in) :: first(3)
        real(dp), intent(in) :: second(3)
        complex(dp) :: product(3)
        product = [ &
            first(2)*second(3) - first(3)*second(2), &
            first(3)*second(1) - first(1)*second(3), &
            first(1)*second(2) - first(2)*second(1)]
    end function cross_product

    pure function cross_product_real(first, second) result(product)
        real(dp), intent(in) :: first(3), second(3)
        real(dp) :: product(3)
        product = [ &
            first(2)*second(3) - first(3)*second(2), &
            first(3)*second(1) - first(1)*second(3), &
            first(1)*second(2) - first(2)*second(1)]
    end function cross_product_real

end module fortfem_maxwell_mfie_rwg_rbc_3d
