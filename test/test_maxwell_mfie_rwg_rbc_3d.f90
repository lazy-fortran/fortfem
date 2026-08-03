program test_maxwell_mfie_rwg_rbc_3d
    use check, only: check_condition, check_summary
    use fortfem_boundary, only: &
        assemble_maxwell_mfie_rwg_rbc_3d, &
        evaluate_maxwell_magnetic_field_rwg_3d
    use fortfem_feec, only: &
        build_maxwell_bc_transformation, &
        evaluate_maxwell_localized_rwg_basis, triangle_duffy_quadrature
    use fortfem_kinds, only: dp
    implicit none

    complex(dp), allocatable :: matrix(:, :), principal_depth1(:, :)
    complex(dp), allocatable :: principal_depth2(:, :), principal_depth3(:, :)
    complex(dp) :: exterior_coarse(6, 6), exterior_fine(6, 6)
    complex(dp) :: exterior_limit(6, 6)
    real(dp) :: error12, error23, trace_error, vertices(3, 4)
    integer :: status, triangles(3, 4)
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

    call assemble_maxwell_mfie_rwg_rbc_3d( &
        vertices, triangles, 1.3_dp, 2, 1.0e-14_dp, 1, principal_depth1, &
        matrix, status)
    call assemble_maxwell_mfie_rwg_rbc_3d( &
        vertices, triangles, 1.3_dp, 2, 1.0e-14_dp, 2, principal_depth2, &
        matrix, status)
    call assemble_maxwell_mfie_rwg_rbc_3d( &
        vertices, triangles, 1.3_dp, 2, 1.0e-14_dp, 3, principal_depth3, &
        matrix, status)
    error12 = sqrt(sum(abs(principal_depth2 - principal_depth1)**2))
    error23 = sqrt(sum(abs(principal_depth3 - principal_depth2)**2))
    call record_condition(status == 0 .and. error23 < 0.65_dp*error12, &
        "adaptive RWG-RBC principal value converges under touching refinement")
    call record_condition(maxval(abs( &
        principal_depth3 - transpose(principal_depth3))) > 1.0e-5_dp, &
        "Maxwell magnetic principal value retains its nonsymmetric dual action")
    call record_condition(all(abs(matrix) < huge(1.0_dp)), &
        "MFIE jump-plus-principal-value matrix remains finite")
    call assemble_exterior_trace(0.12_dp, exterior_coarse, status)
    call assemble_exterior_trace(0.06_dp, exterior_fine, status)
    exterior_limit = 2.0_dp*exterior_fine - exterior_coarse
    trace_error = sqrt(sum(abs(exterior_limit + matrix)**2))/ &
        sqrt(sum(abs(matrix)**2))
    call record_condition(trace_error < 0.25_dp, &
        "MFIE sign and jump agree with an extrapolated exterior trace")

    call check_summary("Three-dimensional Maxwell RWG-RBC MFIE")
    if (.not. all_passed) error stop 1

contains

    subroutine assemble_exterior_trace(offset, trace_matrix, local_status)
        real(dp), intent(in) :: offset
        complex(dp), intent(out) :: trace_matrix(6, 6)
        integer, intent(out) :: local_status

        complex(dp) :: coefficients(6), magnetic_field(3)
        integer, allocatable :: refined_triangles(:, :)
        real(dp), allocatable :: eta(:), refined_vertices(:, :)
        real(dp), allocatable :: transformation(:, :), weights(:), xi(:)
        real(dp) :: bc_value(3), divergence, jacobian, local_value(3)
        real(dp) :: normal(3), observation(3), point(3)
        integer :: local_edge, node, refined_panel, row, test_basis, trial_basis

        trace_matrix = cmplx(0.0_dp, 0.0_dp, dp)
        call build_maxwell_bc_transformation( &
            vertices, triangles, refined_vertices, refined_triangles, &
            transformation, local_status)
        if (local_status /= 0) return
        call triangle_duffy_quadrature( &
            3, xi, eta, weights, local_status)
        if (local_status /= 0) return
        do refined_panel = 1, size(refined_triangles, 2)
            normal = triangle_normal( &
                refined_vertices(:, refined_triangles(:, refined_panel)))
            jacobian = 2.0_dp*triangle_area( &
                refined_vertices(:, refined_triangles(:, refined_panel)))
            do node = 1, size(weights)
                point = triangle_point( &
                    refined_vertices(:, refined_triangles(:, refined_panel)), &
                    xi(node), eta(node))
                observation = point + offset*normal
                do trial_basis = 1, 6
                    coefficients = cmplx(0.0_dp, 0.0_dp, dp)
                    coefficients(trial_basis) = cmplx(1.0_dp, 0.0_dp, dp)
                    call evaluate_maxwell_magnetic_field_rwg_3d( &
                        vertices, triangles, coefficients, observation, &
                        1.3_dp, 10, magnetic_field, local_status)
                    if (local_status /= 0) return
                    do test_basis = 1, 6
                        bc_value = 0.0_dp
                        do local_edge = 1, 3
                            call evaluate_maxwell_localized_rwg_basis( &
                                refined_vertices, refined_triangles, &
                                refined_panel, local_edge, point, local_value, &
                                divergence, local_status)
                            if (local_status /= 0) return
                            row = 3*(refined_panel - 1) + local_edge
                            bc_value = bc_value + &
                                transformation(row, test_basis)*local_value
                        end do
                        trace_matrix(test_basis, trial_basis) = &
                            trace_matrix(test_basis, trial_basis) - &
                            jacobian*weights(node)*sum(bc_value*magnetic_field)
                    end do
                end do
            end do
        end do
    end subroutine assemble_exterior_trace

    pure function triangle_point(points, xi, eta) result(point)
        real(dp), intent(in) :: points(3, 3), xi, eta
        real(dp) :: point(3)

        point = points(:, 1) + xi*(points(:, 2) - points(:, 1)) + &
            eta*(points(:, 3) - points(:, 1))
    end function triangle_point

    pure function triangle_normal(points) result(normal)
        real(dp), intent(in) :: points(3, 3)
        real(dp) :: normal(3)

        normal = cross_product( &
            points(:, 2) - points(:, 1), points(:, 3) - points(:, 1))
        normal = normal/norm2(normal)
    end function triangle_normal

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

end program test_maxwell_mfie_rwg_rbc_3d
