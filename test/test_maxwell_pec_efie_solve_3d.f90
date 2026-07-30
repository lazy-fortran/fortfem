program test_maxwell_pec_efie_solve_3d
    use check, only: check_condition, check_summary
    use fortfem_api, only: &
        assemble_maxwell_efie_rwg_3d, assemble_maxwell_plane_wave_rhs_rwg_3d, &
        build_maxwell_rwg_surface_space, evaluate_maxwell_rwg_basis, &
        evaluate_maxwell_efie_far_field_rwg_3d, &
        evaluate_maxwell_efie_field_rwg_3d, solve_maxwell_pec_efie_rwg_3d
    use fortfem_kinds, only: dp
    implicit none

    complex(dp), allocatable :: density(:), matrix(:, :), right_hand_side(:)
    complex(dp), allocatable :: static_right_hand_side(:)
    complex(dp) :: far_field(3), finite_field(3)
    integer, allocatable :: edge_triangles(:, :), edge_vertices(:, :)
    real(dp) :: area, basis_divergence, basis_value(3), centroid(3)
    real(dp) :: vertices(3, 4)
    integer :: boundary_triangles(3, 4), edge, panel, status
    logical :: all_passed

    all_passed = .true.
    vertices(:, 1) = [0.0_dp, 0.0_dp, 0.0_dp]
    vertices(:, 2) = [1.0_dp, 0.0_dp, 0.0_dp]
    vertices(:, 3) = [0.0_dp, 1.0_dp, 0.0_dp]
    vertices(:, 4) = [0.0_dp, 0.0_dp, 1.0_dp]
    boundary_triangles(:, 1) = [1, 3, 2]
    boundary_triangles(:, 2) = [1, 2, 4]
    boundary_triangles(:, 3) = [1, 4, 3]
    boundary_triangles(:, 4) = [2, 3, 4]
    call build_maxwell_rwg_surface_space( &
        vertices, boundary_triangles, edge_vertices, edge_triangles, status)

    call assemble_maxwell_plane_wave_rhs_rwg_3d( &
        vertices, boundary_triangles, [0.0_dp, 0.0_dp, 1.0_dp], &
        [cmplx(1.0_dp, 0.0_dp, dp), cmplx(0.0_dp, 0.0_dp, dp), &
        cmplx(0.0_dp, 0.0_dp, dp)], 0.0_dp, 6, &
        static_right_hand_side, status)
    allocate(right_hand_side(size(edge_vertices, 2)))
    right_hand_side = cmplx(0.0_dp, 0.0_dp, dp)
    do edge = 1, size(edge_vertices, 2)
        do panel = 1, 2
            centroid = sum(vertices(:, boundary_triangles(:, &
                edge_triangles(panel, edge))), dim=2)/3.0_dp
            area = triangle_area(vertices(:, boundary_triangles(:, &
                edge_triangles(panel, edge))))
            call evaluate_maxwell_rwg_basis( &
                vertices, boundary_triangles, edge_vertices, edge_triangles, &
                edge, edge_triangles(panel, edge), centroid, basis_value, &
                basis_divergence, status)
            right_hand_side(edge) = right_hand_side(edge) - &
                area*basis_value(1)
        end do
    end do
    call record_condition(status == 0 .and. maxval(abs( &
        static_right_hand_side - right_hand_side)) < 2.0e-14_dp, &
        "Plane-wave load has the exact constant-field affine RWG limit")

    call solve_maxwell_pec_efie_rwg_3d( &
        vertices, boundary_triangles, [0.0_dp, 0.0_dp, 1.0_dp], &
        [cmplx(1.0_dp, 0.0_dp, dp), cmplx(0.0_dp, 0.0_dp, dp), &
        cmplx(0.0_dp, 0.0_dp, dp)], 0.8_dp, 1.7_dp, 8, 1.0e-3_dp, 1, &
        density, status)
    call assemble_maxwell_efie_rwg_3d( &
        vertices, boundary_triangles, 0.8_dp, 1.7_dp, 8, 1.0e-3_dp, 1, &
        matrix, status)
    call assemble_maxwell_plane_wave_rhs_rwg_3d( &
        vertices, boundary_triangles, [0.0_dp, 0.0_dp, 1.0_dp], &
        [cmplx(1.0_dp, 0.0_dp, dp), cmplx(0.0_dp, 0.0_dp, dp), &
        cmplx(0.0_dp, 0.0_dp, dp)], 0.8_dp, 8, right_hand_side, status)
    call record_condition(status == 0 .and. sqrt(sum(abs( &
        matmul(matrix, density) - right_hand_side)**2))/ &
        sqrt(sum(abs(right_hand_side)**2)) < 2.0e-12_dp, &
        "PEC EFIE solve satisfies the independently reassembled linear system")
    call evaluate_maxwell_efie_far_field_rwg_3d( &
        vertices, boundary_triangles, density, &
        [1.0_dp, 1.0_dp, 1.0_dp]/sqrt(3.0_dp), 0.8_dp, 1.7_dp, 12, &
        far_field, status)
    call evaluate_maxwell_efie_field_rwg_3d( &
        vertices, boundary_triangles, density, &
        200.0_dp*[1.0_dp, 1.0_dp, 1.0_dp]/sqrt(3.0_dp), &
        0.8_dp, 1.7_dp, 12, finite_field, status)
    call record_condition(status == 0 .and. sqrt(sum(abs( &
        200.0_dp*exp(cmplx(0.0_dp, -160.0_dp, dp))*finite_field - &
        far_field)**2))/sqrt(sum(abs(far_field)**2)) < 1.0e-2_dp, &
        "RWG far field matches the large-radius radiated-field limit")

    call assemble_maxwell_plane_wave_rhs_rwg_3d( &
        vertices, boundary_triangles, [0.0_dp, 0.0_dp, 2.0_dp], &
        [cmplx(1.0_dp, 0.0_dp, dp), cmplx(0.0_dp, 0.0_dp, dp), &
        cmplx(0.0_dp, 0.0_dp, dp)], 0.8_dp, 8, right_hand_side, status)
    call record_condition(status /= 0, &
        "Plane-wave load rejects a non-unit propagation direction")
    call check_summary("Three-dimensional PEC Maxwell EFIE solve")
    if (.not. all_passed) error stop 1

contains

    pure function triangle_area(points) result(value)
        real(dp), intent(in) :: points(3, 3)
        real(dp) :: value

        value = 0.5_dp*norm2(cross_product( &
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

end program test_maxwell_pec_efie_solve_3d
