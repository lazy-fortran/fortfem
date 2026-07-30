program test_adaptive_surface_bem
    use check, only: check_condition, check_summary
    use fortfem_api, only: &
        estimate_helmholtz_p0_two_level_residual_3d, &
        estimate_laplace_p0_two_level_residual_3d, &
        generate_sphere_surface_mesh, mark_bem_dorfler, &
        refine_surface_mesh_marked
    use fortfem_kinds, only: dp
    implicit none

    integer, allocatable :: parent(:), refined_triangles(:, :), triangles(:, :)
    real(dp), allocatable :: density(:), indicators(:), refined_vertices(:, :)
    real(dp), allocatable :: vertices(:, :)
    complex(dp), allocatable :: complex_density(:)
    logical, allocatable :: marked(:)
    real(dp) :: captured, expected_indicator, total
    integer :: edge, first, occurrences, second, status
    logical :: all_passed

    all_passed = .true.
    call generate_sphere_surface_mesh(1.0_dp, 0, vertices, triangles)
    allocate(marked(size(triangles, 2)), source=.false.)
    marked(1) = .true.
    call refine_surface_mesh_marked( &
        vertices, triangles, marked, refined_vertices, refined_triangles, &
        parent, status)
    if (status /= 0) error stop "marked surface refinement failed"
    call record_condition( &
        count(parent == 1) == 4 .and. &
        count(parent /= 1) < size(triangles, 2) + 3, &
        "marked refinement splits its panel and only conforming neighbors")
    call record_condition( &
        all_outward(refined_vertices, refined_triangles), &
        "marked refinement preserves outward surface orientation")
    occurrences = 0
    do first = 1, size(refined_triangles, 2)
        do edge = 1, 3
            occurrences = 0
            do second = 1, size(refined_triangles, 2)
                occurrences = occurrences + count_edge( &
                    refined_triangles(edge, first), &
                    refined_triangles(modulo(edge, 3) + 1, first), &
                    refined_triangles(:, second))
            end do
            if (occurrences /= 2) exit
        end do
        if (occurrences /= 2) exit
    end do
    call record_condition(occurrences == 2, &
        "marked refinement leaves a closed conforming two-manifold")

    allocate(density(size(triangles, 2)))
    density = 0.0_dp
    call estimate_laplace_p0_two_level_residual_3d( &
        vertices, triangles, density, 2.0_dp, 5, indicators, status)
    if (status /= 0) error stop "two-level BEM estimator failed"
    expected_indicator = 2.0_dp*sqrt(triangle_area( &
        vertices(:, triangles(:, 1))))
    call record_condition( &
        maxval(abs(indicators - expected_indicator)) < 1.0e-12_dp, &
        "zero-density indicators equal the exact fine residual norm")
    density = 1.0_dp
    density(1) = 0.0_dp
    call estimate_laplace_p0_two_level_residual_3d( &
        vertices, triangles, density, 1.0_dp, 5, indicators, status)
    if (status /= 0) error stop "perturbed BEM estimator failed"
    call mark_bem_dorfler(indicators, 0.6_dp, marked, status)
    if (status /= 0) error stop "BEM Dorfler marking failed"
    total = sum(indicators**2)
    captured = sum(pack(indicators**2, marked))
    call record_condition( &
        captured >= 0.6_dp*total .and. count(marked) < size(marked), &
        "Dorfler marking captures the requested residual fraction")

    allocate(complex_density(size(triangles, 2)), source=(0.0_dp, 0.0_dp))
    call estimate_helmholtz_p0_two_level_residual_3d( &
        vertices, triangles, complex_density, cmplx(1.0_dp, -2.0_dp, dp), &
        0.7_dp, 5, indicators, status)
    if (status /= 0) error stop "two-level Helmholtz estimator failed"
    expected_indicator = sqrt(5.0_dp*triangle_area( &
        vertices(:, triangles(:, 1))))
    call record_condition( &
        maxval(abs(indicators - expected_indicator)) < 1.0e-12_dp, &
        "zero-density Helmholtz indicators equal the exact residual norm")

    call check_summary("Adaptive triangular-surface BEM")
    if (.not. all_passed) error stop 1

contains

    pure integer function count_edge(first_vertex, second_vertex, cell) &
            result(count_)
        integer, intent(in) :: first_vertex, second_vertex, cell(3)

        count_ = merge(1, 0, any(cell == first_vertex) .and. &
            any(cell == second_vertex))
    end function count_edge

    pure logical function all_outward(points, cells) result(outward)
        real(dp), intent(in) :: points(:, :)
        integer, intent(in) :: cells(:, :)
        real(dp) :: first_edge(3), normal(3), second_edge(3)
        integer :: cell

        outward = .true.
        do cell = 1, size(cells, 2)
            first_edge = &
                points(:, cells(2, cell)) - points(:, cells(1, cell))
            second_edge = &
                points(:, cells(3, cell)) - points(:, cells(1, cell))
            normal = cross_product(first_edge, second_edge)
            if (dot_product(normal, sum(points(:, cells(:, cell)), dim=2)) <= &
                0.0_dp) then
                outward = .false.
                return
            end if
        end do
    end function all_outward

    pure function cross_product(first, second) result(product)
        real(dp), intent(in) :: first(3), second(3)
        real(dp) :: product(3)

        product = [ &
            first(2)*second(3) - first(3)*second(2), &
            first(3)*second(1) - first(1)*second(3), &
            first(1)*second(2) - first(2)*second(1)]
    end function cross_product

    pure real(dp) function triangle_area(points) result(area)
        real(dp), intent(in) :: points(3, 3)

        area = 0.5_dp*norm2(cross_product( &
            points(:, 2) - points(:, 1), points(:, 3) - points(:, 1)))
    end function triangle_area

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_adaptive_surface_bem
