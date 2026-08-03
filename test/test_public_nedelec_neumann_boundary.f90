program test_public_nedelec_neumann_boundary
    use check, only: check_condition, check_summary
    use fortfem_core, only: mesh_t, unit_square_mesh, vector_function_space_t
    use fortfem_feec, only: &
        cell_vector_source, cell_vector_source_t, curl, dx, form_expr_t, &
        inner, operator(*), operator(+), operator(==), solve, &
        vector_bc_t, vector_function, vector_function_space, vector_function_t, &
        vector_test_function, vector_test_function_t, vector_trial_function, &
        vector_trial_function_t
    use fortfem_api_spaces, only: vector_bc_on_edges, &
        vector_neumann_bc_on_edges, vector_neumann_bc_t
    use fortfem_kinds, only: dp
    use fortnum_quadrature, only: gauss_legendre_ab
    use fortfem_triangle_global_dof_map, only: &
        build_triangle_trimmed_dof_map
    use fortfem_triangle_nedelec_arbitrary_order, only: &
        initialize_triangle_nedelec_first_kind, &
        triangle_nedelec_first_kind_t
    use fortfem_triangle_vector_interpolation, only: &
        evaluate_triangle_nedelec_interpolant
    implicit none

    real(dp) :: coarse_error, fine_error
    integer :: order
    logical :: all_passed

    all_passed = .true.
    do order = 1, 4
        coarse_error = solve_mixed_boundary_case(6, order)
        fine_error = solve_mixed_boundary_case(12, order)
        call record_condition(fine_error < 7.0e-2_dp, &
            "Nonzero curl-Neumann data reaches the analytical magnetic field")
        call record_condition(fine_error < 0.7_dp * coarse_error, &
            "Mixed H(curl) boundary solution converges under refinement")
    end do
    call check_summary("Public Nedelec curl-Neumann boundary")
    if (.not. all_passed) error stop 1

contains

    real(dp) function solve_mixed_boundary_case(divisions, order) result(error)
        integer, intent(in) :: divisions, order

        type(cell_vector_source_t) :: source
        type(form_expr_t) :: bilinear_form, linear_form
        type(mesh_t) :: mesh
        type(vector_bc_t) :: essential_condition
        type(vector_function_space_t) :: space
        type(vector_function_t) :: solution
        type(vector_neumann_bc_t) :: natural_condition
        type(vector_test_function_t) :: test_field
        type(vector_trial_function_t) :: trial_field
        logical, allocatable :: essential_edges(:), natural_edges(:)
        real(dp), allocatable :: boundary_values(:), flux_values(:)
        real(dp), allocatable :: source_values(:, :)
        real(dp), allocatable :: nodes(:), weights(:)
        real(dp) :: centroid(2), edge_vector(2), midpoint(2), point(2)
        integer :: degree_of_freedom, edge, moment, node, triangle

        mesh = unit_square_mesh(divisions + 1)
        space = vector_function_space(mesh, "Nedelec", order)
        trial_field = vector_trial_function(space)
        test_field = vector_test_function(space)
        allocate(source_values(2, mesh%data%n_triangles))
        do triangle = 1, mesh%data%n_triangles
            centroid = sum(mesh%data%vertices(:, &
                mesh%data%triangles(:, triangle)), dim=2) / 3.0_dp
            source_values(:, triangle) = exact_field(centroid)
        end do
        source = cell_vector_source(source_values)
        bilinear_form = ( &
            inner(curl(trial_field), curl(test_field)) + &
            inner(trial_field, test_field)) * dx
        linear_form = inner(source, test_field) * dx

        allocate(boundary_values(space%ndof))
        allocate(flux_values(mesh%data%n_edges))
        allocate(essential_edges(mesh%data%n_edges))
        allocate(natural_edges(mesh%data%n_edges))
        boundary_values = 0.0_dp
        flux_values = 0.0_dp
        essential_edges = .false.
        natural_edges = .false.
        allocate(nodes(order + 1), weights(order + 1))
        call gauss_legendre_ab( &
            order + 1, 0.0_dp, 1.0_dp, nodes, weights)
        do edge = 1, mesh%data%n_edges
            if (.not. mesh%data%is_boundary_edge(edge)) cycle
            edge_vector = mesh%data%vertices(:, mesh%data%edges(2, edge)) - &
                mesh%data%vertices(:, mesh%data%edges(1, edge))
            midpoint = 0.5_dp * ( &
                mesh%data%vertices(:, mesh%data%edges(1, edge)) + &
                mesh%data%vertices(:, mesh%data%edges(2, edge)))
            do moment = 1, order
                degree_of_freedom = &
                    mesh%data%edge_to_dof(edge) * order + moment
                do node = 1, size(nodes)
                    point = mesh%data%vertices(:, &
                        mesh%data%edges(1, edge)) + nodes(node) * edge_vector
                    boundary_values(degree_of_freedom) = &
                        boundary_values(degree_of_freedom) + &
                        weights(node) * shifted_legendre( &
                        moment - 1, nodes(node)) * &
                        dot_product(exact_field(point), edge_vector)
                end do
            end do
            if (midpoint(1) < 1.0e-12_dp .or. &
                midpoint(2) < 1.0e-12_dp) then
                essential_edges(edge) = .true.
            else
                natural_edges(edge) = .true.
                flux_values(edge) = 2.0_dp
            end if
        end do
        essential_condition = vector_bc_on_edges( &
            space, boundary_values, essential_edges, "tangential")
        natural_condition = vector_neumann_bc_on_edges( &
            space, flux_values, natural_edges)
        solution = vector_function(space)
        call solve( &
            bilinear_form == linear_form, solution, essential_condition, &
            natural_condition, "direct")
        call measure_error(mesh, order, solution%values(:, 1), error)
    end function solve_mixed_boundary_case

    pure function exact_field(point) result(value)
        real(dp), intent(in) :: point(2)
        real(dp) :: value(2)

        value = [-point(2), point(1)]
    end function exact_field

    subroutine measure_error(mesh, order, coefficients, error)
        type(mesh_t), intent(inout) :: mesh
        integer, intent(in) :: order
        real(dp), intent(in) :: coefficients(:)
        real(dp), intent(out) :: error

        type(triangle_nedelec_first_kind_t) :: basis
        integer, allocatable :: global_dofs(:, :), transforms(:, :)
        real(dp), allocatable :: local_dofs(:)
        real(dp) :: centroid(2), curl_value, value(2), vertices(2, 3)
        integer :: global_count, local_dof, status, triangle

        call build_triangle_trimmed_dof_map( &
            mesh%data, order, global_dofs, transforms, global_count, status)
        if (status /= 0) error stop "curl-Neumann global map failed"
        call initialize_triangle_nedelec_first_kind(order, basis, status)
        if (status /= 0) error stop "curl-Neumann basis failed"
        allocate(local_dofs(size(global_dofs, 1)))
        error = 0.0_dp
        do triangle = 1, mesh%data%n_triangles
            do local_dof = 1, size(local_dofs)
                local_dofs(local_dof) = &
                    real(transforms(local_dof, triangle), dp) * &
                    coefficients(global_dofs(local_dof, triangle))
            end do
            vertices = mesh%data%vertices(:, mesh%data%triangles(:, triangle))
            centroid = sum(vertices, dim=2) / 3.0_dp
            call evaluate_triangle_nedelec_interpolant( &
                vertices, basis, local_dofs, 1.0_dp / 3.0_dp, &
                1.0_dp / 3.0_dp, value, curl_value, status)
            if (status /= 0) error stop "curl-Neumann evaluation failed"
            error = max(error, maxval(abs(value - exact_field(centroid))))
            error = max(error, abs(curl_value - 2.0_dp))
        end do
    end subroutine measure_error

    pure real(dp) function shifted_legendre(degree, parameter) result(value)
        integer, intent(in) :: degree
        real(dp), intent(in) :: parameter

        select case (degree)
        case (0)
            value = 1.0_dp
        case (1)
            value = 2.0_dp * parameter - 1.0_dp
        case (2)
            value = 6.0_dp * parameter**2 - 6.0_dp * parameter + 1.0_dp
        case (3)
            value = 20.0_dp * parameter**3 - 30.0_dp * parameter**2 + &
                12.0_dp * parameter - 1.0_dp
        case default
            value = 0.0_dp
        end select
    end function shifted_legendre

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_public_nedelec_neumann_boundary
