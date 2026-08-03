program test_paper_magnetic_anisotropic_2d
    use check, only: check_condition, check_summary
    use fortfem_core, only: mesh_t, unit_square_mesh, vector_function_space_t
    use fortfem_feec, only: &
        cell_tensor_coefficient, cell_tensor_coefficient_t, &
        cell_vector_source, cell_vector_source_t, curl, dx, form_expr_t, &
        inner, operator(*), operator(+), operator(==), solve, &
        vector_bc_edge_moments, vector_bc_t, &
        vector_function, vector_function_space, &
        vector_function_t, vector_test_function, vector_test_function_t, &
        vector_trial_function, vector_trial_function_t
    use fortfem_kinds, only: dp
    use fortfem_triangle_global_dof_map, only: &
        build_triangle_trimmed_dof_map
    use fortfem_triangle_nedelec_arbitrary_order, only: &
        initialize_triangle_nedelec_first_kind, &
        triangle_nedelec_first_kind_t
    use fortfem_triangle_vector_interpolation, only: &
        evaluate_triangle_nedelec_interpolant
    implicit none

    real(dp) :: coarse_error, fine_error
    logical :: all_passed

    all_passed = .true.
    coarse_error = solve_anisotropic_case(8)
    fine_error = solve_anisotropic_case(16)
    call record_condition(fine_error < 4.0e-2_dp, &
        "Anisotropic transverse magnetic field reaches the analytical field")
    call record_condition(fine_error < 0.65_dp * coarse_error, &
        "Anisotropic transverse magnetic field converges under refinement")
    call check_summary("Paper magnetic anisotropic transverse case")
    if (.not. all_passed) error stop 1

contains

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

    real(dp) function solve_anisotropic_case(divisions) result(error)
        integer, intent(in) :: divisions

        type(cell_tensor_coefficient_t) :: material
        type(cell_vector_source_t) :: source
        type(form_expr_t) :: bilinear_form, linear_form
        type(mesh_t) :: mesh
        type(vector_bc_t) :: boundary_condition
        type(vector_function_space_t) :: space
        type(vector_function_t) :: solution
        type(vector_test_function_t) :: test_field
        type(vector_trial_function_t) :: trial_field
        real(dp), allocatable :: boundary_values(:), source_values(:, :)
        real(dp), allocatable :: tensors(:, :, :)
        real(dp) :: centroid(2), edge_vector(2), midpoint(2)
        integer :: degree_of_freedom, edge, triangle

        mesh = unit_square_mesh(divisions + 1)
        space = vector_function_space(mesh, "Nedelec", 1)
        trial_field = vector_trial_function(space)
        test_field = vector_test_function(space)
        allocate(tensors(2, 2, mesh%data%n_triangles))
        allocate(source_values(2, mesh%data%n_triangles))
        do triangle = 1, mesh%data%n_triangles
            centroid = sum(mesh%data%vertices(:, &
                mesh%data%triangles(:, triangle)), dim=2) / 3.0_dp
            call material_tensor(centroid, tensors(:, :, triangle))
            source_values(:, triangle) = matmul( &
                tensors(:, :, triangle), exact_field(centroid))
        end do
        material = cell_tensor_coefficient(tensors)
        source = cell_vector_source(source_values)
        bilinear_form = ( &
            inner(curl(trial_field), curl(test_field)) + &
            material * inner(trial_field, test_field)) * dx
        linear_form = inner(source, test_field) * dx

        allocate(boundary_values(space%ndof))
        boundary_values = 0.0_dp
        do edge = 1, mesh%data%n_edges
            edge_vector = mesh%data%vertices(:, mesh%data%edges(2, edge)) - &
                mesh%data%vertices(:, mesh%data%edges(1, edge))
            midpoint = 0.5_dp * ( &
                mesh%data%vertices(:, mesh%data%edges(1, edge)) + &
                mesh%data%vertices(:, mesh%data%edges(2, edge)))
            degree_of_freedom = mesh%data%edge_to_dof(edge) + 1
            boundary_values(degree_of_freedom) = &
                dot_product(exact_field(midpoint), edge_vector)
        end do
        boundary_condition = vector_bc_edge_moments( &
            space, boundary_values, "tangential")
        solution = vector_function(space)
        call solve( &
            bilinear_form == linear_form, solution, boundary_condition, &
            "direct")
        call measure_error(mesh, solution%values(:, 1), error)
    end function solve_anisotropic_case

    pure function exact_field(point) result(value)
        real(dp), intent(in) :: point(2)
        real(dp) :: value(2)

        value = [-point(2), point(1)]
    end function exact_field

    pure subroutine material_tensor(point, tensor)
        real(dp), intent(in) :: point(2)
        real(dp), intent(out) :: tensor(2, 2)

        if (point(1) < 0.5_dp) then
            tensor = reshape([2.0_dp, 0.3_dp, 0.3_dp, 1.2_dp], [2, 2])
        else
            tensor = reshape([1.1_dp, -0.2_dp, -0.2_dp, 2.4_dp], [2, 2])
        end if
    end subroutine material_tensor

    subroutine measure_error(mesh, coefficients, error)
        type(mesh_t), intent(inout) :: mesh
        real(dp), intent(in) :: coefficients(:)
        real(dp), intent(out) :: error

        type(triangle_nedelec_first_kind_t) :: basis
        integer, allocatable :: global_dofs(:, :), transforms(:, :)
        real(dp), allocatable :: local_dofs(:)
        real(dp) :: centroid(2), curl_value, value(2), vertices(2, 3)
        integer :: global_count, local_dof, status, triangle

        call build_triangle_trimmed_dof_map( &
            mesh%data, 1, global_dofs, transforms, global_count, status)
        if (status /= 0) error stop "anisotropic global map failed"
        call initialize_triangle_nedelec_first_kind(1, basis, status)
        if (status /= 0) error stop "anisotropic basis failed"
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
            if (status /= 0) error stop "anisotropic field evaluation failed"
            error = max(error, maxval(abs(value - exact_field(centroid))))
        end do
    end subroutine measure_error

end program test_paper_magnetic_anisotropic_2d
