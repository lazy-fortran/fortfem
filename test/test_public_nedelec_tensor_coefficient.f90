program test_public_nedelec_tensor_coefficient
    use check, only: check_condition, check_summary
    use fortfem_core, only: mesh_t, unit_square_mesh
    use fortfem_feec, only: &
        cell_tensor_coefficient, cell_tensor_coefficient_t, &
        cell_vector_source, cell_vector_source_t, curl, dx, form_expr_t, &
        inner, operator(*), operator(+), operator(==), solve, &
        vector_bc, vector_bc_t, vector_function, &
        vector_function_space, vector_function_space_t, vector_function_t, &
        vector_test_function, vector_test_function_t, vector_trial_function, &
        vector_trial_function_t
    use fortfem_kinds, only: dp
    use fortfem_triangle_global_dof_map, only: &
        build_triangle_trimmed_dof_map
    use fortfem_triangle_nedelec_arbitrary_order, only: &
        initialize_triangle_nedelec_first_kind, &
        triangle_nedelec_first_kind_t
    use fortfem_triangle_vector_interpolation, only: &
        evaluate_triangle_nedelec_interpolant
    implicit none

    real(dp) :: error
    integer :: order
    logical :: all_passed

    all_passed = .true.
    do order = 1, 4
        call solve_tensor_case(order, error)
        call record_condition(error < 3.0e-10_dp, &
            "Public tensor-weighted Nedelec solve reproduces the exact field")
    end do
    call check_summary("Public Nedelec cell tensor coefficient")
    if (.not. all_passed) error stop 1

contains

    subroutine solve_tensor_case(order, error)
        integer, intent(in) :: order
        real(dp), intent(out) :: error

        type(cell_tensor_coefficient_t) :: material
        type(cell_vector_source_t) :: source
        type(form_expr_t) :: bilinear_form, linear_form
        type(mesh_t) :: mesh
        type(vector_bc_t) :: boundary_condition
        type(vector_function_space_t) :: space
        type(vector_function_t) :: solution
        type(vector_test_function_t) :: test_field
        type(vector_trial_function_t) :: trial_field
        real(dp), allocatable :: source_values(:, :), tensors(:, :, :)
        real(dp) :: exact_field(2)
        integer :: triangle

        mesh = unit_square_mesh(2)
        space = vector_function_space(mesh, "Nedelec", order)
        trial_field = vector_trial_function(space)
        test_field = vector_test_function(space)
        exact_field = [1.25_dp, -0.75_dp]
        allocate(tensors(2, 2, mesh%data%n_triangles))
        allocate(source_values(2, mesh%data%n_triangles))
        do triangle = 1, mesh%data%n_triangles
            tensors(:, :, triangle) = reshape([ &
                1.5_dp + 0.1_dp * real(triangle, dp), 0.25_dp, &
                0.25_dp, 2.0_dp + 0.2_dp * real(triangle, dp)], [2, 2])
            source_values(:, triangle) = &
                matmul(tensors(:, :, triangle), exact_field)
        end do
        material = cell_tensor_coefficient(tensors)
        source = cell_vector_source(source_values)
        bilinear_form = ( &
            inner(curl(trial_field), curl(test_field)) + &
            material * inner(trial_field, test_field)) * dx
        linear_form = inner(source, test_field) * dx
        boundary_condition = vector_bc(space, exact_field, "tangential")
        solution = vector_function(space)
        call solve( &
            bilinear_form == linear_form, solution, boundary_condition, &
            "direct")
        call measure_constant_field_error( &
            mesh, order, solution%values(:, 1), exact_field, error)
    end subroutine solve_tensor_case

    subroutine measure_constant_field_error( &
            mesh, order, coefficients, exact_value, error)
        type(mesh_t), intent(inout) :: mesh
        integer, intent(in) :: order
        real(dp), intent(in) :: coefficients(:), exact_value(2)
        real(dp), intent(out) :: error

        type(triangle_nedelec_first_kind_t) :: basis
        integer, allocatable :: global_dofs(:, :), transforms(:, :)
        real(dp), allocatable :: local_dofs(:)
        real(dp) :: curl_value, value(2), vertices(2, 3)
        integer :: global_count, local_dof, status, triangle

        call build_triangle_trimmed_dof_map( &
            mesh%data, order, global_dofs, transforms, global_count, status)
        if (status /= 0) error stop "tensor coefficient global map failed"
        call initialize_triangle_nedelec_first_kind(order, basis, status)
        if (status /= 0) error stop "tensor coefficient basis failed"
        allocate(local_dofs(size(global_dofs, 1)))
        error = 0.0_dp
        do triangle = 1, mesh%data%n_triangles
            do local_dof = 1, size(local_dofs)
                local_dofs(local_dof) = &
                    real(transforms(local_dof, triangle), dp) * &
                    coefficients(global_dofs(local_dof, triangle))
            end do
            vertices = mesh%data%vertices(:, mesh%data%triangles(:, triangle))
            call evaluate_triangle_nedelec_interpolant( &
                vertices, basis, local_dofs, 0.23_dp, 0.31_dp, value, &
                curl_value, status)
            if (status /= 0) error stop "tensor coefficient evaluation failed"
            error = max(error, maxval(abs(value - exact_value)))
            error = max(error, abs(curl_value))
        end do
    end subroutine measure_constant_field_error

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_public_nedelec_tensor_coefficient
