program test_public_rt_arbitrary_order_solve
    use check, only: check_condition, check_summary
    use fortfem_api, only: cell_vector_source, cell_vector_source_t, div, &
        dx, form_expr_t, inner, mesh_t, &
        operator(*), operator(+), operator(==), solve, unit_square_mesh, &
        vector_bc, vector_bc_t, vector_function, vector_function_space, &
        vector_function_space_t, vector_function_t, vector_test_function, &
        vector_test_function_t, vector_trial_function, &
        vector_trial_function_t
    use fortfem_kinds, only: dp
    use fortfem_triangle_global_dof_map, only: &
        build_triangle_trimmed_dof_map
    use fortfem_triangle_rt_arbitrary_order, only: &
        initialize_triangle_raviart_thomas, triangle_rt_basis_t
    use fortfem_triangle_vector_interpolation, only: &
        evaluate_triangle_rt_interpolant
    implicit none

    real(dp) :: maximum_error
    integer :: degree

    do degree = 0, 3
        call solve_constant_field(degree, maximum_error)
        call check_condition(maximum_error < 2.0e-10_dp, &
            "Public arbitrary-order RT solve reproduces a constant field")
    end do
    call check_summary("Public arbitrary-order RT solve")

contains

    subroutine solve_constant_field(degree, error)
        integer, intent(in) :: degree
        real(dp), intent(out) :: error

        type(form_expr_t) :: bilinear_form, linear_form
        type(cell_vector_source_t) :: cell_source
        type(mesh_t) :: mesh
        type(vector_bc_t) :: boundary_condition
        type(vector_function_space_t) :: space
        type(vector_function_t) :: source, solution
        type(vector_test_function_t) :: test_field
        type(vector_trial_function_t) :: trial_field
        real(dp), allocatable :: source_values(:, :)
        integer :: order, expected_dof_count

        mesh = unit_square_mesh(3)
        space = vector_function_space(mesh, "RT", degree)
        order = degree + 1
        expected_dof_count = order * mesh%data%n_edges + &
            order * (order - 1) * mesh%data%n_triangles
        if (space%ndof /= expected_dof_count) then
            error stop "public arbitrary-order RT space has wrong dimension"
        end if
        trial_field = vector_trial_function(space)
        test_field = vector_test_function(space)
        bilinear_form = ( &
            inner(div(trial_field), div(test_field)) + &
            inner(trial_field, test_field)) * dx
        source = vector_function(space)
        source%values(:, 1) = 1.25_dp
        source%values(:, 2) = -0.75_dp
        if (mod(degree, 2) == 0) then
            linear_form = inner(source, test_field) * dx
        else
            allocate(source_values(2, mesh%data%n_triangles))
            source_values = spread([1.25_dp, -0.75_dp], 2, &
                mesh%data%n_triangles)
            cell_source = cell_vector_source(source_values)
            linear_form = inner(cell_source, test_field) * dx
        end if
        boundary_condition = vector_bc( &
            space, [1.25_dp, -0.75_dp], "normal")
        solution = vector_function(space)

        call solve( &
            bilinear_form == linear_form, solution, boundary_condition, &
            "direct")
        call measure_constant_field_error( &
            mesh, degree, solution%values(:, 1), &
            [1.25_dp, -0.75_dp], error)
    end subroutine solve_constant_field

    subroutine measure_constant_field_error( &
            public_mesh, degree, coefficients, exact_value, error)
        type(mesh_t), intent(inout) :: public_mesh
        integer, intent(in) :: degree
        real(dp), intent(in) :: coefficients(:), exact_value(2)
        real(dp), intent(out) :: error

        type(triangle_rt_basis_t) :: basis
        integer, allocatable :: global_dofs(:, :), transforms(:, :)
        real(dp), allocatable :: local_dofs(:)
        real(dp) :: divergence_value, value(2), vertices(2, 3)
        integer :: global_count, local_dof, status, triangle

        call build_triangle_trimmed_dof_map( &
            public_mesh%data, degree + 1, global_dofs, transforms, &
            global_count, status)
        if (status /= 0) error stop "arbitrary-order RT global map failed"
        call initialize_triangle_raviart_thomas(degree, basis, status)
        if (status /= 0) error stop "arbitrary-order RT basis failed"
        allocate(local_dofs(size(global_dofs, 1)))
        error = 0.0_dp
        do triangle = 1, public_mesh%data%n_triangles
            do local_dof = 1, size(local_dofs)
                local_dofs(local_dof) = &
                    real(transforms(local_dof, triangle), dp) * &
                    coefficients(global_dofs(local_dof, triangle))
            end do
            vertices = public_mesh%data%vertices(:, &
                public_mesh%data%triangles(:, triangle))
            call evaluate_triangle_rt_interpolant( &
                vertices, basis, local_dofs, 0.23_dp, 0.31_dp, value, &
                divergence_value, status)
            if (status /= 0) error stop "RT field evaluation failed"
            error = max(error, maxval(abs(value - exact_value)))
            error = max(error, abs(divergence_value))
        end do
    end subroutine measure_constant_field_error

end program test_public_rt_arbitrary_order_solve
