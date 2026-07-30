program test_public_nedelec_arbitrary_order_solve
    use check, only: check_condition, check_summary
    use fortfem_api, only: curl, dx, form_expr_t, inner, mesh_t, &
        operator(*), operator(+), operator(==), solve, unit_square_mesh, &
        vector_bc, vector_bc_t, vector_function, vector_function_space, &
        vector_function_space_t, vector_function_t, vector_test_function, &
        vector_test_function_t, vector_trial_function, &
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

    real(dp) :: maximum_error
    integer :: order

    do order = 2, 4
        call solve_constant_field(order, maximum_error)
        call check_condition(maximum_error < 2.0e-10_dp, &
            "Public arbitrary-order curl-curl solve reproduces a constant")
    end do
    call check_summary("Public arbitrary-order Nedelec solve")

contains

    subroutine solve_constant_field(order, error)
        integer, intent(in) :: order
        real(dp), intent(out) :: error

        type(form_expr_t) :: bilinear_form, linear_form
        type(mesh_t) :: mesh
        type(vector_bc_t) :: boundary_condition
        type(vector_function_space_t) :: space
        type(vector_function_t) :: source, solution
        type(vector_test_function_t) :: test_field
        type(vector_trial_function_t) :: trial_field
        integer :: expected_dof_count

        mesh = unit_square_mesh(3)
        space = vector_function_space(mesh, "Nedelec", order)
        expected_dof_count = order * mesh%data%n_edges + &
            order * (order - 1) * mesh%data%n_triangles
        if (space%ndof /= expected_dof_count) then
            error stop "public arbitrary-order space has the wrong dimension"
        end if
        trial_field = vector_trial_function(space)
        test_field = vector_test_function(space)
        bilinear_form = ( &
            inner(curl(trial_field), curl(test_field)) + &
            inner(trial_field, test_field)) * dx
        source = vector_function(space)
        source%values(:, 1) = 1.25_dp
        source%values(:, 2) = -0.75_dp
        linear_form = inner(source, test_field) * dx
        boundary_condition = vector_bc( &
            space, [1.25_dp, -0.75_dp], "tangential")
        solution = vector_function(space)

        call solve( &
            bilinear_form == linear_form, solution, boundary_condition, &
            "direct")
        call measure_constant_field_error( &
            mesh, order, solution%values(:, 1), &
            [1.25_dp, -0.75_dp], error)
    end subroutine solve_constant_field

    subroutine measure_constant_field_error( &
            public_mesh, order, coefficients, exact_value, error)
        type(mesh_t), intent(inout) :: public_mesh
        integer, intent(in) :: order
        real(dp), intent(in) :: coefficients(:), exact_value(2)
        real(dp), intent(out) :: error

        type(triangle_nedelec_first_kind_t) :: basis
        integer, allocatable :: global_dofs(:, :), transforms(:, :)
        real(dp), allocatable :: local_dofs(:)
        real(dp) :: curl_value, value(2), vertices(2, 3)
        integer :: global_count, local_dof, status, triangle

        call build_triangle_trimmed_dof_map( &
            public_mesh%data, order, global_dofs, transforms, &
            global_count, status)
        if (status /= 0) error stop "arbitrary-order global map failed"
        call initialize_triangle_nedelec_first_kind(order, basis, status)
        if (status /= 0) error stop "arbitrary-order basis failed"
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
            call evaluate_triangle_nedelec_interpolant( &
                vertices, basis, local_dofs, 0.23_dp, 0.31_dp, value, &
                curl_value, status)
            if (status /= 0) error stop "order-two field evaluation failed"
            error = max(error, maxval(abs(value - exact_value)))
            error = max(error, abs(curl_value))
        end do
    end subroutine measure_constant_field_error

end program test_public_nedelec_arbitrary_order_solve
