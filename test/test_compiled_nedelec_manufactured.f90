program test_compiled_nedelec_manufactured
    use check, only: check_condition, check_summary
    use fortfem_core, only: mesh_t, rectangle_mesh, vector_function_space_t
    use fortfem_feec, only: cell_vector_source, cell_vector_source_t, curl, &
        dx, form_expr_t, inner, operator(*), operator(+), operator(==), solve, &
        vector_bc, vector_bc_t, vector_function, vector_function_space, &
        vector_function_t, vector_test_function, vector_test_function_t, &
        vector_trial_function, vector_trial_function_t
    use fortfem_gauss_quadrature_2d, only: &
        gauss_quadrature_triangle_t, get_gauss_quadrature_triangle
    use fortfem_kinds, only: dp
    use fortfem_mesh_2d, only: mesh_2d_t
    use fortfem_nedelec_field_2d, only: evaluate_nedelec_field_2d
    implicit none

    real(dp) :: coarse_error, fine_error, first_rate
    real(dp) :: medium_error, second_rate
    logical :: all_passed

    all_passed = .true.
    coarse_error = solve_manufactured_problem(9)
    medium_error = solve_manufactured_problem(17)
    fine_error = solve_manufactured_problem(33)
    first_rate = log(coarse_error / medium_error) / log(2.0_dp)
    second_rate = log(medium_error / fine_error) / log(2.0_dp)

    call record_condition( &
        first_rate > 0.8_dp .and. second_rate > 0.8_dp, &
        "Public Nedelec solve attains first-order field convergence")
    call record_condition(fine_error < 2.0e-2_dp, &
        "Public Nedelec solve reaches the manufactured field")

    call check_summary("Compiled Nedelec manufactured solution")
    if (.not. all_passed) error stop 1

contains

    real(dp) function solve_manufactured_problem(vertex_count) result(error)
        integer, intent(in) :: vertex_count

        type(cell_vector_source_t) :: source
        type(form_expr_t) :: bilinear_form, linear_form
        type(mesh_t) :: mesh
        type(vector_bc_t) :: boundary_condition
        type(vector_function_space_t) :: space
        type(vector_function_t) :: solution
        type(vector_test_function_t) :: test_field
        type(vector_trial_function_t) :: trial_field
        real(dp), allocatable :: source_values(:, :)
        real(dp) :: centroid(2)
        integer :: triangle

        mesh = rectangle_mesh( &
            vertex_count, vertex_count, &
            [0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp])
        space = vector_function_space(mesh, "Nedelec", 1)
        trial_field = vector_trial_function(space)
        test_field = vector_test_function(space)
        bilinear_form = ( &
            inner(curl(trial_field), curl(test_field)) + &
            inner(trial_field, test_field)) * dx

        allocate(source_values(2, mesh%data%n_triangles))
        do triangle = 1, mesh%data%n_triangles
            centroid = sum(mesh%data%vertices(:, &
                mesh%data%triangles(:, triangle)), dim=2) / 3.0_dp
            source_values(:, triangle) = manufactured_source(centroid)
        end do
        source = cell_vector_source(source_values)
        linear_form = inner(source, test_field) * dx
        boundary_condition = vector_bc( &
            space, [0.0_dp, 0.0_dp], "tangential")
        solution = vector_function(space)
        call solve( &
            bilinear_form == linear_form, solution, boundary_condition, &
            "direct")
        error = field_l2_error(mesh%data, solution%values(:, 1))
    end function solve_manufactured_problem

    pure function manufactured_source(point) result(value)
        real(dp), intent(in) :: point(2)
        real(dp) :: value(2)
        real(dp) :: x, y

        x = point(1)
        y = point(2)
        value(1) = (1.0_dp - 2.0_dp * x) * (1.0_dp - 2.0_dp * y)
        value(2) = 2.0_dp * y * (1.0_dp - y) + &
            x * (1.0_dp - x) * y * (1.0_dp - y)
    end function manufactured_source

    pure function exact_field(point) result(value)
        real(dp), intent(in) :: point(2)
        real(dp) :: value(2)

        value(1) = 0.0_dp
        value(2) = point(1) * (1.0_dp - point(1)) * &
            point(2) * (1.0_dp - point(2))
    end function exact_field

    real(dp) function field_l2_error(mesh, dofs) result(error)
        type(mesh_2d_t), intent(in) :: mesh
        real(dp), intent(in) :: dofs(:)

        type(gauss_quadrature_triangle_t) :: quadrature
        complex(dp), allocatable :: complex_dofs(:)
        complex(dp) :: field_value(2)
        real(dp) :: determinant, error_squared, exact_value(2)
        real(dp) :: point(2), vertex_a(2), vertex_b(2), vertex_c(2)
        integer :: point_index, triangle

        quadrature = get_gauss_quadrature_triangle(7)
        allocate(complex_dofs(size(dofs)))
        complex_dofs = cmplx(dofs, 0.0_dp, dp)
        error_squared = 0.0_dp
        do triangle = 1, mesh%n_triangles
            vertex_a = mesh%vertices(:, mesh%triangles(1, triangle))
            vertex_b = mesh%vertices(:, mesh%triangles(2, triangle))
            vertex_c = mesh%vertices(:, mesh%triangles(3, triangle))
            determinant = &
                (vertex_b(1) - vertex_a(1)) * &
                (vertex_c(2) - vertex_a(2)) - &
                (vertex_b(2) - vertex_a(2)) * &
                (vertex_c(1) - vertex_a(1))
            do point_index = 1, quadrature%n_points
                point = vertex_a + &
                    quadrature%xi(point_index) * (vertex_b - vertex_a) + &
                    quadrature%eta(point_index) * (vertex_c - vertex_a)
                call evaluate_nedelec_field_2d( &
                    mesh, triangle, point(1), point(2), complex_dofs, &
                    field_value)
                exact_value = exact_field(point)
                error_squared = error_squared + &
                    determinant * quadrature%weights(point_index) * &
                    sum((real(field_value, dp) - exact_value)**2)
            end do
        end do
        call quadrature%destroy()
        error = sqrt(error_squared)
    end function field_l2_error

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_compiled_nedelec_manufactured
