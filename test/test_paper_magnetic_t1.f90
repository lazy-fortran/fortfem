program test_paper_magnetic_t1
    use check, only: check_condition, check_summary
    use fortfem_feec, only: curl, dx, &
        form_expr_t, init_measures, inner, operator(*), operator(+), &
        operator(==), solve, vector_bc_edge_moments, &
        vector_bc_t, vector_function, vector_function_space, vector_function_t, &
        vector_test_function, &
        vector_test_function_t, vector_trial_function, vector_trial_function_t
    use fortfem_api_spaces, only: cell_coefficient, cell_coefficient_t, &
        vector_function_space_t
    use fortfem_core, only: mesh_t, rectangle_mesh
    use fortfem_assembly_nedelec_2d, only: assemble_nedelec_weighted
    use fortfem_gauss_quadrature_2d, only: &
        gauss_quadrature_triangle_t, get_gauss_quadrature_triangle
    use fortfem_kinds, only: dp
    use fortfem_mesh_2d, only: mesh_2d_t
    use fortfem_nedelec_field_2d, only: evaluate_nedelec_field_2d
    implicit none

    real(dp) :: coarse_error, compiled_coarse_error
    real(dp) :: compiled_fine_error, fine_error
    logical :: all_passed

    all_passed = .true.
    coarse_error = solve_paper_case(10)
    fine_error = solve_paper_case(20)
    compiled_coarse_error = solve_compiled_paper_case(10)
    compiled_fine_error = solve_compiled_paper_case(20)

    call record_condition(fine_error < 7.0e-2_dp, &
        "Paper n=1 transverse field reaches the analytical solution")
    call record_condition(fine_error < 0.65_dp * coarse_error, &
        "Paper n=1 transverse field converges under mesh refinement")
    call record_condition(compiled_fine_error < 9.0e-2_dp, &
        "Compiled paper form reaches the analytical material solution")
    call record_condition( &
        compiled_fine_error < 0.7_dp * compiled_coarse_error, &
        "Compiled paper form converges under mesh refinement")

    call check_summary("Paper magnetic n=1 transverse case")
    if (.not. all_passed) error stop 1

contains

    real(dp) function solve_compiled_paper_case(divisions) &
            result(relative_error)
        integer, intent(in) :: divisions

        type(cell_coefficient_t) :: curl_coefficient, mass_coefficient
        type(form_expr_t) :: bilinear_form, linear_form
        type(mesh_t) :: mesh
        type(vector_bc_t) :: boundary_condition
        type(vector_function_space_t) :: space
        type(vector_function_t) :: solution, source
        type(vector_test_function_t) :: test_field
        type(vector_trial_function_t) :: trial_field
        complex(dp), allocatable :: field_dofs(:)
        real(dp), allocatable :: boundary_values(:)
        real(dp), allocatable :: curl_values(:), mass_values(:)
        real(dp) :: centroid(2)
        integer :: boundary_index, degree_of_freedom, edge, triangle
        integer :: vertex_a, vertex_b

        call init_measures()
        mesh = rectangle_mesh( &
            divisions + 1, divisions + 1, &
            [0.0_dp, 1.0_dp, -0.5_dp, 0.5_dp])
        space = vector_function_space(mesh, "Nedelec", 1)
        trial_field = vector_trial_function(space)
        test_field = vector_test_function(space)
        allocate(curl_values(mesh%data%n_triangles))
        allocate(mass_values(mesh%data%n_triangles))
        do triangle = 1, mesh%data%n_triangles
            centroid = sum(mesh%data%vertices(:, &
                mesh%data%triangles(:, triangle)), dim=2) / 3.0_dp
            curl_values(triangle) = &
                reluctivity(centroid(1), centroid(2)) * centroid(1)
            mass_values(triangle) = &
                reluctivity(centroid(1), centroid(2)) / centroid(1)
        end do
        curl_coefficient = cell_coefficient(curl_values)
        mass_coefficient = cell_coefficient(mass_values)
        bilinear_form = ( &
            curl_coefficient * &
            inner(curl(trial_field), curl(test_field)) + &
            mass_coefficient * inner(trial_field, test_field)) * dx
        source = vector_function(space)
        source%values = 0.0_dp
        linear_form = inner(source, test_field) * dx

        allocate(boundary_values(space%ndof))
        boundary_values = 0.0_dp
        do boundary_index = 1, size(mesh%data%boundary_edges)
            edge = mesh%data%boundary_edges(boundary_index)
            degree_of_freedom = mesh%data%edge_to_dof(edge) + 1
            vertex_a = mesh%data%edges(1, edge)
            vertex_b = mesh%data%edges(2, edge)
            boundary_values(degree_of_freedom) = boundary_moment( &
                mesh%data%vertices(:, vertex_a), &
                mesh%data%vertices(:, vertex_b))
        end do
        boundary_condition = vector_bc_edge_moments( &
            space, boundary_values, "tangential")
        solution = vector_function(space)
        call solve( &
            bilinear_form == linear_form, solution, boundary_condition, &
            "direct")

        allocate(field_dofs(space%ndof))
        field_dofs = cmplx(solution%values(:, 1), 0.0_dp, dp)
        relative_error = magnetic_field_error(mesh%data, field_dofs)
    end function solve_compiled_paper_case

    real(dp) function solve_paper_case(divisions) result(relative_error)
        integer, intent(in) :: divisions

        type(mesh_2d_t) :: mesh
        complex(dp), allocatable :: field_dofs(:)
        real(dp), allocatable :: matrix(:, :), rhs(:), solution(:)
        real(dp), allocatable :: prescribed(:)
        integer :: boundary_dof, boundary_index, edge, interior_count
        integer :: vertex_a, vertex_b

        call mesh%create_rectangular( &
            divisions + 1, divisions + 1, &
            0.0_dp, 1.0_dp, -0.5_dp, 0.5_dp)
        call mesh%build_edge_connectivity()
        call mesh%build_edge_dof_numbering()

        allocate(matrix(mesh%n_edges, mesh%n_edges))
        allocate(prescribed(mesh%n_edges))
        call assemble_nedelec_weighted( &
            mesh, curl_weight, mass_weight, 7, matrix)
        prescribed = 0.0_dp
        do boundary_index = 1, size(mesh%boundary_edges)
            edge = mesh%boundary_edges(boundary_index)
            boundary_dof = mesh%edge_to_dof(edge) + 1
            vertex_a = mesh%edges(1, edge)
            vertex_b = mesh%edges(2, edge)
            prescribed(boundary_dof) = boundary_moment( &
                mesh%vertices(:, vertex_a), mesh%vertices(:, vertex_b))
        end do

        interior_count = mesh%n_interior_dofs
        allocate(rhs(interior_count), solution(interior_count))
        rhs = -matmul( &
            matrix(:interior_count, interior_count + 1:), &
            prescribed(interior_count + 1:))
        call solve_dense(matrix(:interior_count, :interior_count), &
            rhs, solution)

        allocate(field_dofs(mesh%n_edges))
        field_dofs = cmplx(prescribed, 0.0_dp, dp)
        field_dofs(:interior_count) = cmplx(solution, 0.0_dp, dp)
        relative_error = magnetic_field_error(mesh, field_dofs)
    end function solve_paper_case

    pure real(dp) function boundary_moment(first, second) result(moment)
        real(dp), intent(in) :: first(2), second(2)

        real(dp) :: radial_average

        radial_average = 0.5_dp * (first(1) + second(1))
        moment = -radial_average * (second(2) - first(2))
    end function boundary_moment

    real(dp) function magnetic_field_error(mesh, dofs) result(error)
        type(mesh_2d_t), intent(in) :: mesh
        complex(dp), intent(in) :: dofs(:)

        type(gauss_quadrature_triangle_t) :: quadrature
        complex(dp) :: curl_value, field_value(2)
        real(dp) :: det_jacobian, error_squared, exact_squared
        real(dp) :: point(2), vertex_a(2), vertex_b(2), vertex_c(2)
        real(dp) :: weight
        integer :: point_index, triangle

        quadrature = get_gauss_quadrature_triangle(7)
        error_squared = 0.0_dp
        exact_squared = 0.0_dp
        do triangle = 1, mesh%n_triangles
            vertex_a = mesh%vertices(:, mesh%triangles(1, triangle))
            vertex_b = mesh%vertices(:, mesh%triangles(2, triangle))
            vertex_c = mesh%vertices(:, mesh%triangles(3, triangle))
            det_jacobian = &
                (vertex_b(1) - vertex_a(1)) * &
                (vertex_c(2) - vertex_a(2)) - &
                (vertex_b(2) - vertex_a(2)) * &
                (vertex_c(1) - vertex_a(1))
            do point_index = 1, quadrature%n_points
                point = vertex_a + &
                    quadrature%xi(point_index) * (vertex_b - vertex_a) + &
                    quadrature%eta(point_index) * (vertex_c - vertex_a)
                weight = det_jacobian * quadrature%weights(point_index)
                call evaluate_nedelec_field_2d( &
                    mesh, triangle, point(1), point(2), dofs, &
                    field_value, curl_value)
                error_squared = error_squared + weight * ( &
                    (-real(field_value(2), dp) - &
                    analytical_radial_field(point(1)))**2 + &
                    real(field_value(1), dp)**2)
                exact_squared = exact_squared + weight * &
                    analytical_radial_field(point(1))**2
            end do
        end do
        call quadrature%destroy()
        error = sqrt(error_squared / exact_squared)
    end function magnetic_field_error

    pure real(dp) function analytical_radial_field(radial) result(value)
        real(dp), intent(in) :: radial

        if (radial <= 0.4_dp) then
            value = 20000.0_dp * radial / 128927.0_dp
        else if (radial < 0.5_dp) then
            value = 510000.0_dp * radial / 128927.0_dp - &
                78400.0_dp / (128927.0_dp * radial)
        else
            value = 106436.0_dp * radial / 128927.0_dp + &
                22491.0_dp / (128927.0_dp * radial)
        end if
    end function analytical_radial_field

    pure real(dp) function reluctivity(radial, vertical) result(value)
        real(dp), intent(in) :: radial, vertical

        if (radial >= 0.4_dp .and. radial < 0.5_dp) then
            value = 1.0_dp / 50.0_dp
        else
            value = 1.0_dp
        end if
        associate (unused_vertical => vertical)
            if (kind(unused_vertical) /= dp) error stop
        end associate
    end function reluctivity

    pure real(dp) function curl_weight(radial, vertical) result(value)
        real(dp), intent(in) :: radial, vertical

        value = reluctivity(radial, vertical) * radial
    end function curl_weight

    pure real(dp) function mass_weight(radial, vertical) result(value)
        real(dp), intent(in) :: radial, vertical

        value = reluctivity(radial, vertical) / radial
    end function mass_weight

    subroutine solve_dense(matrix, rhs, solution)
        real(dp), intent(in) :: matrix(:, :), rhs(:)
        real(dp), intent(out) :: solution(:)

        real(dp), allocatable :: matrix_work(:, :), rhs_work(:, :)
        integer, allocatable :: pivots(:)
        integer :: info, system_size

        interface
            subroutine dgesv(n, nrhs, a, lda, ipiv, b, ldb, info)
                import :: dp
                integer, intent(in) :: n, nrhs, lda, ldb
                real(dp), intent(inout) :: a(lda, *), b(ldb, *)
                integer, intent(out) :: ipiv(*), info
            end subroutine dgesv
        end interface

        system_size = size(rhs)
        allocate(matrix_work(system_size, system_size))
        allocate(rhs_work(system_size, 1), pivots(system_size))
        matrix_work = matrix
        rhs_work(:, 1) = rhs
        call dgesv( &
            system_size, 1, matrix_work, system_size, pivots, &
            rhs_work, system_size, info)
        if (info /= 0) error stop "Paper magnetic solve failed"
        solution = rhs_work(:, 1)
    end subroutine solve_dense

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_paper_magnetic_t1
