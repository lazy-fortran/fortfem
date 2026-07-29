program test_compiled_vector_form_solve
    use check, only: check_condition, check_summary
    use fortfem_api
    use fortfem_kinds, only: dp
    implicit none

    type(form_expr_t) :: bilinear_form, linear_form
    type(mesh_t) :: mesh
    type(vector_bc_t) :: boundary_condition
    type(vector_function_space_t) :: space
    type(vector_function_t) :: source, solution_one, solution_four
    type(vector_function_t) :: solution_boundary, solution_moments
    type(vector_function_t) :: solution_scaled
    type(vector_test_function_t) :: test_field
    type(vector_trial_function_t) :: trial_field
    real(dp), allocatable :: prescribed_moments(:), reference(:)
    real(dp) :: boundary_error, edge_vector(2), expected_value
    integer :: boundary_index, degree_of_freedom, edge
    logical :: all_passed

    all_passed = .true.
    mesh = unit_square_mesh(4)
    space = vector_function_space(mesh, "Nedelec", 1)
    trial_field = vector_trial_function(space)
    test_field = vector_test_function(space)
    boundary_condition = vector_bc( &
        space, [0.0_dp, 0.0_dp], "tangential")

    bilinear_form = ( &
        inner(curl(trial_field), curl(test_field)) + &
        inner(trial_field, test_field)) * dx
    source = vector_function(space)
    source%values(:, 1) = 1.0_dp
    source%values(:, 2) = 0.0_dp
    linear_form = inner(source, test_field) * dx
    solution_one = vector_function(space)
    call solve( &
        bilinear_form == linear_form, solution_one, boundary_condition, &
        "direct")
    allocate(reference(space%ndof))
    reference = solution_one%values(:, 1)

    source%values(:, 1) = 4.0_dp
    linear_form = inner(source, test_field) * dx
    solution_four = vector_function(space)
    call solve( &
        bilinear_form == linear_form, solution_four, boundary_condition, &
        "direct")
    call record_condition(maxval(abs(reference)) > 1.0e-8_dp .and. &
        maxval(abs(solution_four%values(:, 1) - 4.0_dp * reference)) < &
        3.0e-12_dp, "Compiled vector solve scales with its source field")

    bilinear_form = 2.0_dp * ( &
        inner(curl(trial_field), curl(test_field)) + &
        inner(trial_field, test_field)) * dx
    source%values(:, 1) = 1.0_dp
    linear_form = inner(source, test_field) * dx
    solution_scaled = vector_function(space)
    call solve( &
        bilinear_form == linear_form, solution_scaled, boundary_condition, &
        "direct")
    call record_condition(maxval(abs( &
        solution_scaled%values(:, 1) - 0.5_dp * reference)) < 3.0e-12_dp, &
        "Compiled vector solve preserves its bilinear coefficient")

    source%values = 0.0_dp
    linear_form = inner(source, test_field) * dx
    boundary_condition = vector_bc( &
        space, [1.0_dp, 0.0_dp], "tangential")
    solution_boundary = vector_function(space)
    call solve( &
        bilinear_form == linear_form, solution_boundary, boundary_condition, &
        "direct")
    boundary_error = 0.0_dp
    do boundary_index = 1, size(mesh%data%boundary_edges)
        edge = mesh%data%boundary_edges(boundary_index)
        degree_of_freedom = mesh%data%edge_to_dof(edge) + 1
        edge_vector = mesh%data%vertices(:, mesh%data%edges(2, edge)) - &
            mesh%data%vertices(:, mesh%data%edges(1, edge))
        expected_value = edge_vector(1)
        boundary_error = max(boundary_error, abs( &
            solution_boundary%values(degree_of_freedom, 1) - expected_value))
    end do
    call record_condition(boundary_error < 2.0e-14_dp, &
        "Constant tangential data has the exact boundary edge moments")

    allocate(prescribed_moments(space%ndof))
    prescribed_moments = 0.0_dp
    do boundary_index = 1, size(mesh%data%boundary_edges)
        edge = mesh%data%boundary_edges(boundary_index)
        degree_of_freedom = mesh%data%edge_to_dof(edge) + 1
        edge_vector = mesh%data%vertices(:, mesh%data%edges(2, edge)) - &
            mesh%data%vertices(:, mesh%data%edges(1, edge))
        prescribed_moments(degree_of_freedom) = &
            0.5_dp * sum(mesh%data%vertices(1, &
            mesh%data%edges(:, edge))) * edge_vector(2)
    end do
    boundary_condition = vector_bc_edge_moments( &
        space, prescribed_moments, "tangential")
    prescribed_moments = -100.0_dp
    solution_moments = vector_function(space)
    call solve( &
        bilinear_form == linear_form, solution_moments, boundary_condition, &
        "direct")
    boundary_error = 0.0_dp
    do boundary_index = 1, size(mesh%data%boundary_edges)
        edge = mesh%data%boundary_edges(boundary_index)
        degree_of_freedom = mesh%data%edge_to_dof(edge) + 1
        edge_vector = mesh%data%vertices(:, mesh%data%edges(2, edge)) - &
            mesh%data%vertices(:, mesh%data%edges(1, edge))
        expected_value = 0.5_dp * sum(mesh%data%vertices(1, &
            mesh%data%edges(:, edge))) * edge_vector(2)
        boundary_error = max(boundary_error, abs( &
            solution_moments%values(degree_of_freedom, 1) - expected_value))
    end do
    call record_condition(boundary_error < 2.0e-14_dp, &
        "Owned nonconstant boundary moments are imposed exactly")

    call check_summary("Compiled vector form solve")
    if (.not. all_passed) error stop 1

contains

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_compiled_vector_form_solve
