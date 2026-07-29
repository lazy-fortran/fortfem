program test_compiled_form_solve
    use check, only: check_condition, check_summary
    use fortfem_api
    use fortfem_kinds, only: dp
    implicit none

    type(dirichlet_bc_t) :: boundary_condition
    type(form_expr_t) :: bilinear_form, linear_form
    type(function_space_t) :: space
    type(function_t) :: source, solution_one, solution_four, solution_scaled
    type(mesh_t) :: mesh
    type(solver_options_t) :: options
    type(test_function_t) :: test_field
    type(trial_function_t) :: trial_field
    real(dp) :: center_one
    logical :: all_passed

    all_passed = .true.
    mesh = unit_square_mesh(3)
    space = function_space(mesh, "Lagrange", 1)
    trial_field = trial_function(space)
    test_field = test_function(space)
    boundary_condition = dirichlet_bc(space, 0.0_dp)
    options = solver_options(method="lapack_lu")

    bilinear_form = inner(grad(trial_field), grad(test_field)) * dx
    source = constant(1.0_dp)
    linear_form = source * test_field * dx
    solution_one = function(space)
    call solve( &
        bilinear_form == linear_form, solution_one, boundary_condition, options)
    center_one = solution_one%values(5)

    source = constant(4.0_dp)
    linear_form = source * test_field * dx
    solution_four = function(space)
    call solve( &
        bilinear_form == linear_form, solution_four, boundary_condition, options)
    call record_condition(center_one > 0.0_dp .and. &
        abs(solution_four%values(5) - 4.0_dp * center_one) < 2.0e-13_dp, &
        "Compiled solve scales exactly with the symbolic source")

    bilinear_form = &
        2.0_dp * inner(grad(trial_field), grad(test_field)) * dx
    source = constant(1.0_dp)
    linear_form = source * test_field * dx
    solution_scaled = function(space)
    call solve( &
        bilinear_form == linear_form, solution_scaled, boundary_condition, &
        options)
    call record_condition(abs(solution_scaled%values(5) - &
        0.5_dp * center_one) < 2.0e-13_dp, &
        "Compiled solve preserves the symbolic stiffness coefficient")

    call check_summary("Compiled form solve")
    if (.not. all_passed) error stop 1

contains

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_compiled_form_solve
