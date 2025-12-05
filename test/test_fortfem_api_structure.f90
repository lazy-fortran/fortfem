program test_fortfem_api_structure
    use fortfem_kinds, only: dp
    use fortfem_api_mesh, only: mesh_t, unit_square_mesh
    use fortfem_api_spaces, only: function_space_t, function_space,        &
        trial_function_t, test_function_t, function_t, dirichlet_bc_t,     &
        trial_function, test_function, function, constant, dirichlet_bc
    use fortfem_api_forms, only: form_expr_t, inner, grad, dx,             &
        operator(*), operator(==)
    use fortfem_api_solvers, only: solve
    use fortfem_api, only: mesh_t_api => mesh_t,                           &
        function_space_t_api => function_space_t
    use check, only: check_condition, check_summary
    implicit none

    type(mesh_t) :: mesh_local
    type(mesh_t_api) :: mesh_api
    type(function_space_t) :: V_local
    type(function_space_t_api) :: V_api
    type(trial_function_t) :: u
    type(test_function_t) :: v
    type(function_t) :: f, uh
    type(form_expr_t) :: a, L
    type(dirichlet_bc_t) :: bc
    logical :: same_mesh_size, same_space_size, nontrivial_solution

    write(*,*) "Testing fortfem_api split modules and re-exports..."

    mesh_local = unit_square_mesh(4)
    mesh_api = unit_square_mesh(4)

    V_local = function_space(mesh_local, "Lagrange", 1)
    V_api = function_space(mesh_api, "Lagrange", 1)

    same_mesh_size = mesh_local%data%n_vertices == mesh_api%data%n_vertices
    same_space_size = V_local%ndof == V_api%ndof

    call check_condition(same_mesh_size,                                   &
        "fortfem_api_mesh and fortfem_api provide consistent mesh_t")
    call check_condition(same_space_size,                                   &
        "fortfem_api_spaces and fortfem_api provide consistent spaces")

    u = trial_function(V_local)
    v = test_function(V_local)
    f = constant(1.0_dp)

    a = inner(grad(u), grad(v)) * dx
    L = f * v * dx

    bc = dirichlet_bc(V_local, 0.0_dp)
    uh = function(V_local)

    call solve(a == L, uh, bc)

    nontrivial_solution = any(abs(uh%values) > 0.0_dp)

    call check_condition(nontrivial_solution,                               &
        "solve from split modules produces a nontrivial solution")

    call check_summary("fortfem_api structure")

end program test_fortfem_api_structure
