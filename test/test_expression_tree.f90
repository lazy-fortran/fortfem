program test_expression_tree
    ! Test the forms system with direct dx usage
    use fortfem_api
    use check
    implicit none
    
    type(mesh_t) :: mesh
    type(function_space_t) :: Vh
    type(trial_function_t) :: u
    type(test_function_t) :: v
    type(form_expr_t) :: grad_u, grad_v, a
    character(len=256) :: assembly_code
    
    print *, "Testing forms system with dx"
    print *, "============================"
    
    ! Create mesh and function space
    mesh = unit_square_mesh(4)
    Vh = function_space(mesh, "Lagrange", 1)
    u = trial_function(Vh)
    v = test_function(Vh)
    
    ! Test gradient expressions
    grad_u = grad(u)
    grad_v = grad(v)
    
    call check_condition(len_trim(grad_u%description) > 0, "grad(u) created")
    call check_condition(len_trim(grad_v%description) > 0, "grad(v) created")
    call check_condition(trim(grad_u%form_type) == "trial", "grad(u) is trial function")
    call check_condition(trim(grad_v%form_type) == "test", "grad(v) is test function")
    
    ! Test FEniCS-style syntax with dx
    a = inner(grad_u, grad_v) * dx
    
    call check_condition(len_trim(a%description) > 0, "Form with dx created")
    call check_condition(index(a%description, "dx") > 0, "Form includes dx measure")
    
    print *, "Form expression: ", trim(a%description)
    
    ! Test compilation
    assembly_code = compile_form(a)
    call check_condition(len_trim(assembly_code) > 0, "Form compiles to assembly code")
    
    print *, "Assembly code: ", trim(assembly_code)
    
    ! Clean up
    call grad_u%destroy()
    call grad_v%destroy()
    call a%destroy()
    call mesh%destroy()
    call Vh%destroy()
    
    call check_summary("Forms with dx Test")
    
end program test_expression_tree