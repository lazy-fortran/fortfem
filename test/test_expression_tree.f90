program test_expression_tree
    ! Test the expression tree system for forms
    use fortfem_api
    use check
    implicit none
    
    type(mesh_t) :: mesh
    type(function_space_t) :: Vh
    type(trial_function_t) :: u
    type(test_function_t) :: v
    type(expression_t) :: grad_u, grad_v, bilinear_form
    character(len=256) :: assembly_code
    
    print *, "Testing expression tree system"
    print *, "=============================="
    
    ! Create mesh and function space
    mesh = unit_square_mesh(4)
    Vh = function_space(mesh, "Lagrange", 1)
    u = trial_function(Vh)
    v = test_function(Vh)
    
    ! Test gradient expressions
    grad_u = grad(u)
    grad_v = grad(v)
    
    call check_condition(grad_u%operation%op_type == EXPR_GRAD, "grad(u) is gradient operation")
    call check_condition(grad_v%operation%op_type == EXPR_GRAD, "grad(v) is gradient operation")
    call check_condition(grad_u%tensor_rank == 1, "grad(u) has tensor rank 1")
    call check_condition(trim(grad_u%function_type) == "trial", "grad(u) is trial function")
    call check_condition(trim(grad_v%function_type) == "test", "grad(v) is test function")
    
    ! Test inner product
    bilinear_form = inner(grad_u, grad_v)
    
    call check_condition(bilinear_form%operation%op_type == EXPR_INNER, "inner() creates inner product")
    call check_condition(bilinear_form%tensor_rank == 0, "inner(grad,grad) is scalar")
    call check_condition(trim(bilinear_form%function_type) == "bilinear", "Result is bilinear form")
    
    ! Test expression compilation
    assembly_code = expr_compile(bilinear_form, "bilinear")
    call check_condition(len_trim(assembly_code) > 0, "Expression compiles to assembly code")
    
    print *, "Compiled assembly code:"
    print *, trim(assembly_code)
    
    ! Test operator overloading
    bilinear_form = inner(grad_u, grad_v) * dx
    call check_condition(bilinear_form%operation%op_type == EXPR_INNER, "Operator * preserves form")
    
    ! Clean up
    call grad_u%destroy()
    call grad_v%destroy()
    call bilinear_form%destroy()
    call mesh%destroy()
    call Vh%destroy()
    
    call check_summary("Expression Tree Test")
    
end program test_expression_tree