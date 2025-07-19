program test_forms_api
    ! Test the new forms-based API structure
    use fortfem_api
    use check
    implicit none
    
    type(mesh_t) :: mesh
    type(function_space_t) :: Vh
    type(trial_function_t) :: u
    type(test_function_t) :: v
    type(function_t) :: uh, f
    type(dirichlet_bc_t) :: bc
    
    print *, "Testing forms-based API structure"
    print *, "================================="
    
    ! Test mesh creation
    mesh = unit_square_mesh(4, 4)
    call check_condition(mesh%data%n_vertices > 0, "Unit square mesh created")
    
    ! Test function space creation
    Vh = function_space(mesh, "Lagrange", 1)
    call check_condition(Vh%ndof > 0, "Function space created")
    call check_condition(trim(Vh%element_family) == "Lagrange", "Correct element family")
    call check_condition(Vh%degree == 1, "Correct degree")
    
    ! Test function creation
    u = trial_function(Vh)
    v = test_function(Vh)
    uh = function(Vh)
    f = constant(1.0_dp)
    
    call check_condition(associated(u%space), "Trial function created")
    call check_condition(associated(v%space), "Test function created")
    call check_condition(associated(uh%space), "Function created")
    call check_condition(allocated(uh%values), "Function values allocated")
    
    ! Test boundary condition creation
    bc = dirichlet_bc(Vh, 0.0_dp)
    call check_condition(associated(bc%space), "Boundary condition created")
    call check_condition(bc%value == 0.0_dp, "Correct BC value")
    
    ! Test constant function
    call check_condition(size(f%values) > 0, "Constant function created")
    
    ! Clean up
    call mesh%destroy()
    call Vh%destroy()
    call uh%destroy()
    call f%destroy()
    
    call check_summary("Forms API Structure Test")
    
end program test_forms_api