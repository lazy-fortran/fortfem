program test_debug_projection
    use fortfem_kinds, only: dp
    use fortfem_mesh_2d
    use fortfem_hcurl_space
    use fortfem_l2_projection
    implicit none

    type(mesh_2d_t) :: mesh
    type(hcurl_space_t) :: space
    real(dp), allocatable :: coeff(:), coeff_simple(:)
    integer :: n_dofs
    
    ! Create simple 2x2 mesh
    call mesh%create_rectangular(3, 3, 0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp)
    call mesh%build_edge_connectivity()
    call mesh%build_edge_dof_numbering()
    
    call space%init(mesh)
    n_dofs = space%get_n_dofs()
    
    print *, "Mesh info:"
    print *, "  Triangles:", mesh%n_triangles
    print *, "  Edges:", mesh%n_edges
    print *, "  DOFs:", n_dofs
    
    allocate(coeff(n_dofs))
    allocate(coeff_simple(n_dofs))
    
    ! Test projection with E = [y, x]
    print *, ""
    print *, "Testing projection of E = [y, x]"
    
    ! Try L2 projection
    call project_l2_edge_element(mesh, space, analytical_solution, coeff)
    
    ! Try simple projection  
    call project_edge_dof_analytical(mesh, space, analytical_solution, coeff_simple)
    
    print *, "First 10 coefficients (L2 projection):"
    print *, coeff(1:min(10, n_dofs))
    
    print *, "First 10 coefficients (simple projection):"
    print *, coeff_simple(1:min(10, n_dofs))
    
    ! Test evaluation at a point
    call test_evaluation_at_point(mesh, space, coeff, coeff_simple)
    
    deallocate(coeff, coeff_simple)
    call space%destroy()
    call mesh%destroy()

contains

    subroutine analytical_solution(x, y, E)
        real(dp), intent(in) :: x, y
        real(dp), intent(out) :: E(2)
        E(1) = y
        E(2) = x
    end subroutine

    subroutine test_evaluation_at_point(mesh, space, coeff, coeff_simple)
        type(mesh_2d_t), intent(in) :: mesh
        type(hcurl_space_t), intent(in) :: space
        real(dp), intent(in) :: coeff(:), coeff_simple(:)
        
        real(dp) :: E_computed(2), E_simple(2), E_exact(2)
        real(dp) :: x_test, y_test, xi, eta
        
        ! Test at center of first triangle
        x_test = 0.5_dp
        y_test = 0.5_dp
        xi = 1.0_dp/3.0_dp
        eta = 1.0_dp/3.0_dp
        
        call analytical_solution(x_test, y_test, E_exact)
        call space%evaluate_at_point(1, xi, eta, coeff, E_computed)
        call space%evaluate_at_point(1, xi, eta, coeff_simple, E_simple)
        
        print *, ""
        print *, "Evaluation at triangle 1 center:"
        print *, "  Exact:     ", E_exact
        print *, "  L2 proj:   ", E_computed
        print *, "  Simple:    ", E_simple
        print *, "  L2 error:  ", sqrt(sum((E_computed - E_exact)**2))
        print *, "  Simple err:", sqrt(sum((E_simple - E_exact)**2))
    end subroutine

end program test_debug_projection