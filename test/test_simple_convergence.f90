program test_simple_convergence
    use fortfem_kinds, only: dp
    use fortfem_mesh_2d
    use fortfem_hcurl_space
    use fortfem_l2_projection
    implicit none

    type(mesh_2d_t) :: mesh
    type(hcurl_space_t) :: space
    real(dp), allocatable :: coeff(:)
    integer :: n_dofs, level
    
    print *, "Simple Convergence Test - Constant Field"
    print *, "========================================"
    print *, ""
    print *, "Level    h        DOFs    L2 Error     Max Coeff"
    print *, "----------------------------------------------"
    
    do level = 1, 3
        ! Create mesh with increasing resolution
        call mesh%create_rectangular(2 + level, 2 + level, 0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp)
        call mesh%build_edge_connectivity()
        call mesh%build_edge_dof_numbering()
        
        call space%init(mesh)
        n_dofs = space%get_n_dofs()
        allocate(coeff(n_dofs))
        
        ! Project constant field E = [1, 0]
        call project_edge_dof_analytical(mesh, space, constant_field, coeff)
        
        ! Compute simple L2 error
        call test_projection_accuracy(mesh, space, coeff, level)
        
        deallocate(coeff)
        call space%destroy()
        call mesh%destroy()
    end do

contains

    subroutine constant_field(x, y, E)
        real(dp), intent(in) :: x, y
        real(dp), intent(out) :: E(2)
        ! Simple constant field - should be representable exactly
        E(1) = 1.0_dp
        E(2) = 0.0_dp
    end subroutine

    subroutine test_projection_accuracy(mesh, space, coeff, level)
        type(mesh_2d_t), intent(in) :: mesh
        type(hcurl_space_t), intent(in) :: space
        real(dp), intent(in) :: coeff(:)
        integer, intent(in) :: level
        
        real(dp) :: E_computed(2), E_exact(2), error
        real(dp) :: xi, eta, h, max_coeff
        integer :: t
        
        h = 1.0_dp / real(1 + level, dp)
        max_coeff = maxval(abs(coeff))
        
        error = 0.0_dp
        
        ! Test at center of first triangle
        xi = 1.0_dp/3.0_dp
        eta = 1.0_dp/3.0_dp
        t = 1
        
        call constant_field(0.5_dp, 0.5_dp, E_exact)
        call space%evaluate_at_point(t, xi, eta, coeff, E_computed)
        
        error = sqrt(sum((E_computed - E_exact)**2))
        
        write(*, '(I5, 4(F10.6, 2X))') level, h, real(size(coeff)), error, max_coeff
        
        ! Debug info for first level
        if (level == 1) then
            print *, "Debug info for level 1:"
            print *, "  Exact:    ", E_exact
            print *, "  Computed: ", E_computed
            print *, "  Coeffs:   ", coeff(1:min(8, size(coeff)))
        end if
    end subroutine

end program test_simple_convergence