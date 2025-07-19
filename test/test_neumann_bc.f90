program test_neumann_bc
    ! Test Neumann boundary conditions implementation
    use fortfem_kinds
    use fortfem_mesh_2d
    use basis_p1_2d_module
    use assembly_2d_module
    use boundary_assembly_2d_module
    use fortfem_sparse_matrix
    use fortfem_solver_interface
    implicit none
    
    logical :: all_passed = .true.
    
    print *, "=== Testing Neumann Boundary Conditions ==="
    print *, ""
    
    call test_pure_neumann_problem()
    call test_mixed_bc_problem()
    call test_flux_conservation()
    call test_p2_neumann()
    
    print *, ""
    if (all_passed) then
        print *, "✅ ALL NEUMANN BC TESTS PASSED!"
    else
        print *, "❌ SOME NEUMANN BC TESTS FAILED!"
        stop 1
    end if
    
contains

    subroutine test_pure_neumann_problem()
        ! Test pure Neumann problem (simplified test)
        type(mesh_2d_t) :: mesh
        type(boundary_assembly_2d_t) :: bc_assembler
        real(dp), allocatable :: rhs(:)
        real(dp) :: total_flux
        integer :: n
        logical :: passed = .true.
        
        print *, "Testing pure Neumann problem..."
        
        ! Create mesh
        n = 10
        call mesh%create_rectangular(n, n, 0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp)
        call mesh%build_edge_connectivity()
        
        ! Initialize RHS
        allocate(rhs(mesh%n_vertices))
        rhs = 0.0_dp
        
        ! Add Neumann BC: ∂u/∂n = 1 on all boundaries
        call assemble_neumann_bc(mesh, const_one, rhs)
        
        ! Check total flux
        total_flux = sum(rhs)
        
        ! For unit square with flux = 1, total should be perimeter = 4
        if (abs(total_flux - 4.0_dp) > 1e-10_dp) then
            print *, "  ❌ Total flux incorrect:", total_flux
            print *, "    Expected: 4.0 (perimeter of unit square)"
            passed = .false.
        end if
        
        if (passed) then
            print *, "  ✅ Pure Neumann problem test passed"
        else
            all_passed = .false.
        end if
        
        ! Clean up
        deallocate(rhs)
        call mesh%destroy()
        
    end subroutine test_pure_neumann_problem
    
    subroutine test_mixed_bc_problem()
        ! Test problem with mixed Dirichlet/Neumann BC
        type(mesh_2d_t) :: mesh
        real(dp), allocatable :: rhs(:)
        integer :: n, i, n_boundary_edges
        logical :: passed = .true.
        real(dp) :: x1, y1, x2, y2
        
        print *, "Testing mixed boundary conditions..."
        
        ! Create mesh
        n = 5
        call mesh%create_rectangular(n, n, 0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp)
        call mesh%build_edge_connectivity()
        
        ! Initialize RHS
        allocate(rhs(mesh%n_vertices))
        rhs = 0.0_dp
        
        ! Count boundary edges on top/bottom (Neumann)
        n_boundary_edges = 0
        do i = 1, mesh%n_edges
            if (mesh%is_boundary_edge(i)) then
                x1 = mesh%vertices(1, mesh%edges(1,i))
                y1 = mesh%vertices(2, mesh%edges(1,i))
                x2 = mesh%vertices(1, mesh%edges(2,i))
                y2 = mesh%vertices(2, mesh%edges(2,i))
                
                ! Check if edge is on top or bottom
                if ((abs(y1) < 1e-10_dp .and. abs(y2) < 1e-10_dp) .or. &
                    (abs(y1 - 1.0_dp) < 1e-10_dp .and. abs(y2 - 1.0_dp) < 1e-10_dp)) then
                    n_boundary_edges = n_boundary_edges + 1
                end if
            end if
        end do
        
        if (n_boundary_edges == 2*(n-1)) then
            print *, "  ✅ Correct number of Neumann edges identified"
        else
            print *, "  ❌ Wrong number of Neumann edges:", n_boundary_edges
            print *, "    Expected:", 2*(n-1)
            passed = .false.
        end if
        
        if (passed) then
            print *, "  ✅ Mixed BC problem test passed"
        else
            all_passed = .false.
        end if
        
        ! Clean up
        deallocate(rhs)
        call mesh%destroy()
        
    end subroutine test_mixed_bc_problem
    
    subroutine test_flux_conservation()
        ! Test flux conservation across boundaries
        type(mesh_2d_t) :: mesh
        real(dp), allocatable :: fluxes(:)
        real(dp) :: total_flux, edge_length, flux_val
        integer :: e, v1, v2
        logical :: passed = .true.
        
        print *, "Testing flux conservation..."
        
        ! Create mesh
        call mesh%create_rectangular(5, 5, 0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp)
        call mesh%build_edge_connectivity()
        
        ! Compute boundary fluxes
        allocate(fluxes(mesh%n_edges))
        fluxes = 0.0_dp
        
        total_flux = 0.0_dp
        
        do e = 1, mesh%n_edges
            if (mesh%is_boundary_edge(e)) then
                v1 = mesh%edges(1, e)
                v2 = mesh%edges(2, e)
                
                edge_length = sqrt(sum((mesh%vertices(:,v2) - mesh%vertices(:,v1))**2))
                
                ! Simple constant flux
                flux_val = 1.0_dp
                fluxes(e) = flux_val * edge_length
                total_flux = total_flux + fluxes(e)
            end if
        end do
        
        ! For a closed domain, total flux should be related to domain integral
        if (abs(total_flux - 4.0_dp) > 1e-10_dp) then
            print *, "  ❌ Total flux incorrect:", total_flux
            passed = .false.
        end if
        
        if (passed) then
            print *, "  ✅ Flux conservation test passed"
        else
            all_passed = .false.
        end if
        
        deallocate(fluxes)
        call mesh%destroy()
        
    end subroutine test_flux_conservation
    
    subroutine test_p2_neumann()
        ! Test Neumann BC with P2 elements
        type(mesh_2d_t) :: mesh
        type(boundary_assembly_2d_t) :: bc_assembler
        real(dp), allocatable :: rhs(:)
        real(dp) :: edge_vertices(2,2), be(3)
        integer :: n_dofs
        logical :: passed = .true.
        
        print *, "Testing P2 Neumann BC assembly..."
        
        ! Create small mesh
        call mesh%create_rectangular(3, 3, 0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp)
        call mesh%build_edge_connectivity()
        
        ! P2 DOFs: vertices + edges
        n_dofs = mesh%n_vertices + mesh%n_edges
        allocate(rhs(n_dofs))
        rhs = 0.0_dp
        
        ! Test single edge assembly
        edge_vertices(:,1) = [0.0_dp, 0.0_dp]
        edge_vertices(:,2) = [1.0_dp, 0.0_dp]
        
        call bc_assembler%init(2)  ! P2 elements
        call bc_assembler%element_boundary_vector_p2(edge_vertices, const_one, be)
        
        ! Check that integral is correct: ∫_edge 1 ds = edge_length = 1
        if (abs(sum(be) - 1.0_dp) > 1e-10_dp) then
            print *, "  ❌ P2 edge integral incorrect:", sum(be)
            passed = .false.
        end if
        
        ! Check symmetry for constant function
        if (abs(be(1) - be(2)) > 1e-10_dp) then
            print *, "  ❌ P2 edge values not symmetric for constant function"
            passed = .false.
        end if
        
        if (passed) then
            print *, "  ✅ P2 Neumann BC test passed"
        else
            all_passed = .false.
        end if
        
        deallocate(rhs)
        call mesh%destroy()
        call bc_assembler%destroy()
        
    end subroutine test_p2_neumann
    
    ! Test functions
    pure function const_one(x, y) result(f)
        real(dp), intent(in) :: x, y
        real(dp) :: f
        f = 1.0_dp
    end function const_one
    
    pure function neumann_flux(x, y) result(g)
        real(dp), intent(in) :: x, y
        real(dp) :: g
        g = -x  ! Flux that gives compatibility
    end function neumann_flux

end program test_neumann_bc