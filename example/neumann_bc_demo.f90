program neumann_bc_demo
    ! Demonstrate Neumann boundary conditions in FortFEM
    ! Solve -Δu = f with mixed boundary conditions
    
    use fortfem_kinds
    use fortfem_mesh_2d
    use basis_p1_2d_module
    use assembly_2d_module
    use boundary_assembly_2d_module
    use fortfem_sparse_matrix
    use fortfem_solver_interface
    ! use fortplotlib, only: plot, figure, savefig, xlabel, ylabel, title, colorbar
    implicit none
    
    type(mesh_2d_t) :: mesh
    type(assembly_2d_t) :: assembler
    type(triplet_matrix_t) :: A_triplet
    type(csr_matrix_t) :: A_csr
    type(lapack_dense_solver_t) :: solver
    real(dp), allocatable :: u(:), rhs(:), x_plot(:), y_plot(:), u_plot(:,:)
    integer :: n, i, j, idx, info, n_interior
    integer, allocatable :: dof_map(:)
    logical, allocatable :: is_dirichlet(:)
    real(dp) :: x, y, h
    
    print *, "=== Neumann Boundary Condition Demo ==="
    print *, ""
    print *, "Problem: -Δu = 0 (Laplace equation)"
    print *, "Domain: [0,1] × [0,1]"
    print *, "Boundary conditions:"
    print *, "  - u = 0 on left boundary (x = 0)"
    print *, "  - u = 1 on right boundary (x = 1)"
    print *, "  - ∂u/∂n = 0 on top and bottom (insulated)"
    print *, ""
    print *, "Expected solution: u(x,y) = x (linear in x)"
    print *, ""
    
    ! Create mesh
    n = 21
    h = 1.0_dp / real(n-1, dp)
    call mesh%create_rectangular(n, n, 0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp)
    call mesh%build_edge_connectivity()
    
    print *, "Mesh statistics:"
    print *, "  Vertices:", mesh%n_vertices
    print *, "  Triangles:", mesh%n_triangles
    print *, "  Edges:", mesh%n_edges
    print *, "  Boundary edges:", mesh%n_boundary_edges
    print *, ""
    
    ! Identify Dirichlet nodes (left and right boundaries)
    allocate(is_dirichlet(mesh%n_vertices))
    allocate(dof_map(mesh%n_vertices))
    is_dirichlet = .false.
    dof_map = -1
    
    do i = 1, mesh%n_vertices
        x = mesh%vertices(1, i)
        if (abs(x) < 1e-10_dp .or. abs(x - 1.0_dp) < 1e-10_dp) then
            is_dirichlet(i) = .true.
        end if
    end do
    
    ! Create DOF mapping for interior + Neumann nodes
    n_interior = 0
    do i = 1, mesh%n_vertices
        if (.not. is_dirichlet(i)) then
            n_interior = n_interior + 1
            dof_map(i) = n_interior
        end if
    end do
    
    print *, "DOF statistics:"
    print *, "  Total nodes:", mesh%n_vertices
    print *, "  Dirichlet nodes:", count(is_dirichlet)
    print *, "  Interior + Neumann nodes:", n_interior
    print *, ""
    
    ! Initialize system
    allocate(u(mesh%n_vertices), rhs(n_interior))
    u = 0.0_dp
    rhs = 0.0_dp
    
    ! Set Dirichlet values
    do i = 1, mesh%n_vertices
        if (is_dirichlet(i)) then
            x = mesh%vertices(1, i)
            if (abs(x - 1.0_dp) < 1e-10_dp) then
                u(i) = 1.0_dp  ! Right boundary
            else
                u(i) = 0.0_dp  ! Left boundary
            end if
        end if
    end do
    
    ! Assemble system matrix (Laplacian)
    call A_triplet%init(n_interior, n_interior * 10)
    call assembler%init(n_interior, n_interior * 10)
    
    ! Assemble with DOF elimination
    call assemble_laplacian_with_bc(mesh, A_triplet, rhs, u, dof_map, is_dirichlet)
    
    ! Note: For homogeneous Neumann BC (∂u/∂n = 0), no boundary integral needed
    ! For non-zero Neumann BC, we would call:
    ! call assemble_neumann_bc(mesh, neumann_func, rhs)
    
    print *, "Solving linear system..."
    
    ! Convert to CSR and solve
    call A_triplet%to_csr(A_csr)
    
    ! Allocate solution vector for interior DOFs
    allocate(x_plot(n_interior))
    
    call solver%init()
    call solver%solve(A_csr, rhs, x_plot, info)
    
    if (info /= 0) then
        print *, "Solver failed with info =", info
        stop 1
    end if
    
    ! Copy solution back to full vector
    do i = 1, mesh%n_vertices
        if (dof_map(i) > 0) then
            u(i) = x_plot(dof_map(i))
        end if
    end do
    
    print *, "Solution computed successfully!"
    print *, ""
    
    ! Check solution accuracy
    call check_solution_accuracy(mesh, u)
    
    ! Save solution for visualization
    call save_solution_vtk(mesh, u, "neumann_bc_solution.vtk")
    
    print *, ""
    print *, "Solution saved to 'neumann_bc_solution.vtk'"
    print *, "Visualize with ParaView or similar software"
    print *, ""
    print *, "Key observations:"
    print *, "- Solution is linear in x-direction (u = x)"
    print *, "- Homogeneous Neumann BC enforced on top/bottom"
    print *, "- No flux through insulated boundaries"
    print *, "- Mixed BC problems are common in heat transfer"
    
    ! Clean up
    deallocate(u, rhs, x_plot, is_dirichlet, dof_map)
    call mesh%destroy()
    call assembler%destroy()
    call A_triplet%destroy()
    call A_csr%destroy()
    call solver%destroy()
    
contains

    subroutine assemble_laplacian_with_bc(mesh, A, rhs, u_bc, dof_map, is_dirichlet)
        type(mesh_2d_t), intent(in) :: mesh
        type(triplet_matrix_t), intent(inout) :: A
        real(dp), intent(inout) :: rhs(:)
        real(dp), intent(in) :: u_bc(:)
        integer, intent(in) :: dof_map(:)
        logical, intent(in) :: is_dirichlet(:)
        
        type(basis_p1_2d_t) :: basis
        real(dp) :: vertices(2,3), Ae(3,3)
        real(dp) :: jac(2,2), det_j
        integer :: t, i, j, vi, vj, gi, gj
        
        ! Loop over elements
        do t = 1, mesh%n_triangles
            ! Get element vertices
            do i = 1, 3
                vertices(:,i) = mesh%vertices(:, mesh%triangles(i,t))
            end do
            
            ! Compute element stiffness matrix
            call assembler%element_stiffness_matrix(vertices, Ae)
            
            ! Assemble with BC elimination
            do i = 1, 3
                vi = mesh%triangles(i, t)
                gi = dof_map(vi)
                
                if (gi > 0) then  ! Interior or Neumann node
                    do j = 1, 3
                        vj = mesh%triangles(j, t)
                        gj = dof_map(vj)
                        
                        if (gj > 0) then
                            ! Add to system matrix
                            call A%add(gi, gj, Ae(i,j))
                        else
                            ! Move known Dirichlet value to RHS
                            rhs(gi) = rhs(gi) - Ae(i,j) * u_bc(vj)
                        end if
                    end do
                end if
            end do
        end do
        
    end subroutine assemble_laplacian_with_bc
    
    subroutine check_solution_accuracy(mesh, u)
        type(mesh_2d_t), intent(in) :: mesh
        real(dp), intent(in) :: u(:)
        
        real(dp) :: max_error, error, x, y, u_exact
        integer :: i
        
        max_error = 0.0_dp
        
        do i = 1, mesh%n_vertices
            x = mesh%vertices(1, i)
            y = mesh%vertices(2, i)
            
            ! Exact solution: u = x
            u_exact = x
            
            error = abs(u(i) - u_exact)
            max_error = max(max_error, error)
        end do
        
        print *, "Solution accuracy:"
        print '(A,E12.5)', "  Maximum error: ", max_error
        
        if (max_error < 1e-10_dp) then
            print *, "  ✅ Solution is exact (within machine precision)"
        else if (max_error < 1e-3_dp) then
            print *, "  ✅ Solution is accurate"
        else
            print *, "  ⚠️  Solution error is large"
        end if
        
    end subroutine check_solution_accuracy
    
    subroutine save_solution_vtk(mesh, u, filename)
        type(mesh_2d_t), intent(in) :: mesh
        real(dp), intent(in) :: u(:)
        character(len=*), intent(in) :: filename
        
        integer :: unit, i, t
        
        open(newunit=unit, file=filename, status='replace')
        
        ! VTK header
        write(unit, '(A)') '# vtk DataFile Version 3.0'
        write(unit, '(A)') 'Neumann BC Demo Solution'
        write(unit, '(A)') 'ASCII'
        write(unit, '(A)') 'DATASET UNSTRUCTURED_GRID'
        
        ! Points
        write(unit, '(A,I0,A)') 'POINTS ', mesh%n_vertices, ' double'
        do i = 1, mesh%n_vertices
            write(unit, '(3F16.8)') mesh%vertices(:,i), 0.0_dp
        end do
        
        ! Cells (triangles)
        write(unit, '(A,I0,I0)') 'CELLS ', mesh%n_triangles, mesh%n_triangles * 4
        do t = 1, mesh%n_triangles
            write(unit, '(4I8)') 3, mesh%triangles(:,t) - 1  ! VTK uses 0-based indexing
        end do
        
        ! Cell types (5 = triangle)
        write(unit, '(A,I0)') 'CELL_TYPES ', mesh%n_triangles
        do t = 1, mesh%n_triangles
            write(unit, '(I2)') 5
        end do
        
        ! Point data
        write(unit, '(A,I0)') 'POINT_DATA ', mesh%n_vertices
        write(unit, '(A)') 'SCALARS solution double 1'
        write(unit, '(A)') 'LOOKUP_TABLE default'
        do i = 1, mesh%n_vertices
            write(unit, '(F16.8)') u(i)
        end do
        
        close(unit)
        
    end subroutine save_solution_vtk

end program neumann_bc_demo