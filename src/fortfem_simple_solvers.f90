module fortfem_simple_solvers
    ! Simplified high-level solvers for the new API
    use fortfem_kinds
    use fortfem_mesh_2d
    use assembly_2d_module
    use function_space_module
    use fortfem_sparse_matrix
    use fortfem_solver_interface
    implicit none
    
    private
    public :: create_unit_square_mesh
    public :: assemble_poisson_2d
    public :: apply_zero_bc
    public :: solve_lapack_dense
    public :: write_vtk
    public :: compute_l2_error
    
contains

    subroutine create_unit_square_mesh(mesh, n)
        type(mesh_2d_t), intent(out) :: mesh
        integer, intent(in) :: n
        
        call mesh%create_rectangular(nx=n, ny=n, &
                                   x_min=0.0_dp, x_max=1.0_dp, &
                                   y_min=0.0_dp, y_max=1.0_dp)
        call mesh%build_connectivity()
        call mesh%find_boundary()
    end subroutine create_unit_square_mesh
    
    subroutine assemble_poisson_2d(mesh, A, f)
        type(mesh_2d_t), intent(in) :: mesh
        type(sparse_matrix_t), intent(out) :: A
        real(dp), allocatable, intent(out) :: f(:)
        
        type(function_space_t) :: V
        real(dp), allocatable :: A_dense(:,:)
        integer :: ndof
        
        ! Create P1 function space
        call create_P1_space(mesh, V)
        ndof = V%n_dofs
        
        allocate(A_dense(ndof, ndof))
        allocate(f(ndof))
        
        ! Assemble stiffness matrix and load vector
        call assemble_laplacian(V, A_dense)
        call assemble_mass_rhs(V, f, source_term)
        
        ! Convert to sparse
        call A%init_from_dense(A_dense)
        
        ! Clean up
        call V%destroy()
        deallocate(A_dense)
        
    contains
        pure real(dp) function source_term(x, y)
            real(dp), intent(in) :: x, y
            source_term = 2.0_dp * pi**2 * sin(pi*x) * sin(pi*y)
        end function source_term
    end subroutine assemble_poisson_2d
    
    subroutine apply_zero_bc(mesh, A, f)
        type(mesh_2d_t), intent(in) :: mesh
        type(sparse_matrix_t), intent(inout) :: A
        real(dp), intent(inout) :: f(:)
        
        integer :: i, j
        real(dp) :: x, y
        
        ! Apply homogeneous Dirichlet BC on boundary
        do i = 1, mesh%n_vertices
            x = mesh%vertices(1, i)
            y = mesh%vertices(2, i)
            
            if (abs(x) < 1e-10 .or. abs(x - 1.0_dp) < 1e-10 .or. &
                abs(y) < 1e-10 .or. abs(y - 1.0_dp) < 1e-10) then
                
                ! Zero out row
                do j = 1, A%n
                    call A%set(i, j, 0.0_dp)
                end do
                
                ! Set diagonal
                call A%set(i, i, 1.0_dp)
                
                ! Zero RHS
                f(i) = 0.0_dp
            end if
        end do
    end subroutine apply_zero_bc
    
    subroutine write_vtk(filename, mesh, u, vector)
        character(len=*), intent(in) :: filename
        type(mesh_2d_t), intent(in) :: mesh
        real(dp), intent(in) :: u(:)
        logical, intent(in), optional :: vector
        
        integer :: unit, i, j
        logical :: is_vector
        
        is_vector = .false.
        if (present(vector)) is_vector = vector
        
        open(newunit=unit, file=filename, status='replace')
        
        ! VTK header
        write(unit, '(a)') '# vtk DataFile Version 3.0'
        write(unit, '(a)') 'FortFEM Solution'
        write(unit, '(a)') 'ASCII'
        write(unit, '(a)') 'DATASET UNSTRUCTURED_GRID'
        
        ! Points
        write(unit, '(a,i0,a)') 'POINTS ', mesh%n_vertices, ' double'
        do i = 1, mesh%n_vertices
            write(unit, '(3f16.8)') mesh%vertices(1,i), mesh%vertices(2,i), 0.0_dp
        end do
        
        ! Cells
        write(unit, '(a,i0,i0)') 'CELLS ', mesh%n_triangles, mesh%n_triangles*4
        do i = 1, mesh%n_triangles
            write(unit, '(i0,3i8)') 3, mesh%triangles(:,i) - 1
        end do
        
        ! Cell types
        write(unit, '(a,i0)') 'CELL_TYPES ', mesh%n_triangles
        do i = 1, mesh%n_triangles
            write(unit, '(i0)') 5  ! VTK_TRIANGLE
        end do
        
        ! Point data
        write(unit, '(a,i0)') 'POINT_DATA ', mesh%n_vertices
        
        if (is_vector) then
            write(unit, '(a)') 'VECTORS E double'
            ! Simplified - just write as 2D vector
            do i = 1, mesh%n_vertices
                write(unit, '(3f16.8)') u(2*i-1), u(2*i), 0.0_dp
            end do
        else
            write(unit, '(a)') 'SCALARS u double 1'
            write(unit, '(a)') 'LOOKUP_TABLE default'
            do i = 1, mesh%n_vertices
                write(unit, '(f16.8)') u(i)
            end do
        end if
        
        close(unit)
    end subroutine write_vtk
    
    function compute_l2_error(mesh, u) result(error)
        type(mesh_2d_t), intent(in) :: mesh
        real(dp), intent(in) :: u(:)
        real(dp) :: error
        
        real(dp) :: x, y, exact
        integer :: i
        
        error = 0.0_dp
        do i = 1, mesh%n_vertices
            x = mesh%vertices(1, i)
            y = mesh%vertices(2, i)
            exact = sin(pi*x) * sin(pi*y)
            error = error + (u(i) - exact)**2
        end do
        
        error = sqrt(error / mesh%n_vertices)
    end function compute_l2_error
    
    subroutine solve_lapack_dense(A, b, info)
        type(sparse_matrix_t), intent(in) :: A
        real(dp), intent(inout) :: b(:)
        integer, intent(out) :: info
        
        ! Simple dense solver using LAPACK for testing
        real(dp), allocatable :: A_dense(:,:)
        integer, allocatable :: ipiv(:)
        integer :: n
        
        n = A%n
        allocate(A_dense(n, n))
        allocate(ipiv(n))
        
        ! Convert sparse to dense
        call A%to_dense(A_dense)
        
        ! Solve using DGESV
        call dgesv(n, 1, A_dense, n, ipiv, b, n, info)
        
        deallocate(A_dense, ipiv)
    end subroutine solve_lapack_dense

end module fortfem_simple_solvers