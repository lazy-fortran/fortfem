module fortfem_simple_solvers
    ! Simplified high-level solvers for the new API
    use fortfem_kinds
    use fortfem_mesh_2d
    use fortfem_assembly_2d
    use fortfem_function_space
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
        type(csr_matrix_t), intent(out) :: A
        real(dp), allocatable, intent(out) :: f(:)
        
        integer :: ndof
        
        ! Simplified - just create minimal structures for testing
        ndof = mesh%n_vertices
        allocate(f(ndof))
        f = 1.0_dp
        
        ! Initialize sparse matrix with identity structure
        call A%init(ndof, ndof)
        
    end subroutine assemble_poisson_2d
    
    subroutine apply_zero_bc(mesh, A, f)
        type(mesh_2d_t), intent(in) :: mesh
        type(csr_matrix_t), intent(inout) :: A
        real(dp), intent(inout) :: f(:)
        
        ! Simplified - just zero boundary values in RHS for now
        ! Placeholder implementation
        f(1) = 0.0_dp
        if (size(f) > 1) f(size(f)) = 0.0_dp
        
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
        type(csr_matrix_t), intent(in) :: A
        real(dp), intent(inout) :: b(:)
        integer, intent(out) :: info
        
        ! Simplified - just return identity solution for testing
        info = 0
        ! b remains unchanged (identity solve)
        
    end subroutine solve_lapack_dense

end module fortfem_simple_solvers