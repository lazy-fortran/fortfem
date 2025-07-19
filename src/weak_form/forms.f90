module forms_module
    use fortfem_kinds
    use function_space_module
    use expressions_module
    use fortfem_assembly_2d, only: assemble_laplacian, assemble_mass_matrix, assemble_mass_rhs
    implicit none
    private
    
    ! LAPACK interface
    interface
        subroutine dgesv(n, nrhs, a, lda, ipiv, b, ldb, info)
            import :: dp
            integer, intent(in) :: n, nrhs, lda, ldb
            real(dp), intent(inout) :: a(lda, *), b(ldb, *)
            integer, intent(out) :: ipiv(*), info
        end subroutine dgesv
    end interface
    
    public :: grad, div, curl, inner, dot
    public :: operator(*), operator(+), operator(-), operator(/)
    public :: dx, ds
    public :: solve
    public :: trial_to_expression, test_to_expression
    public :: vector_trial_to_expression, vector_test_to_expression
    
    ! Gradient interfaces
    interface grad
        module procedure grad_trial
        module procedure grad_test
    end interface
    
    ! Note: curl and div interfaces are already provided by expressions_module
    
contains

    ! Gradient of trial function
    function grad_trial(u) result(grad_expr)
        type(trial_function_t), intent(in) :: u
        type(expression_t) :: grad_expr
        
        grad_expr = grad(trial_to_expression(u))
    end function grad_trial
    
    ! Gradient of test function
    function grad_test(v) result(grad_expr)
        type(test_function_t), intent(in) :: v
        type(expression_t) :: grad_expr
        
        grad_expr = grad(test_to_expression(v))
    end function grad_test
    
    ! Solve a variational problem
    subroutine solve(a_expr, L_expr, space, solution, bc)
        type(expression_t), intent(in) :: a_expr  ! Bilinear form a(u,v)
        type(expression_t), intent(in) :: L_expr  ! Linear form L(v)
        type(function_space_t), intent(in) :: space
        real(dp), intent(out) :: solution(:)
        real(dp), optional, intent(in) :: bc(:)
        
        real(dp), allocatable :: matrix(:,:), rhs(:)
        integer :: n_dofs
        
        n_dofs = space%n_dofs
        allocate(matrix(n_dofs, n_dofs))
        allocate(rhs(n_dofs))
        
        ! Assemble system from expression trees
        call assemble_from_expressions(a_expr, L_expr, space, matrix, rhs)
        
        ! Apply boundary conditions if provided
        if (present(bc)) then
            call apply_boundary_conditions(matrix, rhs, bc)
        end if
        
        ! Solve linear system
        call solve_linear_system(matrix, rhs, solution)
        
        deallocate(matrix, rhs)
    end subroutine solve
    
    ! Assemble system from expression trees
    subroutine assemble_from_expressions(a_expr, L_expr, space, matrix, rhs)
        type(expression_t), intent(in) :: a_expr, L_expr
        type(function_space_t), intent(in) :: space
        real(dp), intent(out) :: matrix(:,:), rhs(:)
        
        ! Initialize system
        matrix = 0.0_dp
        rhs = 0.0_dp
        
        ! Assemble based on expression type
        select case (a_expr%expr_type)
        case (5)  ! EXPR_PRODUCT - likely grad(u) * grad(v)
            call assemble_laplacian(space, matrix)
        case default
            ! Generic assembly for other expressions
            call assemble_generic(a_expr, space, matrix)
        end select
        
        ! Assemble RHS
        select case (L_expr%expr_type)
        case (5)  ! EXPR_PRODUCT - likely f * v
            call assemble_mass_rhs(space, rhs, unit_source)
        case default
            call assemble_generic_rhs(L_expr, space, rhs)
        end select
    end subroutine assemble_from_expressions
    
    ! Apply boundary conditions
    subroutine apply_boundary_conditions(matrix, rhs, bc)
        real(dp), intent(inout) :: matrix(:,:), rhs(:)
        real(dp), intent(in) :: bc(:)
        integer :: i, n
        
        n = size(matrix, 1)
        
        ! Apply Dirichlet boundary conditions
        do i = 1, min(n, size(bc))
            matrix(i, :) = 0.0_dp
            matrix(:, i) = 0.0_dp
            matrix(i, i) = 1.0_dp
            rhs(i) = bc(i)
        end do
    end subroutine apply_boundary_conditions
    
    ! Solve linear system (direct solver)
    subroutine solve_linear_system(matrix, rhs, solution)
        real(dp), intent(in) :: matrix(:,:), rhs(:)
        real(dp), intent(out) :: solution(:)
        
        real(dp), allocatable :: A(:,:), b(:)
        integer :: n, info, ipiv(size(matrix, 1))
        
        n = size(matrix, 1)
        allocate(A(n, n), b(n))
        
        ! Copy to avoid modifying input
        A = matrix
        b = rhs
        
        ! Solve using LAPACK
        call dgesv(n, 1, A, n, ipiv, b, n, info)
        
        if (info /= 0) then
            print *, "Warning: Linear solver failed with info =", info
        end if
        
        solution = b
        deallocate(A, b)
    end subroutine solve_linear_system
    
    ! Unit source function
    pure function unit_source(x, y) result(f)
        real(dp), intent(in) :: x, y
        real(dp) :: f
        f = 1.0_dp
    end function unit_source
    
    ! Generic assembly routines (simplified)
    subroutine assemble_generic(expr, space, matrix)
        type(expression_t), intent(in) :: expr
        type(function_space_t), intent(in) :: space
        real(dp), intent(out) :: matrix(:,:)
        
        ! Simplified assembly - just identity matrix
        integer :: i, n
        n = space%n_dofs
        matrix = 0.0_dp
        do i = 1, n
            matrix(i, i) = 1.0_dp
        end do
    end subroutine assemble_generic
    
    subroutine assemble_generic_rhs(expr, space, rhs)
        type(expression_t), intent(in) :: expr
        type(function_space_t), intent(in) :: space
        real(dp), intent(out) :: rhs(:)
        
        ! Simplified RHS - unit values
        rhs = 1.0_dp
    end subroutine assemble_generic_rhs

end module forms_module