module fortfem_api
    ! High-level FEniCS-style API for FortFEM
    use fortfem_kinds
    use fortfem_mesh_2d
    use fortfem_forms_simple
    implicit none
    
    private
    
    ! Public types with _t suffix
    public :: mesh_t
    public :: function_space_t
    public :: function_t
    public :: trial_function_t
    public :: test_function_t
    public :: dirichlet_bc_t
    public :: simple_expression_t
    public :: form_expr_t
    
    ! Public constructors
    public :: unit_square_mesh
    public :: function_space
    public :: function
    public :: trial_function
    public :: test_function
    public :: constant
    public :: dirichlet_bc
    
    ! Public form operations (simplified)
    public :: inner, grad
    public :: dx
    public :: compile_form
    
    ! Mesh type (wrapper around mesh_2d_t)
    type :: mesh_t
        type(mesh_2d_t) :: data
    contains
        procedure :: destroy => mesh_destroy
    end type mesh_t
    
    ! Function space type
    type :: function_space_t
        type(mesh_t), pointer :: mesh => null()
        character(len=32) :: element_family = ""
        integer :: degree = 0
        integer :: ndof = 0
    contains
        procedure :: destroy => function_space_destroy
    end type function_space_t
    
    ! Function type (holds values)
    type :: function_t
        type(function_space_t), pointer :: space => null()
        real(dp), allocatable :: values(:)
    contains
        procedure :: destroy => function_destroy
    end type function_t
    
    ! Trial function type (symbolic)
    type :: trial_function_t
        type(function_space_t), pointer :: space => null()
    end type trial_function_t
    
    ! Test function type (symbolic)
    type :: test_function_t
        type(function_space_t), pointer :: space => null()
    end type test_function_t
    
    ! Boundary condition type
    type :: dirichlet_bc_t
        type(function_space_t), pointer :: space => null()
        real(dp) :: value = 0.0_dp
        logical :: on_boundary = .false.
    end type dirichlet_bc_t
    
    ! Simple expression type for forms
    type :: simple_expression_t
        character(len=64) :: description = ""
    end type simple_expression_t
    
    ! Global measure instance
    type(simple_expression_t), parameter :: dx = simple_expression_t("dx")
    
    ! Operators for expressions
    interface operator(*)
        module procedure expr_times_measure
        module procedure expr_times_expr
    end interface
    
    interface operator(+)
        module procedure expr_plus_expr
    end interface
    
contains

    ! Mesh constructors
    function unit_square_mesh(n) result(mesh)
        integer, intent(in) :: n
        type(mesh_t) :: mesh
        
        call mesh%data%create_rectangular(nx=n, ny=n, &
                                         x_min=0.0_dp, x_max=1.0_dp, &
                                         y_min=0.0_dp, y_max=1.0_dp)
        call mesh%data%build_connectivity()
        call mesh%data%find_boundary()
    end function unit_square_mesh
    
    ! Function space constructor
    function function_space(mesh, family, degree) result(space)
        type(mesh_t), target, intent(in) :: mesh
        character(len=*), intent(in) :: family
        integer, intent(in) :: degree
        type(function_space_t) :: space
        
        space%mesh => mesh
        space%element_family = family
        space%degree = degree
        
        ! Set DOF count based on element type
        select case (trim(family))
        case ("Lagrange", "P")
            if (degree == 1) then
                space%ndof = mesh%data%n_vertices
            end if
        end select
    end function function_space
    
    ! Function constructors
    function function(space) result(f)
        type(function_space_t), target, intent(in) :: space
        type(function_t) :: f
        
        f%space => space
        allocate(f%values(space%ndof))
        f%values = 0.0_dp
    end function function
    
    function trial_function(space) result(u)
        type(function_space_t), target, intent(in) :: space
        type(trial_function_t) :: u
        
        u%space => space
    end function trial_function
    
    function test_function(space) result(v)
        type(function_space_t), target, intent(in) :: space
        type(test_function_t) :: v
        
        v%space => space
    end function test_function
    
    function constant(val) result(f)
        real(dp), intent(in) :: val
        type(function_t) :: f
        
        allocate(f%values(1))
        f%values(1) = val
    end function constant
    
    ! Boundary condition constructor
    function dirichlet_bc(space, value) result(bc)
        type(function_space_t), target, intent(in) :: space
        real(dp), intent(in) :: value
        type(dirichlet_bc_t) :: bc
        
        bc%space => space
        bc%value = value
        bc%on_boundary = .true.
    end function dirichlet_bc
    
    ! Form operations with simple expressions
    function inner(a, b) result(expr)
        class(*), intent(in) :: a, b
        type(form_expr_t) :: expr
        type(form_expr_t) :: expr_a, expr_b
        
        ! Convert inputs to expressions
        select type(a)
        type is (trial_function_t)
            expr_a = create_grad("u", "trial")
        type is (test_function_t)
            expr_a = create_grad("v", "test")
        type is (form_expr_t)
            expr_a = a
        class default
            expr_a%description = "unknown"
            expr_a%form_type = "unknown"
        end select
        
        select type(b)
        type is (trial_function_t)
            expr_b = create_grad("u", "trial")
        type is (test_function_t)
            expr_b = create_grad("v", "test")
        type is (form_expr_t)
            expr_b = b
        class default
            expr_b%description = "unknown"
            expr_b%form_type = "unknown"
        end select
        
        expr = create_inner(expr_a, expr_b)
    end function inner
    
    function grad(u) result(gradu)
        class(*), intent(in) :: u
        type(form_expr_t) :: gradu
        
        select type(u)
        type is (trial_function_t)
            gradu = create_grad("u", "trial")
        type is (test_function_t)
            gradu = create_grad("v", "test")
        class default
            gradu = create_grad("unknown", "unknown")
        end select
    end function grad
    
    ! Operator overloading
    function expr_times_measure(expr, measure) result(form)
        type(form_expr_t), intent(in) :: expr
        type(simple_expression_t), intent(in) :: measure
        type(form_expr_t) :: form
        
        ! For now, return the expression as-is
        ! In a full implementation, this would attach measure metadata
        form = expr
    end function expr_times_measure
    
    function expr_times_expr(a, b) result(product)
        type(form_expr_t), intent(in) :: a, b
        type(form_expr_t) :: product
        
        product%description = "(" // trim(a%description) // " * " // trim(b%description) // ")"
        product%form_type = a%form_type
        product%tensor_rank = a%tensor_rank + b%tensor_rank
    end function expr_times_expr
    
    function expr_plus_expr(a, b) result(sum_expr)
        type(form_expr_t), intent(in) :: a, b
        type(form_expr_t) :: sum_expr
        
        sum_expr%description = "(" // trim(a%description) // " + " // trim(b%description) // ")"
        sum_expr%form_type = a%form_type
        sum_expr%tensor_rank = max(a%tensor_rank, b%tensor_rank)
    end function expr_plus_expr
    
    ! Destructor procedures
    subroutine mesh_destroy(this)
        class(mesh_t), intent(inout) :: this
        call this%data%destroy()
    end subroutine mesh_destroy
    
    subroutine function_space_destroy(this)
        class(function_space_t), intent(inout) :: this
        this%mesh => null()
    end subroutine function_space_destroy
    
    subroutine function_destroy(this)
        class(function_t), intent(inout) :: this
        if (allocated(this%values)) deallocate(this%values)
        this%space => null()
    end subroutine function_destroy

end module fortfem_api