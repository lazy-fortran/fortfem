module fortfem_api
    ! High-level FEniCS-style API for FortFEM
    use fortfem_kinds
    use fortfem_mesh_2d
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
    
    ! Form operations (simplified)
    function inner(a, b) result(expr)
        class(*), intent(in) :: a, b
        type(simple_expression_t) :: expr
        
        expr%description = "inner(grad(u), grad(v))"
    end function inner
    
    function grad(u) result(gradu)
        class(*), intent(in) :: u
        type(simple_expression_t) :: gradu
        
        gradu%description = "grad(u)"
    end function grad
    
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