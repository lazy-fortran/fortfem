module fortfem_api
    ! High-level FEniCS-style API for FortFEM
    use fortfem_kinds
    use fortfem_mesh_2d
    use fortfem_assembly_2d
    use fortfem_solver
    use fortfem_sparse_matrix
    implicit none
    
    private
    
    ! Public types with _t suffix
    public :: mesh_t
    public :: function_space_t
    public :: function_t
    public :: trial_function_t
    public :: test_function_t
    public :: form_t
    public :: bilinear_form_t
    public :: linear_form_t
    public :: dirichlet_bc_t
    public :: measure_t
    
    ! Public constructors
    public :: unit_square_mesh
    public :: rectangle_mesh
    public :: function_space
    public :: function
    public :: trial_function
    public :: test_function
    public :: constant
    public :: dirichlet_bc
    
    ! Public form operations
    public :: inner, dot, cross, outer
    public :: grad, div, curl
    public :: dx, ds, dS
    
    ! Public solver
    public :: solve
    public :: norm
    
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
        integer, allocatable :: dof_map(:,:)
    contains
        procedure :: destroy => function_space_destroy
    end type function_space_t
    
    ! Function type (holds values)
    type :: function_t
        type(function_space_t), pointer :: space => null()
        real(dp), allocatable :: values(:)
    contains
        procedure :: interpolate => function_interpolate
        procedure :: write_vtk => function_write_vtk
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
    
    ! Form types
    type :: form_t
        integer :: rank = 0  ! 0=functional, 1=linear, 2=bilinear
    end type form_t
    
    type, extends(form_t) :: bilinear_form_t
        procedure(bilinear_integrand), pointer, nopass :: integrand => null()
    end type bilinear_form_t
    
    type, extends(form_t) :: linear_form_t
        procedure(linear_integrand), pointer, nopass :: integrand => null()
    end type linear_form_t
    
    ! Measure type
    type :: measure_t
        character(len=8) :: name = ""
    end type measure_t
    
    ! Boundary condition type
    type :: dirichlet_bc_t
        type(function_space_t), pointer :: space => null()
        real(dp) :: value = 0.0_dp
        integer :: boundary_marker = 0
        logical :: on_boundary = .false.
    contains
        procedure :: apply => dirichlet_bc_apply
    end type dirichlet_bc_t
    
    ! Interfaces for form integrands
    abstract interface
        real(dp) function bilinear_integrand(u, v, x, y)
            import :: dp
            real(dp), intent(in) :: u(:), v(:), x, y
        end function bilinear_integrand
        
        real(dp) function linear_integrand(v, x, y)
            import :: dp
            real(dp), intent(in) :: v(:), x, y
        end function linear_integrand
    end interface
    
    ! Operators
    interface operator(*)
        module procedure form_times_measure
    end interface
    
    interface operator(==)
        module procedure forms_equal
    end interface
    
    ! Global measure instances
    type(measure_t), parameter :: dx = measure_t("dx")
    type(measure_t), parameter :: ds = measure_t("ds")
    type(measure_t), parameter :: dS = measure_t("dS")
    
contains

    ! Mesh constructors
    function unit_square_mesh(n, nx, ny) result(mesh)
        integer, intent(in) :: n
        integer, intent(in), optional :: nx, ny
        type(mesh_t) :: mesh
        integer :: nx_, ny_
        
        if (present(nx) .and. present(ny)) then
            nx_ = nx
            ny_ = ny
        else
            nx_ = n
            ny_ = n
        end if
        
        call mesh%data%create_rectangular(nx=nx_, ny=ny_, &
                                         x_min=0.0_dp, x_max=1.0_dp, &
                                         y_min=0.0_dp, y_max=1.0_dp)
        call mesh%data%build_connectivity()
        call mesh%data%find_boundary()
    end function unit_square_mesh
    
    function rectangle_mesh(corner1, corner2, nx, ny) result(mesh)
        real(dp), intent(in) :: corner1(2), corner2(2)
        integer, intent(in) :: nx, ny
        type(mesh_t) :: mesh
        
        call mesh%data%create_rectangular(nx=nx, ny=ny, &
                                         x_min=corner1(1), x_max=corner2(1), &
                                         y_min=corner1(2), y_max=corner2(2))
        call mesh%data%build_connectivity()
        call mesh%data%find_boundary()
    end function rectangle_mesh
    
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
            else if (degree == 2) then
                space%ndof = mesh%data%n_vertices + mesh%data%n_edges
            end if
        case ("Nedelec", "N")
            space%ndof = mesh%data%n_edges
        case ("RT", "Raviart-Thomas")
            space%ndof = mesh%data%n_edges
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
        
        ! Create a constant function (simplified - needs proper implementation)
        allocate(f%values(1))
        f%values(1) = val
    end function constant
    
    ! Boundary condition constructor
    function dirichlet_bc(space, value, marker) result(bc)
        type(function_space_t), target, intent(in) :: space
        real(dp), intent(in) :: value
        integer, intent(in), optional :: marker
        type(dirichlet_bc_t) :: bc
        
        bc%space => space
        bc%value = value
        if (present(marker)) then
            bc%boundary_marker = marker
        else
            bc%on_boundary = .true.
        end if
    end function dirichlet_bc
    
    ! Form operations
    function inner(a, b) result(form)
        class(*), intent(in) :: a, b
        type(bilinear_form_t) :: form
        
        form%rank = 2
        ! Set integrand based on input types
    end function inner
    
    function grad(u) result(gradu)
        class(*), intent(in) :: u
        class(*), allocatable :: gradu
        
        ! Return gradient operator applied to u
    end function grad
    
    function form_times_measure(form, measure) result(integrated_form)
        class(form_t), intent(in) :: form
        type(measure_t), intent(in) :: measure
        class(form_t), allocatable :: integrated_form
        
        ! Attach measure to form
        allocate(integrated_form, source=form)
    end function form_times_measure
    
    function forms_equal(a, L) result(problem)
        type(bilinear_form_t), intent(in) :: a
        type(linear_form_t), intent(in) :: L
        logical :: problem
        
        problem = .true.  ! Placeholder
    end function forms_equal
    
    ! Solver interface
    subroutine solve(a, L, u, bcs)
        type(bilinear_form_t), intent(in) :: a
        type(linear_form_t), intent(in) :: L
        type(function_t), intent(inout) :: u
        type(dirichlet_bc_t), intent(in), optional :: bcs(:)
        
        type(sparse_matrix_t) :: A
        real(dp), allocatable :: b(:)
        
        ! Assemble system
        ! call assemble_system(a, L, A, b)
        
        ! Apply boundary conditions
        if (present(bcs)) then
            ! call apply_bcs(A, b, bcs)
        end if
        
        ! Solve
        ! call solve_sparse(A, b, u%values)
    end subroutine solve
    
    ! Utility functions
    function norm(u, norm_type) result(val)
        type(function_t), intent(in) :: u
        character(len=*), intent(in) :: norm_type
        real(dp) :: val
        
        select case (norm_type)
        case ("L2")
            val = sqrt(sum(u%values**2))
        case default
            val = 0.0_dp
        end select
    end function norm
    
    ! Destructor procedures
    subroutine mesh_destroy(this)
        class(mesh_t), intent(inout) :: this
        call this%data%destroy()
    end subroutine mesh_destroy
    
    subroutine function_space_destroy(this)
        class(function_space_t), intent(inout) :: this
        if (allocated(this%dof_map)) deallocate(this%dof_map)
        this%mesh => null()
    end subroutine function_space_destroy
    
    subroutine function_destroy(this)
        class(function_t), intent(inout) :: this
        if (allocated(this%values)) deallocate(this%values)
        this%space => null()
    end subroutine function_destroy
    
    subroutine function_interpolate(this, f)
        class(function_t), intent(inout) :: this
        interface
            real(dp) function f(x, y)
                import :: dp
                real(dp), intent(in) :: x, y
            end function f
        end interface
        
        integer :: i
        real(dp) :: x, y
        
        ! Simple interpolation at vertices for P1
        do i = 1, this%space%ndof
            x = this%space%mesh%data%vertices(1, i)
            y = this%space%mesh%data%vertices(2, i)
            this%values(i) = f(x, y)
        end do
    end subroutine function_interpolate
    
    subroutine function_write_vtk(this, filename)
        class(function_t), intent(in) :: this
        character(len=*), intent(in) :: filename
        
        ! Write VTK file
        ! call write_vtk_scalar(filename, this%space%mesh%data, this%values)
    end subroutine function_write_vtk
    
    subroutine dirichlet_bc_apply(this, A, b)
        class(dirichlet_bc_t), intent(in) :: this
        type(sparse_matrix_t), intent(inout) :: A
        real(dp), intent(inout) :: b(:)
        
        ! Apply boundary conditions
        ! Implementation needed
    end subroutine dirichlet_bc_apply

end module fortfem_api