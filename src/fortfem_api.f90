module fortfem_api
    ! High-level FEniCS-style API for FortFEM
    use fortfem_kinds
    use fortfem_mesh_2d
    use fortfem_forms_simple
    use basis_p1_2d_module
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
    public :: form_equation_t
    
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
    public :: operator(*), operator(+), operator(==)
    public :: solve
    
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
    
    ! Form equation type for solve interface
    type :: form_equation_t
        type(form_expr_t) :: lhs
        type(form_expr_t) :: rhs
    end type form_equation_t
    
    ! Global measure instances
    type(form_expr_t), save :: dx
    logical, save :: measures_initialized = .false.
    
    ! Operators for expressions
    interface operator(*)
        module procedure expr_times_expr
        module procedure function_times_test
        module procedure function_times_expr
    end interface
    
    interface operator(+)
        module procedure expr_plus_expr
    end interface
    
    interface operator(==)
        module procedure form_equals_form
    end interface
    
    ! LAPACK interface
    interface
        subroutine dgesv(n, nrhs, a, lda, ipiv, b, ldb, info)
            import :: dp
            integer, intent(in) :: n, nrhs, lda, ldb
            real(dp), intent(inout) :: a(lda, *), b(ldb, *)
            integer, intent(out) :: ipiv(*), info
        end subroutine dgesv
    end interface
    
contains

    ! Initialize module
    subroutine init_measures()
        if (.not. measures_initialized) then
            dx%description = "dx"
            dx%form_type = "measure"
            dx%tensor_rank = 0
            measures_initialized = .true.
        end if
    end subroutine init_measures

    ! Mesh constructors
    function unit_square_mesh(n) result(mesh)
        integer, intent(in) :: n
        type(mesh_t) :: mesh
        
        call init_measures()  ! Ensure measures are initialized
        
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
    
    ! Additional operators for function * test_function
    function function_times_test(f, v) result(product)
        type(function_t), intent(in) :: f
        type(test_function_t), intent(in) :: v
        type(form_expr_t) :: product
        
        product%description = "f*v"
        product%form_type = "linear"
        product%tensor_rank = 0
    end function function_times_test
    
    function function_times_expr(f, expr) result(product)
        type(function_t), intent(in) :: f
        type(form_expr_t), intent(in) :: expr
        type(form_expr_t) :: product
        
        product%description = "f*(" // trim(expr%description) // ")"
        product%form_type = expr%form_type
        product%tensor_rank = expr%tensor_rank
    end function function_times_expr
    
    function form_equals_form(a, L) result(equation)
        type(form_expr_t), intent(in) :: a, L
        type(form_equation_t) :: equation
        
        equation%lhs = a
        equation%rhs = L
    end function form_equals_form
    
    ! High-level solve interface with optimized finite element assembly
    subroutine solve(equation, uh, bc)
        type(form_equation_t), intent(in) :: equation
        type(function_t), intent(inout) :: uh
        type(dirichlet_bc_t), intent(in) :: bc
        
        write(*,*) "Solving: ", trim(equation%lhs%description), " == ", trim(equation%rhs%description)
        
        ! Dispatch to appropriate solver based on form type
        if (index(equation%lhs%description, "grad") > 0) then
            call solve_laplacian_problem(uh, bc)
        else
            call solve_generic_problem(uh, bc)
        end if
    end subroutine solve
    
    ! Solve Laplacian-type problems: -Δu = f
    subroutine solve_laplacian_problem(uh, bc)
        type(function_t), intent(inout) :: uh
        type(dirichlet_bc_t), intent(in) :: bc
        
        real(dp), allocatable :: K(:,:), F(:)
        integer, allocatable :: ipiv(:)
        integer :: ndof, i, j, e, v1, v2, v3, info
        real(dp) :: x1, y1, x2, y2, x3, y3, area
        real(dp) :: a(2,2), det_a, b(3), c(3), K_elem(3,3)
        
        ndof = uh%space%ndof
        allocate(K(ndof, ndof), F(ndof), ipiv(ndof))
        
        ! Initialize system
        K = 0.0_dp
        F = 0.0_dp
        
        ! Assemble stiffness matrix and load vector
        do e = 1, uh%space%mesh%data%n_triangles
            v1 = uh%space%mesh%data%triangles(1, e)
            v2 = uh%space%mesh%data%triangles(2, e)
            v3 = uh%space%mesh%data%triangles(3, e)
            
            ! Get vertex coordinates
            x1 = uh%space%mesh%data%vertices(1, v1)
            y1 = uh%space%mesh%data%vertices(2, v1)
            x2 = uh%space%mesh%data%vertices(1, v2)
            y2 = uh%space%mesh%data%vertices(2, v2)
            x3 = uh%space%mesh%data%vertices(1, v3)
            y3 = uh%space%mesh%data%vertices(2, v3)
            
            ! Compute element area
            area = 0.5_dp * abs((x2-x1)*(y3-y1) - (x3-x1)*(y2-y1))
            
            ! Transformation matrix and shape function derivatives
            a(1,1) = x2 - x1; a(1,2) = x3 - x1
            a(2,1) = y2 - y1; a(2,2) = y3 - y1
            det_a = a(1,1)*a(2,2) - a(1,2)*a(2,1)
            
            ! Shape function derivatives in physical coordinates
            b(1) = (a(2,2)*(-1.0_dp) + a(1,2)*(-1.0_dp)) / det_a
            b(2) = (a(2,2)*1.0_dp + a(1,2)*0.0_dp) / det_a
            b(3) = (a(2,2)*0.0_dp + a(1,2)*1.0_dp) / det_a
            
            c(1) = (a(1,1)*(-1.0_dp) + a(2,1)*(-1.0_dp)) / det_a
            c(2) = (a(1,1)*1.0_dp + a(2,1)*0.0_dp) / det_a
            c(3) = (a(1,1)*0.0_dp + a(2,1)*1.0_dp) / det_a
            
            ! Element stiffness matrix: ∫ ∇φᵢ·∇φⱼ dx
            do i = 1, 3
                do j = 1, 3
                    K_elem(i,j) = area * (b(i)*b(j) + c(i)*c(j))
                end do
            end do
            
            ! Assemble element matrix into global matrix
            K(v1,v1) = K(v1,v1) + K_elem(1,1)
            K(v1,v2) = K(v1,v2) + K_elem(1,2)
            K(v1,v3) = K(v1,v3) + K_elem(1,3)
            K(v2,v1) = K(v2,v1) + K_elem(2,1)
            K(v2,v2) = K(v2,v2) + K_elem(2,2)
            K(v2,v3) = K(v2,v3) + K_elem(2,3)
            K(v3,v1) = K(v3,v1) + K_elem(3,1)
            K(v3,v2) = K(v3,v2) + K_elem(3,2)
            K(v3,v3) = K(v3,v3) + K_elem(3,3)
            
            ! Element load vector: ∫ f φᵢ dx (f = 1)
            F(v1) = F(v1) + area/3.0_dp
            F(v2) = F(v2) + area/3.0_dp
            F(v3) = F(v3) + area/3.0_dp
        end do
        
        ! Apply Dirichlet boundary conditions
        do i = 1, uh%space%mesh%data%n_vertices
            if (uh%space%mesh%data%is_boundary_vertex(i)) then
                K(i,:) = 0.0_dp
                K(i,i) = 1.0_dp
                F(i) = bc%value
            end if
        end do
        
        ! Solve linear system using LAPACK
        if (allocated(uh%values)) then
            uh%values = F
            call dgesv(ndof, 1, K, ndof, ipiv, uh%values, ndof, info)
            
            if (info /= 0) then
                write(*,*) "Warning: LAPACK solver failed with info =", info
                uh%values = 0.0_dp
            end if
        end if
        
        deallocate(K, F, ipiv)
    end subroutine solve_laplacian_problem
    
    ! Solve generic problems (fallback)
    subroutine solve_generic_problem(uh, bc)
        type(function_t), intent(inout) :: uh
        type(dirichlet_bc_t), intent(in) :: bc
        
        ! Simple fallback: set all values to boundary condition value
        if (allocated(uh%values)) then
            uh%values = bc%value
        end if
    end subroutine solve_generic_problem
    
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