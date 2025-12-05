module fortfem_api_types
    use fortfem_kinds
    use fortfem_mesh_2d
    implicit none

    private

    public :: mesh_t
    public :: function_space_t
    public :: vector_function_space_t
    public :: function_t
    public :: vector_function_t
    public :: trial_function_t
    public :: test_function_t
    public :: vector_trial_function_t
    public :: vector_test_function_t
    public :: dirichlet_bc_t
    public :: vector_bc_t
    public :: neumann_bc_t
    public :: simple_expression_t

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

    ! Vector function space type (for edge elements)
    type :: vector_function_space_t
        type(mesh_t), pointer :: mesh => null()
        character(len=32) :: element_family = ""
        integer :: degree = 0
        integer :: ndof = 0
        integer :: n_components = 2
    contains
        procedure :: destroy => vector_function_space_destroy
    end type vector_function_space_t

    ! Function type (holds values)
    type :: function_t
        type(function_space_t), pointer :: space => null()
        real(dp), allocatable :: values(:)
    contains
        procedure :: destroy => function_destroy
    end type function_t

    ! Vector function type (holds vector values)
    type :: vector_function_t
        type(vector_function_space_t), pointer :: space => null()
        real(dp), allocatable :: values(:,:)
    contains
        procedure :: destroy => vector_function_destroy
    end type vector_function_t

    ! Trial function type (symbolic)
    type :: trial_function_t
        type(function_space_t), pointer :: space => null()
    end type trial_function_t

    ! Test function type (symbolic)
    type :: test_function_t
        type(function_space_t), pointer :: space => null()
    end type test_function_t

    ! Vector trial function type (symbolic)
    type :: vector_trial_function_t
        type(vector_function_space_t), pointer :: space => null()
    end type vector_trial_function_t

    ! Vector test function type (symbolic)
    type :: vector_test_function_t
        type(vector_function_space_t), pointer :: space => null()
    end type vector_test_function_t

    ! Boundary condition type
    type :: dirichlet_bc_t
        type(function_space_t), pointer :: space => null()
        real(dp) :: value = 0.0_dp
        logical :: on_boundary = .false.
    end type dirichlet_bc_t

    ! Vector boundary condition type
    type :: vector_bc_t
        type(vector_function_space_t), pointer :: space => null()
        real(dp) :: values(2) = [0.0_dp, 0.0_dp]
        character(len=32) :: bc_type = "tangential"
        logical :: on_boundary = .false.
    end type vector_bc_t

    ! Neumann boundary condition type
    type :: neumann_bc_t
        type(function_space_t), pointer :: space => null()
        character(len=32) :: flux_type = "constant"
        real(dp) :: constant_value = 0.0_dp
        character(len=32) :: boundary_part = "all"
        logical :: on_boundary = .true.
    end type neumann_bc_t

    ! Simple expression type for forms
    type :: simple_expression_t
        character(len=64) :: description = ""
    end type simple_expression_t

contains

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

    subroutine vector_function_space_destroy(this)
        class(vector_function_space_t), intent(inout) :: this
        this%mesh => null()
    end subroutine vector_function_space_destroy

    subroutine vector_function_destroy(this)
        class(vector_function_t), intent(inout) :: this
        if (allocated(this%values)) deallocate(this%values)
        this%space => null()
    end subroutine vector_function_destroy

end module fortfem_api_types

