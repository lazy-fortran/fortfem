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
    public :: vector_neumann_bc_t
    public :: neumann_bc_t
    public :: simple_expression_t
    public :: cell_coefficient_t
    public :: cell_tensor_coefficient_t
    public :: cell_vector_source_t

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
        real(dp), allocatable :: edge_values(:)
        logical, allocatable :: edge_mask(:)
        character(len=32) :: bc_type = "tangential"
        logical :: on_boundary = .false.
    contains
        procedure, private :: assign_vector_bc
        generic :: assignment(=) => assign_vector_bc
    end type vector_bc_t

    type :: vector_neumann_bc_t
        real(dp), allocatable :: edge_values(:)
        logical, allocatable :: edge_mask(:)
    contains
        procedure, private :: assign_vector_neumann_bc
        generic :: assignment(=) => assign_vector_neumann_bc
    end type vector_neumann_bc_t

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

    type :: cell_coefficient_t
        real(dp), allocatable :: values(:)
    contains
        procedure, private :: assign_cell_coefficient
        generic :: assignment(=) => assign_cell_coefficient
    end type cell_coefficient_t

    type :: cell_tensor_coefficient_t
        real(dp), allocatable :: values(:, :, :)
    contains
        procedure, private :: assign_cell_tensor_coefficient
        generic :: assignment(=) => assign_cell_tensor_coefficient
    end type cell_tensor_coefficient_t

    type :: cell_vector_source_t
        real(dp), allocatable :: values(:, :)
    contains
        procedure, private :: assign_cell_vector_source
        generic :: assignment(=) => assign_cell_vector_source
    end type cell_vector_source_t

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

    subroutine assign_vector_bc(lhs, rhs)
        class(vector_bc_t), intent(out) :: lhs
        type(vector_bc_t), intent(in) :: rhs

        lhs%space => rhs%space
        lhs%values = rhs%values
        if (allocated(rhs%edge_values)) then
            allocate(lhs%edge_values, source=rhs%edge_values)
        end if
        if (allocated(rhs%edge_mask)) then
            allocate(lhs%edge_mask, source=rhs%edge_mask)
        end if
        lhs%bc_type = rhs%bc_type
        lhs%on_boundary = rhs%on_boundary
    end subroutine assign_vector_bc

    subroutine assign_vector_neumann_bc(lhs, rhs)
        class(vector_neumann_bc_t), intent(out) :: lhs
        type(vector_neumann_bc_t), intent(in) :: rhs

        if (allocated(rhs%edge_values)) then
            allocate(lhs%edge_values, source=rhs%edge_values)
        end if
        if (allocated(rhs%edge_mask)) then
            allocate(lhs%edge_mask, source=rhs%edge_mask)
        end if
    end subroutine assign_vector_neumann_bc

    subroutine assign_cell_coefficient(lhs, rhs)
        class(cell_coefficient_t), intent(out) :: lhs
        type(cell_coefficient_t), intent(in) :: rhs

        if (allocated(rhs%values)) then
            allocate(lhs%values, source=rhs%values)
        end if
    end subroutine assign_cell_coefficient

    subroutine assign_cell_tensor_coefficient(lhs, rhs)
        class(cell_tensor_coefficient_t), intent(out) :: lhs
        type(cell_tensor_coefficient_t), intent(in) :: rhs

        if (allocated(rhs%values)) then
            allocate(lhs%values, source=rhs%values)
        end if
    end subroutine assign_cell_tensor_coefficient

    subroutine assign_cell_vector_source(lhs, rhs)
        class(cell_vector_source_t), intent(out) :: lhs
        type(cell_vector_source_t), intent(in) :: rhs

        if (allocated(rhs%values)) then
            allocate(lhs%values, source=rhs%values)
        end if
    end subroutine assign_cell_vector_source

end module fortfem_api_types
