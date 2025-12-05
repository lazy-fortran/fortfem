module fortfem_api_spaces
    use fortfem_kinds
    use fortfem_api_types, only: mesh_t, function_space_t,                 &
        vector_function_space_t, function_t, vector_function_t,             &
        trial_function_t, test_function_t, vector_trial_function_t,         &
        vector_test_function_t, dirichlet_bc_t, vector_bc_t, neumann_bc_t
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
    public :: function_space
    public :: vector_function_space
    public :: function
    public :: vector_function
    public :: trial_function
    public :: test_function
    public :: vector_trial_function
    public :: vector_test_function
    public :: constant
    public :: dirichlet_bc
    public :: dirichlet_bc_on_boundary
    public :: vector_bc
    public :: neumann_bc_constant
    public :: neumann_bc_on_boundary

contains

    function function_space(mesh, family, degree) result(space)
        type(mesh_t), target, intent(in) :: mesh
        character(len=*), intent(in) :: family
        integer, intent(in) :: degree
        type(function_space_t) :: space

        space%mesh => mesh
        space%element_family = family
        space%degree = degree

        select case (trim(family))
        case ("Lagrange", "P")
            if (degree == 1) then
                space%ndof = mesh%data%n_vertices
            else if (degree == 2) then
                space%ndof = mesh%data%n_vertices + mesh%data%n_edges
            end if
        end select
    end function function_space

    function vector_function_space(mesh, family, degree) result(space)
        type(mesh_t), target, intent(in) :: mesh
        character(len=*), intent(in) :: family
        integer, intent(in) :: degree
        type(vector_function_space_t) :: space

        space%mesh => mesh
        space%element_family = family
        space%degree = degree
        space%n_components = 2

        select case (trim(family))
        case ("Nedelec", "Edge", "RT")
            if (degree == 1) then
                space%ndof = mesh%data%n_edges
            end if
        end select
    end function vector_function_space

    function function(space) result(f)
        type(function_space_t), target, intent(in) :: space
        type(function_t) :: f

        f%space => space
        allocate(f%values(space%ndof))
        f%values = 0.0_dp
    end function function

    function vector_function(space) result(f)
        type(vector_function_space_t), target, intent(in) :: space
        type(vector_function_t) :: f

        f%space => space
        allocate(f%values(space%ndof, space%n_components))
        f%values = 0.0_dp
    end function vector_function

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

    function vector_trial_function(space) result(u)
        type(vector_function_space_t), target, intent(in) :: space
        type(vector_trial_function_t) :: u

        u%space => space
    end function vector_trial_function

    function vector_test_function(space) result(v)
        type(vector_function_space_t), target, intent(in) :: space
        type(vector_test_function_t) :: v

        v%space => space
    end function vector_test_function

    function constant(val) result(f)
        real(dp), intent(in) :: val
        type(function_t) :: f

        allocate(f%values(1))
        f%values(1) = val
    end function constant

    function dirichlet_bc(space, value) result(bc)
        type(function_space_t), target, intent(in) :: space
        real(dp), intent(in) :: value
        type(dirichlet_bc_t) :: bc

        bc%space => space
        bc%value = value
        bc%on_boundary = .true.
    end function dirichlet_bc

    function vector_bc(space, values, bc_type) result(bc)
        type(vector_function_space_t), target, intent(in) :: space
        real(dp), intent(in) :: values(2)
        character(len=*), intent(in), optional :: bc_type
        type(vector_bc_t) :: bc

        bc%space => space
        bc%values = values
        bc%on_boundary = .true.
        if (present(bc_type)) then
            bc%bc_type = bc_type
        else
            bc%bc_type = "tangential"
        end if
    end function vector_bc

    function neumann_bc_constant(space, flux_value) result(bc)
        type(function_space_t), target, intent(in) :: space
        real(dp), intent(in) :: flux_value
        type(neumann_bc_t) :: bc

        bc%space => space
        bc%flux_type = "constant"
        bc%constant_value = flux_value
        bc%boundary_part = "all"
        bc%on_boundary = .true.
    end function neumann_bc_constant

    function neumann_bc_on_boundary(space, flux_value, boundary_part)       &
        result(bc)
        type(function_space_t), target, intent(in) :: space
        real(dp), intent(in) :: flux_value
        character(len=*), intent(in) :: boundary_part
        type(neumann_bc_t) :: bc

        bc%space => space
        bc%flux_type = "constant"
        bc%constant_value = flux_value
        bc%boundary_part = boundary_part
        bc%on_boundary = .true.
    end function neumann_bc_on_boundary

    function dirichlet_bc_on_boundary(space, value, boundary_part)          &
        result(bc)
        type(function_space_t), target, intent(in) :: space
        real(dp), intent(in) :: value
        character(len=*), intent(in) :: boundary_part
        type(dirichlet_bc_t) :: bc

        bc%space => space
        bc%value = value
        bc%on_boundary = .true.
    end function dirichlet_bc_on_boundary

end module fortfem_api_spaces

