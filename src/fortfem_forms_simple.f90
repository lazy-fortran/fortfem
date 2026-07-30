module fortfem_forms_simple
    use fortfem_assembly_nedelec_arbitrary_order_2d, only: &
        assemble_triangle_nedelec_cell_vector_load, &
        assemble_triangle_nedelec_cell_tensor_csc, &
        assemble_triangle_nedelec_curl_mass_csc, &
        assemble_triangle_nedelec_curl_mass_element
    use fortfem_assembly_full_vector_arbitrary_order_2d, only: &
        assemble_triangle_bdm_cell_vector_load, &
        assemble_triangle_bdm_div_mass_csc, &
        assemble_triangle_bdm_div_mass_element, &
        assemble_triangle_nedelec_second_cell_vector_load, &
        assemble_triangle_nedelec_second_curl_mass_csc, &
        assemble_triangle_nedelec_second_curl_mass_element
    use fortfem_assembly_rt_arbitrary_order_2d, only: &
        assemble_triangle_rt_cell_vector_load, &
        assemble_triangle_rt_div_mass_csc, &
        assemble_triangle_rt_div_mass_element
    use fortfem_kinds, only: dp
    use fortfem_mesh_2d, only: mesh_2d_t
    use fortsparse, only: csc_from_triplet, csc_matvec, csc_t, &
        FORTSPARSE_INVALID_MATRIX, fortsparse_status_t, status_set
    implicit none

    private

    integer, parameter :: token_symbol = 1
    integer, parameter :: token_gradient = 2
    integer, parameter :: token_curl = 3
    integer, parameter :: token_inner = 4
    integer, parameter :: token_multiply = 5
    integer, parameter :: token_add = 6
    integer, parameter :: token_scalar = 7
    integer, parameter :: token_measure = 8
    integer, parameter :: token_divergence = 9
    integer, parameter :: token_constant_load = 10
    integer, parameter :: token_cell_coefficient = 11
    integer, parameter :: token_cell_vector_source = 12
    integer, parameter :: token_cell_tensor_coefficient = 13

    integer, parameter :: role_trial = 1
    integer, parameter :: role_test = 2
    integer, parameter :: role_function = 3

    integer, parameter :: derivative_identity = 0
    integer, parameter :: derivative_gradient = 1
    integer, parameter :: derivative_curl = 2
    integer, parameter :: derivative_divergence = 3

    integer, parameter :: item_argument = 1
    integer, parameter :: item_form = 2
    integer, parameter :: item_scalar = 3
    integer, parameter :: item_measure = 4
    integer, parameter :: item_coefficient = 5

    type :: form_token_t
        integer :: token_type = 0
        integer :: role = 0
        integer :: tensor_rank = 0
        real(dp) :: scalar = 0.0_dp
        real(dp) :: vector_value(2) = 0.0_dp
        real(dp), allocatable :: cell_values(:)
        real(dp), allocatable :: cell_vector_values(:, :)
        real(dp), allocatable :: cell_tensor_values(:, :, :)
    contains
        procedure, private :: assign_form_token
        generic :: assignment(=) => assign_form_token
    end type form_token_t

    type, public :: form_expr_t
        character(len=128) :: description = ""
        character(len=32) :: form_type = ""
        integer :: tensor_rank = 0
        type(form_token_t), allocatable :: tokens(:)
    contains
        procedure :: destroy => form_expr_destroy
    end type form_expr_t

    type :: compiler_item_t
        integer :: item_type = 0
        integer :: role = 0
        integer :: derivative = derivative_identity
        integer :: tensor_rank = 0
        integer :: field_rank = 0
        logical :: integrated = .false.
        real(dp) :: scalar = 0.0_dp
        real(dp) :: mass_coefficient = 0.0_dp
        real(dp) :: stiffness_coefficient = 0.0_dp
        real(dp) :: curl_coefficient = 0.0_dp
        real(dp) :: divergence_coefficient = 0.0_dp
        real(dp) :: load_coefficient = 0.0_dp
        real(dp) :: vector_load(2) = 0.0_dp
        real(dp) :: vector_value(2) = 0.0_dp
        logical :: has_vector_load = .false.
        integer :: vector_load_token = 0
        real(dp) :: vector_load_scale = 1.0_dp
        integer :: coefficient_token = 0
        integer :: mass_field_token = 0
        integer :: curl_field_token = 0
        integer :: divergence_field_token = 0
        real(dp) :: mass_field_scale = 0.0_dp
        real(dp) :: curl_field_scale = 0.0_dp
        real(dp) :: divergence_field_scale = 0.0_dp
    end type compiler_item_t

    public :: assignment(=)
    public :: compile_form, compile_form_matrix, compile_form_vector
    public :: compile_vector_form_csc, compile_vector_form_rhs
    public :: compile_vector_form_element
    public :: create_cell_coefficient, create_constant_load
    public :: create_cell_tensor_coefficient
    public :: create_cell_vector_source
    public :: create_curl, create_divergence
    public :: create_grad, create_inner
    public :: create_measure
    public :: create_product, create_scale, create_sum, create_symbol
    public :: create_vector_constant_function

    interface assignment(=)
        module procedure assign_form_expr
    end interface assignment(=)

contains

    impure elemental subroutine assign_form_token(lhs, rhs)
        class(form_token_t), intent(out) :: lhs
        type(form_token_t), intent(in) :: rhs

        lhs%token_type = rhs%token_type
        lhs%role = rhs%role
        lhs%tensor_rank = rhs%tensor_rank
        lhs%scalar = rhs%scalar
        lhs%vector_value = rhs%vector_value
        if (allocated(rhs%cell_values)) then
            allocate(lhs%cell_values, source=rhs%cell_values)
        end if
        if (allocated(rhs%cell_vector_values)) then
            allocate( &
                lhs%cell_vector_values, source=rhs%cell_vector_values)
        end if
        if (allocated(rhs%cell_tensor_values)) then
            allocate( &
                lhs%cell_tensor_values, source=rhs%cell_tensor_values)
        end if
    end subroutine assign_form_token

    subroutine assign_form_expr(lhs, rhs)
        type(form_expr_t), intent(inout) :: lhs
        type(form_expr_t), intent(in) :: rhs

        type(form_token_t), allocatable :: copied_tokens(:)
        character(len=128) :: description
        character(len=32) :: form_type
        integer :: tensor_rank

        description = rhs%description
        form_type = rhs%form_type
        tensor_rank = rhs%tensor_rank
        if (allocated(rhs%tokens)) then
            allocate(copied_tokens(size(rhs%tokens)))
            copied_tokens = rhs%tokens
        end if
        if (allocated(lhs%tokens)) deallocate(lhs%tokens)
        lhs%description = description
        lhs%form_type = form_type
        lhs%tensor_rank = tensor_rank
        if (allocated(copied_tokens)) call move_alloc(copied_tokens, lhs%tokens)
    end subroutine assign_form_expr

    function create_symbol(name, role_name, tensor_rank) result(expr)
        character(len=*), intent(in) :: name, role_name
        integer, intent(in) :: tensor_rank
        type(form_expr_t) :: expr

        allocate(expr%tokens(1))
        expr%description = name
        expr%form_type = role_name
        expr%tensor_rank = tensor_rank
        expr%tokens(1)%token_type = token_symbol
        expr%tokens(1)%role = role_from_name(role_name)
        expr%tokens(1)%tensor_rank = tensor_rank
    end function create_symbol

    function create_grad(func_name, func_type) result(expr)
        character(len=*), intent(in) :: func_name, func_type
        type(form_expr_t) :: expr
        type(form_expr_t) :: symbol

        symbol = create_symbol(func_name, func_type, 0)
        expr = append_unary_token(symbol, token_gradient)
        expr%description = "grad(" // trim(func_name) // ")"
        expr%form_type = func_type
        expr%tensor_rank = 1
    end function create_grad

    function create_curl(func_name, func_type) result(expr)
        character(len=*), intent(in) :: func_name, func_type
        type(form_expr_t) :: expr
        type(form_expr_t) :: symbol

        symbol = create_symbol(func_name, func_type, 1)
        expr = append_unary_token(symbol, token_curl)
        expr%description = "curl(" // trim(func_name) // ")"
        expr%form_type = func_type
        expr%tensor_rank = 0
    end function create_curl

    function create_divergence(func_name, func_type) result(expr)
        character(len=*), intent(in) :: func_name, func_type
        type(form_expr_t) :: expr
        type(form_expr_t) :: symbol

        symbol = create_symbol(func_name, func_type, 1)
        expr = append_unary_token(symbol, token_divergence)
        expr%description = "div(" // trim(func_name) // ")"
        expr%form_type = func_type
        expr%tensor_rank = 0
    end function create_divergence

    function create_inner(a, b) result(expr)
        type(form_expr_t), intent(in) :: a, b
        type(form_expr_t) :: expr

        expr = append_binary_token(a, b, token_inner)
        expr%description = "inner(" // trim(a%description) // ", " // &
            trim(b%description) // ")"
        expr%tensor_rank = 0
        if (has_role(a, role_trial) .and. has_role(b, role_test) .or. &
            has_role(a, role_test) .and. has_role(b, role_trial)) then
            expr%form_type = "bilinear"
        else if (has_role(a, role_test) .or. has_role(b, role_test)) then
            expr%form_type = "linear"
        else
            expr%form_type = "functional"
        end if
    end function create_inner

    function create_product(a, b) result(expr)
        type(form_expr_t), intent(in) :: a, b
        type(form_expr_t) :: expr

        expr = append_binary_token(a, b, token_multiply)
        expr%description = "(" // trim(a%description) // " * " // &
            trim(b%description) // ")"
        expr%form_type = product_form_type(a, b)
        expr%tensor_rank = a%tensor_rank + b%tensor_rank
    end function create_product

    function create_sum(a, b) result(expr)
        type(form_expr_t), intent(in) :: a, b
        type(form_expr_t) :: expr

        expr = append_binary_token(a, b, token_add)
        expr%description = "(" // trim(a%description) // " + " // &
            trim(b%description) // ")"
        if (trim(a%form_type) == trim(b%form_type)) then
            expr%form_type = a%form_type
        else
            expr%form_type = "unknown"
        end if
        expr%tensor_rank = max(a%tensor_rank, b%tensor_rank)
    end function create_sum

    function create_scale(scalar, a) result(expr)
        real(dp), intent(in) :: scalar
        type(form_expr_t), intent(in) :: a
        type(form_expr_t) :: expr
        type(form_expr_t) :: scalar_expr

        allocate(scalar_expr%tokens(1))
        write(scalar_expr%description, '(es16.8)') scalar
        scalar_expr%description = adjustl(scalar_expr%description)
        scalar_expr%form_type = "scalar"
        scalar_expr%tokens(1)%token_type = token_scalar
        scalar_expr%tokens(1)%scalar = scalar
        expr = create_product(scalar_expr, a)
    end function create_scale

    function create_measure() result(expr)
        type(form_expr_t) :: expr

        allocate(expr%tokens(1))
        expr%description = "dx"
        expr%form_type = "measure"
        expr%tensor_rank = 0
        expr%tokens(1)%token_type = token_measure
    end function create_measure

    function create_constant_load(value) result(expr)
        real(dp), intent(in) :: value
        type(form_expr_t) :: expr

        allocate(expr%tokens(1))
        write(expr%description, '("constant(",es16.8,")*v")') value
        expr%description = adjustl(expr%description)
        expr%form_type = "linear"
        expr%tensor_rank = 0
        expr%tokens(1)%token_type = token_constant_load
        expr%tokens(1)%scalar = value
    end function create_constant_load

    function create_cell_coefficient(values) result(expr)
        real(dp), intent(in) :: values(:)
        type(form_expr_t) :: expr

        allocate(expr%tokens(1))
        expr%description = "cell_coefficient"
        expr%form_type = "coefficient"
        expr%tensor_rank = 0
        expr%tokens(1)%token_type = token_cell_coefficient
        allocate(expr%tokens(1)%cell_values, source=values)
    end function create_cell_coefficient

    function create_cell_tensor_coefficient(values) result(expr)
        real(dp), intent(in) :: values(:, :, :)
        type(form_expr_t) :: expr

        allocate(expr%tokens(1))
        expr%description = "cell_tensor_coefficient"
        expr%form_type = "coefficient"
        expr%tensor_rank = 2
        expr%tokens(1)%token_type = token_cell_tensor_coefficient
        allocate(expr%tokens(1)%cell_tensor_values, source=values)
    end function create_cell_tensor_coefficient

    function create_vector_constant_function(value) result(expr)
        real(dp), intent(in) :: value(2)
        type(form_expr_t) :: expr

        allocate(expr%tokens(1))
        expr%description = "constant_vector"
        expr%form_type = "function"
        expr%tensor_rank = 1
        expr%tokens(1)%token_type = token_symbol
        expr%tokens(1)%role = role_function
        expr%tokens(1)%tensor_rank = 1
        expr%tokens(1)%vector_value = value
    end function create_vector_constant_function

    function create_cell_vector_source(values) result(expr)
        real(dp), intent(in) :: values(:, :)
        type(form_expr_t) :: expr

        allocate(expr%tokens(1))
        expr%description = "cell_vector_source"
        expr%form_type = "function"
        expr%tensor_rank = 1
        expr%tokens(1)%token_type = token_cell_vector_source
        expr%tokens(1)%role = role_function
        expr%tokens(1)%tensor_rank = 1
        allocate(expr%tokens(1)%cell_vector_values, source=values)
    end function create_cell_vector_source

    function compile_form(expr) result(assembly_code)
        type(form_expr_t), intent(in) :: expr
        character(len=256) :: assembly_code

        select case (trim(expr%form_type))
        case ("bilinear")
            assembly_code = "assemble_bilinear(" // trim(expr%description) // ")"
        case ("linear")
            assembly_code = "assemble_linear(" // trim(expr%description) // ")"
        case default
            assembly_code = "evaluate(" // trim(expr%description) // ")"
        end select
    end function compile_form

    subroutine compile_form_matrix( &
            expr, vertices, triangles, matrix, status)
        type(form_expr_t), intent(in) :: expr
        real(dp), intent(in) :: vertices(:, :)
        integer, intent(in) :: triangles(:, :)
        real(dp), intent(out) :: matrix(:, :)
        integer, intent(out) :: status

        real(dp) :: mass_coefficient, stiffness_coefficient
        integer :: compiler_status

        matrix = 0.0_dp
        status = 1
        call analyze_scalar_bilinear_form( &
            expr, mass_coefficient, stiffness_coefficient, compiler_status)
        if (compiler_status /= 0) return
        call assemble_scalar_p1_matrix( &
            vertices, triangles, mass_coefficient, stiffness_coefficient, &
            matrix, compiler_status)
        if (compiler_status /= 0) return
        status = 0
    end subroutine compile_form_matrix

    subroutine compile_form_vector( &
            expr, vertices, triangles, vector, status)
        type(form_expr_t), intent(in) :: expr
        real(dp), intent(in) :: vertices(:, :)
        integer, intent(in) :: triangles(:, :)
        real(dp), intent(out) :: vector(:)
        integer, intent(out) :: status

        real(dp) :: load_coefficient
        integer :: compiler_status

        vector = 0.0_dp
        status = 1
        call analyze_scalar_linear_form( &
            expr, load_coefficient, compiler_status)
        if (compiler_status /= 0) return
        call assemble_scalar_p1_load( &
            vertices, triangles, load_coefficient, vector, compiler_status)
        if (compiler_status /= 0) return
        status = 0
    end subroutine compile_form_vector

    subroutine compile_vector_form_element( &
            expr, family, degree, vertices, quadrature_degree, matrix, status)
        type(form_expr_t), intent(in) :: expr
        character(len=*), intent(in) :: family
        integer, intent(in) :: degree
        real(dp), intent(in) :: vertices(2, 3)
        integer, intent(in) :: quadrature_degree
        real(dp), allocatable, intent(out) :: matrix(:, :)
        integer, intent(out) :: status

        real(dp) :: curl_coefficient, divergence_coefficient
        real(dp) :: mass_coefficient
        real(dp) :: curl_field_scale, divergence_field_scale
        real(dp) :: mass_field_scale
        integer :: compiler_status, curl_field_token
        integer :: divergence_field_token, mass_field_token

        status = 1
        call analyze_vector_bilinear_form( &
            expr, mass_coefficient, curl_coefficient, &
            divergence_coefficient, mass_field_token, curl_field_token, &
            divergence_field_token, mass_field_scale, curl_field_scale, &
            divergence_field_scale, compiler_status)
        if (compiler_status /= 0) return
        if (mass_field_token /= 0 .or. curl_field_token /= 0 .or. &
            divergence_field_token /= 0) return
        select case (trim(family))
        case ("Nedelec", "Nedelec1", "Edge")
            if (divergence_coefficient /= 0.0_dp) then
                status = 3
                return
            end if
            call assemble_triangle_nedelec_curl_mass_element( &
                vertices, degree, quadrature_degree, matrix, status, &
                curl_coefficient=curl_coefficient, &
                mass_coefficient=mass_coefficient)
        case ("RT", "Raviart-Thomas")
            if (curl_coefficient /= 0.0_dp) then
                status = 3
                return
            end if
            call assemble_triangle_rt_div_mass_element( &
                vertices, degree, quadrature_degree, matrix, status, &
                divergence_coefficient=divergence_coefficient, &
                mass_coefficient=mass_coefficient)
        case ("Nedelec2", "Nedelec-second")
            if (divergence_coefficient /= 0.0_dp) then
                status = 3
                return
            end if
            call assemble_triangle_nedelec_second_curl_mass_element( &
                vertices, degree, quadrature_degree, matrix, status, &
                curl_coefficient=curl_coefficient, &
                mass_coefficient=mass_coefficient)
        case ("BDM")
            if (curl_coefficient /= 0.0_dp) then
                status = 3
                return
            end if
            call assemble_triangle_bdm_div_mass_element( &
                vertices, degree, quadrature_degree, matrix, status, &
                divergence_coefficient=divergence_coefficient, &
                mass_coefficient=mass_coefficient)
        case default
            status = 3
        end select
    end subroutine compile_vector_form_element

    subroutine compile_vector_form_csc( &
            expr, mesh, family, degree, quadrature_degree, matrix, status)
        type(form_expr_t), intent(in) :: expr
        type(mesh_2d_t), intent(inout) :: mesh
        character(len=*), intent(in) :: family
        integer, intent(in) :: degree, quadrature_degree
        type(csc_t), intent(out) :: matrix
        type(fortsparse_status_t), intent(out) :: status

        real(dp) :: curl_coefficient, divergence_coefficient
        real(dp) :: mass_coefficient
        real(dp), allocatable :: cell_curl(:), cell_mass(:)
        real(dp), allocatable :: cell_tensors(:, :, :)
        real(dp) :: curl_field_scale, divergence_field_scale
        real(dp) :: mass_field_scale
        integer :: compiler_status, curl_field_token
        integer :: divergence_field_token, mass_field_token

        call analyze_vector_bilinear_form( &
            expr, mass_coefficient, curl_coefficient, &
            divergence_coefficient, mass_field_token, curl_field_token, &
            divergence_field_token, mass_field_scale, curl_field_scale, &
            divergence_field_scale, compiler_status)
        if (compiler_status /= 0) then
            call status_set( &
                status, FORTSPARSE_INVALID_MATRIX, &
                "Sparse form compiler requires an integrated vector form")
            return
        end if
        select case (trim(family))
        case ("Nedelec", "Nedelec1", "Edge")
            if (divergence_coefficient /= 0.0_dp) then
                call set_incompatible_family_status(status)
                return
            end if
            if (divergence_field_token /= 0) then
                call set_incompatible_family_status(status)
                return
            end if
            if (is_tensor_token(expr, curl_field_token)) then
                call set_incompatible_family_status(status)
                return
            end if
            if (is_tensor_token(expr, mass_field_token)) then
                call build_cell_values( &
                    expr, mesh%n_triangles, curl_coefficient, &
                    curl_field_token, curl_field_scale, cell_curl, &
                    compiler_status)
                if (compiler_status /= 0) then
                    call set_invalid_cell_coefficient_status(status)
                    return
                end if
                call build_cell_tensors( &
                    expr, mesh%n_triangles, mass_coefficient, &
                    mass_field_token, mass_field_scale, cell_tensors, &
                    compiler_status)
                if (compiler_status /= 0) then
                    call set_invalid_cell_coefficient_status(status)
                    return
                end if
                call assemble_triangle_nedelec_cell_tensor_csc( &
                    mesh, degree, quadrature_degree, cell_curl, cell_tensors, &
                    matrix, status)
            else if (mass_field_token /= 0 .or. curl_field_token /= 0) then
                call build_cell_values( &
                    expr, mesh%n_triangles, curl_coefficient, &
                    curl_field_token, curl_field_scale, cell_curl, &
                    compiler_status)
                if (compiler_status /= 0) then
                    call set_invalid_cell_coefficient_status(status)
                    return
                end if
                call build_cell_values( &
                    expr, mesh%n_triangles, mass_coefficient, &
                    mass_field_token, mass_field_scale, cell_mass, &
                    compiler_status)
                if (compiler_status /= 0) then
                    call set_invalid_cell_coefficient_status(status)
                    return
                end if
                call assemble_cell_weighted_nedelec_csc( &
                    mesh, degree, quadrature_degree, cell_curl, cell_mass, &
                    matrix, status)
            else
                call assemble_triangle_nedelec_curl_mass_csc( &
                    mesh, degree, quadrature_degree, matrix, status, &
                    curl_coefficient, mass_coefficient)
            end if
        case ("RT", "Raviart-Thomas")
            if (curl_coefficient /= 0.0_dp) then
                call set_incompatible_family_status(status)
                return
            end if
            if (mass_field_token /= 0 .or. curl_field_token /= 0 .or. &
                divergence_field_token /= 0) then
                call set_incompatible_family_status(status)
                return
            end if
            call assemble_triangle_rt_div_mass_csc( &
                mesh, degree, quadrature_degree, matrix, status, &
                divergence_coefficient, mass_coefficient)
        case ("Nedelec2", "Nedelec-second")
            if (divergence_coefficient /= 0.0_dp) then
                call set_incompatible_family_status(status)
                return
            end if
            if (mass_field_token /= 0 .or. curl_field_token /= 0 .or. &
                divergence_field_token /= 0) then
                call set_incompatible_family_status(status)
                return
            end if
            call assemble_triangle_nedelec_second_curl_mass_csc( &
                mesh, degree, quadrature_degree, matrix, status, &
                curl_coefficient, mass_coefficient)
        case ("BDM")
            if (curl_coefficient /= 0.0_dp) then
                call set_incompatible_family_status(status)
                return
            end if
            if (mass_field_token /= 0 .or. curl_field_token /= 0 .or. &
                divergence_field_token /= 0) then
                call set_incompatible_family_status(status)
                return
            end if
            call assemble_triangle_bdm_div_mass_csc( &
                mesh, degree, quadrature_degree, matrix, status, &
                divergence_coefficient, mass_coefficient)
        case default
            call status_set( &
                status, FORTSPARSE_INVALID_MATRIX, &
                "Sparse form compiler received an unknown vector family")
        end select
    end subroutine compile_vector_form_csc

    subroutine compile_vector_form_rhs( &
            expr, mesh, family, degree, quadrature_degree, vector, status)
        type(form_expr_t), intent(in) :: expr
        type(mesh_2d_t), intent(inout) :: mesh
        character(len=*), intent(in) :: family
        integer, intent(in) :: degree, quadrature_degree
        real(dp), intent(out) :: vector(:)
        type(fortsparse_status_t), intent(out) :: status

        type(csc_t) :: mass_matrix
        real(dp), allocatable :: cell_values(:, :), source_dofs(:)
        real(dp) :: edge_vector(2), source_scale, source_value(2)
        integer :: compiler_status, dof, edge, source_token

        vector = 0.0_dp
        call analyze_vector_linear_form( &
            expr, source_value, source_token, source_scale, compiler_status)
        if (compiler_status /= 0) then
            call status_set( &
                status, FORTSPARSE_INVALID_MATRIX, &
                "Sparse vector load compiler received an unsupported source")
            return
        end if
        if (trim(family) /= "Nedelec" .and. &
            trim(family) /= "Nedelec1" .and. trim(family) /= "Edge" .and. &
            trim(family) /= "RT" .and. &
            trim(family) /= "Raviart-Thomas" .and. &
            trim(family) /= "BDM" .and. &
            trim(family) /= "Nedelec2" .and. &
            trim(family) /= "Nedelec-second") then
            call status_set( &
                status, FORTSPARSE_INVALID_MATRIX, &
                "Sparse vector load compiler received an unsupported family")
            return
        end if
        if (source_token /= 0) then
            if (.not. allocated( &
                expr%tokens(source_token)%cell_vector_values)) then
                call status_set( &
                    status, FORTSPARSE_INVALID_MATRIX, &
                    "Sparse vector load compiler lost its cell source")
                return
            end if
            if (trim(family) == "Nedelec2" .or. &
                trim(family) == "Nedelec-second") then
                call assemble_triangle_nedelec_second_cell_vector_load( &
                    mesh, degree, quadrature_degree, &
                    expr%tokens(source_token)%cell_vector_values, &
                    vector, status)
            else if (trim(family) == "BDM") then
                call assemble_triangle_bdm_cell_vector_load( &
                    mesh, degree, quadrature_degree, &
                    expr%tokens(source_token)%cell_vector_values, &
                    vector, status)
            else if (trim(family) == "RT" .or. &
                    trim(family) == "Raviart-Thomas") then
                call assemble_triangle_rt_cell_vector_load( &
                    mesh, degree, quadrature_degree, &
                    expr%tokens(source_token)%cell_vector_values, &
                    vector, status)
            else
                call assemble_triangle_nedelec_cell_vector_load( &
                    mesh, degree, quadrature_degree, &
                    expr%tokens(source_token)%cell_vector_values, &
                    vector, status)
            end if
            if (status%code == 0) vector = source_scale * vector
            return
        end if
        if (trim(family) == "Nedelec2" .or. &
            trim(family) == "Nedelec-second") then
            allocate(cell_values(2, mesh%n_triangles))
            cell_values = spread(source_value, 2, mesh%n_triangles)
            call assemble_triangle_nedelec_second_cell_vector_load( &
                mesh, degree, quadrature_degree, cell_values, vector, status)
            return
        end if
        if (trim(family) == "BDM") then
            allocate(cell_values(2, mesh%n_triangles))
            cell_values = spread(source_value, 2, mesh%n_triangles)
            call assemble_triangle_bdm_cell_vector_load( &
                mesh, degree, quadrature_degree, cell_values, vector, status)
            return
        end if
        if (trim(family) == "RT" .or. &
            trim(family) == "Raviart-Thomas") then
            allocate(cell_values(2, mesh%n_triangles))
            cell_values = spread(source_value, 2, mesh%n_triangles)
            call assemble_triangle_rt_cell_vector_load( &
                mesh, degree, quadrature_degree, cell_values, vector, status)
            return
        end if
        if (degree > 1) then
            allocate(cell_values(2, mesh%n_triangles))
            cell_values = spread(source_value, 2, mesh%n_triangles)
            call assemble_triangle_nedelec_cell_vector_load( &
                mesh, degree, quadrature_degree, cell_values, vector, status)
            return
        end if
        if (.not. allocated(mesh%edges)) call mesh%build_edge_connectivity()
        if (.not. allocated(mesh%edge_to_dof)) then
            call mesh%build_edge_dof_numbering()
        end if
        allocate(source_dofs(mesh%n_edges))
        source_dofs = 0.0_dp
        do edge = 1, mesh%n_edges
            edge_vector = mesh%vertices(:, mesh%edges(2, edge)) - &
                mesh%vertices(:, mesh%edges(1, edge))
            dof = mesh%edge_to_dof(edge) + 1
            source_dofs(dof) = dot_product(source_value, edge_vector)
        end do
        call assemble_triangle_nedelec_curl_mass_csc( &
            mesh, degree, quadrature_degree, mass_matrix, status, &
            0.0_dp, 1.0_dp)
        if (status%code /= 0) return
        if (size(vector) /= mass_matrix%nrow) then
            call status_set( &
                status, FORTSPARSE_INVALID_MATRIX, &
                "Sparse vector load output has the wrong dimension")
            return
        end if
        vector = csc_matvec(mass_matrix, source_dofs)
    end subroutine compile_vector_form_rhs

    subroutine set_incompatible_family_status(status)
        type(fortsparse_status_t), intent(out) :: status

        call status_set( &
            status, FORTSPARSE_INVALID_MATRIX, &
            "Vector differential operator is incompatible with the family")
    end subroutine set_incompatible_family_status

    subroutine set_invalid_cell_coefficient_status(status)
        type(fortsparse_status_t), intent(out) :: status

        call status_set( &
            status, FORTSPARSE_INVALID_MATRIX, &
            "Cell coefficient count must equal the triangle count")
    end subroutine set_invalid_cell_coefficient_status

    subroutine build_cell_values( &
            expr, triangle_count, constant_value, field_token, field_scale, &
            values, status)
        type(form_expr_t), intent(in) :: expr
        integer, intent(in) :: triangle_count, field_token
        real(dp), intent(in) :: constant_value, field_scale
        real(dp), allocatable, intent(out) :: values(:)
        integer, intent(out) :: status

        allocate(values(triangle_count))
        values = constant_value
        status = 2
        if (field_token == 0) then
            status = 0
            return
        end if
        if (.not. allocated(expr%tokens)) return
        if (field_token < 1 .or. field_token > size(expr%tokens)) return
        if (.not. allocated(expr%tokens(field_token)%cell_values)) return
        if (size(expr%tokens(field_token)%cell_values) /= triangle_count) return
        values = values + &
            field_scale * expr%tokens(field_token)%cell_values
        status = 0
    end subroutine build_cell_values

    pure logical function is_tensor_token(expr, token_index) result(is_tensor)
        type(form_expr_t), intent(in) :: expr
        integer, intent(in) :: token_index

        is_tensor = .false.
        if (token_index == 0) return
        if (.not. allocated(expr%tokens)) return
        if (token_index < 1 .or. token_index > size(expr%tokens)) return
        is_tensor = allocated(expr%tokens(token_index)%cell_tensor_values)
    end function is_tensor_token

    subroutine build_cell_tensors( &
            expr, triangle_count, constant_value, field_token, field_scale, &
            values, status)
        type(form_expr_t), intent(in) :: expr
        integer, intent(in) :: triangle_count, field_token
        real(dp), intent(in) :: constant_value, field_scale
        real(dp), allocatable, intent(out) :: values(:, :, :)
        integer, intent(out) :: status

        integer :: triangle

        allocate(values(2, 2, triangle_count))
        values = 0.0_dp
        do triangle = 1, triangle_count
            values(1, 1, triangle) = constant_value
            values(2, 2, triangle) = constant_value
        end do
        status = 2
        if (field_token == 0) then
            status = 0
            return
        end if
        if (.not. allocated(expr%tokens)) return
        if (field_token < 1 .or. field_token > size(expr%tokens)) return
        if (.not. allocated( &
            expr%tokens(field_token)%cell_tensor_values)) return
        if (size( &
            expr%tokens(field_token)%cell_tensor_values, 3) /= &
            triangle_count) return
        values = values + field_scale * &
            expr%tokens(field_token)%cell_tensor_values
        status = 0
    end subroutine build_cell_tensors

    subroutine assemble_cell_weighted_nedelec_csc( &
            mesh, degree, quadrature_degree, curl_values, mass_values, &
            matrix, status)
        type(mesh_2d_t), intent(inout) :: mesh
        integer, intent(in) :: degree, quadrature_degree
        real(dp), intent(in) :: curl_values(:), mass_values(:)
        type(csc_t), intent(out) :: matrix
        type(fortsparse_status_t), intent(out) :: status

        integer, allocatable :: columns(:), rows(:)
        real(dp), allocatable :: values(:)
        real(dp), allocatable :: element_matrix(:, :)
        real(dp) :: triangle_vertices(2, 3)
        integer :: edge_dofs(3), edge_orientations(3)
        integer :: entry, first_basis, second_basis
        integer :: assembly_status, triangle

        if (degree /= 1) then
            call status_set( &
                status, FORTSPARSE_INVALID_MATRIX, &
                "Cell-weighted Nedelec assembly supports order one")
            return
        end if
        if (size(curl_values) /= mesh%n_triangles .or. &
            size(mass_values) /= mesh%n_triangles) then
            call set_invalid_cell_coefficient_status(status)
            return
        end if
        if (.not. allocated(mesh%edges)) call mesh%build_edge_connectivity()
        if (.not. allocated(mesh%edge_to_dof)) then
            call mesh%build_edge_dof_numbering()
        end if
        allocate(rows(9 * mesh%n_triangles))
        allocate(columns(9 * mesh%n_triangles))
        allocate(values(9 * mesh%n_triangles))

        entry = 0
        do triangle = 1, mesh%n_triangles
            triangle_vertices = mesh%vertices(:, mesh%triangles(:, triangle))
            call assemble_triangle_nedelec_curl_mass_element( &
                triangle_vertices, degree, quadrature_degree, &
                element_matrix, assembly_status, &
                curl_coefficient=curl_values(triangle), &
                mass_coefficient=mass_values(triangle))
            if (assembly_status /= 0) then
                call status_set( &
                    status, FORTSPARSE_INVALID_MATRIX, &
                    "Cell-weighted Nedelec element assembly failed")
                return
            end if
            call mesh%get_triangle_edge_dofs( &
                triangle, edge_dofs, edge_orientations)
            do second_basis = 1, 3
                do first_basis = 1, 3
                    entry = entry + 1
                    rows(entry) = edge_dofs(first_basis) + 1
                    columns(entry) = edge_dofs(second_basis) + 1
                    values(entry) = real( &
                        edge_orientations(first_basis) * &
                        edge_orientations(second_basis), dp) * &
                        element_matrix(first_basis, second_basis)
                end do
            end do
        end do
        call csc_from_triplet( &
            mesh%n_edges, mesh%n_edges, rows, columns, values, matrix, status)
    end subroutine assemble_cell_weighted_nedelec_csc

    subroutine analyze_scalar_bilinear_form( &
            expr, mass_coefficient, stiffness_coefficient, status)
        type(form_expr_t), intent(in) :: expr
        real(dp), intent(out) :: mass_coefficient, stiffness_coefficient
        integer, intent(out) :: status

        type(compiler_item_t), allocatable :: stack(:)
        integer :: stack_size, token_index

        mass_coefficient = 0.0_dp
        stiffness_coefficient = 0.0_dp
        status = 2
        if (.not. allocated(expr%tokens)) return
        if (size(expr%tokens) < 1) return
        allocate(stack(size(expr%tokens)))
        stack_size = 0
        do token_index = 1, size(expr%tokens)
            call apply_compiler_token( &
                expr%tokens(token_index), token_index, &
                stack, stack_size, status)
            if (status /= 0) return
        end do
        if (stack_size /= 1) then
            status = 2
            return
        end if
        if (stack(1)%item_type /= item_form) then
            status = 2
            return
        end if
        if (.not. stack(1)%integrated) then
            status = 2
            return
        end if
        if (stack(1)%field_rank /= 0) then
            status = 2
            return
        end if
        if (stack(1)%curl_coefficient /= 0.0_dp) then
            status = 2
            return
        end if
        if (stack(1)%divergence_coefficient /= 0.0_dp) then
            status = 2
            return
        end if
        if (stack(1)%load_coefficient /= 0.0_dp) then
            status = 2
            return
        end if
        mass_coefficient = stack(1)%mass_coefficient
        stiffness_coefficient = stack(1)%stiffness_coefficient
        status = 0
    end subroutine analyze_scalar_bilinear_form

    subroutine analyze_scalar_linear_form(expr, load_coefficient, status)
        type(form_expr_t), intent(in) :: expr
        real(dp), intent(out) :: load_coefficient
        integer, intent(out) :: status

        type(compiler_item_t), allocatable :: stack(:)
        integer :: stack_size, token_index

        load_coefficient = 0.0_dp
        status = 2
        if (.not. allocated(expr%tokens)) return
        if (size(expr%tokens) < 1) return
        allocate(stack(size(expr%tokens)))
        stack_size = 0
        do token_index = 1, size(expr%tokens)
            call apply_compiler_token( &
                expr%tokens(token_index), token_index, &
                stack, stack_size, status)
            if (status /= 0) return
        end do
        if (stack_size /= 1) then
            status = 2
            return
        end if
        if (stack(1)%item_type /= item_form) then
            status = 2
            return
        end if
        if (.not. stack(1)%integrated) then
            status = 2
            return
        end if
        if (stack(1)%field_rank /= 0) then
            status = 2
            return
        end if
        if (stack(1)%mass_coefficient /= 0.0_dp .or. &
            stack(1)%stiffness_coefficient /= 0.0_dp) then
            status = 2
            return
        end if
        if (stack(1)%curl_coefficient /= 0.0_dp .or. &
            stack(1)%divergence_coefficient /= 0.0_dp) then
            status = 2
            return
        end if
        load_coefficient = stack(1)%load_coefficient
        status = 0
    end subroutine analyze_scalar_linear_form

    subroutine analyze_vector_bilinear_form( &
            expr, mass_coefficient, curl_coefficient, divergence_coefficient, &
            mass_field_token, curl_field_token, divergence_field_token, &
            mass_field_scale, curl_field_scale, divergence_field_scale, status)
        type(form_expr_t), intent(in) :: expr
        real(dp), intent(out) :: mass_coefficient, curl_coefficient
        real(dp), intent(out) :: divergence_coefficient
        integer, intent(out) :: mass_field_token, curl_field_token
        integer, intent(out) :: divergence_field_token
        real(dp), intent(out) :: mass_field_scale, curl_field_scale
        real(dp), intent(out) :: divergence_field_scale
        integer, intent(out) :: status

        type(compiler_item_t), allocatable :: stack(:)
        integer :: stack_size, token_index

        mass_coefficient = 0.0_dp
        curl_coefficient = 0.0_dp
        divergence_coefficient = 0.0_dp
        mass_field_token = 0
        curl_field_token = 0
        divergence_field_token = 0
        mass_field_scale = 0.0_dp
        curl_field_scale = 0.0_dp
        divergence_field_scale = 0.0_dp
        status = 2
        if (.not. allocated(expr%tokens)) return
        if (size(expr%tokens) < 1) return
        allocate(stack(size(expr%tokens)))
        stack_size = 0
        do token_index = 1, size(expr%tokens)
            call apply_compiler_token( &
                expr%tokens(token_index), token_index, &
                stack, stack_size, status)
            if (status /= 0) return
        end do
        if (stack_size /= 1) then
            status = 2
            return
        end if
        if (stack(1)%item_type /= item_form) then
            status = 2
            return
        end if
        if (.not. stack(1)%integrated) then
            status = 2
            return
        end if
        if (stack(1)%field_rank /= 1) then
            status = 2
            return
        end if
        if (stack(1)%stiffness_coefficient /= 0.0_dp) then
            status = 2
            return
        end if
        if (stack(1)%has_vector_load) then
            status = 2
            return
        end if
        mass_coefficient = stack(1)%mass_coefficient
        curl_coefficient = stack(1)%curl_coefficient
        divergence_coefficient = stack(1)%divergence_coefficient
        mass_field_token = stack(1)%mass_field_token
        curl_field_token = stack(1)%curl_field_token
        divergence_field_token = stack(1)%divergence_field_token
        mass_field_scale = stack(1)%mass_field_scale
        curl_field_scale = stack(1)%curl_field_scale
        divergence_field_scale = stack(1)%divergence_field_scale
        status = 0
    end subroutine analyze_vector_bilinear_form

    subroutine analyze_vector_linear_form( &
            expr, source_value, source_token, source_scale, status)
        type(form_expr_t), intent(in) :: expr
        real(dp), intent(out) :: source_value(2)
        integer, intent(out) :: source_token, status
        real(dp), intent(out) :: source_scale

        type(compiler_item_t), allocatable :: stack(:)
        integer :: stack_size, token_index

        source_value = 0.0_dp
        source_token = 0
        source_scale = 1.0_dp
        status = 2
        if (.not. allocated(expr%tokens)) return
        if (size(expr%tokens) < 1) return
        allocate(stack(size(expr%tokens)))
        stack_size = 0
        do token_index = 1, size(expr%tokens)
            call apply_compiler_token( &
                expr%tokens(token_index), token_index, &
                stack, stack_size, status)
            if (status /= 0) return
        end do
        if (stack_size /= 1) then
            status = 2
            return
        end if
        if (stack(1)%item_type /= item_form) then
            status = 2
            return
        end if
        if (.not. stack(1)%integrated .or. &
            .not. stack(1)%has_vector_load) then
            status = 2
            return
        end if
        if (stack(1)%field_rank /= 1) then
            status = 2
            return
        end if
        if (stack(1)%mass_coefficient /= 0.0_dp .or. &
            stack(1)%stiffness_coefficient /= 0.0_dp) then
            status = 2
            return
        end if
        if (stack(1)%curl_coefficient /= 0.0_dp .or. &
            stack(1)%divergence_coefficient /= 0.0_dp) then
            status = 2
            return
        end if
        source_value = stack(1)%vector_load
        source_token = stack(1)%vector_load_token
        source_scale = stack(1)%vector_load_scale
        if (source_token /= 0 .and. any(source_value /= 0.0_dp)) return
        status = 0
    end subroutine analyze_vector_linear_form

    subroutine apply_compiler_token( &
            token, token_index, stack, stack_size, status)
        type(form_token_t), intent(in) :: token
        integer, intent(in) :: token_index
        type(compiler_item_t), intent(inout) :: stack(:)
        integer, intent(inout) :: stack_size
        integer, intent(out) :: status

        status = 2
        select case (token%token_type)
        case (token_symbol)
            stack_size = stack_size + 1
            stack(stack_size)%item_type = item_argument
            stack(stack_size)%role = token%role
            stack(stack_size)%tensor_rank = token%tensor_rank
            stack(stack_size)%field_rank = token%tensor_rank
            stack(stack_size)%derivative = derivative_identity
            stack(stack_size)%vector_value = token%vector_value
        case (token_cell_vector_source)
            stack_size = stack_size + 1
            stack(stack_size)%item_type = item_argument
            stack(stack_size)%role = role_function
            stack(stack_size)%tensor_rank = 1
            stack(stack_size)%field_rank = 1
            stack(stack_size)%derivative = derivative_identity
            stack(stack_size)%vector_load_token = token_index
        case (token_gradient)
            if (stack_size < 1) return
            if (stack(stack_size)%item_type /= item_argument) return
            if (stack(stack_size)%tensor_rank /= 0) return
            stack(stack_size)%tensor_rank = 1
            stack(stack_size)%derivative = derivative_gradient
        case (token_curl)
            if (stack_size < 1) return
            if (stack(stack_size)%item_type /= item_argument) return
            if (stack(stack_size)%tensor_rank /= 1) return
            stack(stack_size)%tensor_rank = 0
            stack(stack_size)%derivative = derivative_curl
        case (token_divergence)
            if (stack_size < 1) return
            if (stack(stack_size)%item_type /= item_argument) return
            if (stack(stack_size)%tensor_rank /= 1) return
            stack(stack_size)%tensor_rank = 0
            stack(stack_size)%derivative = derivative_divergence
        case (token_inner)
            call apply_inner_token(stack, stack_size, status)
            return
        case (token_scalar)
            stack_size = stack_size + 1
            stack(stack_size)%item_type = item_scalar
            stack(stack_size)%scalar = token%scalar
        case (token_measure)
            stack_size = stack_size + 1
            stack(stack_size)%item_type = item_measure
        case (token_constant_load)
            stack_size = stack_size + 1
            stack(stack_size)%item_type = item_form
            stack(stack_size)%field_rank = 0
            stack(stack_size)%load_coefficient = token%scalar
        case (token_cell_coefficient, token_cell_tensor_coefficient)
            stack_size = stack_size + 1
            stack(stack_size)%item_type = item_coefficient
            stack(stack_size)%coefficient_token = token_index
        case (token_multiply)
            call apply_multiply_token(stack, stack_size, status)
            return
        case (token_add)
            call apply_add_token(stack, stack_size, status)
            return
        case default
            return
        end select
        status = 0
    end subroutine apply_compiler_token

    subroutine apply_inner_token(stack, stack_size, status)
        type(compiler_item_t), intent(inout) :: stack(:)
        integer, intent(inout) :: stack_size
        integer, intent(out) :: status

        type(compiler_item_t) :: first, second

        status = 2
        if (stack_size < 2) return
        first = stack(stack_size - 1)
        second = stack(stack_size)
        if (first%item_type /= item_argument .or. &
            second%item_type /= item_argument) return
        if (function_test_roles(first%role, second%role)) then
            if (first%derivative /= derivative_identity .or. &
                second%derivative /= derivative_identity) return
            if (first%field_rank /= 1 .or. second%field_rank /= 1) return
            stack_size = stack_size - 1
            stack(stack_size) = compiler_item_t()
            stack(stack_size)%item_type = item_form
            stack(stack_size)%field_rank = 1
            stack(stack_size)%has_vector_load = .true.
            if (first%role == role_function) then
                stack(stack_size)%vector_load = first%vector_value
                stack(stack_size)%vector_load_token = first%vector_load_token
            else
                stack(stack_size)%vector_load = second%vector_value
                stack(stack_size)%vector_load_token = second%vector_load_token
            end if
            status = 0
            return
        end if
        if (.not. complementary_roles(first%role, second%role)) return
        if (first%derivative /= second%derivative) return
        if (first%tensor_rank /= second%tensor_rank) return
        stack_size = stack_size - 1
        stack(stack_size) = compiler_item_t()
        stack(stack_size)%item_type = item_form
        stack(stack_size)%field_rank = first%field_rank
        select case (first%derivative)
        case (derivative_identity)
            stack(stack_size)%mass_coefficient = 1.0_dp
        case (derivative_gradient)
            if (first%tensor_rank /= 1) return
            stack(stack_size)%stiffness_coefficient = 1.0_dp
        case (derivative_curl)
            if (first%tensor_rank /= 0) return
            stack(stack_size)%curl_coefficient = 1.0_dp
        case (derivative_divergence)
            if (first%tensor_rank /= 0) return
            stack(stack_size)%divergence_coefficient = 1.0_dp
        case default
            return
        end select
        status = 0
    end subroutine apply_inner_token

    subroutine apply_multiply_token(stack, stack_size, status)
        type(compiler_item_t), intent(inout) :: stack(:)
        integer, intent(inout) :: stack_size
        integer, intent(out) :: status

        type(compiler_item_t) :: first, second
        real(dp) :: factor

        status = 2
        if (stack_size < 2) return
        first = stack(stack_size - 1)
        second = stack(stack_size)
        if (first%item_type == item_scalar .and. &
            second%item_type == item_form) then
            factor = first%scalar
            stack(stack_size - 1) = second
        else if (first%item_type == item_form .and. &
                second%item_type == item_scalar) then
            factor = second%scalar
            stack(stack_size - 1) = first
        else if (first%item_type == item_coefficient .and. &
                second%item_type == item_form) then
            stack(stack_size - 1) = second
            call apply_cell_coefficient( &
                stack(stack_size - 1), first%coefficient_token, status)
            if (status /= 0) return
            stack_size = stack_size - 1
            return
        else if (first%item_type == item_form .and. &
                second%item_type == item_coefficient) then
            stack(stack_size - 1) = first
            call apply_cell_coefficient( &
                stack(stack_size - 1), second%coefficient_token, status)
            if (status /= 0) return
            stack_size = stack_size - 1
            return
        else if (first%item_type == item_measure .and. &
                second%item_type == item_form) then
            if (second%integrated) return
            stack(stack_size - 1) = second
            stack(stack_size - 1)%integrated = .true.
            stack_size = stack_size - 1
            status = 0
            return
        else if (first%item_type == item_form .and. &
                second%item_type == item_measure) then
            if (first%integrated) return
            stack(stack_size - 1) = first
            stack(stack_size - 1)%integrated = .true.
            stack_size = stack_size - 1
            status = 0
            return
        else
            return
        end if
        stack(stack_size - 1)%mass_coefficient = factor * &
            stack(stack_size - 1)%mass_coefficient
        stack(stack_size - 1)%stiffness_coefficient = factor * &
            stack(stack_size - 1)%stiffness_coefficient
        stack(stack_size - 1)%curl_coefficient = factor * &
            stack(stack_size - 1)%curl_coefficient
        stack(stack_size - 1)%divergence_coefficient = factor * &
            stack(stack_size - 1)%divergence_coefficient
        stack(stack_size - 1)%mass_field_scale = factor * &
            stack(stack_size - 1)%mass_field_scale
        stack(stack_size - 1)%curl_field_scale = factor * &
            stack(stack_size - 1)%curl_field_scale
        stack(stack_size - 1)%divergence_field_scale = factor * &
            stack(stack_size - 1)%divergence_field_scale
        stack(stack_size - 1)%load_coefficient = factor * &
            stack(stack_size - 1)%load_coefficient
        stack(stack_size - 1)%vector_load = factor * &
            stack(stack_size - 1)%vector_load
        stack(stack_size - 1)%vector_load_scale = factor * &
            stack(stack_size - 1)%vector_load_scale
        stack_size = stack_size - 1
        status = 0
    end subroutine apply_multiply_token

    subroutine apply_cell_coefficient(item, coefficient_token, status)
        type(compiler_item_t), intent(inout) :: item
        integer, intent(in) :: coefficient_token
        integer, intent(out) :: status

        status = 2
        if (item%mass_coefficient /= 0.0_dp) then
            if (item%mass_field_token /= 0) return
            item%mass_field_token = coefficient_token
            item%mass_field_scale = item%mass_coefficient
            item%mass_coefficient = 0.0_dp
        end if
        if (item%curl_coefficient /= 0.0_dp) then
            if (item%curl_field_token /= 0) return
            item%curl_field_token = coefficient_token
            item%curl_field_scale = item%curl_coefficient
            item%curl_coefficient = 0.0_dp
        end if
        if (item%divergence_coefficient /= 0.0_dp) then
            if (item%divergence_field_token /= 0) return
            item%divergence_field_token = coefficient_token
            item%divergence_field_scale = item%divergence_coefficient
            item%divergence_coefficient = 0.0_dp
        end if
        if (item%mass_field_token == 0 .and. &
            item%curl_field_token == 0 .and. &
            item%divergence_field_token == 0) return
        status = 0
    end subroutine apply_cell_coefficient

    subroutine apply_add_token(stack, stack_size, status)
        type(compiler_item_t), intent(inout) :: stack(:)
        integer, intent(inout) :: stack_size
        integer, intent(out) :: status

        type(compiler_item_t) :: first, second

        status = 2
        if (stack_size < 2) return
        first = stack(stack_size - 1)
        second = stack(stack_size)
        if (first%item_type /= item_form .or. &
            second%item_type /= item_form) return
        if (first%integrated .neqv. second%integrated) return
        if (first%field_rank /= second%field_rank) return
        stack(stack_size - 1) = first
        stack(stack_size - 1)%mass_coefficient = &
            first%mass_coefficient + second%mass_coefficient
        stack(stack_size - 1)%stiffness_coefficient = &
            first%stiffness_coefficient + second%stiffness_coefficient
        stack(stack_size - 1)%curl_coefficient = &
            first%curl_coefficient + second%curl_coefficient
        stack(stack_size - 1)%divergence_coefficient = &
            first%divergence_coefficient + second%divergence_coefficient
        call combine_field_terms( &
            first%mass_field_token, first%mass_field_scale, &
            second%mass_field_token, second%mass_field_scale, &
            stack(stack_size - 1)%mass_field_token, &
            stack(stack_size - 1)%mass_field_scale, status)
        if (status /= 0) return
        call combine_field_terms( &
            first%curl_field_token, first%curl_field_scale, &
            second%curl_field_token, second%curl_field_scale, &
            stack(stack_size - 1)%curl_field_token, &
            stack(stack_size - 1)%curl_field_scale, status)
        if (status /= 0) return
        call combine_field_terms( &
            first%divergence_field_token, first%divergence_field_scale, &
            second%divergence_field_token, second%divergence_field_scale, &
            stack(stack_size - 1)%divergence_field_token, &
            stack(stack_size - 1)%divergence_field_scale, status)
        if (status /= 0) return
        stack(stack_size - 1)%load_coefficient = &
            first%load_coefficient + second%load_coefficient
        stack(stack_size - 1)%vector_load = &
            first%vector_load + second%vector_load
        call combine_field_terms( &
            first%vector_load_token, first%vector_load_scale, &
            second%vector_load_token, second%vector_load_scale, &
            stack(stack_size - 1)%vector_load_token, &
            stack(stack_size - 1)%vector_load_scale, status)
        if (status /= 0) return
        stack(stack_size - 1)%has_vector_load = &
            first%has_vector_load .or. second%has_vector_load
        stack_size = stack_size - 1
        status = 0
    end subroutine apply_add_token

    pure subroutine combine_field_terms( &
            first_token, first_scale, second_token, second_scale, &
            combined_token, combined_scale, status)
        integer, intent(in) :: first_token, second_token
        real(dp), intent(in) :: first_scale, second_scale
        integer, intent(out) :: combined_token, status
        real(dp), intent(out) :: combined_scale

        status = 2
        combined_token = 0
        combined_scale = 0.0_dp
        if (first_token == 0) then
            combined_token = second_token
            combined_scale = second_scale
        else if (second_token == 0) then
            combined_token = first_token
            combined_scale = first_scale
        else if (first_token == second_token) then
            combined_token = first_token
            combined_scale = first_scale + second_scale
        else
            return
        end if
        status = 0
    end subroutine combine_field_terms

    subroutine assemble_scalar_p1_matrix( &
            vertices, triangles, mass_coefficient, stiffness_coefficient, &
            matrix, status)
        real(dp), intent(in) :: vertices(:, :)
        integer, intent(in) :: triangles(:, :)
        real(dp), intent(in) :: mass_coefficient, stiffness_coefficient
        real(dp), intent(out) :: matrix(:, :)
        integer, intent(out) :: status

        real(dp) :: area, determinant, gradients(2, 3), local(3, 3)
        real(dp) :: mass(3, 3), x(3), y(3)
        integer :: first_basis, first_node, second_basis, second_node, triangle

        matrix = 0.0_dp
        status = 3
        if (size(vertices, 1) /= 2 .or. size(vertices, 2) < 3) return
        if (size(triangles, 1) /= 3 .or. size(triangles, 2) < 1) return
        if (size(matrix, 1) /= size(vertices, 2) .or. &
            size(matrix, 2) /= size(vertices, 2)) return
        if (any(triangles < 1) .or. &
            any(triangles > size(vertices, 2))) return
        mass = 1.0_dp
        mass(1, 1) = 2.0_dp
        mass(2, 2) = 2.0_dp
        mass(3, 3) = 2.0_dp
        do triangle = 1, size(triangles, 2)
            x = vertices(1, triangles(:, triangle))
            y = vertices(2, triangles(:, triangle))
            determinant = (x(2) - x(1)) * (y(3) - y(1)) - &
                (x(3) - x(1)) * (y(2) - y(1))
            if (abs(determinant) <= 64.0_dp * epsilon(1.0_dp)) return
            area = 0.5_dp * abs(determinant)
            gradients(1, :) = [y(2) - y(3), y(3) - y(1), y(1) - y(2)] / &
                determinant
            gradients(2, :) = [x(3) - x(2), x(1) - x(3), x(2) - x(1)] / &
                determinant
            local = stiffness_coefficient * area * &
                matmul(transpose(gradients), gradients) + &
                mass_coefficient * area * mass / 12.0_dp
            do second_basis = 1, 3
                second_node = triangles(second_basis, triangle)
                do first_basis = 1, 3
                    first_node = triangles(first_basis, triangle)
                    matrix(first_node, second_node) = &
                        matrix(first_node, second_node) + &
                        local(first_basis, second_basis)
                end do
            end do
        end do
        status = 0
    end subroutine assemble_scalar_p1_matrix

    subroutine assemble_scalar_p1_load( &
            vertices, triangles, load_coefficient, vector, status)
        real(dp), intent(in) :: vertices(:, :)
        integer, intent(in) :: triangles(:, :)
        real(dp), intent(in) :: load_coefficient
        real(dp), intent(out) :: vector(:)
        integer, intent(out) :: status

        real(dp) :: area, determinant, x(3), y(3)
        integer :: basis, node, triangle

        vector = 0.0_dp
        status = 3
        if (size(vertices, 1) /= 2 .or. size(vertices, 2) < 3) return
        if (size(triangles, 1) /= 3 .or. size(triangles, 2) < 1) return
        if (size(vector) /= size(vertices, 2)) return
        if (any(triangles < 1) .or. &
            any(triangles > size(vertices, 2))) return
        do triangle = 1, size(triangles, 2)
            x = vertices(1, triangles(:, triangle))
            y = vertices(2, triangles(:, triangle))
            determinant = (x(2) - x(1)) * (y(3) - y(1)) - &
                (x(3) - x(1)) * (y(2) - y(1))
            if (abs(determinant) <= 64.0_dp * epsilon(1.0_dp)) return
            area = 0.5_dp * abs(determinant)
            do basis = 1, 3
                node = triangles(basis, triangle)
                vector(node) = vector(node) + load_coefficient * area / 3.0_dp
            end do
        end do
        status = 0
    end subroutine assemble_scalar_p1_load

    function append_unary_token(a, token_type) result(expr)
        type(form_expr_t), intent(in) :: a
        integer, intent(in) :: token_type
        type(form_expr_t) :: expr
        integer :: token_count

        token_count = expression_token_count(a)
        allocate(expr%tokens(token_count + 1))
        if (token_count > 0) expr%tokens(1:token_count) = a%tokens
        expr%tokens(token_count + 1)%token_type = token_type
    end function append_unary_token

    function append_binary_token(a, b, token_type) result(expr)
        type(form_expr_t), intent(in) :: a, b
        integer, intent(in) :: token_type
        type(form_expr_t) :: expr
        integer :: first_count, second_count

        first_count = expression_token_count(a)
        second_count = expression_token_count(b)
        allocate(expr%tokens(first_count + second_count + 1))
        if (first_count > 0) expr%tokens(1:first_count) = a%tokens
        if (second_count > 0) then
            expr%tokens(first_count + 1:first_count + second_count) = b%tokens
        end if
        expr%tokens(first_count + second_count + 1)%token_type = token_type
    end function append_binary_token

    pure integer function expression_token_count(expr) result(count)
        type(form_expr_t), intent(in) :: expr

        count = 0
        if (allocated(expr%tokens)) count = size(expr%tokens)
    end function expression_token_count

    pure integer function role_from_name(role_name) result(role)
        character(len=*), intent(in) :: role_name

        select case (trim(role_name))
        case ("trial")
            role = role_trial
        case ("test")
            role = role_test
        case ("function")
            role = role_function
        case default
            role = 0
        end select
    end function role_from_name

    pure logical function has_role(expr, role) result(found)
        type(form_expr_t), intent(in) :: expr
        integer, intent(in) :: role
        integer :: token_index

        found = .false.
        if (.not. allocated(expr%tokens)) return
        do token_index = 1, size(expr%tokens)
            if (expr%tokens(token_index)%token_type == token_symbol .and. &
                expr%tokens(token_index)%role == role) then
                found = .true.
                return
            end if
        end do
    end function has_role

    pure logical function complementary_roles(first, second) result(valid)
        integer, intent(in) :: first, second

        valid = (first == role_trial .and. second == role_test) .or. &
            (first == role_test .and. second == role_trial)
    end function complementary_roles

    pure logical function function_test_roles(first, second) result(valid)
        integer, intent(in) :: first, second

        valid = (first == role_function .and. second == role_test) .or. &
            (first == role_test .and. second == role_function)
    end function function_test_roles

    pure function product_form_type(a, b) result(form_type)
        type(form_expr_t), intent(in) :: a, b
        character(len=32) :: form_type

        if (trim(a%form_type) == "measure") then
            form_type = b%form_type
        else if (trim(b%form_type) == "measure") then
            form_type = a%form_type
        else if (trim(a%form_type) == "scalar") then
            form_type = b%form_type
        else if (trim(b%form_type) == "scalar") then
            form_type = a%form_type
        else if (trim(a%form_type) == "coefficient") then
            form_type = b%form_type
        else if (trim(b%form_type) == "coefficient") then
            form_type = a%form_type
        else
            form_type = a%form_type
        end if
    end function product_form_type

    subroutine form_expr_destroy(this)
        class(form_expr_t), intent(inout) :: this

        if (allocated(this%tokens)) deallocate(this%tokens)
        this%description = ""
        this%form_type = ""
        this%tensor_rank = 0
    end subroutine form_expr_destroy

end module fortfem_forms_simple
