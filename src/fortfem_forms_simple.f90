module fortfem_forms_simple
    use fortfem_assembly_nedelec_arbitrary_order_2d, only: &
        assemble_triangle_nedelec_curl_mass_element
    use fortfem_assembly_full_vector_arbitrary_order_2d, only: &
        assemble_triangle_bdm_div_mass_element, &
        assemble_triangle_nedelec_second_curl_mass_element
    use fortfem_assembly_rt_arbitrary_order_2d, only: &
        assemble_triangle_rt_div_mass_element
    use fortfem_kinds, only: dp
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

    type :: form_token_t
        integer :: token_type = 0
        integer :: role = 0
        integer :: tensor_rank = 0
        real(dp) :: scalar = 0.0_dp
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
    end type compiler_item_t

    public :: assignment(=)
    public :: compile_form, compile_form_matrix, compile_vector_form_element
    public :: create_curl, create_divergence, create_grad, create_inner
    public :: create_measure
    public :: create_product, create_scale, create_sum, create_symbol

    interface assignment(=)
        module procedure assign_form_expr
    end interface assignment(=)

contains

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
        integer :: compiler_status

        status = 1
        call analyze_vector_bilinear_form( &
            expr, mass_coefficient, curl_coefficient, &
            divergence_coefficient, compiler_status)
        if (compiler_status /= 0) return
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
                expr%tokens(token_index), stack, stack_size, status)
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
        mass_coefficient = stack(1)%mass_coefficient
        stiffness_coefficient = stack(1)%stiffness_coefficient
        status = 0
    end subroutine analyze_scalar_bilinear_form

    subroutine analyze_vector_bilinear_form( &
            expr, mass_coefficient, curl_coefficient, divergence_coefficient, &
            status)
        type(form_expr_t), intent(in) :: expr
        real(dp), intent(out) :: mass_coefficient, curl_coefficient
        real(dp), intent(out) :: divergence_coefficient
        integer, intent(out) :: status

        type(compiler_item_t), allocatable :: stack(:)
        integer :: stack_size, token_index

        mass_coefficient = 0.0_dp
        curl_coefficient = 0.0_dp
        divergence_coefficient = 0.0_dp
        status = 2
        if (.not. allocated(expr%tokens)) return
        if (size(expr%tokens) < 1) return
        allocate(stack(size(expr%tokens)))
        stack_size = 0
        do token_index = 1, size(expr%tokens)
            call apply_compiler_token( &
                expr%tokens(token_index), stack, stack_size, status)
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
        mass_coefficient = stack(1)%mass_coefficient
        curl_coefficient = stack(1)%curl_coefficient
        divergence_coefficient = stack(1)%divergence_coefficient
        status = 0
    end subroutine analyze_vector_bilinear_form

    subroutine apply_compiler_token(token, stack, stack_size, status)
        type(form_token_t), intent(in) :: token
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
        stack_size = stack_size - 1
        status = 0
    end subroutine apply_multiply_token

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
        stack_size = stack_size - 1
        status = 0
    end subroutine apply_add_token

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
