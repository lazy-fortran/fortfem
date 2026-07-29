module fortfem_api_forms
    use fortfem_kinds
    use fortfem_forms_simple, only: assignment(=), compile_form, &
        compile_form_matrix, compile_vector_form_element, create_curl, &
        create_divergence, create_grad, create_inner, &
        create_measure, create_product, create_scale, create_sum, &
        create_symbol, form_expr_t
    use fortfem_api_types, only: trial_function_t, test_function_t,         &
        vector_trial_function_t, vector_test_function_t,                    &
        vector_function_t, function_t, simple_expression_t
    implicit none

    private

    public :: form_expr_t
    public :: simple_expression_t
    public :: form_equation_t
    public :: dx
    public :: init_measures
    public :: inner
    public :: grad
    public :: curl
    public :: div
    public :: compile_form
    public :: compile_form_matrix
    public :: compile_vector_form_element
    public :: assignment(=)
    public :: operator(*)
    public :: operator(+)
    public :: operator(==)

    type :: form_equation_t
        type(form_expr_t) :: lhs
        type(form_expr_t) :: rhs
    end type form_equation_t

    type(form_expr_t), save :: dx
    logical, save :: measures_initialized = .false.

    interface operator(*)
        module procedure expr_times_expr
        module procedure function_times_test
        module procedure function_times_expr
        module procedure vector_function_times_vector_test
        module procedure real_times_expr
        module procedure expr_times_real
    end interface operator(*)

    interface operator(+)
        module procedure expr_plus_expr
    end interface operator(+)

    interface operator(==)
        module procedure form_equals_form
    end interface operator(==)

contains

    subroutine init_measures()
        if (.not. measures_initialized) then
            dx = create_measure()
            measures_initialized = .true.
        end if
    end subroutine init_measures

    function inner(a, b) result(expr)
        class(*), intent(in) :: a, b
        type(form_expr_t) :: expr
        type(form_expr_t) :: expr_a, expr_b

        select type(a)
            type is (trial_function_t)
            expr_a = create_symbol("u", "trial", 0)
            type is (test_function_t)
            expr_a = create_symbol("v", "test", 0)
            type is (vector_trial_function_t)
            expr_a = create_symbol("E", "trial", 1)
            type is (vector_test_function_t)
            expr_a = create_symbol("F", "test", 1)
            type is (vector_function_t)
            expr_a = create_symbol("j", "function", 1)
            type is (form_expr_t)
            expr_a = a
        class default
            expr_a%description = "unknown"
            expr_a%form_type = "unknown"
        end select

        select type(b)
            type is (trial_function_t)
            expr_b = create_symbol("u", "trial", 0)
            type is (test_function_t)
            expr_b = create_symbol("v", "test", 0)
            type is (vector_trial_function_t)
            expr_b = create_symbol("E", "trial", 1)
            type is (vector_test_function_t)
            expr_b = create_symbol("F", "test", 1)
            type is (vector_function_t)
            expr_b = create_symbol("j", "function", 1)
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

    function curl(u) result(curlu)
        class(*), intent(in) :: u
        type(form_expr_t) :: curlu

        select type(u)
            type is (vector_trial_function_t)
            curlu = create_curl("u", "trial")
            type is (vector_test_function_t)
            curlu = create_curl("v", "test")
        class default
            curlu = create_curl("unknown", "unknown")
        end select
    end function curl

    function div(u) result(divu)
        class(*), intent(in) :: u
        type(form_expr_t) :: divu

        select type(u)
            type is (vector_trial_function_t)
            divu = create_divergence("u", "trial")
            type is (vector_test_function_t)
            divu = create_divergence("v", "test")
        class default
            divu = create_divergence("unknown", "unknown")
        end select
    end function div

    function expr_times_expr(a, b) result(product)
        type(form_expr_t), intent(in) :: a, b
        type(form_expr_t) :: product

        product = create_product(a, b)
    end function expr_times_expr

    function expr_plus_expr(a, b) result(sum_expr)
        type(form_expr_t), intent(in) :: a, b
        type(form_expr_t) :: sum_expr

        sum_expr = create_sum(a, b)
    end function expr_plus_expr

    function real_times_expr(a, b) result(product)
        real(dp), intent(in) :: a
        type(form_expr_t), intent(in) :: b
        type(form_expr_t) :: product

        product = create_scale(a, b)
    end function real_times_expr

    function expr_times_real(a, b) result(product)
        type(form_expr_t), intent(in) :: a
        real(dp), intent(in) :: b
        type(form_expr_t) :: product

        product = create_scale(b, a)
    end function expr_times_real

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

    function vector_function_times_vector_test(f, v) result(product)
        type(vector_function_t), intent(in) :: f
        type(vector_test_function_t), intent(in) :: v
        type(form_expr_t) :: product

        product%description = "f*v"
        product%form_type = "linear"
        product%tensor_rank = 0
    end function vector_function_times_vector_test

    function form_equals_form(a, L) result(equation)
        type(form_expr_t), intent(in) :: a, L
        type(form_equation_t) :: equation

        equation%lhs = a
        equation%rhs = L
    end function form_equals_form

end module fortfem_api_forms
