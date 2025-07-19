module fortfem_forms_simple
    ! Simplified forms system for mathematical expressions
    use fortfem_kinds
    implicit none
    
    private
    public :: form_expr_t, create_grad, create_inner, compile_form
    
    ! Simple form expression type
    type :: form_expr_t
        character(len=128) :: description = ""
        character(len=32) :: form_type = ""  ! "bilinear", "linear", "functional"
        integer :: tensor_rank = 0
    contains
        procedure :: destroy => form_expr_destroy
    end type form_expr_t
    
contains

    function create_grad(func_name, func_type) result(expr)
        character(len=*), intent(in) :: func_name, func_type
        type(form_expr_t) :: expr
        
        expr%description = "grad(" // trim(func_name) // ")"
        expr%form_type = func_type
        expr%tensor_rank = 1
    end function create_grad
    
    function create_inner(a, b) result(expr)
        type(form_expr_t), intent(in) :: a, b
        type(form_expr_t) :: expr
        
        expr%description = "inner(" // trim(a%description) // ", " // trim(b%description) // ")"
        expr%tensor_rank = 0
        
        ! Determine form type
        if ((trim(a%form_type) == "trial" .or. trim(a%form_type) == "bilinear") .and. &
            (trim(b%form_type) == "test" .or. trim(b%form_type) == "bilinear")) then
            expr%form_type = "bilinear"
        else if (trim(a%form_type) == "test" .or. trim(b%form_type) == "test") then
            expr%form_type = "linear"
        else
            expr%form_type = "functional"
        end if
    end function create_inner
    
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
    
    subroutine form_expr_destroy(this)
        class(form_expr_t), intent(inout) :: this
        ! Nothing to deallocate for simple strings
    end subroutine form_expr_destroy

end module fortfem_forms_simple