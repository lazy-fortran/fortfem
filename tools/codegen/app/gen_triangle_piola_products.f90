program gen_triangle_piola_products
    use fortsym_arena, only: arena_t
    use fortsym_engine, only: engine_result_t
    use fortsym_engine_native, only: make_native_engine, native_engine_t
    use fortsym_expr, only: expr_t, operator(*), operator(+), operator(-), operator(/), sym
    use fortsym_products, only: jvp, vjp
    use fortsym_subs, only: subs_many
    use fortsym_kernel, only: emit_kernel, kernel_spec_t, KERNEL_SUBROUTINE, &
        KERNEL_SNIPPET, CSE_NONE
    use fortsym_string, only: chars, str, str_t
    use fortfem_codegen_provenance, only: fortsym_revision, generated_path
    implicit none
    type(arena_t), target :: arena
    type(native_engine_t) :: engine
    type(expr_t) :: inputs(7), dots(7), bars(3), determinant, mapped(3)
    character(3), parameter :: input_names(7) = &
        ["j11", "j21", "j12", "j22", "v1 ", "v2 ", "s  "]
    type(str_t) :: args(14)
    integer :: i, family, unit
    character(:), allocatable :: family_name

    call arena%init()
    engine = make_native_engine(arena)
    do i = 1, 7
        inputs(i) = sym(arena, trim(input_names(i)))
        dots(i) = sym(arena, trim(input_names(i))//"_dot")
    end do
    do i = 1, 3
        bars(i) = sym(arena, "out"//integer_text(i)//"_bar")
    end do
    determinant = inputs(1)*inputs(4) - inputs(3)*inputs(2)
    open (newunit=unit, &
        file=generated_path("fortfem_triangle_piola_products.f90"), &
        status="replace", action="write")
    do family = 1, 2
        if (family == 1) then
            family_name = "triangle_covariant"
            mapped(1) = (inputs(4)*inputs(5) - inputs(2)*inputs(6))/determinant
            mapped(2) = (inputs(1)*inputs(6) - inputs(3)*inputs(5))/determinant
        else
            family_name = "triangle_contravariant"
            mapped(1) = (inputs(1)*inputs(5) + inputs(3)*inputs(6))/determinant
            mapped(2) = (inputs(2)*inputs(5) + inputs(4)*inputs(6))/determinant
        end if
        mapped(3) = inputs(7)/determinant
        do i = 1, 7
            args(i) = str(trim(input_names(i)))
            args(7 + i) = str(trim(input_names(i))//"_dot")
        end do
        call emit_inline_product(family_name, mapped, family, .false.)
        call emit_inline_product(family_name, jvp(mapped, inputs, dots), &
            family, .true.)
        do i = 1, 3
            args(7 + i) = str("out"//integer_text(i)//"_bar")
        end do
        call emit(family_name//"_vjp", vjp(mapped, inputs, bars), args(:10))
    end do
    close (unit)
contains
    subroutine emit_inline_product(name, expressions, family_id, tangent_mode)
        character(*), intent(in) :: name
        type(expr_t), intent(in) :: expressions(:)
        integer, intent(in) :: family_id
        logical, intent(in) :: tangent_mode
        type(expr_t) :: roots(3), replacements(14), symbols(14)
        type(engine_result_t) :: result
        type(kernel_spec_t) :: spec
        character(:), allocatable :: code, scalar_name, output_suffix, file_suffix
        integer :: i, include_unit, ninput

        scalar_name = "curls"
        if (family_id == 2) scalar_name = "divergences"
        replacements(1) = sym(arena, "jacobian(1,1)")
        replacements(2) = sym(arena, "jacobian(2,1)")
        replacements(3) = sym(arena, "jacobian(1,2)")
        replacements(4) = sym(arena, "jacobian(2,2)")
        replacements(5) = sym(arena, "reference_values(1,basis_dof)")
        replacements(6) = sym(arena, "reference_values(2,basis_dof)")
        replacements(7) = sym(arena, "reference_"//scalar_name//"(basis_dof)")
        symbols(:7) = inputs
        ninput = 7
        output_suffix = ""
        file_suffix = "_primal"
        if (tangent_mode) then
            symbols(8:14) = dots
            replacements(8) = sym(arena, "jacobian_dot(1,1)")
            replacements(9) = sym(arena, "jacobian_dot(2,1)")
            replacements(10) = sym(arena, "jacobian_dot(1,2)")
            replacements(11) = sym(arena, "jacobian_dot(2,2)")
            replacements(12) = sym(arena, "reference_values_dot(1,basis_dof)")
            replacements(13) = sym(arena, "reference_values_dot(2,basis_dof)")
            replacements(14) = sym(arena, "reference_"//scalar_name//"_dot(basis_dof)")
            ninput = 14
            output_suffix = "_dot"
            file_suffix = "_jvp"
        end if
        do i = 1, 3
            result = engine%simplify(expressions(i))
            if (.not. result%ok) error stop "Piola primal simplification failed"
            roots(i) = subs_many(result%value, symbols(:ninput), replacements(:ninput))
        end do
        spec%name = str("generated_"//name//"_inline")
        spec%mode = KERNEL_SNIPPET
        spec%cse_level = CSE_NONE
        spec%generator = str("gen_triangle_piola_products")
        spec%generator_revision = str(fortsym_revision())
        spec%regenerate_command = str("cd tools/codegen && ./generate.sh")
        allocate(spec%args(6), spec%outputs(2), spec%output_references(3))
        spec%args = [str("jacobian"), str("reference_values"), &
            str("reference_"//scalar_name), str("jacobian_dot"), &
            str("reference_values_dot"), str("reference_"//scalar_name//"_dot")]
        spec%outputs = [str("physical_values"//output_suffix), &
            str("physical_"//scalar_name//output_suffix)]
        spec%output_references = [str("physical_values"//output_suffix//"(1,basis_dof)"), &
            str("physical_values"//output_suffix//"(2,basis_dof)"), &
            str("physical_"//scalar_name//output_suffix//"(basis_dof)")]
        code = chars(emit_kernel(roots, spec))
        open(newunit=include_unit, &
            file=generated_path("fortfem_"//name//file_suffix//".inc"), &
            status="replace", action="write")
        write(include_unit, "(a)") code(:len(code) - 1)
        close(include_unit)
    end subroutine emit_inline_product

    subroutine emit(name, expressions, arguments)
        character(*), intent(in) :: name
        type(expr_t), intent(in) :: expressions(:)
        type(str_t), intent(in) :: arguments(:)
        type(expr_t) :: roots(size(expressions))
        type(engine_result_t) :: result
        type(kernel_spec_t) :: spec
        character(:), allocatable :: code
        integer :: root
        do root = 1, size(expressions)
            result = engine%simplify(expressions(root))
            if (.not. result%ok) error stop "Piola simplification failed"
            roots(root) = result%value
        end do
        spec%name = str("generated_"//name)
        spec%module_name = str("fortfem_generated_"//name)
        spec%mode = KERNEL_SUBROUTINE
        spec%generator = str("gen_triangle_piola_products")
        spec%generator_revision = str(fortsym_revision())
        spec%regenerate_command = str("cd tools/codegen && ./generate.sh")
        spec%pure_procedure = .true.
        spec%args = arguments
        allocate(spec%outputs(1), spec%output_shapes(1))
        allocate(spec%output_references(size(roots)))
        spec%outputs = [str("product")]
        spec%output_shapes = [str("("//integer_text(size(roots))//")")]
        do root = 1, size(roots)
            spec%output_references(root) = str("product("//integer_text(root)//")")
        end do
        code = chars(emit_kernel(roots, spec))
        write (unit, "(a)") code(:len(code) - 1)
    end subroutine emit
    function integer_text(value) result(text)
        integer, intent(in) :: value
        character(:), allocatable :: text
        character(32) :: buffer
        write (buffer, "(i0)") value
        text = trim(buffer)
    end function integer_text
end program gen_triangle_piola_products
