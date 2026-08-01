program gen_fci_quadrilateral_area_products
    use fortsym_arena, only: arena_t
    use fortsym_engine, only: engine_result_t
    use fortsym_engine_native, only: make_native_engine, native_engine_t
    use fortsym_expr, only: expr_t, operator(+), operator(-), operator(*), &
        operator(/), real_expr, sym
    use fortsym_kernel, only: emit_kernel, kernel_spec_t, KERNEL_SUBROUTINE
    use fortsym_products, only: jvp, vjp
    use fortsym_string, only: chars, str
    use fortfem_codegen_provenance, only: fortsym_revision, generated_path
    implicit none

    type(arena_t), target :: arena
    type(engine_result_t) :: result
    type(native_engine_t) :: engine
    type(expr_t) :: variables(8), variables_dot(8), area(1), two
    type(expr_t) :: area_jvp(1), area_vjp(8), area_bar(1)
    type(kernel_spec_t) :: spec
    character(:), allocatable :: primal_code, jvp_code, vjp_code
    integer :: index, ios, unit

    call arena%init()
    engine = make_native_engine(arena)
    two = real_expr(arena, 2.0d0)
    variables(1) = sym(arena, "x_1")
    variables(2) = sym(arena, "y_1")
    variables(3) = sym(arena, "x_2")
    variables(4) = sym(arena, "y_2")
    variables(5) = sym(arena, "x_3")
    variables(6) = sym(arena, "y_3")
    variables(7) = sym(arena, "x_4")
    variables(8) = sym(arena, "y_4")
    area(1) = (variables(1)*variables(4) + &
        variables(3)*variables(6) + variables(5)*variables(8) + &
        variables(7)*variables(2) - variables(2)*variables(3) - &
        variables(4)*variables(5) - variables(6)*variables(7) - &
        variables(8)*variables(1))/two
    call simplify_all(area)

    call initialize_spec(spec, "generated_fci_quadrilateral_cell_area", &
        "fortfem_generated_fci_quadrilateral_area", 8, 1)
    spec%args(1) = str("x_1")
    spec%args(2) = str("y_1")
    spec%args(3) = str("x_2")
    spec%args(4) = str("y_2")
    spec%args(5) = str("x_3")
    spec%args(6) = str("y_3")
    spec%args(7) = str("x_4")
    spec%args(8) = str("y_4")
    spec%outputs(1) = str("area")
    primal_code = chars(emit_kernel(area, spec))

    variables_dot(1) = sym(arena, "x_1_dot")
    variables_dot(2) = sym(arena, "y_1_dot")
    variables_dot(3) = sym(arena, "x_2_dot")
    variables_dot(4) = sym(arena, "y_2_dot")
    variables_dot(5) = sym(arena, "x_3_dot")
    variables_dot(6) = sym(arena, "y_3_dot")
    variables_dot(7) = sym(arena, "x_4_dot")
    variables_dot(8) = sym(arena, "y_4_dot")
    area_bar(1) = sym(arena, "area_bar")
    area_jvp = jvp(area, variables, variables_dot)
    area_vjp = vjp(area, variables, area_bar)
    call simplify_all(area_jvp)
    call simplify_all(area_vjp)

    call initialize_spec(spec, "generated_fci_quadrilateral_cell_area_jvp", &
        "fortfem_generated_fci_quadrilateral_area_jvp", 16, 1)
    spec%args(1) = str("x_1")
    spec%args(2) = str("y_1")
    spec%args(3) = str("x_2")
    spec%args(4) = str("y_2")
    spec%args(5) = str("x_3")
    spec%args(6) = str("y_3")
    spec%args(7) = str("x_4")
    spec%args(8) = str("y_4")
    spec%args(9) = str("x_1_dot")
    spec%args(10) = str("y_1_dot")
    spec%args(11) = str("x_2_dot")
    spec%args(12) = str("y_2_dot")
    spec%args(13) = str("x_3_dot")
    spec%args(14) = str("y_3_dot")
    spec%args(15) = str("x_4_dot")
    spec%args(16) = str("y_4_dot")
    spec%outputs(1) = str("area_dot")
    jvp_code = chars(emit_kernel(area_jvp, spec))

    call initialize_spec(spec, "generated_fci_quadrilateral_cell_area_vjp", &
        "fortfem_generated_fci_quadrilateral_area_vjp", 9, 8)
    spec%args(1) = str("x_1")
    spec%args(2) = str("y_1")
    spec%args(3) = str("x_2")
    spec%args(4) = str("y_2")
    spec%args(5) = str("x_3")
    spec%args(6) = str("y_3")
    spec%args(7) = str("x_4")
    spec%args(8) = str("y_4")
    spec%args(9) = str("area_bar")
    spec%outputs(1) = str("x_1_bar")
    spec%outputs(2) = str("y_1_bar")
    spec%outputs(3) = str("x_2_bar")
    spec%outputs(4) = str("y_2_bar")
    spec%outputs(5) = str("x_3_bar")
    spec%outputs(6) = str("y_3_bar")
    spec%outputs(7) = str("x_4_bar")
    spec%outputs(8) = str("y_4_bar")
    vjp_code = chars(emit_kernel(area_vjp, spec))

    open (newunit=unit, file=generated_path( &
        "fortfem_fci_quadrilateral_area_products.f90"), status="replace", &
        action="write", iostat=ios)
    if (ios /= 0) error stop "cannot write FCI quadrilateral area products"
    write (unit, "(a)") primal_code(:len(primal_code) - 1)
    write (unit, "(a)") jvp_code(:len(jvp_code) - 1)
    write (unit, "(a)") vjp_code(:len(vjp_code) - 1)
    close (unit)

contains

    subroutine initialize_spec(kernel_spec, name, module_name, argument_count, &
            output_count)
        type(kernel_spec_t), intent(out) :: kernel_spec
        character(*), intent(in) :: name, module_name
        integer, intent(in) :: argument_count, output_count

        kernel_spec%name = str(name)
        kernel_spec%module_name = str(module_name)
        kernel_spec%mode = KERNEL_SUBROUTINE
        kernel_spec%generator = str("gen_fci_quadrilateral_area_products")
        kernel_spec%generator_revision = str(fortsym_revision())
        kernel_spec%regenerate_command = str( &
            "cd tools/codegen && ./generate.sh")
        kernel_spec%pure_procedure = .true.
        allocate(kernel_spec%args(argument_count), kernel_spec%outputs(output_count))
    end subroutine initialize_spec

    subroutine simplify_all(expressions)
        type(expr_t), intent(inout) :: expressions(:)

        do index = 1, size(expressions)
            result = engine%simplify(expressions(index))
            if (result%ok) expressions(index) = result%value
        end do
    end subroutine simplify_all

end program gen_fci_quadrilateral_area_products
