program gen_fci_curved_quadrilateral_area_products
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
    type(expr_t) :: variables(16), variables_dot(16), area(1), area_bar(1)
    type(expr_t) :: area_jvp(1), area_vjp(16)
    type(expr_t) :: start_x, start_y, control_x, control_y
    type(expr_t) :: end_x, end_y, coefficient_a_x, coefficient_a_y
    type(expr_t) :: coefficient_b_x, coefficient_b_y
    type(expr_t) :: half, two, three
    type(kernel_spec_t) :: spec
    character(:), allocatable :: primal_code, jvp_code, vjp_code
    integer :: edge, next_edge, index, ios, unit

    call arena%init()
    engine = make_native_engine(arena)
    half = real_expr(arena, 0.5d0)
    two = real_expr(arena, 2.0d0)
    three = real_expr(arena, 3.0d0)
    variables = [ &
        sym(arena, "x_1"), sym(arena, "y_1"), &
        sym(arena, "x_2"), sym(arena, "y_2"), &
        sym(arena, "x_3"), sym(arena, "y_3"), &
        sym(arena, "x_4"), sym(arena, "y_4"), &
        sym(arena, "control_x_1"), sym(arena, "control_y_1"), &
        sym(arena, "control_x_2"), sym(arena, "control_y_2"), &
        sym(arena, "control_x_3"), sym(arena, "control_y_3"), &
        sym(arena, "control_x_4"), sym(arena, "control_y_4")]
    variables_dot = [ &
        sym(arena, "x_1_dot"), sym(arena, "y_1_dot"), &
        sym(arena, "x_2_dot"), sym(arena, "y_2_dot"), &
        sym(arena, "x_3_dot"), sym(arena, "y_3_dot"), &
        sym(arena, "x_4_dot"), sym(arena, "y_4_dot"), &
        sym(arena, "control_x_1_dot"), sym(arena, "control_y_1_dot"), &
        sym(arena, "control_x_2_dot"), sym(arena, "control_y_2_dot"), &
        sym(arena, "control_x_3_dot"), sym(arena, "control_y_3_dot"), &
        sym(arena, "control_x_4_dot"), sym(arena, "control_y_4_dot")]
    area(1) = real_expr(arena, 0.0d0)
    do edge = 1, 4
        next_edge = mod(edge, 4) + 1
        start_x = variables(2*edge - 1)
        start_y = variables(2*edge)
        end_x = variables(2*next_edge - 1)
        end_y = variables(2*next_edge)
        control_x = variables(8 + 2*edge - 1)
        control_y = variables(8 + 2*edge)
        coefficient_a_x = start_x - two*control_x + end_x
        coefficient_a_y = start_y - two*control_y + end_y
        coefficient_b_x = two*(control_x - start_x)
        coefficient_b_y = two*(control_y - start_y)
        area(1) = area(1) + half*( &
            (coefficient_b_x*coefficient_a_y - &
            coefficient_b_y*coefficient_a_x)/three + &
            start_x*coefficient_a_y - start_y*coefficient_a_x + &
            start_x*coefficient_b_y - start_y*coefficient_b_x)
    end do
    call simplify_all(area)

    call initialize_spec(spec, "generated_fci_curved_quadrilateral_cell_area", &
        "fortfem_generated_fci_curved_quadrilateral_area", 16, 1)
    spec%args = primal_arguments()
    spec%outputs(1) = str("area")
    primal_code = chars(emit_kernel(area, spec))

    area_bar(1) = sym(arena, "area_bar")
    area_jvp = jvp(area, variables, variables_dot)
    area_vjp = vjp(area, variables, area_bar)
    call simplify_all(area_jvp)
    call simplify_all(area_vjp)

    call initialize_spec(spec, &
        "generated_fci_curved_quadrilateral_cell_area_jvp", &
        "fortfem_generated_fci_curved_quadrilateral_area_jvp", 32, 1)
    spec%args = jvp_arguments()
    spec%outputs(1) = str("area_dot")
    jvp_code = chars(emit_kernel(area_jvp, spec))

    call initialize_spec(spec, &
        "generated_fci_curved_quadrilateral_cell_area_vjp", &
        "fortfem_generated_fci_curved_quadrilateral_area_vjp", 17, 16)
    spec%args = vjp_arguments()
    spec%outputs = vjp_outputs()
    vjp_code = chars(emit_kernel(area_vjp, spec))

    open (newunit=unit, file=generated_path( &
        "fortfem_fci_curved_quadrilateral_area_products.f90"), &
        status="replace", action="write", iostat=ios)
    if (ios /= 0) error stop "cannot write FCI curved quadrilateral products"
    write (unit, "(a)") primal_code(:len(primal_code) - 1)
    write (unit, "(a)") jvp_code(:len(jvp_code) - 1)
    write (unit, "(a)") vjp_code(:len(vjp_code) - 1)
    close (unit)

contains

    function primal_arguments() result(arguments)
        use fortsym_string, only: str_t
        type(str_t), allocatable :: arguments(:)

        allocate(arguments(16))
        arguments = [ &
            str("x_1"), str("y_1"), str("x_2"), str("y_2"), &
            str("x_3"), str("y_3"), str("x_4"), str("y_4"), &
            str("control_x_1"), str("control_y_1"), &
            str("control_x_2"), str("control_y_2"), &
            str("control_x_3"), str("control_y_3"), &
            str("control_x_4"), str("control_y_4")]
    end function primal_arguments

    function jvp_arguments() result(arguments)
        use fortsym_string, only: str_t
        type(str_t), allocatable :: arguments(:)

        allocate(arguments(32))
        arguments(1:16) = primal_arguments()
        arguments(17:32) = [ &
            str("x_1_dot"), str("y_1_dot"), str("x_2_dot"), &
            str("y_2_dot"), str("x_3_dot"), str("y_3_dot"), &
            str("x_4_dot"), str("y_4_dot"), str("control_x_1_dot"), &
            str("control_y_1_dot"), str("control_x_2_dot"), &
            str("control_y_2_dot"), str("control_x_3_dot"), &
            str("control_y_3_dot"), str("control_x_4_dot"), &
            str("control_y_4_dot")]
    end function jvp_arguments

    function vjp_arguments() result(arguments)
        use fortsym_string, only: str_t
        type(str_t), allocatable :: arguments(:)

        allocate(arguments(17))
        arguments(1:16) = primal_arguments()
        arguments(17) = str("area_bar")
    end function vjp_arguments

    function vjp_outputs() result(outputs)
        use fortsym_string, only: str_t
        type(str_t), allocatable :: outputs(:)

        allocate(outputs(16))
        outputs = [ &
            str("x_1_bar"), str("y_1_bar"), str("x_2_bar"), str("y_2_bar"), &
            str("x_3_bar"), str("y_3_bar"), str("x_4_bar"), str("y_4_bar"), &
            str("control_x_1_bar"), str("control_y_1_bar"), &
            str("control_x_2_bar"), str("control_y_2_bar"), &
            str("control_x_3_bar"), str("control_y_3_bar"), &
            str("control_x_4_bar"), str("control_y_4_bar")]
    end function vjp_outputs

    subroutine initialize_spec(kernel_spec, name, module_name, argument_count, &
            output_count)
        type(kernel_spec_t), intent(out) :: kernel_spec
        character(*), intent(in) :: name, module_name
        integer, intent(in) :: argument_count, output_count

        kernel_spec%name = str(name)
        kernel_spec%module_name = str(module_name)
        kernel_spec%mode = KERNEL_SUBROUTINE
        kernel_spec%generator = str("gen_fci_curved_quadrilateral_area_products")
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

end program gen_fci_curved_quadrilateral_area_products
