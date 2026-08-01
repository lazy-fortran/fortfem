program gen_fci_polygon_edge_area_products
    use fortsym_arena, only: arena_t
    use fortsym_engine, only: engine_result_t
    use fortsym_engine_native, only: make_native_engine, native_engine_t
    use fortsym_expr, only: expr_t, operator(-), operator(*), &
        operator(/), real_expr, sym
    use fortsym_kernel, only: emit_kernel, kernel_spec_t, KERNEL_SUBROUTINE
    use fortsym_products, only: jvp, vjp
    use fortsym_string, only: chars, str
    use fortfem_codegen_provenance, only: fortsym_revision, generated_path
    implicit none

    type(arena_t), target :: arena
    type(engine_result_t) :: result
    type(native_engine_t) :: engine
    type(expr_t) :: variables(4), variables_dot(4), area(1), area_bar(1)
    type(expr_t) :: area_jvp(1), area_vjp(4), two
    type(kernel_spec_t) :: spec
    character(:), allocatable :: primal_code, jvp_code, vjp_code
    integer :: index, ios, unit

    call arena%init()
    engine = make_native_engine(arena)
    two = real_expr(arena, 2.0d0)
    variables = [ &
        sym(arena, "x_start"), sym(arena, "y_start"), &
        sym(arena, "x_end"), sym(arena, "y_end")]
    area(1) = (variables(1)*variables(4) - variables(2)*variables(3))/two
    call simplify_all(area)

    call initialize_spec(spec, "generated_fci_polygon_edge_area", &
        "fortfem_generated_fci_polygon_edge_area", 4, 1)
    spec%args = primal_arguments()
    spec%outputs(1) = str("edge_area")
    primal_code = chars(emit_kernel(area, spec))

    variables_dot = [ &
        sym(arena, "x_start_dot"), sym(arena, "y_start_dot"), &
        sym(arena, "x_end_dot"), sym(arena, "y_end_dot")]
    area_bar(1) = sym(arena, "edge_area_bar")
    area_jvp = jvp(area, variables, variables_dot)
    area_vjp = vjp(area, variables, area_bar)
    call simplify_all(area_jvp)
    call simplify_all(area_vjp)

    call initialize_spec(spec, "generated_fci_polygon_edge_area_jvp", &
        "fortfem_generated_fci_polygon_edge_area_jvp", 8, 1)
    spec%args = jvp_arguments()
    spec%outputs(1) = str("edge_area_dot")
    jvp_code = chars(emit_kernel(area_jvp, spec))

    call initialize_spec(spec, "generated_fci_polygon_edge_area_vjp", &
        "fortfem_generated_fci_polygon_edge_area_vjp", 5, 4)
    spec%args = vjp_arguments()
    spec%outputs = vjp_outputs()
    vjp_code = chars(emit_kernel(area_vjp, spec))

    open (newunit=unit, file=generated_path( &
        "fortfem_fci_polygon_edge_area_products.f90"), status="replace", &
        action="write", iostat=ios)
    if (ios /= 0) error stop "cannot write FCI polygon edge area products"
    write (unit, "(a)") primal_code(:len(primal_code) - 1)
    write (unit, "(a)") jvp_code(:len(jvp_code) - 1)
    write (unit, "(a)") vjp_code(:len(vjp_code) - 1)
    close (unit)

contains

    function primal_arguments() result(arguments)
        use fortsym_string, only: str_t
        type(str_t), allocatable :: arguments(:)

        allocate(arguments(4))
        arguments = [str("x_start"), str("y_start"), str("x_end"), &
            str("y_end")]
    end function primal_arguments

    function jvp_arguments() result(arguments)
        use fortsym_string, only: str_t
        type(str_t), allocatable :: arguments(:)

        allocate(arguments(8))
        arguments(1:4) = primal_arguments()
        arguments(5:8) = [str("x_start_dot"), str("y_start_dot"), &
            str("x_end_dot"), str("y_end_dot")]
    end function jvp_arguments

    function vjp_arguments() result(arguments)
        use fortsym_string, only: str_t
        type(str_t), allocatable :: arguments(:)

        allocate(arguments(5))
        arguments(1:4) = primal_arguments()
        arguments(5) = str("edge_area_bar")
    end function vjp_arguments

    function vjp_outputs() result(outputs)
        use fortsym_string, only: str_t
        type(str_t), allocatable :: outputs(:)

        allocate(outputs(4))
        outputs = [str("x_start_bar"), str("y_start_bar"), &
            str("x_end_bar"), str("y_end_bar")]
    end function vjp_outputs

    subroutine initialize_spec(kernel_spec, name, module_name, argument_count, &
            output_count)
        type(kernel_spec_t), intent(out) :: kernel_spec
        character(*), intent(in) :: name, module_name
        integer, intent(in) :: argument_count, output_count

        kernel_spec%name = str(name)
        kernel_spec%module_name = str(module_name)
        kernel_spec%mode = KERNEL_SUBROUTINE
        kernel_spec%generator = str("gen_fci_polygon_edge_area_products")
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

end program gen_fci_polygon_edge_area_products
