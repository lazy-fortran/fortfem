program gen_cgl_pressure_tensor_products
    use fortsym_arena, only: arena_t
    use fortsym_engine, only: engine_result_t
    use fortsym_engine_native, only: make_native_engine, native_engine_t
    use fortsym_expr, only: expr_t, operator(*), operator(+), operator(-), sym
    use fortsym_kernel, only: emit_kernel, kernel_spec_t, KERNEL_SUBROUTINE
    use fortsym_products, only: jvp, vjp
    use fortsym_string, only: chars, str
    use fortfem_codegen_provenance, only: fortsym_revision, generated_path
    implicit none

    type(arena_t), target :: arena
    type(engine_result_t) :: result
    type(native_engine_t) :: engine
    type(expr_t) :: variables(5), variables_dot(5)
    type(expr_t) :: tensor(6), tensor_jvp(6), tensor_vjp(5)
    type(expr_t) :: tensor_bar(6)
    type(kernel_spec_t) :: spec
    character(:), allocatable :: primal_code, jvp_code, vjp_code
    integer :: ios, root, unit
    type(expr_t) :: delta

    call arena%init()
    engine = make_native_engine(arena)
    variables = [ &
        sym(arena, "p_parallel"), sym(arena, "p_perpendicular"), &
        sym(arena, "direction_1"), sym(arena, "direction_2"), &
        sym(arena, "direction_3")]
    variables_dot = [ &
        sym(arena, "p_parallel_dot"), sym(arena, "p_perpendicular_dot"), &
        sym(arena, "direction_1_dot"), sym(arena, "direction_2_dot"), &
        sym(arena, "direction_3_dot")]
    tensor_bar = [ &
        sym(arena, "pressure_11_bar"), sym(arena, "pressure_22_bar"), &
        sym(arena, "pressure_33_bar"), sym(arena, "pressure_12_bar"), &
        sym(arena, "pressure_13_bar"), sym(arena, "pressure_23_bar")]

    delta = variables(1) - variables(2)
    tensor = [ &
        variables(2) + delta*variables(3)*variables(3), &
        variables(2) + delta*variables(4)*variables(4), &
        variables(2) + delta*variables(5)*variables(5), &
        delta*variables(3)*variables(4), &
        delta*variables(3)*variables(5), &
        delta*variables(4)*variables(5)]
    tensor_jvp = jvp(tensor, variables, variables_dot)
    tensor_vjp = vjp(tensor, variables, tensor_bar)
    call simplify_all(tensor)
    call simplify_all(tensor_jvp)
    call simplify_all(tensor_vjp)

    call initialize_spec( &
        spec, "generated_cgl_pressure_tensor", &
        "fortfem_generated_cgl_pressure_tensor", 5, 6)
    spec%args = [ &
        str("p_parallel"), str("p_perpendicular"), str("direction_1"), &
        str("direction_2"), str("direction_3")]
    spec%outputs = [ &
        str("pressure_11"), str("pressure_22"), str("pressure_33"), &
        str("pressure_12"), str("pressure_13"), str("pressure_23")]
    primal_code = chars(emit_kernel(tensor, spec))

    call initialize_spec( &
        spec, "generated_cgl_pressure_tensor_jvp", &
        "fortfem_generated_cgl_pressure_tensor_jvp", 10, 6)
    spec%args = [ &
        str("p_parallel"), str("p_perpendicular"), str("direction_1"), &
        str("direction_2"), str("direction_3"), str("p_parallel_dot"), &
        str("p_perpendicular_dot"), str("direction_1_dot"), &
        str("direction_2_dot"), str("direction_3_dot")]
    spec%outputs = [ &
        str("pressure_11_dot"), str("pressure_22_dot"), &
        str("pressure_33_dot"), str("pressure_12_dot"), &
        str("pressure_13_dot"), str("pressure_23_dot")]
    jvp_code = chars(emit_kernel(tensor_jvp, spec))

    call initialize_spec( &
        spec, "generated_cgl_pressure_tensor_vjp", &
        "fortfem_generated_cgl_pressure_tensor_vjp", 11, 5)
    spec%args = [ &
        str("p_parallel"), str("p_perpendicular"), str("direction_1"), &
        str("direction_2"), str("direction_3"), str("pressure_11_bar"), &
        str("pressure_22_bar"), str("pressure_33_bar"), &
        str("pressure_12_bar"), str("pressure_13_bar"), &
        str("pressure_23_bar")]
    spec%outputs = [ &
        str("p_parallel_bar"), str("p_perpendicular_bar"), &
        str("direction_1_bar"), str("direction_2_bar"), &
        str("direction_3_bar")]
    vjp_code = chars(emit_kernel(tensor_vjp, spec))

    open (newunit=unit, &
        file=generated_path("fortfem_cgl_pressure_tensor_products.f90"), &
        status="replace", action="write", iostat=ios)
    if (ios /= 0) error stop "cannot write CGL pressure tensor products"
    write (unit, "(a)") primal_code(:len(primal_code) - 1)
    write (unit, "(a)") jvp_code(:len(jvp_code) - 1)
    write (unit, "(a)") vjp_code(:len(vjp_code) - 1)
    close (unit)

contains

    subroutine initialize_spec( &
            kernel_spec, name, module_name, argument_count, output_count)
        type(kernel_spec_t), intent(out) :: kernel_spec
        character(*), intent(in) :: name, module_name
        integer, intent(in) :: argument_count, output_count

        kernel_spec%name = str(name)
        kernel_spec%module_name = str(module_name)
        kernel_spec%mode = KERNEL_SUBROUTINE
        kernel_spec%generator = str("gen_cgl_pressure_tensor_products")
        kernel_spec%generator_revision = str(fortsym_revision())
        kernel_spec%regenerate_command = str( &
            "cd tools/codegen && ./generate.sh")
        kernel_spec%pure_procedure = .true.
        allocate ( &
            kernel_spec%args(argument_count), &
            kernel_spec%outputs(output_count))
    end subroutine initialize_spec

    subroutine simplify_all(expressions)
        type(expr_t), intent(inout) :: expressions(:)

        do root = 1, size(expressions)
            result = engine%simplify(expressions(root))
            if (result%ok) expressions(root) = result%value
        end do
    end subroutine simplify_all

end program gen_cgl_pressure_tensor_products
