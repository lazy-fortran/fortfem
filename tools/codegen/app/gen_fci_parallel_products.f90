program gen_fci_parallel_products
    use fortsym_arena, only: arena_t
    use fortsym_engine, only: engine_result_t
    use fortsym_engine_native, only: make_native_engine, native_engine_t
    use fortsym_expr, only: expr_t, operator(*), operator(+), operator(-), &
        operator(/), sym
    use fortsym_kernel, only: emit_kernel, kernel_spec_t, KERNEL_SUBROUTINE
    use fortsym_products, only: jvp, vjp
    use fortsym_string, only: chars, str
    use fortfem_codegen_provenance, only: fortsym_revision, generated_path
    implicit none

    type(arena_t), target :: arena
    type(engine_result_t) :: result
    type(native_engine_t) :: engine
    type(expr_t) :: variables(5), variables_dot(5), gradient(1)
    type(expr_t) :: gradient_jvp(1), gradient_vjp(5), gradient_bar(1)
    type(kernel_spec_t) :: spec
    character(:), allocatable :: primal_code, jvp_code, vjp_code
    integer :: ios, root, unit

    call arena%init()
    engine = make_native_engine(arena)
    variables = [ &
        sym(arena, "forward_value"), sym(arena, "upper_field"), &
        sym(arena, "backward_value"), sym(arena, "lower_field"), &
        sym(arena, "line_length")]
    variables_dot = [ &
        sym(arena, "forward_value_dot"), sym(arena, "upper_field_dot"), &
        sym(arena, "backward_value_dot"), sym(arena, "lower_field_dot"), &
        sym(arena, "line_length_dot")]
    gradient_bar(1) = sym(arena, "gradient_bar")
    gradient(1) = (variables(1)*variables(2) - &
        variables(3)*variables(4))/variables(5)
    gradient_jvp = jvp(gradient, variables, variables_dot)
    gradient_vjp = vjp(gradient, variables, gradient_bar)
    call simplify_all(gradient)
    call simplify_all(gradient_jvp)
    call simplify_all(gradient_vjp)

    call initialize_spec( &
        spec, "generated_fci_parallel_gradient", &
        "fortfem_generated_fci_parallel_gradient", 5, 1)
    spec%args = [ &
        str("forward_value"), str("upper_field"), str("backward_value"), &
        str("lower_field"), str("line_length")]
    spec%outputs = [str("gradient")]
    primal_code = chars(emit_kernel(gradient, spec))

    call initialize_spec( &
        spec, "generated_fci_parallel_gradient_jvp", &
        "fortfem_generated_fci_parallel_gradient_jvp", 10, 1)
    spec%args = [ &
        str("forward_value"), str("upper_field"), str("backward_value"), &
        str("lower_field"), str("line_length"), str("forward_value_dot"), &
        str("upper_field_dot"), str("backward_value_dot"), &
        str("lower_field_dot"), str("line_length_dot")]
    spec%outputs = [str("gradient_dot")]
    jvp_code = chars(emit_kernel(gradient_jvp, spec))

    call initialize_spec( &
        spec, "generated_fci_parallel_gradient_vjp", &
        "fortfem_generated_fci_parallel_gradient_vjp", 6, 5)
    spec%args = [ &
        str("forward_value"), str("upper_field"), str("backward_value"), &
        str("lower_field"), str("line_length"), str("gradient_bar")]
    spec%outputs = [ &
        str("forward_value_bar"), str("upper_field_bar"), &
        str("backward_value_bar"), str("lower_field_bar"), &
        str("line_length_bar")]
    vjp_code = chars(emit_kernel(gradient_vjp, spec))

    open (newunit=unit, &
        file=generated_path("fortfem_fci_parallel_products.f90"), &
        status="replace", action="write", iostat=ios)
    if (ios /= 0) error stop "cannot write FCI parallel products"
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
        kernel_spec%generator = str("gen_fci_parallel_products")
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

end program gen_fci_parallel_products
