program gen_laplace_bem_kernel_products
    use fortsym_arena, only: arena_t
    use fortsym_engine, only: engine_result_t
    use fortsym_engine_native, only: make_native_engine, native_engine_t
    use fortsym_expr, only: expr_t, operator(+), operator(-), operator(*), &
        operator(/), sqrt, sym
    use fortsym_kernel, only: emit_kernel, kernel_spec_t, KERNEL_SUBROUTINE
    use fortsym_products, only: jvp, vjp
    use fortsym_string, only: chars, str
    use fortfem_codegen_provenance, only: fortsym_revision, generated_path
    implicit none

    type(arena_t), target :: arena
    type(native_engine_t) :: engine
    type(engine_result_t) :: result
    type(expr_t) :: variables(8), variables_dot(8), displacement(3)
    type(expr_t) :: kernel_scale, outputs(1), output_bar(1)
    type(expr_t) :: jvp_roots(1), vjp_roots(8), radius
    type(kernel_spec_t) :: spec
    character(:), allocatable :: primal_code, jvp_code, vjp_code
    integer :: component, ios, root, unit

    call arena%init()
    engine = make_native_engine(arena)
    variables = [ &
        sym(arena, "first_point_1"), sym(arena, "first_point_2"), &
        sym(arena, "first_point_3"), sym(arena, "second_point_1"), &
        sym(arena, "second_point_2"), sym(arena, "second_point_3"), &
        sym(arena, "first_jacobian"), sym(arena, "second_jacobian")]
    variables_dot = [ &
        sym(arena, "first_point_1_dot"), sym(arena, "first_point_2_dot"), &
        sym(arena, "first_point_3_dot"), sym(arena, "second_point_1_dot"), &
        sym(arena, "second_point_2_dot"), sym(arena, "second_point_3_dot"), &
        sym(arena, "first_jacobian_dot"), &
        sym(arena, "second_jacobian_dot")]
    kernel_scale = sym(arena, "kernel_scale")
    output_bar(1) = sym(arena, "value_bar")
    do component = 1, 3
        displacement(component) = variables(component) - &
            variables(component + 3)
    end do
    radius = sqrt( &
        displacement(1)*displacement(1) + &
        displacement(2)*displacement(2) + &
        displacement(3)*displacement(3))
    outputs(1) = kernel_scale*variables(7)*variables(8)/radius

    jvp_roots = jvp(outputs, variables, variables_dot)
    vjp_roots = vjp(outputs, variables, output_bar)
    call simplify_all(jvp_roots)
    call simplify_all(vjp_roots)

    call initialize_spec( &
        spec, "generated_laplace_single_layer_integrand", &
        "fortfem_generated_laplace_single_layer_integrand", 9, 1)
    spec%args = [ &
        str("first_point_1"), str("first_point_2"), str("first_point_3"), &
        str("second_point_1"), str("second_point_2"), &
        str("second_point_3"), str("first_jacobian"), &
        str("second_jacobian"), str("kernel_scale")]
    spec%outputs = [str("value")]
    primal_code = chars(emit_kernel(outputs, spec))

    call initialize_spec( &
        spec, "generated_laplace_single_layer_integrand_jvp", &
        "fortfem_generated_laplace_single_layer_integrand_jvp", 17, 1)
    spec%args = [ &
        str("first_point_1"), str("first_point_2"), str("first_point_3"), &
        str("second_point_1"), str("second_point_2"), &
        str("second_point_3"), str("first_jacobian"), &
        str("second_jacobian"), str("kernel_scale"), &
        str("first_point_1_dot"), str("first_point_2_dot"), &
        str("first_point_3_dot"), str("second_point_1_dot"), &
        str("second_point_2_dot"), str("second_point_3_dot"), &
        str("first_jacobian_dot"), str("second_jacobian_dot")]
    spec%outputs = [str("value_dot")]
    jvp_code = chars(emit_kernel(jvp_roots, spec))

    call initialize_spec( &
        spec, "generated_laplace_single_layer_integrand_vjp", &
        "fortfem_generated_laplace_single_layer_integrand_vjp", 10, 8)
    spec%args = [ &
        str("first_point_1"), str("first_point_2"), str("first_point_3"), &
        str("second_point_1"), str("second_point_2"), &
        str("second_point_3"), str("first_jacobian"), &
        str("second_jacobian"), str("kernel_scale"), str("value_bar")]
    spec%outputs = [ &
        str("first_point_1_bar"), str("first_point_2_bar"), &
        str("first_point_3_bar"), str("second_point_1_bar"), &
        str("second_point_2_bar"), str("second_point_3_bar"), &
        str("first_jacobian_bar"), str("second_jacobian_bar")]
    vjp_code = chars(emit_kernel(vjp_roots, spec))

    open (newunit=unit, &
        file=generated_path("fortfem_laplace_bem_kernel_products.f90"), &
        status="replace", action="write", iostat=ios)
    if (ios /= 0) error stop "cannot write Laplace BEM kernel products"
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
        kernel_spec%generator = str("gen_laplace_bem_kernel_products")
        kernel_spec%generator_revision = str(fortsym_revision())
        kernel_spec%regenerate_command = str("cd tools/codegen && ./generate.sh")
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

end program gen_laplace_bem_kernel_products
