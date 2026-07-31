program gen_field_aligned_flux_products
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
    type(expr_t) :: variables(8), variables_dot(8), flux(3)
    type(expr_t) :: flux_jvp(3), flux_vjp(8), flux_bar(3)
    type(kernel_spec_t) :: spec
    character(:), allocatable :: primal_code, jvp_code, vjp_code
    integer :: ios, root, unit
    type(expr_t) :: delta, projection

    call arena%init()
    engine = make_native_engine(arena)
    variables = [ &
        sym(arena, "parallel_coefficient"), &
        sym(arena, "perpendicular_coefficient"), &
        sym(arena, "direction_1"), sym(arena, "direction_2"), &
        sym(arena, "direction_3"), sym(arena, "gradient_1"), &
        sym(arena, "gradient_2"), sym(arena, "gradient_3")]
    variables_dot = [ &
        sym(arena, "parallel_coefficient_dot"), &
        sym(arena, "perpendicular_coefficient_dot"), &
        sym(arena, "direction_1_dot"), sym(arena, "direction_2_dot"), &
        sym(arena, "direction_3_dot"), sym(arena, "gradient_1_dot"), &
        sym(arena, "gradient_2_dot"), sym(arena, "gradient_3_dot")]
    flux_bar = [ &
        sym(arena, "flux_1_bar"), sym(arena, "flux_2_bar"), &
        sym(arena, "flux_3_bar")]
    delta = variables(1) - variables(2)
    projection = variables(3)*variables(6) + variables(4)*variables(7) + &
        variables(5)*variables(8)
    flux = [ &
        variables(2)*variables(6) + delta*variables(3)*projection, &
        variables(2)*variables(7) + delta*variables(4)*projection, &
        variables(2)*variables(8) + delta*variables(5)*projection]
    flux_jvp = jvp(flux, variables, variables_dot)
    flux_vjp = vjp(flux, variables, flux_bar)
    call simplify_all(flux)
    call simplify_all(flux_jvp)
    call simplify_all(flux_vjp)

    call initialize_spec( &
        spec, "generated_field_aligned_flux", &
        "fortfem_generated_field_aligned_flux", 8, 3)
    spec%args = [ &
        str("parallel_coefficient"), str("perpendicular_coefficient"), &
        str("direction_1"), str("direction_2"), str("direction_3"), &
        str("gradient_1"), str("gradient_2"), str("gradient_3")]
    spec%outputs = [str("flux_1"), str("flux_2"), str("flux_3")]
    primal_code = chars(emit_kernel(flux, spec))

    call initialize_spec( &
        spec, "generated_field_aligned_flux_jvp", &
        "fortfem_generated_field_aligned_flux_jvp", 16, 3)
    spec%args = [ &
        str("parallel_coefficient"), str("perpendicular_coefficient"), &
        str("direction_1"), str("direction_2"), str("direction_3"), &
        str("gradient_1"), str("gradient_2"), str("gradient_3"), &
        str("parallel_coefficient_dot"), &
        str("perpendicular_coefficient_dot"), str("direction_1_dot"), &
        str("direction_2_dot"), str("direction_3_dot"), str("gradient_1_dot"), &
        str("gradient_2_dot"), str("gradient_3_dot")]
    spec%outputs = [str("flux_1_dot"), str("flux_2_dot"), str("flux_3_dot")]
    jvp_code = chars(emit_kernel(flux_jvp, spec))

    call initialize_spec( &
        spec, "generated_field_aligned_flux_vjp", &
        "fortfem_generated_field_aligned_flux_vjp", 11, 8)
    spec%args = [ &
        str("parallel_coefficient"), str("perpendicular_coefficient"), &
        str("direction_1"), str("direction_2"), str("direction_3"), &
        str("gradient_1"), str("gradient_2"), str("gradient_3"), &
        str("flux_1_bar"), str("flux_2_bar"), str("flux_3_bar")]
    spec%outputs = [ &
        str("parallel_coefficient_bar"), &
        str("perpendicular_coefficient_bar"), str("direction_1_bar"), &
        str("direction_2_bar"), str("direction_3_bar"), str("gradient_1_bar"), &
        str("gradient_2_bar"), str("gradient_3_bar")]
    vjp_code = chars(emit_kernel(flux_vjp, spec))

    open (newunit=unit, &
        file=generated_path("fortfem_field_aligned_flux_products.f90"), &
        status="replace", action="write", iostat=ios)
    if (ios /= 0) error stop "cannot write field-aligned flux products"
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
        kernel_spec%generator = str("gen_field_aligned_flux_products")
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

end program gen_field_aligned_flux_products
