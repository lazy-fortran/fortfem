program gen_regularized_surface_current_products
    !! Emit a normalized Gaussian sheet-current sample and tangent products.
    use fortsym_arena, only: arena_t
    use fortsym_engine, only: engine_result_t
    use fortsym_engine_native, only: make_native_engine, native_engine_t
    use fortsym_expr, only: exp, expr_t, operator(*), operator(-), &
        operator(/), operator(**), sym
    use fortsym_kernel, only: emit_kernel, kernel_spec_t, KERNEL_SUBROUTINE
    use fortsym_products, only: jvp, vjp
    use fortsym_string, only: chars, str
    use fortfem_codegen_provenance, only: fortsym_revision, generated_path
    implicit none

    type(arena_t), target :: arena
    type(native_engine_t) :: engine
    type(engine_result_t) :: simplified
    type(expr_t) :: variables(3), variables_dot(3), value(1), value_dot(1)
    type(expr_t) :: value_vjp(3), value_bar(1), inverse_sqrt_pi, profile
    type(kernel_spec_t) :: spec
    character(:), allocatable :: primal_code, jvp_code, vjp_code
    integer :: root, unit, ios

    call arena%init()
    engine = make_native_engine(arena)
    variables = [ &
        sym(arena, "signed_distance"), sym(arena, "sheet_current"), &
        sym(arena, "thickness")]
    variables_dot = [ &
        sym(arena, "signed_distance_dot"), sym(arena, "sheet_current_dot"), &
        sym(arena, "thickness_dot")]
    inverse_sqrt_pi = sym(arena, "inverse_sqrt_pi")
    value_bar(1) = sym(arena, "volume_current_bar")
    profile = inverse_sqrt_pi/variables(3)* &
        exp(-(variables(1)/variables(3))**2)
    value(1) = profile*variables(2)
    value_dot = jvp(value, variables, variables_dot)
    value_vjp = vjp(value, variables, value_bar)
    call simplify_all(value)
    call simplify_all(value_dot)
    call simplify_all(value_vjp)

    call initialize_spec( &
        spec, "generated_regularized_surface_current", &
        "fortfem_generated_regularized_surface_current", 4, 1)
    spec%args = [ &
        str("signed_distance"), str("sheet_current"), str("thickness"), &
        str("inverse_sqrt_pi")]
    spec%outputs = [str("volume_current")]
    primal_code = chars(emit_kernel(value, spec))

    call initialize_spec( &
        spec, "generated_regularized_surface_current_jvp", &
        "fortfem_generated_regularized_surface_current_jvp", 7, 1)
    spec%args = [ &
        str("signed_distance"), str("sheet_current"), str("thickness"), &
        str("inverse_sqrt_pi"), str("signed_distance_dot"), &
        str("sheet_current_dot"), str("thickness_dot")]
    spec%outputs = [str("volume_current_dot")]
    jvp_code = chars(emit_kernel(value_dot, spec))

    call initialize_spec( &
        spec, "generated_regularized_surface_current_vjp", &
        "fortfem_generated_regularized_surface_current_vjp", 5, 3)
    spec%args = [ &
        str("signed_distance"), str("sheet_current"), str("thickness"), &
        str("inverse_sqrt_pi"), str("volume_current_bar")]
    spec%outputs = [ &
        str("signed_distance_bar"), str("sheet_current_bar"), &
        str("thickness_bar")]
    vjp_code = chars(emit_kernel(value_vjp, spec))

    open (newunit=unit, file=generated_path( &
        "fortfem_regularized_surface_current_products.f90"), &
        status="replace", action="write", iostat=ios)
    if (ios /= 0) error stop "cannot write regularized surface-current products"
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
        kernel_spec%generator = str("gen_regularized_surface_current_products")
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
            simplified = engine%simplify(expressions(root))
            if (simplified%ok) expressions(root) = simplified%value
        end do
    end subroutine simplify_all

end program gen_regularized_surface_current_products
