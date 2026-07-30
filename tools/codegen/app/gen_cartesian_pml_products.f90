program gen_cartesian_pml_products
    use fortsym_arena, only: arena_t
    use fortsym_engine, only: engine_result_t
    use fortsym_engine_native, only: make_native_engine, native_engine_t
    use fortsym_expr, only: expr_t, operator(*), operator(/), sym
    use fortsym_kernel, only: emit_kernel, kernel_spec_t, KERNEL_SUBROUTINE
    use fortsym_products, only: jvp, vjp
    use fortsym_string, only: chars, str
    use fortfem_codegen_provenance, only: fortsym_revision, generated_path
    implicit none

    type(arena_t), target :: arena
    type(native_engine_t) :: engine
    type(engine_result_t) :: result
    type(expr_t) :: stretch(3), stretch_dot(3)
    type(expr_t) :: scalar_outputs(4), scalar_outputs_bar(4)
    type(expr_t) :: curl_outputs(6), curl_outputs_bar(6)
    type(expr_t) :: scalar_jvp(4), scalar_vjp(3)
    type(expr_t) :: curl_jvp(6), curl_vjp(3)
    type(kernel_spec_t) :: spec
    character(:), allocatable :: scalar_jvp_code, scalar_vjp_code
    character(:), allocatable :: curl_jvp_code, curl_vjp_code
    integer :: ios, root, unit

    call arena%init()
    engine = make_native_engine(arena)
    stretch = [sym(arena, "stretch_1"), sym(arena, "stretch_2"), &
        sym(arena, "stretch_3")]
    stretch_dot = [sym(arena, "stretch_1_dot"), &
        sym(arena, "stretch_2_dot"), sym(arena, "stretch_3_dot")]
    scalar_outputs = [ &
        stretch(2)*stretch(3)/stretch(1), &
        stretch(1)*stretch(3)/stretch(2), &
        stretch(1)*stretch(2)/stretch(3), &
        stretch(1)*stretch(2)*stretch(3)]
    scalar_outputs_bar = [ &
        sym(arena, "gradient_1_bar"), sym(arena, "gradient_2_bar"), &
        sym(arena, "gradient_3_bar"), sym(arena, "mass_bar")]
    curl_outputs = [ &
        stretch(1)/(stretch(2)*stretch(3)), &
        stretch(2)/(stretch(1)*stretch(3)), &
        stretch(3)/(stretch(1)*stretch(2)), scalar_outputs(1:3)]
    curl_outputs_bar = [ &
        sym(arena, "curl_1_bar"), sym(arena, "curl_2_bar"), &
        sym(arena, "curl_3_bar"), sym(arena, "mass_1_bar"), &
        sym(arena, "mass_2_bar"), sym(arena, "mass_3_bar")]

    scalar_jvp = jvp(scalar_outputs, stretch, stretch_dot)
    scalar_vjp = vjp(scalar_outputs, stretch, scalar_outputs_bar)
    curl_jvp = jvp(curl_outputs, stretch, stretch_dot)
    curl_vjp = vjp(curl_outputs, stretch, curl_outputs_bar)
    call simplify_all(scalar_jvp)
    call simplify_all(scalar_vjp)
    call simplify_all(curl_jvp)
    call simplify_all(curl_vjp)

    call initialize_spec( &
        spec, "generated_cartesian_scalar_pml_jvp", &
        "fortfem_generated_cartesian_scalar_pml_jvp", 6, 4)
    spec%args = [ &
        str("stretch_1"), str("stretch_2"), str("stretch_3"), &
        str("stretch_1_dot"), str("stretch_2_dot"), str("stretch_3_dot")]
    spec%outputs = [ &
        str("gradient_1_dot"), str("gradient_2_dot"), &
        str("gradient_3_dot"), str("mass_dot")]
    scalar_jvp_code = emit_text(scalar_jvp, spec)

    call initialize_spec( &
        spec, "generated_cartesian_scalar_pml_vjp", &
        "fortfem_generated_cartesian_scalar_pml_vjp", 7, 3)
    spec%args = [ &
        str("stretch_1"), str("stretch_2"), str("stretch_3"), &
        str("gradient_1_bar"), str("gradient_2_bar"), &
        str("gradient_3_bar"), str("mass_bar")]
    spec%outputs = [ &
        str("stretch_1_bar"), str("stretch_2_bar"), str("stretch_3_bar")]
    scalar_vjp_code = emit_text(scalar_vjp, spec)

    call initialize_spec( &
        spec, "generated_cartesian_curl_curl_pml_jvp", &
        "fortfem_generated_cartesian_curl_curl_pml_jvp", 6, 6)
    spec%args = [ &
        str("stretch_1"), str("stretch_2"), str("stretch_3"), &
        str("stretch_1_dot"), str("stretch_2_dot"), str("stretch_3_dot")]
    spec%outputs = [ &
        str("curl_1_dot"), str("curl_2_dot"), str("curl_3_dot"), &
        str("mass_1_dot"), str("mass_2_dot"), str("mass_3_dot")]
    curl_jvp_code = emit_text(curl_jvp, spec)

    call initialize_spec( &
        spec, "generated_cartesian_curl_curl_pml_vjp", &
        "fortfem_generated_cartesian_curl_curl_pml_vjp", 9, 3)
    spec%args = [ &
        str("stretch_1"), str("stretch_2"), str("stretch_3"), &
        str("curl_1_bar"), str("curl_2_bar"), str("curl_3_bar"), &
        str("mass_1_bar"), str("mass_2_bar"), str("mass_3_bar")]
    spec%outputs = [ &
        str("stretch_1_bar"), str("stretch_2_bar"), str("stretch_3_bar")]
    curl_vjp_code = emit_text(curl_vjp, spec)

    open (newunit=unit, &
        file=generated_path("fortfem_cartesian_pml_products.f90"), &
        status="replace", action="write", iostat=ios)
    if (ios /= 0) error stop "cannot write Cartesian PML products"
    write (unit, "(a)") scalar_jvp_code(:len(scalar_jvp_code) - 1)
    write (unit, "(a)") scalar_vjp_code(:len(scalar_vjp_code) - 1)
    write (unit, "(a)") curl_jvp_code(:len(curl_jvp_code) - 1)
    write (unit, "(a)") curl_vjp_code(:len(curl_vjp_code) - 1)
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
        kernel_spec%generator = str("gen_cartesian_pml_products")
        kernel_spec%generator_revision = str(fortsym_revision())
        kernel_spec%regenerate_command = str("cd tools/codegen && ./generate.sh")
        kernel_spec%scalar_type = str("complex(dp)")
        kernel_spec%pure_procedure = .true.
        allocate ( &
            kernel_spec%args(argument_count), &
            kernel_spec%outputs(output_count))
    end subroutine initialize_spec

    function emit_text(roots, kernel_spec) result(text)
        type(expr_t), intent(in) :: roots(:)
        type(kernel_spec_t), intent(in) :: kernel_spec
        character(:), allocatable :: text

        text = chars(emit_kernel(roots, kernel_spec))
    end function emit_text

    subroutine simplify_all(expressions)
        type(expr_t), intent(inout) :: expressions(:)

        do root = 1, size(expressions)
            result = engine%simplify(expressions(root))
            if (result%ok) expressions(root) = result%value
        end do
    end subroutine simplify_all

end program gen_cartesian_pml_products
