program gen_fci_support_volume_products
    use fortsym_arena, only: arena_t
    use fortsym_engine, only: engine_result_t
    use fortsym_engine_native, only: make_native_engine, native_engine_t
    use fortsym_expr, only: expr_t, operator(+), operator(*), sym
    use fortsym_kernel, only: emit_kernel, kernel_spec_t, KERNEL_SUBROUTINE
    use fortsym_products, only: jvp, vjp
    use fortsym_string, only: chars, str
    use fortfem_codegen_provenance, only: fortsym_revision, generated_path
    implicit none

    type(arena_t), target :: arena
    type(engine_result_t) :: result
    type(native_engine_t) :: engine
    type(expr_t) :: variables(4), variables_dot(4), volume(1)
    type(expr_t) :: volume_jvp(1), volume_vjp(4), volume_bar(1)
    type(kernel_spec_t) :: spec
    character(:), allocatable :: primal_code, jvp_code, vjp_code
    integer :: ios, root, unit

    call arena%init()
    engine = make_native_engine(arena)
    variables = [ &
        sym(arena, "forward_flux_expansion"), &
        sym(arena, "backward_flux_expansion"), &
        sym(arena, "base_cell_area"), sym(arena, "toroidal_field")]
    variables_dot = [ &
        sym(arena, "forward_flux_expansion_dot"), &
        sym(arena, "backward_flux_expansion_dot"), &
        sym(arena, "base_cell_area_dot"), sym(arena, "toroidal_field_dot")]
    volume_bar(1) = sym(arena, "staggered_volume_bar")
    volume(1) = (variables(1) + variables(2))*variables(3)*variables(4)
    volume_jvp = jvp(volume, variables, variables_dot)
    volume_vjp = vjp(volume, variables, volume_bar)
    call simplify_all(volume)
    call simplify_all(volume_jvp)
    call simplify_all(volume_vjp)

    call initialize_spec( &
        spec, "generated_fci_staggered_flux_box_volume", &
        "fortfem_generated_fci_support_volume", 4, 1)
    spec%args = [ &
        str("forward_flux_expansion"), str("backward_flux_expansion"), &
        str("base_cell_area"), str("toroidal_field")]
    spec%outputs = [str("staggered_volume")]
    primal_code = chars(emit_kernel(volume, spec))

    call initialize_spec( &
        spec, "generated_fci_staggered_flux_box_volume_jvp", &
        "fortfem_generated_fci_support_volume_jvp", 8, 1)
    spec%args = [ &
        str("forward_flux_expansion"), str("backward_flux_expansion"), &
        str("base_cell_area"), str("toroidal_field"), &
        str("forward_flux_expansion_dot"), &
        str("backward_flux_expansion_dot"), str("base_cell_area_dot"), &
        str("toroidal_field_dot")]
    spec%outputs = [str("staggered_volume_dot")]
    jvp_code = chars(emit_kernel(volume_jvp, spec))

    call initialize_spec( &
        spec, "generated_fci_staggered_flux_box_volume_vjp", &
        "fortfem_generated_fci_support_volume_vjp", 5, 4)
    spec%args = [ &
        str("forward_flux_expansion"), str("backward_flux_expansion"), &
        str("base_cell_area"), str("toroidal_field"), &
        str("staggered_volume_bar")]
    spec%outputs = [ &
        str("forward_flux_expansion_bar"), &
        str("backward_flux_expansion_bar"), str("base_cell_area_bar"), &
        str("toroidal_field_bar")]
    vjp_code = chars(emit_kernel(volume_vjp, spec))

    open (newunit=unit, &
        file=generated_path("fortfem_fci_support_volume_products.f90"), &
        status="replace", action="write", iostat=ios)
    if (ios /= 0) error stop "cannot write FCI support-volume products"
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
        kernel_spec%generator = str("gen_fci_support_volume_products")
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

end program gen_fci_support_volume_products
