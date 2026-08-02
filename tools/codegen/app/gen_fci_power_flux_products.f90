program gen_fci_power_flux_products
    !! Emit pointwise products used by the conservative FCI power ledger.
    !!
    !! The runtime wrapper owns array traversal and fixed topology.  FortSym
    !! owns the scalar product, product-rule tangent, and reverse products so
    !! the parallel flux/power signs cannot diverge between value and adjoint
    !! paths.
    use fortsym_arena, only: arena_t
    use fortsym_engine, only: engine_result_t
    use fortsym_engine_native, only: make_native_engine, native_engine_t
    use fortsym_expr, only: expr_t, operator(*), operator(-), sym
    use fortsym_kernel, only: emit_kernel, kernel_spec_t, KERNEL_SUBROUTINE
    use fortsym_products, only: jvp, vjp
    use fortsym_string, only: chars, str
    use fortfem_codegen_provenance, only: fortsym_revision, generated_path
    implicit none

    type(arena_t), target :: arena
    type(native_engine_t) :: engine
    type(engine_result_t) :: result
    type(expr_t) :: parallel_variables(3), parallel_variables_dot(3)
    type(expr_t) :: parallel_outputs(2), parallel_outputs_jvp(2)
    type(expr_t) :: parallel_output_bar(2), parallel_outputs_vjp(3)
    type(expr_t) :: perpendicular_variables(3), perpendicular_variables_dot(3)
    type(expr_t) :: perpendicular_output(1), perpendicular_output_jvp(1)
    type(expr_t) :: perpendicular_output_bar(1), perpendicular_output_vjp(3)
    type(kernel_spec_t) :: spec
    character(:), allocatable :: parallel_code, parallel_jvp_code, parallel_vjp_code
    character(:), allocatable :: perpendicular_code, perpendicular_jvp_code
    character(:), allocatable :: perpendicular_vjp_code
    integer :: ios, root, unit

    call arena%init()
    engine = make_native_engine(arena)

    parallel_variables = [ &
        sym(arena, "gradient"), sym(arena, "coefficient"), &
        sym(arena, "staggered_volume")]
    parallel_variables_dot = [ &
        sym(arena, "gradient_dot"), sym(arena, "coefficient_dot"), &
        sym(arena, "staggered_volume_dot")]
    parallel_outputs = [ &
        -parallel_variables(2)*parallel_variables(1), &
        -parallel_variables(3)*parallel_variables(2)*parallel_variables(1)* &
        parallel_variables(1)]
    parallel_outputs_jvp = jvp( &
        parallel_outputs, parallel_variables, parallel_variables_dot)
    parallel_output_bar = [ &
        sym(arena, "parallel_flux_bar"), sym(arena, "parallel_power_bar")]
    parallel_outputs_vjp = vjp( &
        parallel_outputs, parallel_variables, parallel_output_bar)
    call simplify_all(parallel_outputs)
    call simplify_all(parallel_outputs_jvp)
    call simplify_all(parallel_outputs_vjp)

    call initialize_spec(spec, "generated_fci_parallel_flux_power", &
        "fortfem_generated_fci_parallel_flux_power", 3, 2)
    spec%args = [str("gradient"), str("coefficient"), str("staggered_volume")]
    spec%outputs = [str("parallel_flux"), str("parallel_power")]
    parallel_code = chars(emit_kernel(parallel_outputs, spec))

    call initialize_spec(spec, "generated_fci_parallel_flux_power_jvp", &
        "fortfem_generated_fci_parallel_flux_power_jvp", 6, 2)
    spec%args = [ &
        str("gradient"), str("coefficient"), str("staggered_volume"), &
        str("gradient_dot"), str("coefficient_dot"), str("staggered_volume_dot")]
    spec%outputs = [str("parallel_flux_dot"), str("parallel_power_dot")]
    parallel_jvp_code = chars(emit_kernel(parallel_outputs_jvp, spec))

    call initialize_spec(spec, "generated_fci_parallel_flux_power_vjp", &
        "fortfem_generated_fci_parallel_flux_power_vjp", 5, 3)
    spec%args = [ &
        str("gradient"), str("coefficient"), str("staggered_volume"), &
        str("parallel_flux_bar"), str("parallel_power_bar")]
    spec%outputs = [ &
        str("gradient_bar"), str("coefficient_bar"), str("staggered_volume_bar")]
    parallel_vjp_code = chars(emit_kernel(parallel_outputs_vjp, spec))

    perpendicular_variables = [ &
        sym(arena, "field"), sym(arena, "action"), sym(arena, "canonical_volume")]
    perpendicular_variables_dot = [ &
        sym(arena, "field_dot"), sym(arena, "action_dot"), &
        sym(arena, "canonical_volume_dot")]
    perpendicular_output(1) = perpendicular_variables(1)*perpendicular_variables(2)* &
        perpendicular_variables(3)
    perpendicular_output_jvp = jvp( &
        perpendicular_output, perpendicular_variables, perpendicular_variables_dot)
    perpendicular_output_bar(1) = sym(arena, "perpendicular_power_bar")
    perpendicular_output_vjp = vjp( &
        perpendicular_output, perpendicular_variables, perpendicular_output_bar)
    call simplify_all(perpendicular_output)
    call simplify_all(perpendicular_output_jvp)
    call simplify_all(perpendicular_output_vjp)

    call initialize_spec(spec, "generated_fci_perpendicular_power", &
        "fortfem_generated_fci_perpendicular_power", 3, 1)
    spec%args = [str("field"), str("action"), str("canonical_volume")]
    spec%outputs = [str("perpendicular_power")]
    perpendicular_code = chars(emit_kernel(perpendicular_output, spec))

    call initialize_spec(spec, "generated_fci_perpendicular_power_jvp", &
        "fortfem_generated_fci_perpendicular_power_jvp", 6, 1)
    spec%args = [ &
        str("field"), str("action"), str("canonical_volume"), str("field_dot"), &
        str("action_dot"), str("canonical_volume_dot")]
    spec%outputs = [str("perpendicular_power_dot")]
    perpendicular_jvp_code = chars(emit_kernel(perpendicular_output_jvp, spec))

    call initialize_spec(spec, "generated_fci_perpendicular_power_vjp", &
        "fortfem_generated_fci_perpendicular_power_vjp", 4, 3)
    spec%args = [ &
        str("field"), str("action"), str("canonical_volume"), &
        str("perpendicular_power_bar")]
    spec%outputs = [ &
        str("field_bar"), str("action_bar"), str("canonical_volume_bar")]
    perpendicular_vjp_code = chars(emit_kernel(perpendicular_output_vjp, spec))

    open (newunit=unit, file=generated_path( &
        "fortfem_fci_power_flux_products.f90"), status="replace", &
        action="write", iostat=ios)
    if (ios /= 0) error stop "cannot write FCI power/flux products"
    write (unit, "(a)") parallel_code(:len(parallel_code) - 1)
    write (unit, "(a)") parallel_jvp_code(:len(parallel_jvp_code) - 1)
    write (unit, "(a)") parallel_vjp_code(:len(parallel_vjp_code) - 1)
    write (unit, "(a)") perpendicular_code(:len(perpendicular_code) - 1)
    write (unit, "(a)") perpendicular_jvp_code(:len(perpendicular_jvp_code) - 1)
    write (unit, "(a)") perpendicular_vjp_code(:len(perpendicular_vjp_code) - 1)
    close (unit)

contains

    subroutine initialize_spec(kernel_spec, name, module_name, argument_count, output_count)
        type(kernel_spec_t), intent(out) :: kernel_spec
        character(*), intent(in) :: name, module_name
        integer, intent(in) :: argument_count, output_count

        kernel_spec%name = str(name)
        kernel_spec%module_name = str(module_name)
        kernel_spec%mode = KERNEL_SUBROUTINE
        kernel_spec%generator = str("gen_fci_power_flux_products")
        kernel_spec%generator_revision = str(fortsym_revision())
        kernel_spec%regenerate_command = str("cd tools/codegen && ./generate.sh")
        kernel_spec%pure_procedure = .true.
        allocate (kernel_spec%args(argument_count), kernel_spec%outputs(output_count))
    end subroutine initialize_spec

    subroutine simplify_all(expressions)
        type(expr_t), intent(inout) :: expressions(:)

        do root = 1, size(expressions)
            result = engine%simplify(expressions(root))
            if (result%ok) expressions(root) = result%value
        end do
    end subroutine simplify_all

end program gen_fci_power_flux_products
