program gen_surface_shape_objective_products
    !! Emit pointwise weighted squared shape-mismatch products.
    !!
    !! The geometry wrapper owns the fixed-topology reductions and validation;
    !! FortSym owns the scalar mismatch, tangent, and reverse products.
    use fortsym_arena, only: arena_t
    use fortsym_engine, only: engine_result_t
    use fortsym_engine_native, only: make_native_engine, native_engine_t
    use fortsym_expr, only: expr_t, operator(-), operator(*), real_expr, sym
    use fortsym_kernel, only: emit_kernel, kernel_spec_t, KERNEL_SUBROUTINE
    use fortsym_products, only: jvp, vjp
    use fortsym_string, only: chars, str
    use fortfem_codegen_provenance, only: fortsym_revision, generated_path
    implicit none

    type(arena_t), target :: arena
    type(native_engine_t) :: engine
    type(engine_result_t) :: result
    type(expr_t) :: variables(3), variables_dot(3), outputs(1)
    type(expr_t) :: outputs_jvp(1), output_bar(1), outputs_vjp(3)
    type(expr_t) :: half
    type(kernel_spec_t) :: spec
    character(:), allocatable :: value_code, jvp_code, vjp_code
    integer :: unit, ios

    call arena%init()
    engine = make_native_engine(arena)
    half = real_expr(arena, 0.5d0)
    variables = [sym(arena, "candidate"), sym(arena, "target"), &
        sym(arena, "weight")]
    variables_dot = [sym(arena, "candidate_dot"), sym(arena, "target_dot"), &
        sym(arena, "weight_dot")]
    outputs(1) = half*variables(3)*(variables(1) - variables(2))* &
        (variables(1) - variables(2))
    outputs_jvp = jvp(outputs, variables, variables_dot)
    output_bar(1) = sym(arena, "contribution_bar")
    outputs_vjp = vjp(outputs, variables, output_bar)
    call simplify(outputs(1))
    call simplify(outputs_jvp(1))
    call simplify(outputs_vjp(1))
    call simplify(outputs_vjp(2))
    call simplify(outputs_vjp(3))

    call initialize_spec(spec, "generated_surface_shape_objective_contribution", &
        "fortfem_generated_surface_shape_objective_contribution", 3, 1)
    spec%args = [str("candidate"), str("target"), str("weight")]
    spec%outputs = [str("contribution")]
    value_code = chars(emit_kernel(outputs, spec))

    call initialize_spec(spec, "generated_surface_shape_objective_contribution_jvp", &
        "fortfem_generated_surface_shape_objective_contribution_jvp", 6, 1)
    spec%args = [str("candidate"), str("target"), str("weight"), &
        str("candidate_dot"), str("target_dot"), str("weight_dot")]
    spec%outputs = [str("contribution_dot")]
    jvp_code = chars(emit_kernel(outputs_jvp, spec))

    call initialize_spec(spec, "generated_surface_shape_objective_contribution_vjp", &
        "fortfem_generated_surface_shape_objective_contribution_vjp", 4, 3)
    spec%args = [str("candidate"), str("target"), str("weight"), &
        str("contribution_bar")]
    spec%outputs = [str("candidate_bar"), str("target_bar"), str("weight_bar")]
    vjp_code = chars(emit_kernel(outputs_vjp, spec))

    open (newunit=unit, file=generated_path( &
        "fortfem_surface_shape_objective_products.f90"), status="replace", &
        action="write", iostat=ios)
    if (ios /= 0) error stop "cannot write surface-shape objective products"
    write (unit, "(a)") value_code(:len(value_code) - 1)
    write (unit, "(a)") jvp_code(:len(jvp_code) - 1)
    write (unit, "(a)") vjp_code(:len(vjp_code) - 1)
    close (unit)

contains

    subroutine initialize_spec(kernel_spec, name, module_name, argument_count, &
            output_count)
        type(kernel_spec_t), intent(out) :: kernel_spec
        character(*), intent(in) :: name, module_name
        integer, intent(in) :: argument_count, output_count

        kernel_spec%name = str(name)
        kernel_spec%module_name = str(module_name)
        kernel_spec%mode = KERNEL_SUBROUTINE
        kernel_spec%generator = str("gen_surface_shape_objective_products")
        kernel_spec%generator_revision = str(fortsym_revision())
        kernel_spec%regenerate_command = str( &
            "cd tools/codegen && ./generate.sh")
        kernel_spec%pure_procedure = .true.
        allocate(kernel_spec%args(argument_count), kernel_spec%outputs(output_count))
    end subroutine initialize_spec

    subroutine simplify(expression)
        type(expr_t), intent(inout) :: expression

        result = engine%simplify(expression)
        if (result%ok) expression = result%value
    end subroutine simplify

end program gen_surface_shape_objective_products
