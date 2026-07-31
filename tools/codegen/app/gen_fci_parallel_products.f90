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
    type(expr_t) :: diffusion_variables(8), diffusion_variables_dot(8)
    type(expr_t) :: diffusion_outputs(2), diffusion_jvp(2)
    type(expr_t) :: diffusion_vjp(8), diffusion_output_bar(2)
    type(expr_t) :: diagonal_variables(5), diagonal_output(1)
    type(kernel_spec_t) :: spec
    character(:), allocatable :: primal_code, jvp_code, vjp_code
    character(:), allocatable :: diffusion_code, diffusion_jvp_code
    character(:), allocatable :: diffusion_vjp_code
    character(:), allocatable :: diagonal_code
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

    diffusion_variables = [ &
        sym(arena, "forward_value"), sym(arena, "upper_field"), &
        sym(arena, "backward_value"), sym(arena, "lower_field"), &
        sym(arena, "line_length"), sym(arena, "parallel_coefficient"), &
        sym(arena, "canonical_volume"), sym(arena, "staggered_volume")]
    diffusion_variables_dot = [ &
        sym(arena, "forward_value_dot"), sym(arena, "upper_field_dot"), &
        sym(arena, "backward_value_dot"), sym(arena, "lower_field_dot"), &
        sym(arena, "line_length_dot"), &
        sym(arena, "parallel_coefficient_dot"), &
        sym(arena, "canonical_volume_dot"), &
        sym(arena, "staggered_volume_dot")]
    diffusion_output_bar = [ &
        sym(arena, "lower_contribution_bar"), &
        sym(arena, "upper_contribution_bar")]
    diffusion_outputs(1) = diffusion_variables(3)*diffusion_variables(8)* &
        diffusion_variables(6)*gradient(1)/diffusion_variables(5)/ &
        diffusion_variables(7)
    diffusion_outputs(2) = -diffusion_variables(1)*diffusion_variables(8)* &
        diffusion_variables(6)*gradient(1)/diffusion_variables(5)/ &
        diffusion_variables(7)
    diffusion_jvp = jvp(diffusion_outputs, diffusion_variables, &
        diffusion_variables_dot)
    diffusion_vjp = vjp(diffusion_outputs, diffusion_variables, &
        diffusion_output_bar)
    call simplify_all(diffusion_outputs)
    call simplify_all(diffusion_jvp)
    call simplify_all(diffusion_vjp)

    call initialize_spec( &
        spec, "generated_fci_parallel_diffusion", &
        "fortfem_generated_fci_parallel_diffusion", 8, 2)
    spec%args = [ &
        str("forward_value"), str("upper_field"), str("backward_value"), &
        str("lower_field"), str("line_length"), &
        str("parallel_coefficient"), str("canonical_volume"), &
        str("staggered_volume")]
    spec%outputs = [str("lower_contribution"), str("upper_contribution")]
    diffusion_code = chars(emit_kernel(diffusion_outputs, spec))

    call initialize_spec( &
        spec, "generated_fci_parallel_diffusion_jvp", &
        "fortfem_generated_fci_parallel_diffusion_jvp", 16, 2)
    spec%args = [ &
        str("forward_value"), str("upper_field"), str("backward_value"), &
        str("lower_field"), str("line_length"), &
        str("parallel_coefficient"), str("canonical_volume"), &
        str("staggered_volume"), str("forward_value_dot"), &
        str("upper_field_dot"), str("backward_value_dot"), &
        str("lower_field_dot"), str("line_length_dot"), &
        str("parallel_coefficient_dot"), str("canonical_volume_dot"), &
        str("staggered_volume_dot")]
    spec%outputs = [str("lower_contribution_dot"), &
        str("upper_contribution_dot")]
    diffusion_jvp_code = chars(emit_kernel(diffusion_jvp, spec))

    call initialize_spec( &
        spec, "generated_fci_parallel_diffusion_vjp", &
        "fortfem_generated_fci_parallel_diffusion_vjp", 10, 8)
    spec%args = [ &
        str("forward_value"), str("upper_field"), str("backward_value"), &
        str("lower_field"), str("line_length"), &
        str("parallel_coefficient"), str("canonical_volume"), &
        str("staggered_volume"), str("lower_contribution_bar"), &
        str("upper_contribution_bar")]
    spec%outputs = [ &
        str("forward_value_bar"), str("upper_field_bar"), &
        str("backward_value_bar"), str("lower_field_bar"), &
        str("line_length_bar"), str("parallel_coefficient_bar"), &
        str("canonical_volume_bar"), str("staggered_volume_bar")]
    diffusion_vjp_code = chars(emit_kernel(diffusion_vjp, spec))

    diagonal_variables = [ &
        sym(arena, "map_value"), sym(arena, "line_length"), &
        sym(arena, "parallel_coefficient"), sym(arena, "canonical_volume"), &
        sym(arena, "staggered_volume")]
    diagonal_output(1) = diagonal_variables(1)*diagonal_variables(1)* &
        diagonal_variables(3)*diagonal_variables(5)/ &
        (diagonal_variables(2)*diagonal_variables(2)*diagonal_variables(4))
    call simplify_all(diagonal_output)

    call initialize_spec( &
        spec, "generated_fci_parallel_diffusion_diagonal", &
        "fortfem_generated_fci_parallel_diagonal", 5, 1)
    spec%args = [ &
        str("map_value"), str("line_length"), str("parallel_coefficient"), &
        str("canonical_volume"), str("staggered_volume")]
    spec%outputs = [str("diagonal_contribution")]
    diagonal_code = chars(emit_kernel(diagonal_output, spec))

    open (newunit=unit, &
        file=generated_path("fortfem_fci_parallel_products.f90"), &
        status="replace", action="write", iostat=ios)
    if (ios /= 0) error stop "cannot write FCI parallel products"
    write (unit, "(a)") primal_code(:len(primal_code) - 1)
    write (unit, "(a)") jvp_code(:len(jvp_code) - 1)
    write (unit, "(a)") vjp_code(:len(vjp_code) - 1)
    write (unit, "(a)") diffusion_code(:len(diffusion_code) - 1)
    write (unit, "(a)") diffusion_jvp_code(:len(diffusion_jvp_code) - 1)
    write (unit, "(a)") diffusion_vjp_code(:len(diffusion_vjp_code) - 1)
    write (unit, "(a)") diagonal_code(:len(diagonal_code) - 1)
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
