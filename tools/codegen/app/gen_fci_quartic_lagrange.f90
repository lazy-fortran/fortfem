program gen_fci_quartic_lagrange
    use fortsym_arena, only: arena_t
    use fortsym_engine, only: engine_result_t
    use fortsym_engine_native, only: make_native_engine, native_engine_t
    use fortsym_expr, only: expr_t, operator(*), operator(-), operator(/), sym
    use fortsym_kernel, only: emit_kernel, kernel_spec_t, KERNEL_SUBROUTINE
    use fortsym_products, only: jvp, vjp
    use fortsym_string, only: chars, str
    use fortfem_codegen_provenance, only: fortsym_revision, generated_path
    implicit none

    type(arena_t), target :: arena
    type(engine_result_t) :: result
    type(native_engine_t) :: engine
    type(expr_t) :: variables(6), variables_dot(6), outputs(5)
    type(expr_t) :: outputs_jvp(5), outputs_vjp(6), output_bar(5)
    type(kernel_spec_t) :: spec
    character(:), allocatable :: primal_code, jvp_code, vjp_code
    integer :: index, ios, unit

    call arena%init()
    engine = make_native_engine(arena)
    variables = [ &
        sym(arena, "target"), sym(arena, "node_1"), sym(arena, "node_2"), &
        sym(arena, "node_3"), sym(arena, "node_4"), sym(arena, "node_5")]
    outputs(1) = (variables(1) - variables(3))* &
        (variables(1) - variables(4))* &
        (variables(1) - variables(5))* &
        (variables(1) - variables(6))/ &
        ((variables(2) - variables(3))* &
        (variables(2) - variables(4))* &
        (variables(2) - variables(5))* &
        (variables(2) - variables(6)))
    outputs(2) = (variables(1) - variables(2))* &
        (variables(1) - variables(4))* &
        (variables(1) - variables(5))* &
        (variables(1) - variables(6))/ &
        ((variables(3) - variables(2))* &
        (variables(3) - variables(4))* &
        (variables(3) - variables(5))* &
        (variables(3) - variables(6)))
    outputs(3) = (variables(1) - variables(2))* &
        (variables(1) - variables(3))* &
        (variables(1) - variables(5))* &
        (variables(1) - variables(6))/ &
        ((variables(4) - variables(2))* &
        (variables(4) - variables(3))* &
        (variables(4) - variables(5))* &
        (variables(4) - variables(6)))
    outputs(4) = (variables(1) - variables(2))* &
        (variables(1) - variables(3))* &
        (variables(1) - variables(4))* &
        (variables(1) - variables(6))/ &
        ((variables(5) - variables(2))* &
        (variables(5) - variables(3))* &
        (variables(5) - variables(4))* &
        (variables(5) - variables(6)))
    outputs(5) = (variables(1) - variables(2))* &
        (variables(1) - variables(3))* &
        (variables(1) - variables(4))* &
        (variables(1) - variables(5))/ &
        ((variables(6) - variables(2))* &
        (variables(6) - variables(3))* &
        (variables(6) - variables(4))* &
        (variables(6) - variables(5)))
    do index = 1, size(outputs)
        result = engine%simplify(outputs(index))
        if (result%ok) outputs(index) = result%value
    end do

    call initialize_spec(spec, "generated_fci_quartic_lagrange_weights", &
        "fortfem_generated_fci_quartic_lagrange", 6, 5)
    spec%args = [str("target"), str("node_1"), str("node_2"), str("node_3"), &
        str("node_4"), str("node_5")]
    spec%outputs = [str("weight_1"), str("weight_2"), str("weight_3"), &
        str("weight_4"), str("weight_5")]
    primal_code = chars(emit_kernel(outputs, spec))

    variables_dot = [ &
        sym(arena, "target_dot"), sym(arena, "node_1_dot"), &
        sym(arena, "node_2_dot"), sym(arena, "node_3_dot"), &
        sym(arena, "node_4_dot"), sym(arena, "node_5_dot")]
    output_bar = [sym(arena, "weight_1_bar"), sym(arena, "weight_2_bar"), &
        sym(arena, "weight_3_bar"), sym(arena, "weight_4_bar"), &
        sym(arena, "weight_5_bar")]
    outputs_jvp = jvp(outputs, variables, variables_dot)
    outputs_vjp = vjp(outputs, variables, output_bar)
    do index = 1, size(outputs)
        result = engine%simplify(outputs_jvp(index))
        if (result%ok) outputs_jvp(index) = result%value
    end do
    do index = 1, size(outputs_vjp)
        result = engine%simplify(outputs_vjp(index))
        if (result%ok) outputs_vjp(index) = result%value
    end do

    call initialize_spec(spec, "generated_fci_quartic_lagrange_weights_jvp", &
        "fortfem_generated_fci_quartic_lagrange_jvp", 12, 5)
    spec%args = [str("target"), str("node_1"), str("node_2"), str("node_3"), &
        str("node_4"), str("node_5"), str("target_dot"), str("node_1_dot"), &
        str("node_2_dot"), str("node_3_dot"), str("node_4_dot"), &
        str("node_5_dot")]
    spec%outputs = [str("weight_1_dot"), str("weight_2_dot"), &
        str("weight_3_dot"), str("weight_4_dot"), str("weight_5_dot")]
    jvp_code = chars(emit_kernel(outputs_jvp, spec))

    call initialize_spec(spec, "generated_fci_quartic_lagrange_weights_vjp", &
        "fortfem_generated_fci_quartic_lagrange_vjp", 11, 6)
    spec%args = [str("target"), str("node_1"), str("node_2"), str("node_3"), &
        str("node_4"), str("node_5"), str("weight_1_bar"), str("weight_2_bar"), &
        str("weight_3_bar"), str("weight_4_bar"), str("weight_5_bar")]
    spec%outputs = [str("target_bar"), str("node_1_bar"), str("node_2_bar"), &
        str("node_3_bar"), str("node_4_bar"), str("node_5_bar")]
    vjp_code = chars(emit_kernel(outputs_vjp, spec))

    open (newunit=unit, file=generated_path( &
        "fortfem_fci_quartic_lagrange.f90"), status="replace", &
        action="write", iostat=ios)
    if (ios /= 0) error stop "cannot write FCI quartic Lagrange products"
    write (unit, "(a)") primal_code(:len(primal_code) - 1)
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
        kernel_spec%generator = str("gen_fci_quartic_lagrange")
        kernel_spec%generator_revision = str(fortsym_revision())
        kernel_spec%regenerate_command = str( &
            "cd tools/codegen && ./generate.sh")
        kernel_spec%pure_procedure = .true.
        allocate(kernel_spec%args(argument_count), kernel_spec%outputs(output_count))
    end subroutine initialize_spec

end program gen_fci_quartic_lagrange
