program gen_fci_quintic_lagrange
    use fortsym_arena, only: arena_t
    use fortsym_engine, only: engine_result_t
    use fortsym_engine_native, only: make_native_engine, native_engine_t
    use fortsym_expr, only: expr_t, operator(*), operator(-), operator(/), sym
    use fortsym_kernel, only: emit_kernel, kernel_spec_t, KERNEL_SUBROUTINE
    use fortsym_products, only: jvp, vjp
    use fortsym_string, only: chars, str
    use fortfem_codegen_provenance, only: fortsym_revision, generated_path
    implicit none

    integer, parameter :: node_count = 6
    type(arena_t), target :: arena
    type(engine_result_t) :: result
    type(native_engine_t) :: engine
    type(expr_t) :: variables(node_count + 1), variables_dot(node_count + 1)
    type(expr_t) :: outputs(node_count), outputs_jvp(node_count)
    type(expr_t) :: outputs_vjp(node_count + 1), output_bar(node_count)
    type(expr_t) :: denominator
    type(kernel_spec_t) :: spec
    character(:), allocatable :: primal_code, jvp_code, vjp_code
    integer :: first_other, index, ios, node, other, unit

    call arena%init()
    engine = make_native_engine(arena)
    variables(1) = sym(arena, "target")
    variables(2) = sym(arena, "node_1")
    variables(3) = sym(arena, "node_2")
    variables(4) = sym(arena, "node_3")
    variables(5) = sym(arena, "node_4")
    variables(6) = sym(arena, "node_5")
    variables(7) = sym(arena, "node_6")

    do node = 1, node_count
        first_other = 1
        if (first_other == node) first_other = 2
        outputs(node) = variables(1) - variables(first_other + 1)
        denominator = variables(node + 1) - variables(first_other + 1)
        do other = 1, node_count
            if (other == node .or. other == first_other) cycle
            outputs(node) = outputs(node)* &
                (variables(1) - variables(other + 1))
            denominator = denominator* &
                (variables(node + 1) - variables(other + 1))
        end do
        outputs(node) = outputs(node)/denominator
    end do
    do index = 1, size(outputs)
        result = engine%simplify(outputs(index))
        if (result%ok) outputs(index) = result%value
    end do

    call initialize_spec(spec, "generated_fci_quintic_lagrange_weights", &
        "fortfem_generated_fci_quintic_lagrange", 7, 6)
    spec%args(1) = str("target")
    spec%args(2) = str("node_1")
    spec%args(3) = str("node_2")
    spec%args(4) = str("node_3")
    spec%args(5) = str("node_4")
    spec%args(6) = str("node_5")
    spec%args(7) = str("node_6")
    spec%outputs(1) = str("weight_1")
    spec%outputs(2) = str("weight_2")
    spec%outputs(3) = str("weight_3")
    spec%outputs(4) = str("weight_4")
    spec%outputs(5) = str("weight_5")
    spec%outputs(6) = str("weight_6")
    primal_code = chars(emit_kernel(outputs, spec))

    variables_dot(1) = sym(arena, "target_dot")
    variables_dot(2) = sym(arena, "node_1_dot")
    variables_dot(3) = sym(arena, "node_2_dot")
    variables_dot(4) = sym(arena, "node_3_dot")
    variables_dot(5) = sym(arena, "node_4_dot")
    variables_dot(6) = sym(arena, "node_5_dot")
    variables_dot(7) = sym(arena, "node_6_dot")
    output_bar(1) = sym(arena, "weight_1_bar")
    output_bar(2) = sym(arena, "weight_2_bar")
    output_bar(3) = sym(arena, "weight_3_bar")
    output_bar(4) = sym(arena, "weight_4_bar")
    output_bar(5) = sym(arena, "weight_5_bar")
    output_bar(6) = sym(arena, "weight_6_bar")
    outputs_jvp = jvp(outputs, variables, variables_dot)
    outputs_vjp = vjp(outputs, variables, output_bar)
    do index = 1, size(outputs_jvp)
        result = engine%simplify(outputs_jvp(index))
        if (result%ok) outputs_jvp(index) = result%value
    end do
    do index = 1, size(outputs_vjp)
        result = engine%simplify(outputs_vjp(index))
        if (result%ok) outputs_vjp(index) = result%value
    end do

    call initialize_spec(spec, "generated_fci_quintic_lagrange_weights_jvp", &
        "fortfem_generated_fci_quintic_lagrange_jvp", 14, 6)
    spec%args(1) = str("target")
    spec%args(2) = str("node_1")
    spec%args(3) = str("node_2")
    spec%args(4) = str("node_3")
    spec%args(5) = str("node_4")
    spec%args(6) = str("node_5")
    spec%args(7) = str("node_6")
    spec%args(8) = str("target_dot")
    spec%args(9) = str("node_1_dot")
    spec%args(10) = str("node_2_dot")
    spec%args(11) = str("node_3_dot")
    spec%args(12) = str("node_4_dot")
    spec%args(13) = str("node_5_dot")
    spec%args(14) = str("node_6_dot")
    spec%outputs(1) = str("weight_1_dot")
    spec%outputs(2) = str("weight_2_dot")
    spec%outputs(3) = str("weight_3_dot")
    spec%outputs(4) = str("weight_4_dot")
    spec%outputs(5) = str("weight_5_dot")
    spec%outputs(6) = str("weight_6_dot")
    jvp_code = chars(emit_kernel(outputs_jvp, spec))

    call initialize_spec(spec, "generated_fci_quintic_lagrange_weights_vjp", &
        "fortfem_generated_fci_quintic_lagrange_vjp", 13, 7)
    spec%args(1) = str("target")
    spec%args(2) = str("node_1")
    spec%args(3) = str("node_2")
    spec%args(4) = str("node_3")
    spec%args(5) = str("node_4")
    spec%args(6) = str("node_5")
    spec%args(7) = str("node_6")
    spec%args(8) = str("weight_1_bar")
    spec%args(9) = str("weight_2_bar")
    spec%args(10) = str("weight_3_bar")
    spec%args(11) = str("weight_4_bar")
    spec%args(12) = str("weight_5_bar")
    spec%args(13) = str("weight_6_bar")
    spec%outputs(1) = str("target_bar")
    spec%outputs(2) = str("node_1_bar")
    spec%outputs(3) = str("node_2_bar")
    spec%outputs(4) = str("node_3_bar")
    spec%outputs(5) = str("node_4_bar")
    spec%outputs(6) = str("node_5_bar")
    spec%outputs(7) = str("node_6_bar")
    vjp_code = chars(emit_kernel(outputs_vjp, spec))

    open (newunit=unit, file=generated_path( &
        "fortfem_fci_quintic_lagrange.f90"), status="replace", &
        action="write", iostat=ios)
    if (ios /= 0) error stop "cannot write FCI quintic Lagrange products"
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
        kernel_spec%generator = str("gen_fci_quintic_lagrange")
        kernel_spec%generator_revision = str(fortsym_revision())
        kernel_spec%regenerate_command = str( &
            "cd tools/codegen && ./generate.sh")
        kernel_spec%pure_procedure = .true.
        allocate(kernel_spec%args(argument_count), kernel_spec%outputs(output_count))
    end subroutine initialize_spec

end program gen_fci_quintic_lagrange
