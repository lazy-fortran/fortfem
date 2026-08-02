program gen_block_graph_products
    !! Emit the scalar product used by the packed block-graph action.
    !!
    !! The graph topology and its loops stay runtime-owned.  The local product
    !! rule is deliberately generated so the value, tangent, and reverse paths
    !! cannot drift as more field blocks are composed.
    use fortsym_arena, only: arena_t
    use fortsym_engine, only: engine_result_t
    use fortsym_engine_native, only: make_native_engine, native_engine_t
    use fortsym_expr, only: expr_t, operator(*), sym
    use fortsym_kernel, only: emit_kernel, kernel_spec_t, KERNEL_SUBROUTINE
    use fortsym_products, only: jvp, vjp
    use fortsym_string, only: chars, str
    use fortfem_codegen_provenance, only: fortsym_revision, generated_path
    implicit none

    type(arena_t), target :: arena
    type(native_engine_t) :: engine
    type(engine_result_t) :: simplified
    type(expr_t) :: block_value, state_value, block_value_dot, state_value_dot
    type(expr_t) :: product(1), product_dot(1), product_vjp(2), product_bar
    type(kernel_spec_t) :: spec
    character(:), allocatable :: primal_code, jvp_code, vjp_code
    integer :: root, unit, ios

    call arena%init()
    engine = make_native_engine(arena)
    block_value = sym(arena, "block_value")
    state_value = sym(arena, "state_value")
    block_value_dot = sym(arena, "block_value_dot")
    state_value_dot = sym(arena, "state_value_dot")
    product_bar = sym(arena, "product_bar")
    product = [block_value*state_value]
    product_dot = jvp(product, [block_value, state_value], &
        [block_value_dot, state_value_dot])
    product_vjp = vjp(product, [block_value, state_value], [product_bar])
    call simplify_all(product)
    call simplify_all(product_dot)
    call simplify_all(product_vjp)

    call initialize_spec(spec, "generated_block_graph_product", &
        "fortfem_generated_block_graph", 2, 1)
    spec%args = [str("block_value"), str("state_value")]
    spec%outputs = [str("contribution")]
    primal_code = chars(emit_kernel(product, spec))

    call initialize_spec(spec, "generated_block_graph_product_jvp", &
        "fortfem_generated_block_graph_product_jvp", 4, 1)
    spec%args = [str("block_value"), str("state_value"), &
        str("block_value_dot"), str("state_value_dot")]
    spec%outputs = [str("contribution_dot")]
    jvp_code = chars(emit_kernel(product_dot, spec))

    call initialize_spec(spec, "generated_block_graph_product_vjp", &
        "fortfem_generated_block_graph_product_vjp", 3, 2)
    spec%args = [str("block_value"), str("state_value"), str("product_bar")]
    spec%outputs = [str("block_value_bar"), str("state_value_bar")]
    vjp_code = chars(emit_kernel(product_vjp, spec))

    open (newunit=unit, file=generated_path( &
        "fortfem_block_graph_products.f90"), status="replace", action="write", iostat=ios)
    if (ios /= 0) error stop "cannot write block graph products"
    write (unit, "(a)") primal_code(:len(primal_code) - 1)
    write (unit, "(a)") jvp_code(:len(jvp_code) - 1)
    write (unit, "(a)") vjp_code(:len(vjp_code) - 1)
    close (unit)

contains

    subroutine initialize_spec(kernel_spec, name, module_name, argument_count, output_count)
        type(kernel_spec_t), intent(out) :: kernel_spec
        character(*), intent(in) :: name, module_name
        integer, intent(in) :: argument_count, output_count

        kernel_spec%name = str(name)
        kernel_spec%module_name = str(module_name)
        kernel_spec%mode = KERNEL_SUBROUTINE
        kernel_spec%generator = str("gen_block_graph_products")
        kernel_spec%generator_revision = str(fortsym_revision())
        kernel_spec%regenerate_command = str("cd tools/codegen && ./generate.sh")
        kernel_spec%pure_procedure = .true.
        allocate (kernel_spec%args(argument_count), kernel_spec%outputs(output_count))
    end subroutine initialize_spec

    subroutine simplify_all(expressions)
        type(expr_t), intent(inout) :: expressions(:)

        do root = 1, size(expressions)
            simplified = engine%simplify(expressions(root))
            if (simplified%ok) expressions(root) = simplified%value
        end do
    end subroutine simplify_all

end program gen_block_graph_products
