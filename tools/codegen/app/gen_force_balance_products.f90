program gen_force_balance_products
    !! Emit the scalar weak force-product and its tangent/reverse products.
    !!
    !! Assembly topology remains runtime-owned.  This local contraction is
    !! generated once so volume, boundary, and sheet terms share one product
    !! rule and cannot acquire divergent derivative formulas.
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
    type(expr_t) :: variables(3), variables_dot(3), product(1), product_dot(1)
    type(expr_t) :: product_vjp(3), product_bar
    type(kernel_spec_t) :: spec
    character(:), allocatable :: primal_code, jvp_code, vjp_code
    integer :: root, unit, ios

    call arena%init()
    engine = make_native_engine(arena)
    variables = [ &
        sym(arena, "weight"), sym(arena, "test_value"), &
        sym(arena, "force_value")]
    variables_dot = [ &
        sym(arena, "weight_dot"), sym(arena, "test_value_dot"), &
        sym(arena, "force_value_dot")]
    product_bar = sym(arena, "contribution_bar")
    product = [variables(1)*variables(2)*variables(3)]
    product_dot = jvp(product, variables, variables_dot)
    product_vjp = vjp(product, variables, [product_bar])
    call simplify_all(product)
    call simplify_all(product_dot)
    call simplify_all(product_vjp)

    call initialize_spec(spec, "generated_force_balance_product", &
        "fortfem_generated_force_balance_product", 3, 1)
    spec%args = [str("weight"), str("test_value"), str("force_value")]
    spec%outputs = [str("contribution")]
    primal_code = chars(emit_kernel(product, spec))

    call initialize_spec(spec, "generated_force_balance_product_jvp", &
        "fortfem_generated_force_balance_product_jvp", 6, 1)
    spec%args = [ &
        str("weight"), str("test_value"), str("force_value"), &
        str("weight_dot"), str("test_value_dot"), str("force_value_dot")]
    spec%outputs = [str("contribution_dot")]
    jvp_code = chars(emit_kernel(product_dot, spec))

    call initialize_spec(spec, "generated_force_balance_product_vjp", &
        "fortfem_generated_force_balance_product_vjp", 4, 3)
    spec%args = [ &
        str("weight"), str("test_value"), str("force_value"), &
        str("contribution_bar")]
    spec%outputs = [ &
        str("weight_bar"), str("test_value_bar"), str("force_value_bar")]
    vjp_code = chars(emit_kernel(product_vjp, spec))

    open (newunit=unit, file=generated_path( &
        "fortfem_force_balance_products.f90"), status="replace", &
        action="write", iostat=ios)
    if (ios /= 0) error stop "cannot write force-balance products"
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
        kernel_spec%generator = str("gen_force_balance_products")
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

end program gen_force_balance_products
