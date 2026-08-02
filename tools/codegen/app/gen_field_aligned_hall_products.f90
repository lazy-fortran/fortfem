program gen_field_aligned_hall_products
    !! Emit the three independent products in k_H [b]_x and their derivatives.
    !!
    !! The symmetric projector remains owned by the existing generated CGL
    !! kernel.  Generating only q_i = k_H b_i avoids duplicating it while the
    !! runtime wrapper owns the fixed skew-matrix layout.
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
    type(expr_t) :: variables(4), variables_dot(4)
    type(expr_t) :: products(3), products_dot(3), products_bar(3)
    type(expr_t) :: products_vjp(4)
    type(kernel_spec_t) :: spec
    character(:), allocatable :: primal_code, jvp_code, vjp_code
    integer :: root, unit, ios

    call arena%init()
    engine = make_native_engine(arena)
    variables = [ &
        sym(arena, "hall_coefficient"), sym(arena, "direction_1"), &
        sym(arena, "direction_2"), sym(arena, "direction_3")]
    variables_dot = [ &
        sym(arena, "hall_coefficient_dot"), sym(arena, "direction_1_dot"), &
        sym(arena, "direction_2_dot"), sym(arena, "direction_3_dot")]
    products_bar = [ &
        sym(arena, "hall_direction_1_bar"), &
        sym(arena, "hall_direction_2_bar"), &
        sym(arena, "hall_direction_3_bar")]
    products = [ &
        variables(1)*variables(2), variables(1)*variables(3), &
        variables(1)*variables(4)]
    products_dot = jvp(products, variables, variables_dot)
    products_vjp = vjp(products, variables, products_bar)
    call simplify_all(products)
    call simplify_all(products_dot)
    call simplify_all(products_vjp)

    call initialize_spec(spec, "generated_field_aligned_hall", &
        "fortfem_generated_field_aligned_hall", 4, 3)
    spec%args = [ &
        str("hall_coefficient"), str("direction_1"), str("direction_2"), &
        str("direction_3")]
    spec%outputs = [ &
        str("hall_direction_1"), str("hall_direction_2"), &
        str("hall_direction_3")]
    primal_code = chars(emit_kernel(products, spec))

    call initialize_spec(spec, "generated_field_aligned_hall_jvp", &
        "fortfem_generated_field_aligned_hall_jvp", 8, 3)
    spec%args = [ &
        str("hall_coefficient"), str("direction_1"), str("direction_2"), &
        str("direction_3"), str("hall_coefficient_dot"), &
        str("direction_1_dot"), str("direction_2_dot"), &
        str("direction_3_dot")]
    spec%outputs = [ &
        str("hall_direction_1_dot"), str("hall_direction_2_dot"), &
        str("hall_direction_3_dot")]
    jvp_code = chars(emit_kernel(products_dot, spec))

    call initialize_spec(spec, "generated_field_aligned_hall_vjp", &
        "fortfem_generated_field_aligned_hall_vjp", 7, 4)
    spec%args = [ &
        str("hall_coefficient"), str("direction_1"), str("direction_2"), &
        str("direction_3"), str("hall_direction_1_bar"), &
        str("hall_direction_2_bar"), str("hall_direction_3_bar")]
    spec%outputs = [ &
        str("hall_coefficient_bar"), str("direction_1_bar"), &
        str("direction_2_bar"), str("direction_3_bar")]
    vjp_code = chars(emit_kernel(products_vjp, spec))

    open (newunit=unit, file=generated_path( &
        "fortfem_field_aligned_hall_products.f90"), status="replace", &
        action="write", iostat=ios)
    if (ios /= 0) error stop "cannot write field-aligned Hall products"
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
        kernel_spec%generator = str("gen_field_aligned_hall_products")
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

end program gen_field_aligned_hall_products
