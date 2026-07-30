program gen_nurbs_geometry_products
    use fortsym_arena, only: arena_t
    use fortsym_engine, only: engine_result_t
    use fortsym_engine_native, only: make_native_engine, native_engine_t
    use fortsym_expr, only: expr_t, operator(+), operator(-), operator(*), &
        operator(/), sym
    use fortsym_kernel, only: emit_kernel, kernel_spec_t, KERNEL_SUBROUTINE
    use fortsym_products, only: jvp, vjp
    use fortsym_string, only: chars, str
    use fortfem_codegen_provenance, only: fortsym_revision, generated_path
    implicit none

    type(arena_t), target :: arena
    type(native_engine_t) :: engine
    type(engine_result_t) :: result
    type(expr_t) :: outputs(3), variables(6), directions(6)
    type(expr_t) :: cotangents(3), jvp_roots(3), vjp_roots(6)
    type(kernel_spec_t) :: jvp_spec, vjp_spec
    character(:), allocatable :: jvp_code, vjp_code
    integer :: ios, root, unit

    call arena%init()
    engine = make_native_engine(arena)
    variables = [ &
        sym(arena, "numerator"), &
        sym(arena, "derivative_numerator_x"), &
        sym(arena, "derivative_numerator_y"), &
        sym(arena, "denominator"), &
        sym(arena, "derivative_denominator_x"), &
        sym(arena, "derivative_denominator_y")]
    directions = [ &
        sym(arena, "numerator_dot"), &
        sym(arena, "derivative_numerator_x_dot"), &
        sym(arena, "derivative_numerator_y_dot"), &
        sym(arena, "denominator_dot"), &
        sym(arena, "derivative_denominator_x_dot"), &
        sym(arena, "derivative_denominator_y_dot")]
    cotangents = [ &
        sym(arena, "point_bar"), &
        sym(arena, "jacobian_x_bar"), &
        sym(arena, "jacobian_y_bar")]

    outputs(1) = variables(1)/variables(4)
    outputs(2) = (variables(2) - &
        outputs(1)*variables(5))/variables(4)
    outputs(3) = (variables(3) - &
        outputs(1)*variables(6))/variables(4)
    jvp_roots = jvp(outputs, variables, directions)
    vjp_roots = vjp(outputs, variables, cotangents)
    call simplify_all(jvp_roots)
    call simplify_all(vjp_roots)

    jvp_spec%name = str("generated_nurbs_geometry_quotient_jvp")
    jvp_spec%module_name = str("fortfem_generated_nurbs_geometry_jvp")
    jvp_spec%mode = KERNEL_SUBROUTINE
    jvp_spec%generator = str("gen_nurbs_geometry_products")
    jvp_spec%generator_revision = str(fortsym_revision())
    jvp_spec%regenerate_command = str("cd tools/codegen && ./generate.sh")
    jvp_spec%pure_procedure = .true.
    allocate (jvp_spec%args(12), jvp_spec%outputs(3))
    jvp_spec%args = [ &
        str("numerator"), str("derivative_numerator_x"), &
        str("derivative_numerator_y"), str("denominator"), &
        str("derivative_denominator_x"), &
        str("derivative_denominator_y"), str("numerator_dot"), &
        str("derivative_numerator_x_dot"), &
        str("derivative_numerator_y_dot"), str("denominator_dot"), &
        str("derivative_denominator_x_dot"), &
        str("derivative_denominator_y_dot")]
    jvp_spec%outputs = [ &
        str("point_dot"), str("jacobian_x_dot"), str("jacobian_y_dot")]

    vjp_spec%name = str("generated_nurbs_geometry_quotient_vjp")
    vjp_spec%module_name = str("fortfem_generated_nurbs_geometry_vjp")
    vjp_spec%mode = KERNEL_SUBROUTINE
    vjp_spec%generator = str("gen_nurbs_geometry_products")
    vjp_spec%generator_revision = str(fortsym_revision())
    vjp_spec%regenerate_command = str("cd tools/codegen && ./generate.sh")
    vjp_spec%pure_procedure = .true.
    allocate (vjp_spec%args(9), vjp_spec%outputs(6))
    vjp_spec%args = [ &
        str("numerator"), str("derivative_numerator_x"), &
        str("derivative_numerator_y"), str("denominator"), &
        str("derivative_denominator_x"), &
        str("derivative_denominator_y"), str("point_bar"), &
        str("jacobian_x_bar"), str("jacobian_y_bar")]
    vjp_spec%outputs = [ &
        str("numerator_bar"), str("derivative_numerator_x_bar"), &
        str("derivative_numerator_y_bar"), str("denominator_bar"), &
        str("derivative_denominator_x_bar"), &
        str("derivative_denominator_y_bar")]

    jvp_code = chars(emit_kernel(jvp_roots, jvp_spec))
    vjp_code = chars(emit_kernel(vjp_roots, vjp_spec))
    open (newunit=unit, &
        file=generated_path("fortfem_nurbs_geometry_products.f90"), &
        status="replace", action="write", iostat=ios)
    if (ios /= 0) error stop "cannot write NURBS geometry products"
    write (unit, "(a)") jvp_code(:len(jvp_code) - 1)
    write (unit, "(a)") vjp_code(:len(vjp_code) - 1)
    close (unit)

contains

    subroutine simplify_all(expressions)
        type(expr_t), intent(inout) :: expressions(:)

        do root = 1, size(expressions)
            result = engine%simplify(expressions(root))
            if (result%ok) expressions(root) = result%value
        end do
    end subroutine simplify_all

end program gen_nurbs_geometry_products
