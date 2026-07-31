program gen_laplace_singular_edge_products
    use fortsym_arena, only: arena_t
    use fortsym_engine, only: engine_result_t
    use fortsym_engine_native, only: make_native_engine, native_engine_t
    use fortsym_expr, only: expr_t, log, operator(+), operator(-), &
        operator(*), operator(/), sqrt, sym
    use fortsym_kernel, only: emit_kernel, kernel_spec_t, KERNEL_SUBROUTINE
    use fortsym_products, only: jvp, vjp
    use fortsym_string, only: chars, str
    use fortfem_codegen_provenance, only: fortsym_revision, generated_path
    implicit none

    type(arena_t), target :: arena
    type(native_engine_t) :: engine
    type(engine_result_t) :: result
    type(expr_t) :: variables(9), variables_dot(9), first(3), second(3)
    type(expr_t) :: edge_vector(3), area_vector(3), output(1), output_bar(1)
    type(expr_t) :: first_norm, second_norm, length, height
    type(expr_t) :: jvp_roots(1), vjp_roots(9)
    type(kernel_spec_t) :: spec
    character(:), allocatable :: primal_code, jvp_code, vjp_code
    integer :: component, ios, root, unit

    call arena%init()
    engine = make_native_engine(arena)
    variables = [ &
        sym(arena, "point_1"), sym(arena, "point_2"), &
        sym(arena, "point_3"), sym(arena, "first_vertex_1"), &
        sym(arena, "first_vertex_2"), sym(arena, "first_vertex_3"), &
        sym(arena, "second_vertex_1"), sym(arena, "second_vertex_2"), &
        sym(arena, "second_vertex_3")]
    variables_dot = [ &
        sym(arena, "point_1_dot"), sym(arena, "point_2_dot"), &
        sym(arena, "point_3_dot"), sym(arena, "first_vertex_1_dot"), &
        sym(arena, "first_vertex_2_dot"), sym(arena, "first_vertex_3_dot"), &
        sym(arena, "second_vertex_1_dot"), &
        sym(arena, "second_vertex_2_dot"), &
        sym(arena, "second_vertex_3_dot")]
    output_bar(1) = sym(arena, "value_bar")
    do component = 1, 3
        first(component) = variables(component + 3) - variables(component)
        second(component) = variables(component + 6) - variables(component)
        edge_vector(component) = second(component) - first(component)
    end do
    area_vector = [ &
        first(2)*second(3) - first(3)*second(2), &
        first(3)*second(1) - first(1)*second(3), &
        first(1)*second(2) - first(2)*second(1)]
    first_norm = vector_norm(first)
    second_norm = vector_norm(second)
    length = vector_norm(edge_vector)
    height = vector_norm(area_vector)/length
    output(1) = height*log( &
        (first_norm + second_norm + length)/ &
        (first_norm + second_norm - length))

    jvp_roots = jvp(output, variables, variables_dot)
    vjp_roots = vjp(output, variables, output_bar)
    call simplify_all(jvp_roots)
    call simplify_all(vjp_roots)

    call initialize_spec( &
        spec, "generated_laplace_singular_edge_potential", &
        "fortfem_generated_laplace_singular_edge_potential", 9, 1)
    spec%args = primal_arguments()
    spec%outputs = [str("value")]
    primal_code = chars(emit_kernel(output, spec))

    call initialize_spec( &
        spec, "generated_laplace_singular_edge_potential_jvp", &
        "fortfem_generated_laplace_singular_edge_potential_jvp", 18, 1)
    spec%args(1:9) = primal_arguments()
    spec%args(10:18) = [ &
        str("point_1_dot"), str("point_2_dot"), str("point_3_dot"), &
        str("first_vertex_1_dot"), str("first_vertex_2_dot"), &
        str("first_vertex_3_dot"), str("second_vertex_1_dot"), &
        str("second_vertex_2_dot"), str("second_vertex_3_dot")]
    spec%outputs = [str("value_dot")]
    jvp_code = chars(emit_kernel(jvp_roots, spec))

    call initialize_spec( &
        spec, "generated_laplace_singular_edge_potential_vjp", &
        "fortfem_generated_laplace_singular_edge_potential_vjp", 10, 9)
    spec%args(1:9) = primal_arguments()
    spec%args(10) = str("value_bar")
    spec%outputs = [ &
        str("point_1_bar"), str("point_2_bar"), str("point_3_bar"), &
        str("first_vertex_1_bar"), str("first_vertex_2_bar"), &
        str("first_vertex_3_bar"), str("second_vertex_1_bar"), &
        str("second_vertex_2_bar"), str("second_vertex_3_bar")]
    vjp_code = chars(emit_kernel(vjp_roots, spec))

    open (newunit=unit, &
        file=generated_path("fortfem_laplace_singular_edge_products.f90"), &
        status="replace", action="write", iostat=ios)
    if (ios /= 0) error stop "cannot write Laplace singular edge products"
    write (unit, "(a)") primal_code(:len(primal_code) - 1)
    write (unit, "(a)") jvp_code(:len(jvp_code) - 1)
    write (unit, "(a)") vjp_code(:len(vjp_code) - 1)
    close (unit)

contains

    function vector_norm(vector) result(norm)
        type(expr_t), intent(in) :: vector(3)
        type(expr_t) :: norm

        norm = sqrt( &
            vector(1)*vector(1) + vector(2)*vector(2) + &
            vector(3)*vector(3))
    end function vector_norm

    function primal_arguments() result(arguments)
        use fortsym_string, only: str_t
        type(str_t), allocatable :: arguments(:)

        allocate(arguments(9))
        arguments = [ &
            str("point_1"), str("point_2"), str("point_3"), &
            str("first_vertex_1"), str("first_vertex_2"), &
            str("first_vertex_3"), str("second_vertex_1"), &
            str("second_vertex_2"), str("second_vertex_3")]
    end function primal_arguments

    subroutine initialize_spec( &
            kernel_spec, name, module_name, argument_count, output_count)
        type(kernel_spec_t), intent(out) :: kernel_spec
        character(*), intent(in) :: name, module_name
        integer, intent(in) :: argument_count, output_count

        kernel_spec%name = str(name)
        kernel_spec%module_name = str(module_name)
        kernel_spec%mode = KERNEL_SUBROUTINE
        kernel_spec%generator = str("gen_laplace_singular_edge_products")
        kernel_spec%generator_revision = str(fortsym_revision())
        kernel_spec%regenerate_command = str("cd tools/codegen && ./generate.sh")
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

end program gen_laplace_singular_edge_products
