program gen_surface_triangle_geometry_products
    use fortsym_arena, only: arena_t
    use fortsym_engine, only: engine_result_t
    use fortsym_engine_native, only: make_native_engine, native_engine_t
    use fortsym_expr, only: expr_t, operator(+), operator(-), operator(*), &
        operator(/), sqrt, sym
    use fortsym_kernel, only: emit_kernel, kernel_spec_t, KERNEL_SUBROUTINE
    use fortsym_products, only: jvp, vjp
    use fortsym_string, only: chars, str
    use fortfem_codegen_provenance, only: fortsym_revision, generated_path
    implicit none

    type(arena_t), target :: arena
    type(native_engine_t) :: engine
    type(engine_result_t) :: result
    type(expr_t) :: variables(11), variables_dot(11), outputs(7)
    type(expr_t) :: output_bar(7), jvp_roots(7), vjp_roots(11)
    type(expr_t) :: first(3), second(3), area_vector(3), jacobian
    type(kernel_spec_t) :: spec
    character(:), allocatable :: primal_code, jvp_code, vjp_code
    integer :: component, ios, root, unit

    call arena%init()
    engine = make_native_engine(arena)
    variables = [ &
        sym(arena, "vertex_11"), sym(arena, "vertex_21"), &
        sym(arena, "vertex_31"), sym(arena, "vertex_12"), &
        sym(arena, "vertex_22"), sym(arena, "vertex_32"), &
        sym(arena, "vertex_13"), sym(arena, "vertex_23"), &
        sym(arena, "vertex_33"), sym(arena, "xi"), sym(arena, "eta")]
    variables_dot = [ &
        sym(arena, "vertex_11_dot"), sym(arena, "vertex_21_dot"), &
        sym(arena, "vertex_31_dot"), sym(arena, "vertex_12_dot"), &
        sym(arena, "vertex_22_dot"), sym(arena, "vertex_32_dot"), &
        sym(arena, "vertex_13_dot"), sym(arena, "vertex_23_dot"), &
        sym(arena, "vertex_33_dot"), sym(arena, "xi_dot"), &
        sym(arena, "eta_dot")]
    output_bar = [ &
        sym(arena, "point_1_bar"), sym(arena, "point_2_bar"), &
        sym(arena, "point_3_bar"), sym(arena, "jacobian_bar"), &
        sym(arena, "normal_1_bar"), sym(arena, "normal_2_bar"), &
        sym(arena, "normal_3_bar")]

    do component = 1, 3
        first(component) = variables(component + 3) - variables(component)
        second(component) = variables(component + 6) - variables(component)
    end do
    area_vector = [ &
        first(2)*second(3) - first(3)*second(2), &
        first(3)*second(1) - first(1)*second(3), &
        first(1)*second(2) - first(2)*second(1)]
    jacobian = sqrt( &
        area_vector(1)*area_vector(1) + &
        area_vector(2)*area_vector(2) + &
        area_vector(3)*area_vector(3))
    do component = 1, 3
        outputs(component) = variables(component) + &
            variables(10)*first(component) + &
            variables(11)*second(component)
    end do
    outputs(4) = jacobian
    do component = 1, 3
        outputs(component + 4) = area_vector(component)/jacobian
    end do

    jvp_roots = jvp(outputs, variables, variables_dot)
    vjp_roots = vjp(outputs, variables, output_bar)
    call simplify_all(jvp_roots)
    call simplify_all(vjp_roots)

    call initialize_spec( &
        spec, "generated_surface_triangle_geometry_3d", &
        "fortfem_generated_surface_triangle_geometry_3d", 11, 7)
    spec%args = [ &
        str("vertex_11"), str("vertex_21"), str("vertex_31"), &
        str("vertex_12"), str("vertex_22"), str("vertex_32"), &
        str("vertex_13"), str("vertex_23"), str("vertex_33"), &
        str("xi"), str("eta")]
    spec%outputs = [ &
        str("point_1"), str("point_2"), str("point_3"), &
        str("jacobian"), str("normal_1"), str("normal_2"), str("normal_3")]
    primal_code = chars(emit_kernel(outputs, spec))

    call initialize_spec( &
        spec, "generated_surface_triangle_geometry_3d_jvp", &
        "fortfem_generated_surface_triangle_geometry_3d_jvp", 22, 7)
    spec%args = [ &
        str("vertex_11"), str("vertex_21"), str("vertex_31"), &
        str("vertex_12"), str("vertex_22"), str("vertex_32"), &
        str("vertex_13"), str("vertex_23"), str("vertex_33"), &
        str("xi"), str("eta"), &
        str("vertex_11_dot"), str("vertex_21_dot"), &
        str("vertex_31_dot"), str("vertex_12_dot"), &
        str("vertex_22_dot"), str("vertex_32_dot"), &
        str("vertex_13_dot"), str("vertex_23_dot"), &
        str("vertex_33_dot"), str("xi_dot"), str("eta_dot")]
    spec%outputs = [ &
        str("point_1_dot"), str("point_2_dot"), str("point_3_dot"), &
        str("jacobian_dot"), str("normal_1_dot"), str("normal_2_dot"), &
        str("normal_3_dot")]
    jvp_code = chars(emit_kernel(jvp_roots, spec))

    call initialize_spec( &
        spec, "generated_surface_triangle_geometry_3d_vjp", &
        "fortfem_generated_surface_triangle_geometry_3d_vjp", 18, 11)
    spec%args = [ &
        str("vertex_11"), str("vertex_21"), str("vertex_31"), &
        str("vertex_12"), str("vertex_22"), str("vertex_32"), &
        str("vertex_13"), str("vertex_23"), str("vertex_33"), &
        str("xi"), str("eta"), str("point_1_bar"), str("point_2_bar"), &
        str("point_3_bar"), str("jacobian_bar"), str("normal_1_bar"), &
        str("normal_2_bar"), str("normal_3_bar")]
    spec%outputs = [ &
        str("vertex_11_bar"), str("vertex_21_bar"), &
        str("vertex_31_bar"), str("vertex_12_bar"), &
        str("vertex_22_bar"), str("vertex_32_bar"), &
        str("vertex_13_bar"), str("vertex_23_bar"), &
        str("vertex_33_bar"), str("xi_bar"), str("eta_bar")]
    vjp_code = chars(emit_kernel(vjp_roots, spec))

    open (newunit=unit, &
        file=generated_path("fortfem_surface_triangle_geometry_products.f90"), &
        status="replace", action="write", iostat=ios)
    if (ios /= 0) error stop "cannot write surface triangle geometry products"
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
        kernel_spec%generator = str("gen_surface_triangle_geometry_products")
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

end program gen_surface_triangle_geometry_products
