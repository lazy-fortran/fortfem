program gen_sphere_curved_panel_products
    use fortsym_arena, only: arena_t
    use fortsym_engine, only: engine_result_t
    use fortsym_engine_native, only: make_native_engine, native_engine_t
    use fortsym_expr, only: expr_t, operator(+), operator(-), operator(*), &
        operator(/), sqrt, sym
    use fortsym_kernel, only: emit_kernel, kernel_spec_t, KERNEL_SUBROUTINE
    use fortsym_products, only: jvp, vjp
    use fortsym_string, only: chars, str, str_t
    use fortfem_codegen_provenance, only: fortsym_revision, generated_path
    implicit none

    type(arena_t), target :: arena
    type(engine_result_t) :: result
    type(native_engine_t) :: engine
    type(expr_t) :: variables(12), variables_dot(12), outputs(10)
    type(expr_t) :: output_bar(10), jvp_roots(10), vjp_roots(12)
    type(expr_t) :: affine_point(3), derivative_xi(3), derivative_eta(3)
    type(expr_t) :: tangent_xi(3), tangent_eta(3), cross_product(3)
    type(expr_t) :: norm_affine, dot_xi, dot_eta, norm_cubed
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
        sym(arena, "vertex_33"), sym(arena, "radius"), &
        sym(arena, "xi"), sym(arena, "eta")]
    variables_dot = [ &
        sym(arena, "vertex_11_dot"), sym(arena, "vertex_21_dot"), &
        sym(arena, "vertex_31_dot"), sym(arena, "vertex_12_dot"), &
        sym(arena, "vertex_22_dot"), sym(arena, "vertex_32_dot"), &
        sym(arena, "vertex_13_dot"), sym(arena, "vertex_23_dot"), &
        sym(arena, "vertex_33_dot"), sym(arena, "radius_dot"), &
        sym(arena, "xi_dot"), sym(arena, "eta_dot")]
    output_bar = [ &
        sym(arena, "point_1_bar"), sym(arena, "point_2_bar"), &
        sym(arena, "point_3_bar"), sym(arena, "tangent_xi_1_bar"), &
        sym(arena, "tangent_xi_2_bar"), sym(arena, "tangent_xi_3_bar"), &
        sym(arena, "tangent_eta_1_bar"), sym(arena, "tangent_eta_2_bar"), &
        sym(arena, "tangent_eta_3_bar"), sym(arena, "surface_jacobian_bar")]

    do component = 1, 3
        affine_point(component) = variables(component) + &
            variables(11)*(variables(component + 3) - variables(component)) + &
            variables(12)*(variables(component + 6) - variables(component))
        derivative_xi(component) = variables(component + 3) - variables(component)
        derivative_eta(component) = variables(component + 6) - variables(component)
    end do
    norm_affine = sqrt( &
        affine_point(1)*affine_point(1) + &
        affine_point(2)*affine_point(2) + &
        affine_point(3)*affine_point(3))
    dot_xi = affine_point(1)*derivative_xi(1) + &
        affine_point(2)*derivative_xi(2) + affine_point(3)*derivative_xi(3)
    dot_eta = affine_point(1)*derivative_eta(1) + &
        affine_point(2)*derivative_eta(2) + affine_point(3)*derivative_eta(3)
    norm_cubed = norm_affine*norm_affine*norm_affine
    do component = 1, 3
        outputs(component) = variables(10)*affine_point(component)/norm_affine
        tangent_xi(component) = variables(10)*( &
            derivative_xi(component)/norm_affine - &
            affine_point(component)*dot_xi/norm_cubed)
        tangent_eta(component) = variables(10)*( &
            derivative_eta(component)/norm_affine - &
            affine_point(component)*dot_eta/norm_cubed)
    end do
    outputs(4:6) = tangent_xi
    outputs(7:9) = tangent_eta
    cross_product = [ &
        tangent_xi(2)*tangent_eta(3) - tangent_xi(3)*tangent_eta(2), &
        tangent_xi(3)*tangent_eta(1) - tangent_xi(1)*tangent_eta(3), &
        tangent_xi(1)*tangent_eta(2) - tangent_xi(2)*tangent_eta(1)]
    outputs(10) = sqrt( &
        cross_product(1)*cross_product(1) + &
        cross_product(2)*cross_product(2) + &
        cross_product(3)*cross_product(3))

    jvp_roots = jvp(outputs, variables, variables_dot)
    vjp_roots = vjp(outputs, variables, output_bar)
    call simplify_all(jvp_roots)
    call simplify_all(vjp_roots)

    call initialize_spec( &
        spec, "generated_sphere_curved_panel", &
        "fortfem_generated_sphere_curved_panel", 12, 10)
    spec%args = primal_arguments()
    spec%outputs = primal_outputs()
    primal_code = chars(emit_kernel(outputs, spec))

    call initialize_spec( &
        spec, "generated_sphere_curved_panel_jvp", &
        "fortfem_generated_sphere_curved_panel_jvp", 24, 10)
    spec%args(1:12) = primal_arguments()
    spec%args(13:24) = direction_arguments()
    spec%outputs = jvp_outputs()
    jvp_code = chars(emit_kernel(jvp_roots, spec))

    call initialize_spec( &
        spec, "generated_sphere_curved_panel_vjp", &
        "fortfem_generated_sphere_curved_panel_vjp", 22, 12)
    spec%args(1:12) = primal_arguments()
    spec%args(13:22) = cotangent_arguments()
    spec%outputs = variable_bar_outputs()
    vjp_code = chars(emit_kernel(vjp_roots, spec))

    open (newunit=unit, &
        file=generated_path("fortfem_sphere_curved_panel_products.f90"), &
        status="replace", action="write", iostat=ios)
    if (ios /= 0) error stop "cannot write sphere curved panel products"
    write (unit, "(a)") primal_code(:len(primal_code) - 1)
    write (unit, "(a)") jvp_code(:len(jvp_code) - 1)
    write (unit, "(a)") vjp_code(:len(vjp_code) - 1)
    close (unit)

contains

    function primal_arguments() result(arguments)
        type(str_t), allocatable :: arguments(:)

        allocate (arguments(12))
        arguments = [ &
            str("vertex_11"), str("vertex_21"), str("vertex_31"), &
            str("vertex_12"), str("vertex_22"), str("vertex_32"), &
            str("vertex_13"), str("vertex_23"), str("vertex_33"), &
            str("radius"), str("xi"), str("eta")]
    end function primal_arguments

    function direction_arguments() result(arguments)
        type(str_t), allocatable :: arguments(:)

        allocate (arguments(12))
        arguments = [ &
            str("vertex_11_dot"), str("vertex_21_dot"), &
            str("vertex_31_dot"), str("vertex_12_dot"), &
            str("vertex_22_dot"), str("vertex_32_dot"), &
            str("vertex_13_dot"), str("vertex_23_dot"), &
            str("vertex_33_dot"), str("radius_dot"), &
            str("xi_dot"), str("eta_dot")]
    end function direction_arguments

    function cotangent_arguments() result(arguments)
        type(str_t), allocatable :: arguments(:)

        allocate (arguments(10))
        arguments = [ &
            str("point_1_bar"), str("point_2_bar"), str("point_3_bar"), &
            str("tangent_xi_1_bar"), str("tangent_xi_2_bar"), &
            str("tangent_xi_3_bar"), str("tangent_eta_1_bar"), &
            str("tangent_eta_2_bar"), str("tangent_eta_3_bar"), &
            str("surface_jacobian_bar")]
    end function cotangent_arguments

    function primal_outputs() result(names)
        type(str_t), allocatable :: names(:)

        allocate (names(10))
        names = [ &
            str("point_1"), str("point_2"), str("point_3"), &
            str("tangent_xi_1"), str("tangent_xi_2"), &
            str("tangent_xi_3"), str("tangent_eta_1"), &
            str("tangent_eta_2"), str("tangent_eta_3"), &
            str("surface_jacobian")]
    end function primal_outputs

    function jvp_outputs() result(names)
        type(str_t), allocatable :: names(:)

        allocate (names(10))
        names = [ &
            str("point_1_dot"), str("point_2_dot"), str("point_3_dot"), &
            str("tangent_xi_1_dot"), str("tangent_xi_2_dot"), &
            str("tangent_xi_3_dot"), str("tangent_eta_1_dot"), &
            str("tangent_eta_2_dot"), str("tangent_eta_3_dot"), &
            str("surface_jacobian_dot")]
    end function jvp_outputs

    function variable_bar_outputs() result(names)
        type(str_t), allocatable :: names(:)

        allocate (names(12))
        names = [ &
            str("vertex_11_bar"), str("vertex_21_bar"), str("vertex_31_bar"), &
            str("vertex_12_bar"), str("vertex_22_bar"), str("vertex_32_bar"), &
            str("vertex_13_bar"), str("vertex_23_bar"), str("vertex_33_bar"), &
            str("radius_bar"), str("xi_bar"), str("eta_bar")]
    end function variable_bar_outputs

    subroutine initialize_spec( &
            kernel_spec, name, module_name, argument_count, output_count)
        type(kernel_spec_t), intent(out) :: kernel_spec
        character(*), intent(in) :: name, module_name
        integer, intent(in) :: argument_count, output_count

        kernel_spec%name = str(name)
        kernel_spec%module_name = str(module_name)
        kernel_spec%mode = KERNEL_SUBROUTINE
        kernel_spec%generator = str("gen_sphere_curved_panel_products")
        kernel_spec%generator_revision = str(fortsym_revision())
        kernel_spec%regenerate_command = str("cd tools/codegen && ./generate.sh")
        kernel_spec%pure_procedure = .true.
        allocate ( &
            kernel_spec%args(argument_count), kernel_spec%outputs(output_count))
    end subroutine initialize_spec

    subroutine simplify_all(expressions)
        type(expr_t), intent(inout) :: expressions(:)

        do root = 1, size(expressions)
            result = engine%simplify(expressions(root))
            if (result%ok) expressions(root) = result%value
        end do
    end subroutine simplify_all

end program gen_sphere_curved_panel_products
