program gen_torus_curved_panel_products
    use fortsym_arena, only: arena_t
    use fortsym_engine, only: engine_result_t
    use fortsym_engine_native, only: make_native_engine, native_engine_t
    use fortsym_expr, only: cos, expr_t, operator(+), operator(-), operator(*), &
        operator(/), real_expr, sin, sqrt, sym
    use fortsym_kernel, only: emit_kernel, kernel_spec_t, KERNEL_SUBROUTINE
    use fortsym_products, only: jvp, vjp
    use fortsym_string, only: chars, str, str_t
    use fortfem_codegen_provenance, only: fortsym_revision, generated_path
    implicit none

    type(arena_t), target :: arena
    type(engine_result_t) :: result
    type(native_engine_t) :: engine
    type(expr_t) :: variables(10), variables_dot(10), outputs(10)
    type(expr_t) :: output_bar(10), jvp_roots(10), vjp_roots(10)
    type(expr_t) :: theta, phi, radial, derivative_theta(3), derivative_phi(3)
    type(expr_t) :: tangent_xi(3), tangent_eta(3), cross_product(3)
    type(kernel_spec_t) :: spec
    character(:), allocatable :: primal_code, jvp_code, vjp_code
    integer :: component, ios, root, unit

    call arena%init()
    engine = make_native_engine(arena)
    variables = [ &
        sym(arena, "parameter_11"), sym(arena, "parameter_21"), &
        sym(arena, "parameter_12"), sym(arena, "parameter_22"), &
        sym(arena, "parameter_13"), sym(arena, "parameter_23"), &
        sym(arena, "major_radius"), sym(arena, "minor_radius"), &
        sym(arena, "xi"), sym(arena, "eta")]
    variables_dot = [ &
        sym(arena, "parameter_11_dot"), sym(arena, "parameter_21_dot"), &
        sym(arena, "parameter_12_dot"), sym(arena, "parameter_22_dot"), &
        sym(arena, "parameter_13_dot"), sym(arena, "parameter_23_dot"), &
        sym(arena, "major_radius_dot"), sym(arena, "minor_radius_dot"), &
        sym(arena, "xi_dot"), sym(arena, "eta_dot")]
    output_bar = [ &
        sym(arena, "point_1_bar"), sym(arena, "point_2_bar"), &
        sym(arena, "point_3_bar"), sym(arena, "tangent_xi_1_bar"), &
        sym(arena, "tangent_xi_2_bar"), sym(arena, "tangent_xi_3_bar"), &
        sym(arena, "tangent_eta_1_bar"), sym(arena, "tangent_eta_2_bar"), &
        sym(arena, "tangent_eta_3_bar"), sym(arena, "surface_jacobian_bar")]

    theta = variables(1) + variables(9)*(variables(3) - variables(1)) + &
        variables(10)*(variables(5) - variables(1))
    phi = variables(2) + variables(9)*(variables(4) - variables(2)) + &
        variables(10)*(variables(6) - variables(2))
    radial = variables(7) + variables(8)*cos(theta)
    outputs(1:3) = [radial*cos(phi), radial*sin(phi), variables(8)*sin(theta)]
    derivative_theta = [ &
        -variables(8)*sin(theta)*cos(phi), &
        -variables(8)*sin(theta)*sin(phi), variables(8)*cos(theta)]
    derivative_phi = [&
        -radial*sin(phi), radial*cos(phi), zero_expr(arena)]
    do component = 1, 3
        tangent_xi(component) = (variables(3) - variables(1))* &
            derivative_theta(component) + (variables(4) - variables(2))* &
            derivative_phi(component)
        tangent_eta(component) = (variables(5) - variables(1))* &
            derivative_theta(component) + (variables(6) - variables(2))* &
            derivative_phi(component)
    end do
    cross_product = [ &
        tangent_xi(2)*tangent_eta(3) - tangent_xi(3)*tangent_eta(2), &
        tangent_xi(3)*tangent_eta(1) - tangent_xi(1)*tangent_eta(3), &
        tangent_xi(1)*tangent_eta(2) - tangent_xi(2)*tangent_eta(1)]
    outputs(4:6) = tangent_xi
    outputs(7:9) = tangent_eta
    outputs(10) = sqrt( &
        cross_product(1)*cross_product(1) + &
        cross_product(2)*cross_product(2) + &
        cross_product(3)*cross_product(3))

    jvp_roots = jvp(outputs, variables, variables_dot)
    vjp_roots = vjp(outputs, variables, output_bar)
    call simplify_all(jvp_roots)
    call simplify_all(vjp_roots)

    call initialize_spec( &
        spec, "generated_torus_curved_panel", &
        "fortfem_generated_torus_curved_panel", 10, 10)
    spec%args = primal_arguments()
    spec%outputs = primal_outputs()
    primal_code = chars(emit_kernel(outputs, spec))

    call initialize_spec( &
        spec, "generated_torus_curved_panel_jvp", &
        "fortfem_generated_torus_curved_panel_jvp", 20, 10)
    spec%args(1:10) = primal_arguments()
    spec%args(11:20) = direction_arguments()
    spec%outputs = jvp_outputs()
    jvp_code = chars(emit_kernel(jvp_roots, spec))

    call initialize_spec( &
        spec, "generated_torus_curved_panel_vjp", &
        "fortfem_generated_torus_curved_panel_vjp", 20, 10)
    spec%args(1:10) = primal_arguments()
    spec%args(11:20) = cotangent_arguments()
    spec%outputs = variable_bar_outputs()
    vjp_code = chars(emit_kernel(vjp_roots, spec))

    open (newunit=unit, &
        file=generated_path("fortfem_torus_curved_panel_products.f90"), &
        status="replace", action="write", iostat=ios)
    if (ios /= 0) error stop "cannot write torus curved panel products"
    write (unit, "(a)") primal_code(:len(primal_code) - 1)
    write (unit, "(a)") jvp_code(:len(jvp_code) - 1)
    write (unit, "(a)") vjp_code(:len(vjp_code) - 1)
    close (unit)

contains

    function primal_arguments() result(arguments)
        type(str_t), allocatable :: arguments(:)

        allocate (arguments(10))
        arguments = [ &
            str("parameter_11"), str("parameter_21"), &
            str("parameter_12"), str("parameter_22"), &
            str("parameter_13"), str("parameter_23"), &
            str("major_radius"), str("minor_radius"), &
            str("xi"), str("eta")]
    end function primal_arguments

    function direction_arguments() result(arguments)
        type(str_t), allocatable :: arguments(:)

        allocate (arguments(10))
        arguments = [ &
            str("parameter_11_dot"), str("parameter_21_dot"), &
            str("parameter_12_dot"), str("parameter_22_dot"), &
            str("parameter_13_dot"), str("parameter_23_dot"), &
            str("major_radius_dot"), str("minor_radius_dot"), &
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

        allocate (names(10))
        names = [ &
            str("parameter_11_bar"), str("parameter_21_bar"), &
            str("parameter_12_bar"), str("parameter_22_bar"), &
            str("parameter_13_bar"), str("parameter_23_bar"), &
            str("major_radius_bar"), str("minor_radius_bar"), &
            str("xi_bar"), str("eta_bar")]
    end function variable_bar_outputs

    subroutine initialize_spec( &
            kernel_spec, name, module_name, argument_count, output_count)
        type(kernel_spec_t), intent(out) :: kernel_spec
        character(*), intent(in) :: name, module_name
        integer, intent(in) :: argument_count, output_count

        kernel_spec%name = str(name)
        kernel_spec%module_name = str(module_name)
        kernel_spec%mode = KERNEL_SUBROUTINE
        kernel_spec%generator = str("gen_torus_curved_panel_products")
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

    function zero_expr(arena) result(value)
        type(arena_t), intent(inout), target :: arena
        type(expr_t) :: value

        value = real_expr(arena, 0.0d0)
    end function zero_expr

end program gen_torus_curved_panel_products
