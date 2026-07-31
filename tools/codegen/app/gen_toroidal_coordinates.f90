program gen_toroidal_coordinates
    use fortsym_arena, only: arena_t
    use fortsym_chart, only: chart_t
    use fortsym_expr, only: atan2, cos, cosh, expr_t, log, operator(*), &
        operator(+), operator(-), operator(/), operator(**), real_expr, sin, &
        sinh, sqrt, sym
    use fortsym_kernel, only: emit_kernel, kernel_spec_t, KERNEL_SUBROUTINE
    use fortsym_products, only: jvp, vjp
    use fortsym_string, only: chars, str
    use fortsym_toroidal, only: make_toroidal_chart
    use fortfem_codegen_provenance, only: fortsym_revision, generated_path
    implicit none

    type(arena_t), target :: arena
    type(chart_t) :: chart
    type(expr_t) :: forward_variables(4), forward_directions(4)
    type(expr_t) :: forward_outputs(3), forward_output_bars(3)
    type(expr_t) :: inverse_variables(4), inverse_directions(4)
    type(expr_t) :: inverse_outputs(3), inverse_output_bars(3)
    type(expr_t) :: forward_jvp(3), forward_vjp(4)
    type(expr_t) :: inverse_jvp(3), inverse_vjp(4)
    type(expr_t) :: vector_variables(6), vector_directions(6)
    type(expr_t) :: vector_outputs(3), vector_output_bars(3)
    type(expr_t) :: vector_jvp(3), vector_vjp(6)
    type(expr_t) :: eta, phi, scale, theta, cylindrical_radius
    type(expr_t) :: point_x, point_y, point_z, eta_expression
    type(expr_t) :: theta_expression, phi_expression
    type(expr_t) :: denominator, eta_radial, eta_vertical
    type(expr_t) :: theta_radial, theta_vertical
    type(expr_t) :: component_1, component_2, component_3
    type(expr_t) :: half, two
    type(kernel_spec_t) :: spec
    character(:), allocatable :: forward_code, forward_jvp_code
    character(:), allocatable :: forward_vjp_code, inverse_code
    character(:), allocatable :: inverse_jvp_code, inverse_vjp_code
    character(:), allocatable :: vector_code, vector_jvp_code, vector_vjp_code
    integer :: ios, unit

    call arena%init()
    half = real_expr(arena, 0.5d0)
    two = real_expr(arena, 2.0d0)
    eta = sym(arena, "eta")
    theta = sym(arena, "theta")
    phi = sym(arena, "phi")
    scale = sym(arena, "scale")
    forward_variables = [scale, eta, theta, phi]
    forward_directions = [ &
        sym(arena, "scale_dot"), sym(arena, "eta_dot"), &
        sym(arena, "theta_dot"), sym(arena, "phi_dot")]
    chart = make_toroidal_chart(arena, eta, theta, phi, scale)
    forward_outputs = chart%x
    forward_output_bars = [ &
        sym(arena, "point_1_bar"), sym(arena, "point_2_bar"), &
        sym(arena, "point_3_bar")]
    forward_jvp = jvp(forward_outputs, forward_variables, forward_directions)
    forward_vjp = vjp(forward_outputs, forward_variables, forward_output_bars)

    point_x = sym(arena, "point_1")
    point_y = sym(arena, "point_2")
    point_z = sym(arena, "point_3")
    inverse_variables = [point_x, point_y, point_z, scale]
    inverse_directions = [ &
        sym(arena, "point_1_dot"), sym(arena, "point_2_dot"), &
        sym(arena, "point_3_dot"), sym(arena, "scale_dot")]
    cylindrical_radius = sqrt(point_x*point_x + point_y*point_y)
    eta_expression = log( &
        ((cylindrical_radius + scale)**2 + point_z*point_z)/ &
        ((cylindrical_radius - scale)**2 + point_z*point_z))*half
    theta_expression = atan2( &
        two*scale*point_z, cylindrical_radius*cylindrical_radius + &
        point_z*point_z - scale*scale)
    phi_expression = atan2(point_y, point_x)
    inverse_outputs = [eta_expression, theta_expression, phi_expression]
    inverse_output_bars = [ &
        sym(arena, "eta_bar"), sym(arena, "theta_bar"), sym(arena, "phi_bar")]
    inverse_jvp = jvp(inverse_outputs, inverse_variables, inverse_directions)
    inverse_vjp = vjp(inverse_outputs, inverse_variables, inverse_output_bars)

    component_1 = sym(arena, "component_1")
    component_2 = sym(arena, "component_2")
    component_3 = sym(arena, "component_3")
    vector_variables = [eta, theta, phi, component_1, component_2, component_3]
    vector_directions = [ &
        sym(arena, "eta_dot"), sym(arena, "theta_dot"), sym(arena, "phi_dot"), &
        sym(arena, "component_1_dot"), sym(arena, "component_2_dot"), &
        sym(arena, "component_3_dot")]
    denominator = cosh(eta) - cos(theta)
    eta_radial = (real_expr(arena, 1.0d0) - cosh(eta)*cos(theta))/denominator
    eta_vertical = -sinh(eta)*sin(theta)/denominator
    theta_radial = eta_vertical
    theta_vertical = -eta_radial
    vector_outputs = [ &
        component_1*eta_radial*cos(phi) + component_2*theta_radial*cos(phi) - &
        component_3*sin(phi), &
        component_1*eta_radial*sin(phi) + component_2*theta_radial*sin(phi) + &
        component_3*cos(phi), &
        component_1*eta_vertical + component_2*theta_vertical]
    vector_output_bars = [ &
        sym(arena, "cartesian_1_bar"), sym(arena, "cartesian_2_bar"), &
        sym(arena, "cartesian_3_bar")]
    vector_jvp = jvp(vector_outputs, vector_variables, vector_directions)
    vector_vjp = vjp(vector_outputs, vector_variables, vector_output_bars)

    call initialize_spec( &
        spec, "generated_toroidal_point_to_cartesian", &
        "fortfem_generated_toroidal_coordinates", 4, 3)
    spec%args = [str("scale"), str("eta"), str("theta"), str("phi")]
    spec%outputs = [str("point")]
    spec%output_shapes = [str("(3)")]
    spec%output_references = [ &
        str("point(1)"), str("point(2)"), str("point(3)")]
    forward_code = chars(emit_kernel(forward_outputs, spec))

    call initialize_spec( &
        spec, "generated_toroidal_point_to_cartesian_jvp", &
        "fortfem_generated_toroidal_coordinates_jvp", 8, 3)
    spec%args = [ &
        str("scale"), str("eta"), str("theta"), str("phi"), &
        str("scale_dot"), str("eta_dot"), str("theta_dot"), str("phi_dot")]
    spec%outputs = [str("point_dot")]
    spec%output_shapes = [str("(3)")]
    spec%output_references = [ &
        str("point_dot(1)"), str("point_dot(2)"), str("point_dot(3)")]
    forward_jvp_code = chars(emit_kernel(forward_jvp, spec))

    call initialize_spec( &
        spec, "generated_toroidal_point_to_cartesian_vjp", &
        "fortfem_generated_toroidal_coordinates_vjp", 7, 4)
    spec%args = [ &
        str("scale"), str("eta"), str("theta"), str("phi"), &
        str("point_1_bar"), str("point_2_bar"), str("point_3_bar")]
    spec%outputs = [ &
        str("scale_bar"), str("eta_bar"), str("theta_bar"), str("phi_bar")]
    forward_vjp_code = chars(emit_kernel(forward_vjp, spec))

    call initialize_spec( &
        spec, "generated_cartesian_to_toroidal", &
        "fortfem_generated_cartesian_to_toroidal", 4, 3)
    spec%args = [str("point_1"), str("point_2"), str("point_3"), str("scale")]
    spec%outputs = [str("eta"), str("theta"), str("phi")]
    inverse_code = chars(emit_kernel(inverse_outputs, spec))

    call initialize_spec( &
        spec, "generated_cartesian_to_toroidal_jvp", &
        "fortfem_generated_cartesian_to_toroidal_jvp", 8, 3)
    spec%args = [ &
        str("point_1"), str("point_2"), str("point_3"), str("scale"), &
        str("point_1_dot"), str("point_2_dot"), str("point_3_dot"), &
        str("scale_dot")]
    spec%outputs = [str("eta_dot"), str("theta_dot"), str("phi_dot")]
    inverse_jvp_code = chars(emit_kernel(inverse_jvp, spec))

    call initialize_spec( &
        spec, "generated_cartesian_to_toroidal_vjp", &
        "fortfem_generated_cartesian_to_toroidal_vjp", 7, 4)
    spec%args = [ &
        str("point_1"), str("point_2"), str("point_3"), str("scale"), &
        str("eta_bar"), str("theta_bar"), str("phi_bar")]
    spec%outputs = [ &
        str("point_1_bar"), str("point_2_bar"), str("point_3_bar"), &
        str("scale_bar")]
    inverse_vjp_code = chars(emit_kernel(inverse_vjp, spec))

    call initialize_spec( &
        spec, "generated_toroidal_vector_to_cartesian", &
        "fortfem_generated_toroidal_vector_to_cartesian", 6, 3)
    spec%args = [ &
        str("eta"), str("theta"), str("phi"), str("component_1"), &
        str("component_2"), str("component_3")]
    spec%outputs = [str("cartesian_1"), str("cartesian_2"), str("cartesian_3")]
    vector_code = chars(emit_kernel(vector_outputs, spec))

    call initialize_spec( &
        spec, "generated_toroidal_vector_to_cartesian_jvp", &
        "fortfem_generated_toroidal_vector_to_cartesian_jvp", 12, 3)
    spec%args = [ &
        str("eta"), str("theta"), str("phi"), str("component_1"), &
        str("component_2"), str("component_3"), str("eta_dot"), &
        str("theta_dot"), str("phi_dot"), str("component_1_dot"), &
        str("component_2_dot"), str("component_3_dot")]
    spec%outputs = [ &
        str("cartesian_1_dot"), str("cartesian_2_dot"), &
        str("cartesian_3_dot")]
    vector_jvp_code = chars(emit_kernel(vector_jvp, spec))

    call initialize_spec( &
        spec, "generated_toroidal_vector_to_cartesian_vjp", &
        "fortfem_generated_toroidal_vector_to_cartesian_vjp", 9, 6)
    spec%args = [ &
        str("eta"), str("theta"), str("phi"), str("component_1"), &
        str("component_2"), str("component_3"), str("cartesian_1_bar"), &
        str("cartesian_2_bar"), str("cartesian_3_bar")]
    spec%outputs = [ &
        str("eta_bar"), str("theta_bar"), str("phi_bar"), &
        str("component_1_bar"), str("component_2_bar"), &
        str("component_3_bar")]
    vector_vjp_code = chars(emit_kernel(vector_vjp, spec))

    open (newunit=unit, &
        file=generated_path("fortfem_toroidal_coordinates.f90"), &
        status="replace", action="write", iostat=ios)
    if (ios /= 0) error stop "cannot write toroidal coordinate kernels"
    write (unit, "(a)") forward_code(:len(forward_code) - 1)
    write (unit, "(a)") forward_jvp_code(:len(forward_jvp_code) - 1)
    write (unit, "(a)") forward_vjp_code(:len(forward_vjp_code) - 1)
    write (unit, "(a)") inverse_code(:len(inverse_code) - 1)
    write (unit, "(a)") inverse_jvp_code(:len(inverse_jvp_code) - 1)
    write (unit, "(a)") inverse_vjp_code(:len(inverse_vjp_code) - 1)
    write (unit, "(a)") vector_code(:len(vector_code) - 1)
    write (unit, "(a)") vector_jvp_code(:len(vector_jvp_code) - 1)
    write (unit, "(a)") vector_vjp_code(:len(vector_vjp_code) - 1)
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
        kernel_spec%generator = str("gen_toroidal_coordinates")
        kernel_spec%generator_revision = str(fortsym_revision())
        kernel_spec%regenerate_command = str("cd tools/codegen && ./generate.sh")
        kernel_spec%pure_procedure = .true.
        allocate( &
            kernel_spec%args(argument_count), kernel_spec%outputs(output_count))
    end subroutine initialize_spec

end program gen_toroidal_coordinates
