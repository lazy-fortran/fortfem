program gen_toroidal_poisson_products
    use fortsym_arena, only: arena_t
    use fortsym_engine, only: engine_result_t
    use fortsym_engine_native, only: make_native_engine, native_engine_t
    use fortsym_expr, only: cos, cosh, expr_t, operator(*), operator(+), &
        operator(-), operator(/), real_expr, sin, sinh, sqrt, sym
    use fortsym_kernel, only: emit_kernel, kernel_spec_t, KERNEL_SUBROUTINE
    use fortsym_products, only: jvp, vjp
    use fortsym_string, only: chars, str
    use fortfem_codegen_provenance, only: fortsym_revision, generated_path
    implicit none

    type(arena_t), target :: arena
    type(engine_result_t) :: result
    type(native_engine_t) :: engine
    type(expr_t) :: variables(8), variables_dot(8), outputs(6), output_bar(6)
    type(expr_t) :: jvp_roots(6), vjp_roots(8)
    type(expr_t) :: denominator, root_denominator, harmonic
    type(expr_t) :: eta_derivative, theta_derivative, phi_derivative
    type(expr_t) :: dtn_value, normal_derivative
    type(kernel_spec_t) :: spec
    character(:), allocatable :: primal_code, jvp_code, vjp_code
    integer :: ios, root, unit

    call arena%init()
    engine = make_native_engine(arena)
    variables = [ &
        sym(arena, "degree"), sym(arena, "order"), sym(arena, "scale"), &
        sym(arena, "eta"), sym(arena, "theta"), sym(arena, "phi"), &
        sym(arena, "radial"), sym(arena, "radial_derivative")]
    variables_dot = [ &
        sym(arena, "degree_dot"), sym(arena, "order_dot"), &
        sym(arena, "scale_dot"), sym(arena, "eta_dot"), &
        sym(arena, "theta_dot"), sym(arena, "phi_dot"), &
        sym(arena, "radial_dot"), sym(arena, "radial_derivative_dot")]
    output_bar = [ &
        sym(arena, "harmonic_bar"), sym(arena, "field_1_bar"), &
        sym(arena, "field_2_bar"), sym(arena, "field_3_bar"), &
        sym(arena, "dtn_bar"), sym(arena, "normal_derivative_bar")]

    denominator = cosh(variables(4)) - cos(variables(5))
    root_denominator = sqrt(denominator)
    harmonic = root_denominator*variables(7)* &
        cos(variables(1)*variables(5))*cos(variables(2)*variables(6))
    eta_derivative = ( &
        sinh(variables(4))*variables(7)/(2*root_denominator) + &
        root_denominator*variables(8)*sinh(variables(4)))* &
        cos(variables(1)*variables(5))*cos(variables(2)*variables(6))
    theta_derivative = ( &
        sin(variables(5))*variables(7)*cos(variables(1)*variables(5))/ &
        (2*root_denominator) - root_denominator*variables(7)* &
        variables(1)*sin(variables(1)*variables(5)))* &
        cos(variables(2)*variables(6))
    phi_derivative = -root_denominator*variables(7)*variables(2)* &
        cos(variables(1)*variables(5))*sin(variables(2)*variables(6))
    dtn_value = -denominator/variables(3)*( &
        sinh(variables(4))/(2*denominator) + &
        sinh(variables(4))*variables(8)/variables(7))
    normal_derivative = dtn_value*harmonic
    outputs = [harmonic, -denominator*eta_derivative/variables(3), &
        -denominator*theta_derivative/variables(3), &
        -denominator*phi_derivative/(variables(3)*sinh(variables(4))), &
        dtn_value, normal_derivative]

    jvp_roots = jvp(outputs, variables, variables_dot)
    vjp_roots = vjp(outputs, variables, output_bar)
    call simplify_all(jvp_roots)
    call simplify_all(vjp_roots)

    call initialize_spec( &
        spec, "generated_toroidal_poisson_products", &
        "fortfem_generated_toroidal_poisson_products", 8, 6)
    spec%args = primal_arguments()
    spec%outputs = primal_outputs()
    primal_code = chars(emit_kernel(outputs, spec))

    call initialize_spec( &
        spec, "generated_toroidal_poisson_products_jvp", &
        "fortfem_generated_toroidal_poisson_products_jvp", 16, 6)
    spec%args = [primal_arguments(), direction_arguments()]
    spec%outputs = jvp_outputs()
    jvp_code = chars(emit_kernel(jvp_roots, spec))

    call initialize_spec( &
        spec, "generated_toroidal_poisson_products_vjp", &
        "fortfem_generated_toroidal_poisson_products_vjp", 14, 8)
    spec%args = [primal_arguments(), cotangent_arguments()]
    spec%outputs = variable_bar_outputs()
    vjp_code = chars(emit_kernel(vjp_roots, spec))

    open (newunit=unit, &
        file=generated_path("fortfem_toroidal_poisson_products.f90"), &
        status="replace", action="write", iostat=ios)
    if (ios /= 0) error stop "cannot write toroidal Poisson products"
    write (unit, "(a)") primal_code(:len(primal_code) - 1)
    write (unit, "(a)") jvp_code(:len(jvp_code) - 1)
    write (unit, "(a)") vjp_code(:len(vjp_code) - 1)
    close (unit)

contains

    function primal_arguments() result(arguments)
        use fortsym_string, only: str_t
        type(str_t), allocatable :: arguments(:)

        allocate(arguments(8))
        arguments = [ &
            str("degree"), str("order"), str("scale"), str("eta"), &
            str("theta"), str("phi"), str("radial"), &
            str("radial_derivative")]
    end function primal_arguments

    function direction_arguments() result(arguments)
        use fortsym_string, only: str_t
        type(str_t), allocatable :: arguments(:)

        allocate(arguments(8))
        arguments = [ &
            str("degree_dot"), str("order_dot"), str("scale_dot"), &
            str("eta_dot"), str("theta_dot"), str("phi_dot"), &
            str("radial_dot"), str("radial_derivative_dot")]
    end function direction_arguments

    function cotangent_arguments() result(arguments)
        use fortsym_string, only: str_t
        type(str_t), allocatable :: arguments(:)

        allocate(arguments(6))
        arguments = [ &
            str("harmonic_bar"), str("field_1_bar"), str("field_2_bar"), &
            str("field_3_bar"), str("dtn_bar"), &
            str("normal_derivative_bar")]
    end function cotangent_arguments

    function primal_outputs() result(names)
        use fortsym_string, only: str_t
        type(str_t), allocatable :: names(:)

        allocate(names(6))
        names = [str("harmonic"), str("field_1"), str("field_2"), &
            str("field_3"), str("dtn_value"), str("normal_derivative")]
    end function primal_outputs

    function jvp_outputs() result(names)
        use fortsym_string, only: str_t
        type(str_t), allocatable :: names(:)

        allocate(names(6))
        names = [str("harmonic_dot"), str("field_1_dot"), &
            str("field_2_dot"), str("field_3_dot"), str("dtn_value_dot"), &
            str("normal_derivative_dot")]
    end function jvp_outputs

    function variable_bar_outputs() result(names)
        use fortsym_string, only: str_t
        type(str_t), allocatable :: names(:)

        allocate(names(8))
        names = [str("degree_bar"), str("order_bar"), str("scale_bar"), &
            str("eta_bar"), str("theta_bar"), str("phi_bar"), &
            str("radial_bar"), str("radial_derivative_bar")]
    end function variable_bar_outputs

    subroutine initialize_spec( &
            kernel_spec, name, module_name, argument_count, output_count)
        use fortsym_string, only: str_t
        type(kernel_spec_t), intent(out) :: kernel_spec
        character(*), intent(in) :: name, module_name
        integer, intent(in) :: argument_count, output_count

        kernel_spec%name = str(name)
        kernel_spec%module_name = str(module_name)
        kernel_spec%mode = KERNEL_SUBROUTINE
        kernel_spec%generator = str("gen_toroidal_poisson_products")
        kernel_spec%generator_revision = str(fortsym_revision())
        kernel_spec%regenerate_command = str("cd tools/codegen && ./generate.sh")
        kernel_spec%pure_procedure = .true.
        allocate(kernel_spec%args(argument_count), kernel_spec%outputs(output_count))
    end subroutine initialize_spec

    subroutine simplify_all(expressions)
        type(expr_t), intent(inout) :: expressions(:)

        do root = 1, size(expressions)
            result = engine%simplify(expressions(root))
            if (result%ok) expressions(root) = result%value
        end do
    end subroutine simplify_all

end program gen_toroidal_poisson_products
