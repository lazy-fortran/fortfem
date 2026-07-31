program gen_helmholtz_bem_smooth_products
    use fortsym_arena, only: arena_t
    use fortsym_engine, only: engine_result_t
    use fortsym_engine_native, only: make_native_engine, native_engine_t
    use fortsym_expr, only: cos, expr_t, num, operator(+), operator(-), &
        operator(*), operator(/), sin, sqrt, sym
    use fortsym_kernel, only: emit_kernel, kernel_spec_t, KERNEL_SUBROUTINE
    use fortsym_products, only: jvp, vjp
    use fortsym_string, only: chars, str, str_t
    use fortfem_codegen_provenance, only: fortsym_revision, generated_path
    implicit none

    type(arena_t), target :: arena
    type(native_engine_t) :: engine
    type(engine_result_t) :: result
    type(expr_t) :: variables(9), variables_dot(9), displacement(3)
    type(expr_t) :: amplitude, kernel_scale, outputs(2), output_bar(2), phase
    type(expr_t) :: jvp_roots(2), vjp_roots(9), radius
    type(kernel_spec_t) :: spec
    character(:), allocatable :: primal_code, jvp_code, vjp_code
    integer :: component, ios, root, unit

    call arena%init()
    engine = make_native_engine(arena)
    variables = [ &
        sym(arena, "first_point_1"), sym(arena, "first_point_2"), &
        sym(arena, "first_point_3"), sym(arena, "second_point_1"), &
        sym(arena, "second_point_2"), sym(arena, "second_point_3"), &
        sym(arena, "first_jacobian"), sym(arena, "second_jacobian"), &
        sym(arena, "wave_number")]
    variables_dot = [ &
        sym(arena, "first_point_1_dot"), sym(arena, "first_point_2_dot"), &
        sym(arena, "first_point_3_dot"), sym(arena, "second_point_1_dot"), &
        sym(arena, "second_point_2_dot"), sym(arena, "second_point_3_dot"), &
        sym(arena, "first_jacobian_dot"), &
        sym(arena, "second_jacobian_dot"), sym(arena, "wave_number_dot")]
    kernel_scale = sym(arena, "kernel_scale")
    output_bar = [sym(arena, "value_real_bar"), sym(arena, "value_imag_bar")]
    do component = 1, 3
        displacement(component) = variables(component) - &
            variables(component + 3)
    end do
    radius = sqrt( &
        displacement(1)*displacement(1) + &
        displacement(2)*displacement(2) + &
        displacement(3)*displacement(3))
    amplitude = kernel_scale*variables(7)*variables(8)/radius
    phase = variables(9)*radius
    outputs = [ &
        amplitude*(cos(phase) - num(arena, 1)), amplitude*sin(phase)]

    jvp_roots = jvp(outputs, variables, variables_dot)
    vjp_roots = vjp(outputs, variables, output_bar)
    call simplify_all(jvp_roots)
    call simplify_all(vjp_roots)

    call initialize_spec( &
        spec, "generated_helmholtz_single_layer_smooth_integrand", &
        "fortfem_generated_helmholtz_single_layer_smooth_integrand", 10, 2)
    spec%args = primal_arguments()
    spec%outputs = [str("value_real"), str("value_imag")]
    primal_code = chars(emit_kernel(outputs, spec))

    call initialize_spec( &
        spec, "generated_helmholtz_single_layer_smooth_integrand_jvp", &
        "fortfem_generated_helmholtz_single_layer_smooth_integrand_jvp", &
        19, 2)
    spec%args(1:10) = primal_arguments()
    spec%args(11:19) = tangent_arguments()
    spec%outputs = [str("value_real_dot"), str("value_imag_dot")]
    jvp_code = chars(emit_kernel(jvp_roots, spec))

    call initialize_spec( &
        spec, "generated_helmholtz_single_layer_smooth_integrand_vjp", &
        "fortfem_generated_helmholtz_single_layer_smooth_integrand_vjp", &
        12, 9)
    spec%args(1:10) = primal_arguments()
    spec%args(11:12) = [str("value_real_bar"), str("value_imag_bar")]
    spec%outputs = [ &
        str("first_point_1_bar"), str("first_point_2_bar"), &
        str("first_point_3_bar"), str("second_point_1_bar"), &
        str("second_point_2_bar"), str("second_point_3_bar"), &
        str("first_jacobian_bar"), str("second_jacobian_bar"), &
        str("wave_number_bar")]
    vjp_code = chars(emit_kernel(vjp_roots, spec))

    open (newunit=unit, &
        file=generated_path("fortfem_helmholtz_bem_smooth_products.f90"), &
        status="replace", action="write", iostat=ios)
    if (ios /= 0) error stop "cannot write Helmholtz BEM smooth products"
    write (unit, "(a)") primal_code(:len(primal_code) - 1)
    write (unit, "(a)") jvp_code(:len(jvp_code) - 1)
    write (unit, "(a)") vjp_code(:len(vjp_code) - 1)
    close (unit)

contains

    function primal_arguments() result(arguments)
        type(str_t), allocatable :: arguments(:)

        allocate(arguments(10))
        arguments = [ &
            str("first_point_1"), str("first_point_2"), &
            str("first_point_3"), str("second_point_1"), &
            str("second_point_2"), str("second_point_3"), &
            str("first_jacobian"), str("second_jacobian"), &
            str("wave_number"), str("kernel_scale")]
    end function primal_arguments

    function tangent_arguments() result(arguments)
        type(str_t), allocatable :: arguments(:)

        allocate(arguments(9))
        arguments = [ &
            str("first_point_1_dot"), str("first_point_2_dot"), &
            str("first_point_3_dot"), str("second_point_1_dot"), &
            str("second_point_2_dot"), str("second_point_3_dot"), &
            str("first_jacobian_dot"), str("second_jacobian_dot"), &
            str("wave_number_dot")]
    end function tangent_arguments

    subroutine initialize_spec( &
            kernel_spec, name, module_name, argument_count, output_count)
        type(kernel_spec_t), intent(out) :: kernel_spec
        character(*), intent(in) :: name, module_name
        integer, intent(in) :: argument_count, output_count

        kernel_spec%name = str(name)
        kernel_spec%module_name = str(module_name)
        kernel_spec%mode = KERNEL_SUBROUTINE
        kernel_spec%generator = str("gen_helmholtz_bem_smooth_products")
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

end program gen_helmholtz_bem_smooth_products
