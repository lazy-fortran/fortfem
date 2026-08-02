program gen_fci_quintic_bezier_edge_area_products
    use, intrinsic :: iso_fortran_env, only: real64
    use fortsym_arena, only: arena_t
    use fortsym_diff, only: diff
    use fortsym_engine, only: engine_result_t
    use fortsym_engine_native, only: make_native_engine, native_engine_t
    use fortsym_expr, only: expr_t, operator(+), operator(-), operator(*), &
        operator(/), operator(**), real_expr, sym
    use fortsym_kernel, only: emit_kernel, kernel_spec_t, KERNEL_SUBROUTINE
    use fortsym_products, only: jvp, vjp
    use fortsym_string, only: chars, str
    use fortfem_codegen_provenance, only: fortsym_revision, generated_path
    implicit none

    type(arena_t), target :: arena
    type(engine_result_t) :: result
    type(native_engine_t) :: engine
    type(expr_t) :: variables(12), variables_dot(12), edge_area(1)
    type(expr_t) :: edge_area_bar(1), edge_area_jvp(1), edge_area_vjp(12)
    type(expr_t) :: x_coeff(0:5), y_coeff(0:5), curve_parameter
    type(expr_t) :: half, two, three, four, five, six, ten
    type(expr_t) :: edge_x, edge_y, derivative_x, derivative_y
    type(kernel_spec_t) :: spec
    character(:), allocatable :: primal_code, jvp_code, vjp_code
    integer :: index, degree, ios, unit

    call arena%init()
    engine = make_native_engine(arena)
    half = real_expr(arena, 0.5d0)
    two = real_expr(arena, 2.0d0)
    three = real_expr(arena, 3.0d0)
    four = real_expr(arena, 4.0d0)
    five = real_expr(arena, 5.0d0)
    six = real_expr(arena, 6.0d0)
    ten = real_expr(arena, 10.0d0)
    variables = [ &
        sym(arena, "x_start"), sym(arena, "y_start"), &
        sym(arena, "control_1_x"), sym(arena, "control_1_y"), &
        sym(arena, "control_2_x"), sym(arena, "control_2_y"), &
        sym(arena, "control_3_x"), sym(arena, "control_3_y"), &
        sym(arena, "control_4_x"), sym(arena, "control_4_y"), &
        sym(arena, "x_end"), sym(arena, "y_end")]
    curve_parameter = sym(arena, "parameter")
    x_coeff(0) = variables(1)
    x_coeff(1) = five*(variables(3) - variables(1))
    x_coeff(2) = ten*(variables(1) - two*variables(3) + variables(5))
    x_coeff(3) = ten*(-variables(1) + three*variables(3) - &
        three*variables(5) + variables(7))
    x_coeff(4) = five*(variables(1) - four*variables(3) + &
        six*variables(5) - four*variables(7) + variables(9))
    x_coeff(5) = variables(11) - variables(1) + five*variables(3) - &
        ten*variables(5) + ten*variables(7) - five*variables(9)
    y_coeff(0) = variables(2)
    y_coeff(1) = five*(variables(4) - variables(2))
    y_coeff(2) = ten*(variables(2) - two*variables(4) + variables(6))
    y_coeff(3) = ten*(-variables(2) + three*variables(4) - &
        three*variables(6) + variables(8))
    y_coeff(4) = five*(variables(2) - four*variables(4) + &
        six*variables(6) - four*variables(8) + variables(10))
    y_coeff(5) = variables(12) - variables(2) + five*variables(4) - &
        ten*variables(6) + ten*variables(8) - five*variables(10)
    edge_x = x_coeff(0)
    edge_y = y_coeff(0)
    do degree = 1, 5
        edge_x = edge_x + x_coeff(degree)*curve_parameter**degree
        edge_y = edge_y + y_coeff(degree)*curve_parameter**degree
    end do
    derivative_x = diff(edge_x, curve_parameter)
    derivative_y = diff(edge_y, curve_parameter)
    edge_area(1) = real_expr(arena, 0.0d0)
    do index = 0, 5
        do degree = 1, 5
            edge_area(1) = edge_area(1) + half*real_expr(arena, &
                real(degree, real64))*(x_coeff(index)*y_coeff(degree) - &
                y_coeff(index)*x_coeff(degree))/real_expr(arena, &
                real(index + degree, real64))
        end do
    end do
    if (edge_x%node_count() < 1 .or. derivative_x%node_count() < 1 .or. &
        edge_y%node_count() < 1 .or. derivative_y%node_count() < 1) then
        error stop "quintic Bezier generator produced an empty curve"
    end if
    call simplify_all(edge_area)

    call initialize_spec(spec, "generated_fci_quintic_bezier_edge_area", &
        "fortfem_generated_fci_quintic_bezier_edge_area", 12, 1)
    spec%args = primal_arguments()
    spec%outputs(1) = str("edge_area")
    primal_code = chars(emit_kernel(edge_area, spec))

    variables_dot = [ &
        sym(arena, "x_start_dot"), sym(arena, "y_start_dot"), &
        sym(arena, "control_1_x_dot"), sym(arena, "control_1_y_dot"), &
        sym(arena, "control_2_x_dot"), sym(arena, "control_2_y_dot"), &
        sym(arena, "control_3_x_dot"), sym(arena, "control_3_y_dot"), &
        sym(arena, "control_4_x_dot"), sym(arena, "control_4_y_dot"), &
        sym(arena, "x_end_dot"), sym(arena, "y_end_dot")]
    edge_area_bar(1) = sym(arena, "edge_area_bar")
    edge_area_jvp = jvp(edge_area, variables, variables_dot)
    edge_area_vjp = vjp(edge_area, variables, edge_area_bar)
    call simplify_all(edge_area_jvp)
    call simplify_all(edge_area_vjp)

    call initialize_spec(spec, "generated_fci_quintic_bezier_edge_area_jvp", &
        "fortfem_generated_fci_quintic_bezier_edge_area_jvp", 24, 1)
    spec%args = jvp_arguments()
    spec%outputs(1) = str("edge_area_dot")
    jvp_code = chars(emit_kernel(edge_area_jvp, spec))

    call initialize_spec(spec, "generated_fci_quintic_bezier_edge_area_vjp", &
        "fortfem_generated_fci_quintic_bezier_edge_area_vjp", 13, 12)
    spec%args = vjp_arguments()
    spec%outputs = vjp_outputs()
    vjp_code = chars(emit_kernel(edge_area_vjp, spec))

    open (newunit=unit, file=generated_path( &
        "fortfem_fci_quintic_bezier_edge_area_products.f90"), &
        status="replace", action="write", iostat=ios)
    if (ios /= 0) error stop "cannot write FCI quintic Bezier edge products"
    write (unit, "(a)") primal_code(:len(primal_code) - 1)
    write (unit, "(a)") jvp_code(:len(jvp_code) - 1)
    write (unit, "(a)") vjp_code(:len(vjp_code) - 1)
    close (unit)

contains

    function primal_arguments() result(arguments)
        use fortsym_string, only: str_t
        type(str_t), allocatable :: arguments(:)

        allocate(arguments(12))
        arguments = [str("x_start"), str("y_start"), str("control_1_x"), &
            str("control_1_y"), str("control_2_x"), str("control_2_y"), &
            str("control_3_x"), str("control_3_y"), str("control_4_x"), &
            str("control_4_y"), str("x_end"), str("y_end")]
    end function primal_arguments

    function jvp_arguments() result(arguments)
        use fortsym_string, only: str_t
        type(str_t), allocatable :: arguments(:)

        allocate(arguments(24))
        arguments(1:12) = primal_arguments()
        arguments(13:24) = [str("x_start_dot"), str("y_start_dot"), &
            str("control_1_x_dot"), str("control_1_y_dot"), &
            str("control_2_x_dot"), str("control_2_y_dot"), &
            str("control_3_x_dot"), str("control_3_y_dot"), &
            str("control_4_x_dot"), str("control_4_y_dot"), &
            str("x_end_dot"), str("y_end_dot")]
    end function jvp_arguments

    function vjp_arguments() result(arguments)
        use fortsym_string, only: str_t
        type(str_t), allocatable :: arguments(:)

        allocate(arguments(13))
        arguments(1:12) = primal_arguments()
        arguments(13) = str("edge_area_bar")
    end function vjp_arguments

    function vjp_outputs() result(outputs)
        use fortsym_string, only: str_t
        type(str_t), allocatable :: outputs(:)

        allocate(outputs(12))
        outputs = [str("x_start_bar"), str("y_start_bar"), &
            str("control_1_x_bar"), str("control_1_y_bar"), &
            str("control_2_x_bar"), str("control_2_y_bar"), &
            str("control_3_x_bar"), str("control_3_y_bar"), &
            str("control_4_x_bar"), str("control_4_y_bar"), &
            str("x_end_bar"), str("y_end_bar")]
    end function vjp_outputs

    subroutine initialize_spec(kernel_spec, name, module_name, argument_count, &
            output_count)
        type(kernel_spec_t), intent(out) :: kernel_spec
        character(*), intent(in) :: name, module_name
        integer, intent(in) :: argument_count, output_count

        kernel_spec%name = str(name)
        kernel_spec%module_name = str(module_name)
        kernel_spec%mode = KERNEL_SUBROUTINE
        kernel_spec%generator = str("gen_fci_quintic_bezier_edge_area_products")
        kernel_spec%generator_revision = str(fortsym_revision())
        kernel_spec%regenerate_command = str( &
            "cd tools/codegen && ./generate.sh")
        kernel_spec%pure_procedure = .true.
        allocate(kernel_spec%args(argument_count), kernel_spec%outputs(output_count))
    end subroutine initialize_spec

    subroutine simplify_all(expressions)
        type(expr_t), intent(inout) :: expressions(:)

        do index = 1, size(expressions)
            result = engine%simplify(expressions(index))
            if (result%ok) expressions(index) = result%value
        end do
    end subroutine simplify_all

end program gen_fci_quintic_bezier_edge_area_products
