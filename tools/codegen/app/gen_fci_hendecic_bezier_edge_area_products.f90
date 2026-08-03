program gen_fci_hendecic_bezier_edge_area_products
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

    integer, parameter :: polynomial_degree = 11
    integer, parameter :: variable_count = 2*(polynomial_degree + 1)
    type(arena_t), target :: arena
    type(engine_result_t) :: result
    type(native_engine_t) :: engine
    type(expr_t) :: variables(variable_count), variables_dot(variable_count)
    type(expr_t) :: control_x(0:polynomial_degree), control_y(0:polynomial_degree)
    type(expr_t) :: x_coeff(0:polynomial_degree), y_coeff(0:polynomial_degree)
    type(expr_t) :: edge_area(1), edge_area_bar(1), edge_area_jvp(1)
    type(expr_t) :: edge_area_vjp(variable_count), curve_parameter
    type(expr_t) :: edge_x, edge_y, derivative_x, derivative_y, coefficient
    type(kernel_spec_t) :: spec
    character(:), allocatable :: primal_code, jvp_code, vjp_code
    integer :: index, degree, control, ios, unit

    call arena%init()
    engine = make_native_engine(arena)
    variables = [ &
        sym(arena, "x_start"), sym(arena, "y_start"), &
        sym(arena, "control_1_x"), sym(arena, "control_1_y"), &
        sym(arena, "control_2_x"), sym(arena, "control_2_y"), &
        sym(arena, "control_3_x"), sym(arena, "control_3_y"), &
        sym(arena, "control_4_x"), sym(arena, "control_4_y"), &
        sym(arena, "control_5_x"), sym(arena, "control_5_y"), &
        sym(arena, "control_6_x"), sym(arena, "control_6_y"), &
        sym(arena, "control_7_x"), sym(arena, "control_7_y"), &
        sym(arena, "control_8_x"), sym(arena, "control_8_y"), &
        sym(arena, "control_9_x"), sym(arena, "control_9_y"), &
        sym(arena, "control_10_x"), sym(arena, "control_10_y"), &
        sym(arena, "x_end"), sym(arena, "y_end")]
    do control = 0, polynomial_degree
        control_x(control) = variables(2*control + 1)
        control_y(control) = variables(2*control + 2)
    end do

    do degree = 0, polynomial_degree
        x_coeff(degree) = real_expr(arena, 0.0d0)
        y_coeff(degree) = real_expr(arena, 0.0d0)
        do control = 0, degree
            coefficient = real_expr(arena, real( &
                binomial(polynomial_degree, degree)*binomial(degree, control)* &
                (-1)**(degree - control), real64))
            x_coeff(degree) = x_coeff(degree) + coefficient*control_x(control)
            y_coeff(degree) = y_coeff(degree) + coefficient*control_y(control)
        end do
    end do

    curve_parameter = sym(arena, "parameter")
    edge_x = x_coeff(0)
    edge_y = y_coeff(0)
    do degree = 1, polynomial_degree
        edge_x = edge_x + x_coeff(degree)*curve_parameter**degree
        edge_y = edge_y + y_coeff(degree)*curve_parameter**degree
    end do
    derivative_x = diff(edge_x, curve_parameter)
    derivative_y = diff(edge_y, curve_parameter)
    edge_area(1) = real_expr(arena, 0.0d0)
    do index = 0, polynomial_degree
        do degree = 1, polynomial_degree
            edge_area(1) = edge_area(1) + real_expr(arena, 0.5d0)* &
                real_expr(arena, real(degree, real64))* &
                (x_coeff(index)*y_coeff(degree) - &
                y_coeff(index)*x_coeff(degree))/ &
                real_expr(arena, real(index + degree, real64))
        end do
    end do
    if (edge_x%node_count() < 1 .or. derivative_x%node_count() < 1 .or. &
        edge_y%node_count() < 1 .or. derivative_y%node_count() < 1) then
        error stop "hendecic Bezier generator produced an empty curve"
    end if
    call simplify_all(edge_area)

    call initialize_spec(spec, "generated_fci_hendecic_bezier_edge_area", &
        "fortfem_generated_fci_hendecic_bezier_edge_area", variable_count, 1)
    spec%args = primal_arguments()
    spec%outputs(1) = str("edge_area")
    primal_code = chars(emit_kernel(edge_area, spec))

    variables_dot = [ &
        sym(arena, "x_start_dot"), sym(arena, "y_start_dot"), &
        sym(arena, "control_1_x_dot"), sym(arena, "control_1_y_dot"), &
        sym(arena, "control_2_x_dot"), sym(arena, "control_2_y_dot"), &
        sym(arena, "control_3_x_dot"), sym(arena, "control_3_y_dot"), &
        sym(arena, "control_4_x_dot"), sym(arena, "control_4_y_dot"), &
        sym(arena, "control_5_x_dot"), sym(arena, "control_5_y_dot"), &
        sym(arena, "control_6_x_dot"), sym(arena, "control_6_y_dot"), &
        sym(arena, "control_7_x_dot"), sym(arena, "control_7_y_dot"), &
        sym(arena, "control_8_x_dot"), sym(arena, "control_8_y_dot"), &
        sym(arena, "control_9_x_dot"), sym(arena, "control_9_y_dot"), &
        sym(arena, "control_10_x_dot"), sym(arena, "control_10_y_dot"), &
        sym(arena, "x_end_dot"), sym(arena, "y_end_dot")]
    edge_area_bar(1) = sym(arena, "edge_area_bar")
    edge_area_jvp = jvp(edge_area, variables, variables_dot)
    edge_area_vjp = vjp(edge_area, variables, edge_area_bar)
    call simplify_all(edge_area_jvp)
    call simplify_all(edge_area_vjp)

    call initialize_spec(spec, "generated_fci_hendecic_bezier_edge_area_jvp", &
        "fortfem_generated_fci_hendecic_bezier_edge_area_jvp", &
        2*variable_count, 1)
    spec%args = jvp_arguments()
    spec%outputs(1) = str("edge_area_dot")
    jvp_code = chars(emit_kernel(edge_area_jvp, spec))

    call initialize_spec(spec, "generated_fci_hendecic_bezier_edge_area_vjp", &
        "fortfem_generated_fci_hendecic_bezier_edge_area_vjp", &
        variable_count + 1, variable_count)
    spec%args = vjp_arguments()
    spec%outputs = vjp_outputs()
    vjp_code = chars(emit_kernel(edge_area_vjp, spec))

    open (newunit=unit, file=generated_path( &
        "fortfem_fci_hendecic_bezier_edge_area_products.f90"), &
        status="replace", action="write", iostat=ios)
    if (ios /= 0) error stop "cannot write FCI hendecic Bezier edge products"
    write (unit, "(a)") primal_code(:len(primal_code) - 1)
    write (unit, "(a)") jvp_code(:len(jvp_code) - 1)
    write (unit, "(a)") vjp_code(:len(vjp_code) - 1)
    close (unit)

contains

    pure integer function binomial(n, k) result(value)
        integer, intent(in) :: n, k
        integer :: factor

        value = 1
        do factor = 1, k
            value = value*(n - factor + 1)/factor
        end do
    end function binomial

    function primal_arguments() result(arguments)
        use fortsym_string, only: str_t
        type(str_t), allocatable :: arguments(:)

        allocate(arguments(variable_count))
        arguments = [str("x_start"), str("y_start"), &
            str("control_1_x"), str("control_1_y"), &
            str("control_2_x"), str("control_2_y"), &
            str("control_3_x"), str("control_3_y"), &
            str("control_4_x"), str("control_4_y"), &
            str("control_5_x"), str("control_5_y"), &
            str("control_6_x"), str("control_6_y"), &
            str("control_7_x"), str("control_7_y"), &
            str("control_8_x"), str("control_8_y"), &
            str("control_9_x"), str("control_9_y"), &
            str("control_10_x"), str("control_10_y"), &
            str("x_end"), str("y_end")]
    end function primal_arguments

    function jvp_arguments() result(arguments)
        use fortsym_string, only: str_t
        type(str_t), allocatable :: arguments(:)

        allocate(arguments(2*variable_count))
        arguments(1:variable_count) = primal_arguments()
        arguments(variable_count + 1:) = [ &
            str("x_start_dot"), str("y_start_dot"), &
            str("control_1_x_dot"), str("control_1_y_dot"), &
            str("control_2_x_dot"), str("control_2_y_dot"), &
            str("control_3_x_dot"), str("control_3_y_dot"), &
            str("control_4_x_dot"), str("control_4_y_dot"), &
            str("control_5_x_dot"), str("control_5_y_dot"), &
            str("control_6_x_dot"), str("control_6_y_dot"), &
            str("control_7_x_dot"), str("control_7_y_dot"), &
            str("control_8_x_dot"), str("control_8_y_dot"), &
            str("control_9_x_dot"), str("control_9_y_dot"), &
            str("control_10_x_dot"), str("control_10_y_dot"), &
            str("x_end_dot"), str("y_end_dot")]
    end function jvp_arguments

    function vjp_arguments() result(arguments)
        use fortsym_string, only: str_t
        type(str_t), allocatable :: arguments(:)

        allocate(arguments(variable_count + 1))
        arguments(1:variable_count) = primal_arguments()
        arguments(variable_count + 1) = str("edge_area_bar")
    end function vjp_arguments

    function vjp_outputs() result(outputs)
        use fortsym_string, only: str_t
        type(str_t), allocatable :: outputs(:)

        allocate(outputs(variable_count))
        outputs = [str("x_start_bar"), str("y_start_bar"), &
            str("control_1_x_bar"), str("control_1_y_bar"), &
            str("control_2_x_bar"), str("control_2_y_bar"), &
            str("control_3_x_bar"), str("control_3_y_bar"), &
            str("control_4_x_bar"), str("control_4_y_bar"), &
            str("control_5_x_bar"), str("control_5_y_bar"), &
            str("control_6_x_bar"), str("control_6_y_bar"), &
            str("control_7_x_bar"), str("control_7_y_bar"), &
            str("control_8_x_bar"), str("control_8_y_bar"), &
            str("control_9_x_bar"), str("control_9_y_bar"), &
            str("control_10_x_bar"), str("control_10_y_bar"), &
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
        kernel_spec%generator = str("gen_fci_hendecic_bezier_edge_area_products")
        kernel_spec%generator_revision = str(fortsym_revision())
        kernel_spec%regenerate_command = str( &
            "cd tools/codegen && ./generate.sh")
        kernel_spec%pure_procedure = .true.
        allocate(kernel_spec%args(argument_count), &
            kernel_spec%outputs(output_count))
    end subroutine initialize_spec

    subroutine simplify_all(expressions)
        type(expr_t), intent(inout) :: expressions(:)

        do index = 1, size(expressions)
            result = engine%simplify(expressions(index))
            if (result%ok) expressions(index) = result%value
        end do
    end subroutine simplify_all

end program gen_fci_hendecic_bezier_edge_area_products
