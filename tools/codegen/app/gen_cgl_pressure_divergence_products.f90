program gen_cgl_pressure_divergence_products
    use fortsym_arena, only: arena_t
    use fortsym_engine, only: engine_result_t
    use fortsym_engine_native, only: make_native_engine, native_engine_t
    use fortsym_expr, only: expr_t, operator(*), operator(+), operator(-), sym
    use fortsym_kernel, only: emit_kernel, kernel_spec_t, KERNEL_SUBROUTINE
    use fortsym_products, only: jvp, vjp
    use fortsym_string, only: chars, str
    use fortfem_codegen_provenance, only: fortsym_revision, generated_path
    implicit none

    type(arena_t), target :: arena
    type(engine_result_t) :: result
    type(native_engine_t) :: engine
    type(expr_t) :: variables(20), variables_dot(20), force(3)
    type(expr_t) :: force_jvp(3), force_vjp(20), force_bar(3)
    type(kernel_spec_t) :: spec
    character(:), allocatable :: primal_code, jvp_code, vjp_code
    integer :: i, j, index, ios, root, unit
    type(expr_t) :: delta_pressure, pressure_gradient

    call arena%init()
    engine = make_native_engine(arena)
    variables(1) = sym(arena, "p_parallel")
    variables(2) = sym(arena, "p_perpendicular")
    variables(3:5) = [ &
        sym(arena, "direction_1"), sym(arena, "direction_2"), &
        sym(arena, "direction_3")]
    variables(6:8) = [ &
        sym(arena, "parallel_gradient_1"), sym(arena, "parallel_gradient_2"), &
        sym(arena, "parallel_gradient_3")]
    variables(9:11) = [ &
        sym(arena, "perpendicular_gradient_1"), &
        sym(arena, "perpendicular_gradient_2"), &
        sym(arena, "perpendicular_gradient_3")]
    index = 12
    do i = 1, 3
        do j = 1, 3
            variables(index) = sym(arena, &
                "direction_gradient_"//itoa(i)//itoa(j))
            index = index + 1
        end do
    end do

    variables_dot(1) = sym(arena, "p_parallel_dot")
    variables_dot(2) = sym(arena, "p_perpendicular_dot")
    variables_dot(3:5) = [ &
        sym(arena, "direction_1_dot"), sym(arena, "direction_2_dot"), &
        sym(arena, "direction_3_dot")]
    variables_dot(6:8) = [ &
        sym(arena, "parallel_gradient_1_dot"), &
        sym(arena, "parallel_gradient_2_dot"), &
        sym(arena, "parallel_gradient_3_dot")]
    variables_dot(9:11) = [ &
        sym(arena, "perpendicular_gradient_1_dot"), &
        sym(arena, "perpendicular_gradient_2_dot"), &
        sym(arena, "perpendicular_gradient_3_dot")]
    index = 12
    do i = 1, 3
        do j = 1, 3
            variables_dot(index) = sym(arena, &
                "direction_gradient_"//itoa(i)//itoa(j)//"_dot")
            index = index + 1
        end do
    end do

    force_bar = [ &
        sym(arena, "force_1_bar"), sym(arena, "force_2_bar"), &
        sym(arena, "force_3_bar")]
    delta_pressure = variables(1) - variables(2)
    do i = 1, 3
        force(i) = variables(8 + i)
        do j = 1, 3
            pressure_gradient = variables(5 + j) - variables(8 + j)
            index = 11 + 3*(i - 1) + j
            force(i) = force(i) + variables(2 + i)*variables(2 + j)* &
                pressure_gradient + delta_pressure*( &
                variables(index)*variables(2 + j) + &
                variables(2 + i)*variables(11 + 3*(j - 1) + j))
        end do
    end do
    force_jvp = jvp(force, variables, variables_dot)
    force_vjp = vjp(force, variables, force_bar)
    call simplify_all(force)
    call simplify_all(force_jvp)
    call simplify_all(force_vjp)

    call initialize_spec( &
        spec, "generated_cgl_pressure_divergence", &
        "fortfem_generated_cgl_pressure_divergence", 20, 3)
    call set_primal_names(spec)
    spec%outputs = [str("force_1"), str("force_2"), str("force_3")]
    primal_code = chars(emit_kernel(force, spec))

    call initialize_spec( &
        spec, "generated_cgl_pressure_divergence_jvp", &
        "fortfem_generated_cgl_pressure_divergence_jvp", 40, 3)
    call set_primal_names(spec)
    call set_tangent_names(spec)
    spec%outputs = [ &
        str("force_1_dot"), str("force_2_dot"), str("force_3_dot")]
    jvp_code = chars(emit_kernel(force_jvp, spec))

    call initialize_spec( &
        spec, "generated_cgl_pressure_divergence_vjp", &
        "fortfem_generated_cgl_pressure_divergence_vjp", 23, 20)
    call set_primal_names(spec)
    spec%args(21:23) = [ &
        str("force_1_bar"), str("force_2_bar"), str("force_3_bar")]
    call set_vjp_names(spec)
    vjp_code = chars(emit_kernel(force_vjp, spec))

    open (newunit=unit, &
        file=generated_path("fortfem_cgl_pressure_divergence_products.f90"), &
        status="replace", action="write", iostat=ios)
    if (ios /= 0) error stop "cannot write CGL pressure divergence products"
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
        kernel_spec%generator = str("gen_cgl_pressure_divergence_products")
        kernel_spec%generator_revision = str(fortsym_revision())
        kernel_spec%regenerate_command = str( &
            "cd tools/codegen && ./generate.sh")
        kernel_spec%pure_procedure = .true.
        allocate ( &
            kernel_spec%args(argument_count), &
            kernel_spec%outputs(output_count))
    end subroutine initialize_spec

    subroutine set_primal_names(kernel_spec)
        type(kernel_spec_t), intent(inout) :: kernel_spec

        kernel_spec%args(1:20) = [ &
            str("p_parallel"), str("p_perpendicular"), str("direction_1"), &
            str("direction_2"), str("direction_3"), &
            str("parallel_gradient_1"), str("parallel_gradient_2"), &
            str("parallel_gradient_3"), str("perpendicular_gradient_1"), &
            str("perpendicular_gradient_2"), str("perpendicular_gradient_3"), &
            str("direction_gradient_11"), str("direction_gradient_12"), &
            str("direction_gradient_13"), str("direction_gradient_21"), &
            str("direction_gradient_22"), str("direction_gradient_23"), &
            str("direction_gradient_31"), str("direction_gradient_32"), &
            str("direction_gradient_33")]
    end subroutine set_primal_names

    subroutine set_tangent_names(kernel_spec)
        type(kernel_spec_t), intent(inout) :: kernel_spec

        kernel_spec%args(21:40) = [ &
            str("p_parallel_dot"), str("p_perpendicular_dot"), &
            str("direction_1_dot"), str("direction_2_dot"), &
            str("direction_3_dot"), str("parallel_gradient_1_dot"), &
            str("parallel_gradient_2_dot"), str("parallel_gradient_3_dot"), &
            str("perpendicular_gradient_1_dot"), &
            str("perpendicular_gradient_2_dot"), &
            str("perpendicular_gradient_3_dot"), str("direction_gradient_11_dot"), &
            str("direction_gradient_12_dot"), str("direction_gradient_13_dot"), &
            str("direction_gradient_21_dot"), str("direction_gradient_22_dot"), &
            str("direction_gradient_23_dot"), str("direction_gradient_31_dot"), &
            str("direction_gradient_32_dot"), str("direction_gradient_33_dot")]
    end subroutine set_tangent_names

    subroutine set_vjp_names(kernel_spec)
        type(kernel_spec_t), intent(inout) :: kernel_spec

        kernel_spec%outputs = [ &
            str("p_parallel_bar"), str("p_perpendicular_bar"), &
            str("direction_1_bar"), str("direction_2_bar"), &
            str("direction_3_bar"), str("parallel_gradient_1_bar"), &
            str("parallel_gradient_2_bar"), str("parallel_gradient_3_bar"), &
            str("perpendicular_gradient_1_bar"), &
            str("perpendicular_gradient_2_bar"), &
            str("perpendicular_gradient_3_bar"), str("direction_gradient_11_bar"), &
            str("direction_gradient_12_bar"), str("direction_gradient_13_bar"), &
            str("direction_gradient_21_bar"), str("direction_gradient_22_bar"), &
            str("direction_gradient_23_bar"), str("direction_gradient_31_bar"), &
            str("direction_gradient_32_bar"), str("direction_gradient_33_bar")]
    end subroutine set_vjp_names

    function itoa(value) result(text)
        integer, intent(in) :: value
        character(:), allocatable :: text

        character(2) :: buffer
        write (buffer, "(i0)") value
        text = trim(buffer)
    end function itoa

    subroutine simplify_all(expressions)
        type(expr_t), intent(inout) :: expressions(:)

        do root = 1, size(expressions)
            result = engine%simplify(expressions(root))
            if (result%ok) expressions(root) = result%value
        end do
    end subroutine simplify_all

end program gen_cgl_pressure_divergence_products
