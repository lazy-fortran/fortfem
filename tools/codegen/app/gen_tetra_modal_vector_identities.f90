program gen_tetra_modal_vector_identities
    use fortsym_arena, only: arena_t
    use fortsym_diff, only: diff
    use fortsym_engine, only: engine_result_t
    use fortsym_engine_native, only: make_native_engine, native_engine_t
    use fortsym_expr, only: expr_t, num, operator(*), operator(+), &
        operator(-), sym
    use fortsym_kernel, only: emit_kernel, kernel_spec_t, KERNEL_SUBROUTINE
    use fortsym_string, only: chars, str
    use fortfem_codegen_provenance, only: fortsym_revision, generated_path
    implicit none

    type(arena_t), target :: arena
    type(native_engine_t) :: engine
    type(expr_t) :: component_curls(3, 3), cross_curls(3, 3)
    type(expr_t) :: cross_curls_dot(3, 3), cross_values(3, 3)
    type(expr_t) :: cross_values_dot(3, 3), roots(27), roots_dot(18)
    type(expr_t) :: dx, dy, dz, phi, x, y, z, zero
    type(expr_t) :: dx_dot, dy_dot, dz_dot, phi_dot
    type(expr_t) :: x_dot, y_dot, z_dot
    type(kernel_spec_t) :: spec
    character(:), allocatable :: code, filename
    integer :: column, component, ios, root, unit

    call arena%init()
    engine = make_native_engine(arena)
    x = sym(arena, "x")
    y = sym(arena, "y")
    z = sym(arena, "z")
    phi = sym(arena, "phi")
    dx = sym(arena, "dx")
    dy = sym(arena, "dy")
    dz = sym(arena, "dz")
    x_dot = sym(arena, "x_dot")
    y_dot = sym(arena, "y_dot")
    z_dot = sym(arena, "z_dot")
    phi_dot = sym(arena, "phi_dot")
    dx_dot = sym(arena, "dx_dot")
    dy_dot = sym(arena, "dy_dot")
    dz_dot = sym(arena, "dz_dot")
    zero = num(arena, 0)

    component_curls = zero
    component_curls(:, 1) = [zero, dz, -dy]
    component_curls(:, 2) = [-dz, zero, dx]
    component_curls(:, 3) = [dy, -dx, zero]

    cross_values(:, 1) = [-y*phi, x*phi, zero]
    cross_values(:, 2) = [-z*phi, zero, x*phi]
    cross_values(:, 3) = [zero, -z*phi, y*phi]
    cross_curls(:, 1) = [ &
        -x*dz, -y*dz, 2*phi + x*dx + y*dy]
    cross_curls(:, 2) = [ &
        x*dy, -2*phi - x*dx - z*dz, z*dy]
    cross_curls(:, 3) = [ &
        2*phi + y*dy + z*dz, -y*dx, -z*dx]

    root = 0
    call append_matrix(component_curls)
    call append_matrix(cross_values)
    call append_matrix(cross_curls)

    spec%name = str("evaluate_tetra_modal_vector_identities")
    spec%module_name = str("fortfem_generated_tetra_modal_vector_identities")
    spec%mode = KERNEL_SUBROUTINE
    spec%temp_prefix = str("t")
    spec%generator = str("gen_tetra_modal_vector_identities")
    spec%generator_revision = str(fortsym_revision())
    spec%regenerate_command = str("cd tools/codegen && ./generate.sh")
    spec%pure_procedure = .true.
    allocate(spec%args(7), spec%outputs(3), spec%output_shapes(3))
    allocate(spec%output_references(27))
    spec%args = [ &
        str("x"), str("y"), str("z"), str("phi"), &
        str("dx"), str("dy"), str("dz")]
    spec%outputs = [ &
        str("component_curls"), str("cross_values"), str("cross_curls")]
    spec%output_shapes = [str("(3,3)"), str("(3,3)"), str("(3,3)")]
    root = 0
    call append_references("component_curls")
    call append_references("cross_values")
    call append_references("cross_curls")

    filename = generated_path( &
        "fortfem_tetra_modal_vector_identities.f90")
    open( &
        newunit=unit, file=filename, status="replace", action="write", &
        iostat=ios)
    if (ios /= 0) error stop "cannot write modal vector identity kernel"
    code = chars(emit_kernel(roots, spec))
    write(unit, "(a)") code(:len(code) - 1)
    close(unit)

    do column = 1, 3
        do component = 1, 3
            cross_values_dot(component, column) = &
                directional_derivative(cross_values(component, column))
            cross_curls_dot(component, column) = &
                directional_derivative(cross_curls(component, column))
        end do
    end do
    root = 0
    do column = 1, 3
        do component = 1, 3
            root = root + 1
            roots_dot(root) = cross_values_dot(component, column)
        end do
    end do
    do column = 1, 3
        do component = 1, 3
            root = root + 1
            roots_dot(root) = cross_curls_dot(component, column)
        end do
    end do
    call simplify_all(roots_dot)
    spec = kernel_spec_t()
    spec%name = str("evaluate_tetra_modal_vector_identities_jvp")
    spec%module_name = str( &
        "fortfem_generated_tetra_modal_vector_identities_jvp")
    spec%mode = KERNEL_SUBROUTINE
    spec%temp_prefix = str("t")
    spec%generator = str("gen_tetra_modal_vector_identities")
    spec%generator_revision = str(fortsym_revision())
    spec%regenerate_command = str("cd tools/codegen && ./generate.sh")
    spec%pure_procedure = .true.
    allocate(spec%args(14), spec%outputs(2), spec%output_shapes(2))
    allocate(spec%output_references(18))
    spec%args = [ &
        str("x"), str("y"), str("z"), str("phi"), &
        str("dx"), str("dy"), str("dz"), str("x_dot"), str("y_dot"), &
        str("z_dot"), str("phi_dot"), str("dx_dot"), str("dy_dot"), &
        str("dz_dot")]
    spec%outputs = [str("cross_values_dot"), str("cross_curls_dot")]
    spec%output_shapes = [str("(3,3)"), str("(3,3)")]
    root = 0
    call append_references("cross_values_dot")
    call append_references("cross_curls_dot")
    filename = generated_path( &
        "fortfem_tetra_modal_vector_identities_jvp.f90")
    open( &
        newunit=unit, file=filename, status="replace", action="write", &
        iostat=ios)
    if (ios /= 0) error stop "cannot write modal vector identity JVP kernel"
    code = chars(emit_kernel(roots_dot, spec))
    write(unit, "(a)") code(:len(code) - 1)
    close(unit)

contains

    subroutine append_matrix(matrix)
        type(expr_t), intent(in) :: matrix(3, 3)

        do column = 1, 3
            do component = 1, 3
                root = root + 1
                roots(root) = matrix(component, column)
            end do
        end do
    end subroutine append_matrix

    subroutine append_references(name)
        character(*), intent(in) :: name

        do column = 1, 3
            do component = 1, 3
                root = root + 1
                spec%output_references(root) = str( &
                    name//"("//integer_text(component)//","// &
                    integer_text(column)//")")
            end do
        end do
    end subroutine append_references

    function directional_derivative(expression) result(value)
        type(expr_t), intent(in) :: expression
        type(expr_t) :: value

        value = diff(expression, x)*x_dot + diff(expression, y)*y_dot + &
            diff(expression, z)*z_dot + diff(expression, phi)*phi_dot + &
            diff(expression, dx)*dx_dot + diff(expression, dy)*dy_dot + &
            diff(expression, dz)*dz_dot
    end function directional_derivative

    subroutine simplify_all(expressions)
        type(expr_t), intent(inout) :: expressions(:)

        type(engine_result_t) :: result
        integer :: expression

        do expression = 1, size(expressions)
            result = engine%simplify(expressions(expression))
            if (result%ok) expressions(expression) = result%value
        end do
    end subroutine simplify_all

    function integer_text(value) result(text)
        integer, intent(in) :: value
        character(:), allocatable :: text
        character(16) :: buffer

        write(buffer, "(i0)") value
        text = trim(buffer)
    end function integer_text

end program gen_tetra_modal_vector_identities
