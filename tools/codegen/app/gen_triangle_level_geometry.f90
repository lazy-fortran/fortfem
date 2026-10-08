program gen_triangle_level_geometry
    use fortsym_arena, only: arena_t
    use fortsym_expr, only: expr_t, sym, num, operator(+), operator(-), &
        operator(*), operator(/), sin, pi_expr
    use fortsym_diff, only: diff
    use fortsym_subs, only: subs
    use fortsym_engine, only: engine_result_t
    use fortsym_engine_native, only: native_engine_t, make_native_engine
    use fortsym_kernel, only: kernel_spec_t, emit_kernel, KERNEL_SUBROUTINE
    use fortsym_string, only: str_t, str, chars
    use fortfem_codegen_provenance, only: fortsym_revision, generated_path
    implicit none
    type(arena_t), target :: arena
    type(native_engine_t) :: engine
    type(expr_t) :: c(6), x, y, q, edges(6), discriminant(3), inner(3)
    type(expr_t) :: roots(9), slice(4), scalar(2), interval(2), dc, qb, qa
    type(expr_t) :: lower, upper, unit_node, unit_weight, zero, angle
    character(:), allocatable :: filename
    character(8) :: number
    integer :: i, unit, ios

    call arena%init()
    engine = make_native_engine(arena)
    x = sym(arena, 'xi')
    y = sym(arena, 'eta')
    zero = num(arena, 0)
    do i = 1, 6
        write(number, '(i0)') i
        c(i) = sym(arena, 'coefficients('//trim(number)//')')
    end do
    ! Native reference quadratic: constant, xi, eta, xi^2, xi*eta, eta^2.
    q = c(1) + c(2)*x + c(3)*y + c(4)*x*x + c(5)*x*y + c(6)*y*y
    inner = polynomial_coefficients(q, y)
    edges(1:3) = polynomial_coefficients(subs(q, y, zero), x)
    edges(4:6) = polynomial_coefficients(subs(q, y, 1 - x), x)
    dc = inner(1)
    qb = inner(2)
    qa = inner(3)
    discriminant = polynomial_coefficients(qb*qb - 4*qa*dc, x)
    roots = [edges, discriminant]

    filename = generated_path('fortfem_triangle_level_geometry.f90')
    open(newunit=unit, file=filename, status='replace', action='write', iostat=ios)
    if (ios /= 0) error stop 'cannot write triangle level geometry products'
    call write_product(roots, 'events', [str('coefficients')], [str('(6)')], &
        [str('edge_coefficients'), str('discriminant_coefficients')], &
        [str('(3,2)'), str('(3)')], &
        [str('edge_coefficients(1,1)'), str('edge_coefficients(2,1)'), &
        str('edge_coefficients(3,1)'), str('edge_coefficients(1,2)'), &
        str('edge_coefficients(2,2)'), str('edge_coefficients(3,2)'), &
        str('discriminant_coefficients(1)'), str('discriminant_coefficients(2)'), &
        str('discriminant_coefficients(3)')])
    slice = [inner, 1 - x]
    call write_product(slice, 'slice', [str('coefficients'), str('xi')], &
        [str('(6)'), str('')], [str('inner_coefficients'), str('upper_eta')], &
        [str('(3)'), str('')], &
        [str('inner_coefficients(1)'), str('inner_coefficients(2)'), &
        str('inner_coefficients(3)'), str('upper_eta')])

    q = c(1) + c(2)*x + c(3)*x*x
    scalar = [q, diff(q, x)]
    call write_product(scalar, 'polynomial', [str('coefficients'), str('xi')], &
        [str('(3)'), str('')], [str('value'), str('slope')], &
        [str(''), str('')], [str('value'), str('slope')])

    lower = sym(arena, 'lower')
    upper = sym(arena, 'upper')
    unit_node = sym(arena, 'unit_node')
    unit_weight = sym(arena, 'unit_weight')
    interval(1) = lower + (upper - lower)*unit_node
    interval(2) = diff(interval(1), unit_node)*unit_weight
    call write_product(interval, 'interval', &
        [str('lower'), str('upper'), str('unit_node'), str('unit_weight')], &
        [str(''), str(''), str(''), str('')], &
        [str('node'), str('weight')], [str(''), str('')], &
        [str('node'), str('weight')])
    ! Sine-squared coordinates regularize quadratic-root tangencies at either end.
    ! The positive measure is the native derivative of this exact node definition.
    angle = pi_expr(arena)*unit_node/2
    interval(1) = lower + (upper - lower)*sin(angle)*sin(angle)
    interval(2) = diff(interval(1), unit_node)*unit_weight
    call write_product(interval, 'tangent_interval', &
        [str('lower'), str('upper'), str('unit_node'), str('unit_weight')], &
        [str(''), str(''), str(''), str('')], &
        [str('node'), str('weight')], [str(''), str('')], &
        [str('node'), str('weight')])
    close(unit)
contains
    function polynomial_coefficients(expression, variable) result(coefficients)
        type(expr_t), intent(in) :: expression, variable
        type(expr_t) :: coefficients(3)
        coefficients(1) = subs(expression, variable, zero)
        coefficients(2) = subs(diff(expression, variable), variable, zero)
        coefficients(3) = subs(diff(diff(expression, variable), variable)/2, &
            variable, zero)
    end function polynomial_coefficients

    subroutine write_product(expressions, name, args, arg_shapes, outputs, &
            output_shapes, references)
        type(expr_t), intent(inout) :: expressions(:)
        character(*), intent(in) :: name
        type(str_t), intent(in) :: args(:), arg_shapes(:), outputs(:)
        type(str_t), intent(in) :: output_shapes(:), references(:)
        type(kernel_spec_t) :: spec
        type(engine_result_t) :: result
        character(:), allocatable :: code, message
        logical :: ok
        integer :: root
        do root = 1, size(expressions)
            result = engine%simplify(expressions(root))
            if (.not. result%ok) error stop 'triangle level simplification failed'
            expressions(root) = result%value
        end do
        spec%name = str('generated_triangle_level_'//name)
        spec%module_name = str('fortfem_generated_triangle_level_'//name)
        spec%mode = KERNEL_SUBROUTINE
        spec%generator = str('gen_triangle_level_geometry')
        spec%generator_revision = str(fortsym_revision())
        spec%regenerate_command = str('cd tools/codegen && ./generate.sh')
        spec%pure_procedure = .true.
        spec%args = args
        spec%arg_shapes = arg_shapes
        spec%outputs = outputs
        spec%output_shapes = output_shapes
        spec%output_references = references
        code = chars(emit_kernel(expressions, spec, ok=ok, message=message))
        if (.not. ok) then
            print *, message
            error stop 'triangle level product emission failed'
        end if
        write(unit, '(a)') code(:len(code) - 1)
    end subroutine write_product
end program gen_triangle_level_geometry
