program gen_harmonic_field_products
    use fortsym_arena, only: arena_t
    use fortsym_expr, only: expr_t, sym, operator(+), operator(-), &
        operator(*), operator(/)
    use fortsym_diff, only: diff
    use fortsym_kernel, only: kernel_spec_t, emit_kernel, KERNEL_SUBROUTINE
    use fortsym_products, only: jvp, vjp
    use fortsym_engine_native, only: native_engine_t, make_native_engine
    use fortsym_engine, only: engine_result_t, VERDICT_TRUE
    use fortsym_string, only: chars, str, str_t
    use fortfem_codegen_provenance, only: generated_path, fortsym_revision
    implicit none
    type(arena_t), target :: arena
    type(expr_t) :: v(7), product(1), x, y, rr, hx, hy
    type(str_t) :: names(7)
    type(native_engine_t) :: engine
    type(engine_result_t) :: identity
    integer :: i, unit
    call arena%init()
    engine = make_native_engine(arena)
    x = sym(arena, 'x')
    y = sym(arena, 'y')
    rr = x*x+y*y
    hx = -y/rr
    hy = x/rr
    ! Domain x*x+y*y>0. Check harmonicity of grad(azimuth).
    identity = engine%zero_test(diff(hx, x)+diff(hy, y))
    if (identity%verdict /= VERDICT_TRUE) error stop 'toroidal divergence'
    identity = engine%zero_test(diff(hy, x)-diff(hx, y))
    if (identity%verdict /= VERDICT_TRUE) error stop 'toroidal curl'
    identity = engine%zero_test(hx*hx+hy*hy-1/rr)
    if (identity%verdict /= VERDICT_TRUE) error stop 'toroidal energy density'
    names = [str('weight'), str('a1'), str('a2'), str('a3'), &
        str('b1'), str('b2'), str('b3')]
    do i = 1, 7
        v(i) = sym(arena, chars(names(i)))
    end do
    product = [v(1)*(v(2)*v(5)+v(3)*v(6)+v(4)*v(7))]
    open(newunit=unit, &
        file=generated_path('fortfem_harmonic_field_products.f90'), &
        status='replace', action='write')
    call emit('field_inner_product', names, v, product, [str('value')])
    close(unit)
contains
    subroutine emit(name, names, variables, roots, outputs)
        character(*), intent(in) :: name
        type(str_t), intent(in) :: names(:), outputs(:)
        type(expr_t), intent(in) :: variables(:), roots(:)
        type(kernel_spec_t) :: spec
        type(expr_t) :: directions(size(variables)), derivatives(size(roots))
        type(str_t) :: dot_names(size(names)), dot_outputs(size(outputs))
        character(:), allocatable :: code
        type(expr_t) :: cotangents(size(roots)), adjoints(size(variables))
        type(str_t) :: bar_names(size(outputs)), bar_outputs(size(names))
        integer :: k
        spec%name = str(name)
        spec%module_name = str('fortfem_generated_'//name)
        spec%generator = str('gen_harmonic_field_products')
        spec%generator_revision = str(fortsym_revision())
        spec%regenerate_command = str('cd tools/codegen && ./generate.sh')
        spec%mode = KERNEL_SUBROUTINE
        spec%pure_procedure = .true.
        spec%args = names
        spec%outputs = outputs
        code = chars(emit_kernel(roots, spec))
        write(unit, '(a)') code(:len(code)-1)
        do k = 1, size(names)
            dot_names(k) = str(chars(names(k))//'_dot')
            directions(k) = sym(arena, chars(dot_names(k)))
        end do
        do k = 1, size(outputs)
            dot_outputs(k) = str(chars(outputs(k))//'_dot')
        end do
        derivatives = jvp(roots, variables, directions)
        spec%name = str(name//'_jvp')
        spec%module_name = str('fortfem_generated_'//name//'_jvp')
        spec%args = [names, dot_names]
        spec%outputs = dot_outputs
        code = chars(emit_kernel(derivatives, spec))
        write(unit, '(a)') code(:len(code)-1)
        do k = 1, size(outputs)
            bar_names(k) = str(chars(outputs(k))//'_bar')
            cotangents(k) = sym(arena, chars(bar_names(k)))
        end do
        do k = 1, size(names)
            bar_outputs(k) = str(chars(names(k))//'_bar')
        end do
        adjoints = vjp(roots, variables, cotangents)
        spec%name = str(name//'_vjp')
        spec%module_name = str('fortfem_generated_'//name//'_vjp')
        spec%args = [names, bar_names]
        spec%outputs = bar_outputs
        code = chars(emit_kernel(adjoints, spec))
        write(unit, '(a)') code(:len(code)-1)
    end subroutine
end program
