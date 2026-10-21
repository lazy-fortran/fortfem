program gen_bem_high_order_products
    use fortsym_arena, only: arena_t
    use fortsym_expr, only: expr_t, sym, operator(+), operator(-), operator(*)
    use fortsym_diff, only: diff
    use fortsym_kernel, only: kernel_spec_t, emit_kernel, KERNEL_SUBROUTINE
    use fortsym_products, only: jvp, vjp
    use fortsym_engine_native, only: native_engine_t, make_native_engine
    use fortsym_engine, only: engine_result_t, VERDICT_TRUE
    use fortsym_string, only: chars, str, str_t
    use fortfem_codegen_provenance, only: generated_path, fortsym_revision
    implicit none
    type(arena_t), target :: arena
    type(expr_t) :: v(7), product(1), x, y, phi(6), l(3), total
    type(str_t) :: names(7)
    type(native_engine_t) :: engine
    type(engine_result_t) :: identity
    integer :: i, unit
    call arena%init()
    engine = make_native_engine(arena)
    x = sym(arena, 'x'); y = sym(arena, 'y')
    l = [1-x-y,x,y]
    phi = [l(1)*(2*l(1)-1),l(2)*(2*l(2)-1),l(3)*(2*l(3)-1), &
        4*l(1)*l(2),4*l(2)*l(3),4*l(3)*l(1)]
    total=phi(1)
    do i=2,6
        total=total+phi(i)
    end do
    call prove(total-1,'P2 partition of unity')
    call prove(l(1)+l(2)+l(3)-1,'P1 partition of unity')
    call prove(diff(total,x),'P2 constant xi derivative')
    call prove(diff(total,y),'P2 constant eta derivative')
    names = [str('kernel'),str('jx'),str('jy'),str('test'),str('trial'), &
        str('unused1'),str('unused2')]
    do i=1,7
        v(i)=sym(arena,chars(names(i)))
    end do
    open(newunit=unit,file=generated_path('fortfem_bem_high_order_products.f90'), &
        status='replace',action='write')
    product=[v(1)*v(2)*v(3)*v(4)*v(5)]
    call emit('bem_trace_product',names(:5),v(:5),product,[str('value')])
    names=[str('green'),str('tx'),str('ty'),str('tz'), &
        str('sx'),str('sy'),str('sz')]
    do i=1,7
        v(i)=sym(arena,chars(names(i)))
    end do
    product=[v(1)*(v(2)*v(5)+v(3)*v(6)+v(4)*v(7))]
    call emit('bem_hyper_product',names,v,product,[str('value')])
    close(unit)
contains
    subroutine prove(expression,label)
        type(expr_t), intent(in) :: expression
        character(*), intent(in) :: label
        identity=engine%zero_test(expression)
        if (identity%verdict /= VERDICT_TRUE) then
            write(*,'(a)') label
            error stop 'BEM trace identity failed'
        end if
    end subroutine
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
        spec%generator = str('gen_bem_high_order_products')
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
