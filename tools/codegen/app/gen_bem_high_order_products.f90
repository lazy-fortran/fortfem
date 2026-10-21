program gen_bem_high_order_products
    use fortsym_arena, only: arena_t
    use fortsym_expr, only: expr_t, sym, operator(+), operator(-), operator(*), &
        operator(/), operator(**)
    use fortsym_diff, only: diff
    use fortsym_kernel, only: kernel_spec_t, emit_kernel, KERNEL_SUBROUTINE
    use fortsym_products, only: jvp, vjp
    use fortsym_engine_native, only: native_engine_t, make_native_engine
    use fortsym_engine, only: engine_result_t, VERDICT_TRUE
    use fortsym_string, only: chars, str, str_t
    use fortfem_codegen_provenance, only: generated_path, fortsym_revision
    implicit none
    type(arena_t), target :: arena
    type(expr_t) :: v(7), product(1), x, y, z, phi(6), l(3), total
    type(expr_t) :: poly, lap, euler, rr, primitive, legendre
    type(expr_t) :: fan(9), mapped(3), signed_jacobian
    type(str_t) :: fan_names(9)
    type(str_t) :: names(7)
    type(native_engine_t) :: engine
    type(engine_result_t) :: identity
    integer :: i, unit, degree, exponent
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
    ! For harmonic homogeneous P_l, u=P_l/r**(2*l+1) is exterior harmonic.
    ! Product/chain rule: r**(2-s) Laplacian(r**s P) =
    ! r**2 Laplacian(P)+2*s*x.grad(P)+s*(s+1)*P. Domain r>0.
    z=sym(arena,'z'); rr=x*x+y*y+z*z
    do degree=1,2
        poly=z
        if (degree==2) poly=(2*z*z-x*x-y*y)/2
        lap=diff(diff(poly,x),x)+diff(diff(poly,y),y)+diff(diff(poly,z),z)
        euler=x*diff(poly,x)+y*diff(poly,y)+z*diff(poly,z)
        call prove(lap,'solid spherical harmonic Laplacian')
        call prove(euler-degree*poly,'solid harmonic Euler identity')
        exponent=-(2*degree+1)
        call prove(rr*lap+2*exponent*euler+exponent*(exponent+1)*poly, &
            'exterior spherical harmonic Laplacian')
        call prove(euler+exponent*poly+(degree+1)*poly, &
            'sphere exterior DtN eigenvalue')
        legendre=z; primitive=z**3/3
        if (degree==2) then
            legendre=(3*z*z-1)/2
            primitive=(9*z**5/5-2*z**3+z)/4
        end if
        call prove(diff(primitive,z)-legendre*legendre, &
            'sphere harmonic squared-norm antiderivative')
    end do
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
    fan_names=[str('apex1'),str('apex2'),str('e1u'),str('e1v'), &
        str('e2u'),str('e2v'),str('radius'),str('angle'),str('determinant')]
    do i=1,9
        fan(i)=sym(arena,chars(fan_names(i)))
    end do
    mapped(1)=fan(1)+fan(7)**2*((1-fan(8))*fan(3)+fan(8)*fan(5))
    mapped(2)=fan(2)+fan(7)**2*((1-fan(8))*fan(4)+fan(8)*fan(6))
    mapped(3)=2*fan(7)**3*fan(9)
    signed_jacobian=diff(mapped(1),fan(7))*diff(mapped(2),fan(8))- &
        diff(mapped(1),fan(8))*diff(mapped(2),fan(7))
    call prove(signed_jacobian-2*fan(7)**3*(fan(3)*fan(6)-fan(4)*fan(5)), &
        'squared radial fan Jacobian')
    call emit('bem_fan_map',fan_names,fan,mapped, &
        [str('uv1'),str('uv2'),str('weight')])
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
