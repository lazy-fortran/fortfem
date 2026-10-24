program gen_bem_high_order_products
    use fortsym_arena, only: arena_t
    use fortsym_expr, only: expr_t, sym, operator(+), operator(-), operator(*), &
        operator(/), operator(**), sqrt, asinh, sinh, cosh, log
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
    type(expr_t) :: normal_gap,radial_distance,radial_norm
    type(str_t) :: fan_names(9)
    type(str_t) :: names(7), pair_names(12)
    type(expr_t) :: pair(12), pair_map(5), pair_det
    type(str_t) :: angle_names(11)
    type(expr_t) :: angle(11), angle_map(2), first(3), delta_ray(3), cross_ray(3)
    type(expr_t) :: center, width, normsq, lower, upper, angle_t
    type(native_engine_t) :: engine
    type(engine_result_t) :: identity
    integer :: i, j, unit, degree, exponent
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
    ! Six sectors of the reference difference hexagon. Each sector has
    ! affine target/source anchors and a common residual triangle coordinate.
    ! The measure is r*(1-r)**2 dz dr ds; the weak 1/r singularity cancels.
    pair_names=[str('z1'),str('z2'),str('r'),str('s'),str('ta1'),str('ta2'), &
        str('tb1'),str('tb2'),str('sa1'),str('sa2'),str('sb1'),str('sb2')]
    do i=1,12
        pair(i)=sym(arena,chars(pair_names(i)))
    end do
    do j=1,2
        pair_map(j)=(1-pair(3))*pair(j)+pair(3)* &
            ((1-pair(4))*pair(j+4)+pair(4)*pair(j+6))
        pair_map(j+2)=(1-pair(3))*pair(j)+pair(3)* &
            ((1-pair(4))*pair(j+8)+pair(4)*pair(j+10))
    end do
    pair_map(5)=pair(3)*(1-pair(3))**2
    ! Difference-sector Jacobian for the canonical consecutive directions
    ! (1,0),(0,1); all six sectors have this determinant by construction.
    pair_det=diff(pair(3)*(1-pair(4)),pair(3))* &
        diff(pair(3)*pair(4),pair(4))- &
        diff(pair(3)*(1-pair(4)),pair(4))* &
        diff(pair(3)*pair(4),pair(3))
    call prove(pair_det*(1-pair(3))**2-pair_map(5), &
        'coincident difference-sector Jacobian')
    call emit('bem_coincident_map',pair_names,pair,pair_map, &
        [str('x1'),str('x2'),str('y1'),str('y2'),str('weight')], .false.)
    ! A tangent-metric sinh map removes the narrow angular peak of skinny
    ! panels: |first+s*delta| = |delta|*sqrt(width**2+(s-center)**2).
    angle_names=[str('ax'),str('ay'),str('az'),str('bx'),str('by'),str('bz'), &
        str('du1'),str('dv1'),str('du2'),str('dv2'),str('q')]
    do i=1,11
        angle(i)=sym(arena,chars(angle_names(i)))
    end do
    do j=1,3
        first(j)=angle(j)*angle(7)+angle(j+3)*angle(8)
        delta_ray(j)=angle(j)*(angle(9)-angle(7))+ &
            angle(j+3)*(angle(10)-angle(8))
    end do
    normsq=delta_ray(1)**2+delta_ray(2)**2+delta_ray(3)**2
    center=-(first(1)*delta_ray(1)+first(2)*delta_ray(2)+ &
        first(3)*delta_ray(3))/normsq
    cross_ray=[first(2)*delta_ray(3)-first(3)*delta_ray(2), &
        first(3)*delta_ray(1)-first(1)*delta_ray(3), &
        first(1)*delta_ray(2)-first(2)*delta_ray(1)]
    width=sqrt(cross_ray(1)**2+cross_ray(2)**2+cross_ray(3)**2)/normsq
    lower=asinh(-center/width); upper=asinh((1-center)/width)
    angle_t=(1-angle(11))*lower+angle(11)*upper
    angle_map=[center+width*sinh(angle_t),width*cosh(angle_t)*(upper-lower)]
    call prove(diff(angle_map(1),angle(11))-angle_map(2), &
        'metric angle Jacobian')
    call emit('bem_self_angle',angle_names,angle,angle_map, &
        [str('s'),str('weight')], .false.)
    ! For normal_gap>0 this primitive has endpoint difference one on
    ! [0,infinity); angular integration gives the double-layer jump 1/2.
    normal_gap=sym(arena,'normal_gap')
    radial_distance=sym(arena,'radial_distance')
    radial_norm=sqrt(normal_gap*normal_gap+radial_distance*radial_distance)
    call prove(diff(-normal_gap/radial_norm,radial_distance)- &
        normal_gap*radial_distance/(radial_norm**3), &
        'near-normal peak radial antiderivative')
    call emit_shared_edge_maps()
    close(unit)
contains
    subroutine emit_shared_edge_maps()
        ! Seven charts cover T x T up to zero-measure boundaries. Their
        ! positive Jacobians are parametric measures, never 1/distance.
        type(str_t) :: names(5), outputs(5)
        type(expr_t) :: vars(5), rho, e, z, q, h, a, b, xray, yray
        type(expr_t) :: eta, lo, hi, ray, limit, zz, jac, delta, t, u, ss, vv
        type(expr_t) :: mapped(5), angular, hh, abstract_x, abstract_y
        character(1) :: suffix
        integer :: chart, sector
        names=[str('rho'),str('edge'),str('inner'),str('angle'),str('aspect')]
        outputs=[str('t'),str('u'),str('s'),str('v'),str('weight')]
        do chart=1,5
            vars(chart)=sym(arena,chars(names(chart)))
        end do
        rho=vars(1); e=vars(2); z=vars(3); q=vars(4); h=vars(5)
        ! First prove the algebraic denominator identity on independent
        ! ray coordinates. On the open ray triangle, d=2-x-y>0 and
        ! a=x+y-1, hence (1-a)*(1-y)/d=1-y (or 1-x for N2) by field
        ! cancellation. Compose the proved identity with the hyperbolic map
        ! only afterwards, rather than expanding reciprocals of sinh/log.
        abstract_x=sym(arena,'ray_x'); abstract_y=sym(arena,'ray_y')
        call prove(1-(abstract_x+abstract_y-1)-(2-abstract_x-abstract_y), &
            'ray inverse denominator')
        do chart=1,7
            select case(chart)
            case(1)
                sector=1; a=q; b=z; angular=1+0*rho
            case(2,3)
                sector=3; eta=q*asinh(h); xray=sinh(eta)/h
                if (chart==2) then
                    a=(1-xray)*z
                else
                    xray=-xray; a=-xray+(1+xray)*z
                end if
                b=1-a-xray
                angular=(1-sinh(eta)/h)*cosh(eta)*asinh(h)/h
            case(4:7)
                sector=2
                if (chart>=6) sector=4
                hh=h*h+1
                lo=-asinh(1/h); hi=log(h)
                if (mod(chart,2)==1) then
                    lo=log(h); hi=asinh(h)
                end if
                eta=(1-q)*lo+q*hi
                zz=(1+h*sinh(eta))/hh
                limit=1/(1-zz)
                if (mod(chart,2)==1) limit=1/zz
                ray=1+(limit-1)*z
                xray=ray*zz; yray=ray*(1-zz)
                a=xray+yray-1
                b=(1-yray)/(2-xray-yray)
                if (sector==4) b=(1-xray)/(2-xray-yray)
                angular=ray*(limit-1)*h*cosh(eta)*(hi-lo)/hh
            end select
            delta=rho*(1-a); u=rho*a; vv=rho*b
            if (sector==2 .or. sector==4) then
                delta=rho*(1-a)*b; vv=rho
            end if
            t=(1-rho)*e+delta; ss=(1-rho)*e
            if (sector>=3) then
                ss=t; t=(1-rho)*e
                xray=u; u=vv; vv=xray
            end if
            jac=rho*rho*(1-rho)*angular
            mapped=[t,u,ss,vv,jac]
            write(suffix,'(i1)') chart
            call emit('bem_edge_chart_'//suffix,names,vars,mapped,outputs,.false.)
        end do
    end subroutine
    subroutine prove(expression,label)
        type(expr_t), intent(in) :: expression
        character(*), intent(in) :: label
        identity=engine%zero_test(expression)
        if (identity%verdict /= VERDICT_TRUE) then
            write(*,'(a)') label
            write(*,'(a,i0)') 'FortSym verdict=',identity%verdict
            error stop 'BEM trace identity failed'
        end if
    end subroutine
    subroutine emit(name, names, variables, roots, outputs, derivatives_needed)
        character(*), intent(in) :: name
        logical, intent(in), optional :: derivatives_needed
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
        if (present(derivatives_needed)) then
            if (.not. derivatives_needed) return
        end if
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
