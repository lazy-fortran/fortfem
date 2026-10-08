program gen_polynomial_candidate_jets
    use fortsym_arena,only:arena_t
    use fortsym_expr,only:expr_t,sym,num,operator(+),operator(-),operator(*),operator(/),operator(**)
    use fortsym_diff,only:diff
    use fortsym_subs,only:subs_many
    use fortsym_engine,only:engine_result_t
    use fortsym_engine_native,only:native_engine_t,make_native_engine
    use fortsym_kernel,only:kernel_spec_t,emit_kernel,KERNEL_SUBROUTINE,KERNEL_SNIPPET, CSE_NONE
    use fortsym_string,only:str,str_t,chars
    use fortfem_codegen_provenance,only:generated_path,fortsym_revision
    implicit none
    type(arena_t),target::arena
    type(native_engine_t)::engine
    integer::unit,dimension,order,kind
    call arena%init();engine=make_native_engine(arena)
    open(newunit=unit,file=generated_path('fortfem_polynomial_candidate_jets.f90'),status='replace',action='write')
    call power_jet()
    do dimension=2,3
        do order=1,2
            call monomial_jet(dimension,order)
        end do
    end do
    do kind=1,3
        do order=1,2
            call candidate_jet(kind,order)
        end do
    end do
    do kind=1,3
        do order=1,2
            call tetra_component_jet(kind,order)
        end do
    end do
    close(unit)
    call multiply_power()
contains
    function text(i)result(s)
        integer,intent(in)::i
        character(:),allocatable::s
        character(12)::buffer
        write(buffer,'(i0)')i;s=trim(buffer)
    end function
    function simplify(e)result(v)
        type(expr_t),intent(in)::e
        type(expr_t)::v
        type(engine_result_t)::r
        r=engine%simplify(e)
        if(.not.r%ok)error stop 'Native polynomial simplification failed'
        v=r%value
    end function
    subroutine initialize(spec,name)
        type(kernel_spec_t),intent(out)::spec
        character(*),intent(in)::name
        spec%name=str(name);spec%module_name=str('fortfem_generated_'//name)
        spec%mode=KERNEL_SUBROUTINE;spec%pure_procedure=.true.
        spec%generator=str('gen_polynomial_candidate_jets');spec%generator_revision=str(fortsym_revision())
        spec%regenerate_command=str('cd tools/codegen && ./generate.sh')
    end subroutine
    subroutine emit(r,spec)
        type(expr_t),intent(inout)::r(:)
        type(kernel_spec_t),intent(in)::spec
        character(:),allocatable::code
        character(:),allocatable::message
        logical::ok
        integer::i
        do i=1,size(r)
            r(i)=simplify(r(i))
        end do
        code=chars(emit_kernel(r,spec,ok=ok,message=message))
        if(.not.ok)then
            print *,message
            error stop 'Polynomial jet kernel emission failed'
        end if
        write(unit,'(a)')code(:len(code)-1)
    end subroutine
    subroutine emit_inline(roots,spec,filename)
        type(expr_t),intent(inout)::roots(:)
        type(kernel_spec_t),intent(inout)::spec
        character(*),intent(in)::filename
        integer::saved_unit
        saved_unit=unit
        spec%mode=KERNEL_SNIPPET
        spec%cse_level=CSE_NONE
        open(newunit=unit,file=generated_path(filename),status='replace',action='write')
        call emit(roots,spec)
        close(unit)
        unit=saved_unit
    end subroutine emit_inline
    subroutine multiply_power()
        type(expr_t)::r(1)
        type(kernel_spec_t)::spec
        call initialize(spec,'multiply_power')
        spec%mode=KERNEL_SNIPPET
        spec%cse_level=CSE_NONE
        spec%args=[str('previous_power'),str('coordinate')]
        spec%outputs=[str('power')]
        r(1)=sym(arena,'previous_power')*sym(arena,'coordinate')
        open(newunit=unit,file=generated_path('fortfem_power_multiply.inc'),status='replace',action='write')
        call emit(r,spec)
        close(unit)
    end subroutine
    subroutine power_jet()
        type(expr_t)::x,p,r(3),from(3),to(3)
        type(kernel_spec_t)::spec
        integer::i
        x=sym(arena,'x');p=sym(arena,'degree')
        r(1)=x**p;r(2)=diff(r(1),x);r(3)=diff(r(2),x)
        from=[x**p,x**(p-1),x**(p-2)]
        to=[sym(arena,'power'),sym(arena,'previous_power'),sym(arena,'second_previous_power')]
        do i=1,3
            from(i)=simplify(from(i))
        end do
        do i=1,3
            r(i)=subs_many(simplify(r(i)),from,to)
        end do
        call initialize(spec,'cached_power_jet')
        spec%args=[str('power'),str('previous_power'),str('second_previous_power'),str('degree')]
        spec%outputs=[str('value'),str('gradient'),str('hessian')]
        call emit(r,spec)
        call emit_inline(r,spec,'fortfem_cached_power_jet.inc')
    end subroutine
    subroutine monomial_jet(d,o)
        integer,intent(in)::d,o
        type(expr_t)::t(d),zero(d),p,f
        type(expr_t),allocatable::r(:)
        type(kernel_spec_t)::spec
        integer::i,j,k
        character(:),allocatable::shape
        do i=1,d
            t(i)=sym(arena,'t'//text(i));zero(i)=num(arena,0)
        end do
        p=num(arena,1)
        do i=1,d
            f=sym(arena,'factor_values('//text(i)//')')+sym(arena,'factor_gradients('//text(i)//')')*t(i)
            if(o==2)f=f+sym(arena,'factor_hessians('//text(i)//')')*t(i)*t(i)/2
            p=p*f
        end do
        k=1+d
        if(o==2)k=k+d*d
        allocate(r(k));r(1)=subs_many(p,t,zero)
        do i=1,d
            r(1+i)=subs_many(diff(p,t(i)),t,zero)
        end do
        if(o==2)then
            do j=1,d
                do i=1,d
                    r(1+d+(j-1)*d+i)=subs_many(diff(diff(p,t(i)),t(j)),t,zero)
                end do
            end do
        end if
        call initialize(spec,'monomial_jet'//text(d)//'_order'//text(o))
        shape='('//text(d)//')'
        spec%args=[str('factor_values'),str('factor_gradients')];spec%arg_shapes=[str(shape),str(shape)]
        spec%outputs=[str('value'),str('gradient')];spec%output_shapes=[str(''),str(shape)]
        if(o==2)then
            spec%args=[spec%args,str('factor_hessians')];spec%arg_shapes=[spec%arg_shapes,str(shape)]
            spec%outputs=[spec%outputs,str('hessian')];spec%output_shapes=[spec%output_shapes,str('('//text(d)//','//text(d)//')')]
        end if
        allocate(spec%output_references(k));spec%output_references(1)=str('value')
        do i=1,d
            spec%output_references(1+i)=str('gradient('//text(i)//')')
        end do
        if(o==2)then
            do j=1,d
                do i=1,d
                    spec%output_references(1+d+(j-1)*d+i)=str('hessian('//text(i)//','//text(j)//')')
                end do
            end do
        end if
        call emit(r,spec)
        call emit_inline(r,spec,'fortfem_monomial_jet'//text(d)//'_order'//text(o)//'.inc')
    end subroutine
    subroutine tetra_component_jet(component,o)
        integer,intent(in)::component,o
        type(expr_t)::t(3),zero(3),m,f(3),curl(3),r(3),direction(3),base(3)
        type(kernel_spec_t)::spec
        integer::i,j,k,other(2)
        do i=1,3
            t(i)=sym(arena,'t'//text(i))
            direction(i)=sym(arena,'point_dot('//text(i)//')')
        end do
        zero=num(arena,0)
        m=sym(arena,'value')
        do i=1,3
            m=m+sym(arena,'gradient('//text(i)//')')*t(i)
            if(o/=2)cycle
            m=m+sym(arena,'hessian('//text(i)//','//text(i)//')')*t(i)*t(i)/2
            do j=i+1,3
                m=m+sym(arena,'hessian('//text(i)//','//text(j)//')')*t(i)*t(j)
            end do
        end do
        f=num(arena,0);f(component)=sym(arena,'coefficient')*m
        curl=[diff(f(3),t(2))-diff(f(2),t(3)), &
            diff(f(1),t(3))-diff(f(3),t(1)),diff(f(2),t(1))-diff(f(1),t(2))]
        k=0
        do i=1,3
            if(i==component)cycle
            k=k+1;other(k)=i
        end do
        base=[f(component),curl(other(1)),curl(other(2))];r=base
        if(o==2)then
            do i=1,3
                r(i)=num(arena,0)
                do j=1,3
                    r(i)=r(i)+diff(base(i),t(j))*direction(j)
                end do
            end do
        end if
        do i=1,3
            r(i)=subs_many(r(i),t,zero)
        end do
        call initialize(spec,'tetra_component'//text(component)//'_order'//text(o))
        spec%args=[str('value'),str('gradient'),str('coefficient'), &
            str('value_accumulator'),str('curl_accumulators')]
        spec%arg_shapes=[str(''),str('(3)'),str(''),str(''),str('(3)')]
        if(o==2)then
            spec%args=[spec%args,str('hessian'),str('point_dot')]
            spec%arg_shapes=[spec%arg_shapes,str('(3,3)'),str('(3)')]
        end if
        if(o==1)then
            spec%outputs=[str('values'),str('curls')]
        else
            spec%outputs=[str('values_dot'),str('curls_dot')]
        end if
        spec%output_shapes=[str('(3,:)'),str('(3,:)')]
        allocate(spec%output_references(3))
        spec%output_references(1)=str(chars(spec%outputs(1))//'('//text(component)//',candidate)')
        r(1)=r(1)+sym(arena,'value_accumulator')
        do i=1,2
            spec%output_references(i+1)=str(chars(spec%outputs(2))//'('//text(other(i))//',candidate)')
            r(i+1)=r(i+1)+sym(arena,'curl_accumulators('//text(other(i))//')')
        end do
        call emit_inline(r,spec,'fortfem_tetra_component'//text(component)//'_order'//text(o)//'.inc')
    end subroutine tetra_component_jet
    subroutine candidate_jet(kind,o)
        integer,intent(in)::kind,o
        type(expr_t)::t(2),zero(2),m,f(2),curl,r(3),base(3),direction(2)
        type(kernel_spec_t)::spec
        integer::i
        t=[sym(arena,'tx'),sym(arena,'ty')];zero=num(arena,0)
        m=sym(arena,'value')+sym(arena,'gradient(1)')*t(1)+sym(arena,'gradient(2)')*t(2)
        if(o==2)m=m+sym(arena,'hessian(1,1)')*t(1)*t(1)/2+ &
            sym(arena,'hessian(1,2)')*t(1)*t(2)+sym(arena,'hessian(2,2)')*t(2)*t(2)/2
        select case(kind)
        case(1)
            f=[m,num(arena,0)]
        case(2)
            f=[num(arena,0),m]
        case(3)
            f=[-(sym(arena,'point(2)')+t(2))*m,(sym(arena,'point(1)')+t(1))*m]
        end select
        curl=diff(f(2),t(1))-diff(f(1),t(2));base=[f,curl];r=base
        if(o==2)then
            direction=[sym(arena,'direction(1)'),sym(arena,'direction(2)')]
            do i=1,3
                r(i)=diff(base(i),t(1))*direction(1)+diff(base(i),t(2))*direction(2)
            end do
        end if
        do i=1,3
            r(i)=subs_many(r(i),t,zero)
        end do
        call initialize(spec,'triangle_candidate'//text(kind)//'_order'//text(o))
        spec%args=[str('point'),str('value'),str('gradient')];spec%arg_shapes=[str('(2)'),str(''),str('(2)')]
        if(o==2)then
            spec%args=[spec%args,str('hessian'),str('direction')];spec%arg_shapes=[spec%arg_shapes,str('(2,2)'),str('(2)')]
        end if
        spec%outputs=[str('field'),str('curl')];spec%output_shapes=[str('(2)'),str('')]
        spec%output_references=[str('field(1)'),str('field(2)'),str('curl')]
        call emit(r,spec)
        do i=1,3
            r(i)=subs_many(r(i),[sym(arena,'point(1)'),sym(arena,'point(2)'), &
                sym(arena,'value')],[sym(arena,'xi'),sym(arena,'eta'),sym(arena,'monomial')])
            if(o==2)r(i)=subs_many(r(i),[sym(arena,'direction(1)'),sym(arena,'direction(2)')], &
                [sym(arena,'xi_dot'),sym(arena,'eta_dot')])
        end do
        spec%args=[str('xi'),str('eta'),str('monomial'),str('gradient')]
        spec%arg_shapes=[str(''),str(''),str(''),str('(2)')]
        if(o==2)then
            spec%args=[spec%args,str('hessian'),str('xi_dot'),str('eta_dot')]
            spec%arg_shapes=[spec%arg_shapes,str('(2,2)'),str(''),str('')]
            spec%outputs=[str('values_dot'),str('curls_dot')]
            spec%output_references=[str('values_dot(1,candidate)'),str('values_dot(2,candidate)'), &
                str('curls_dot(candidate)')]
        else
            spec%outputs=[str('values'),str('curls')]
            spec%output_references=[str('values(1,candidate)'),str('values(2,candidate)'), &
                str('curls(candidate)')]
        end if
        spec%output_shapes=[str('(2,:)'),str('(:)')]
        call emit_inline(r,spec,'fortfem_triangle_candidate'//text(kind)//'_order'//text(o)//'.inc')
    end subroutine
end program
