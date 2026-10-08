program gen_reference_basis_products
    use fortsym_arena, only: arena_t
    use fortsym_engine, only: engine_result_t, VERDICT_TRUE
    use fortsym_engine_native, only: make_native_engine, native_engine_t
    use fortsym_expr, only: expr_t, sym, num, operator(+), operator(-), &
        operator(*), operator(/), operator(**)
    use fortsym_diff, only: diff
    use fortsym_subs, only: subs_many
    use fortsym_kernel, only: emit_kernel, kernel_spec_t, KERNEL_SUBROUTINE
    use fortsym_rigorous_emit, only: rigorous_kernel_spec_t, interval_runtime, &
        emit_rigorous_kernel
    use fortsym_string, only: chars, str
    use fortfem_codegen_provenance, only: fortsym_revision, generated_path
    implicit none
    type(arena_t), target :: arena
    type(native_engine_t) :: engine
    type(engine_result_t) :: result
    type(expr_t) :: x, y, l(3), p1(3), p2(6), q1(4), line(2)
    type(expr_t) :: edge(2,3), rt(2,3), roots(6)
    type(kernel_spec_t) :: spec
    character(:), allocatable :: code, name
    character(4) :: index_string
    integer :: unit, i, j

    call arena%init()
    engine = make_native_engine(arena)
    x = sym(arena, 'xi'); y = sym(arena, 'eta')
    l(1)=1-x-y; l(2)=x; l(3)=y
    do i=1,3
        p1(i)=l(i)
    end do
    do i=1,3
        p2(i)=l(i)*(2*l(i)-1)
    end do
    p2(4)=4*l(1)*l(2); p2(5)=4*l(2)*l(3); p2(6)=4*l(3)*l(1)
    q1(1)=(1-x)*(1-y); q1(2)=(1+x)*(1-y)
    q1(3)=(1+x)*(1+y); q1(4)=(1-x)*(1+y)
    ! Exact rational normalization on the square reference element.
    do i=1,4
        q1(i)=q1(i)/4
    end do
    line(1)=1-x; line(2)=x
    edge(1,1)=1-y; edge(2,1)=x
    edge(1,2)=-y; edge(2,2)=x
    edge(1,3)=-y; edge(2,3)=x-1
    rt(1,1)=x; rt(2,1)=y-1
    rt(1,2)=x; rt(2,2)=y
    rt(1,3)=x-1; rt(2,3)=y
    open(newunit=unit, file=generated_path('fortfem_reference_basis_products.f90'), &
        status='replace', action='write')
    call scalar_family('p1',p1)
    call scalar_family('p2',p2)
    call scalar_family('q1',q1)
    call scalar_family('line',line)
    call vector_family('edge',edge)
    call vector_family('rt',rt)
    call affine_geometry()
    call nodal_scalar_coefficients("p1", p1)
    call nodal_scalar_coefficients("p2", p2)
    close(unit)
    call interval_affine_p2()
contains
    subroutine interval_affine_p2()
        type(expr_t) :: r, z, vertex_r(3), vertex_z(3), nodes(6), field
        type(expr_t) :: determinant, xi, eta, mapped, roots_interval(4)
        type(expr_t) :: locations_x(6), locations_y(6), oracle, recovered
        type(expr_t) :: centered(5), differences(5), coefficients(5), symbols(5)
        type(expr_t) :: polynomial, origin(2), coefficient_roots(5)
        type(rigorous_kernel_spec_t) :: interval_spec
        character(1) :: digit
        integer :: node, degree_x, degree_y

        locations_x = [num(arena,0), num(arena,1), num(arena,0), &
            num(arena,1)/2, num(arena,1)/2, num(arena,0)]
        locations_y = [num(arena,0), num(arena,0), num(arena,1), &
            num(arena,0), num(arena,1)/2, num(arena,1)/2]
        do degree_x = 0, 2
            do degree_y = 0, 2-degree_x
                oracle = x**degree_x*y**degree_y
                recovered = num(arena,0)
                do node = 1, 6
                    recovered = recovered + p2(node)*subs_many(oracle, &
                        [x,y], [locations_x(node),locations_y(node)])
                end do
                result = engine%zero_test(recovered-oracle)
                if (result%verdict /= VERDICT_TRUE) &
                    error stop 'P2 polynomial reproduction unproved'
            end do
        end do
        r = sym(arena,'r'); z = sym(arena,'z')
        do node = 1, 3
            write(digit,'(i0)') node
            vertex_r(node) = sym(arena,'r'//digit)
            vertex_z(node) = sym(arena,'z'//digit)
        end do
        do node = 1, 6
            write(digit,'(i0)') node
            nodes(node) = sym(arena,'u'//digit)
        end do
        determinant = (vertex_r(2)-vertex_r(1))*(vertex_z(3)-vertex_z(1)) &
            -(vertex_r(3)-vertex_r(1))*(vertex_z(2)-vertex_z(1))
        xi = ((r-vertex_r(1))*(vertex_z(3)-vertex_z(1)) &
            -(z-vertex_z(1))*(vertex_r(3)-vertex_r(1)))/determinant
        eta = ((vertex_r(2)-vertex_r(1))*(z-vertex_z(1)) &
            -(vertex_z(2)-vertex_z(1))*(r-vertex_r(1)))/determinant
        result = engine%zero_test(vertex_r(1)+(vertex_r(2)-vertex_r(1))*xi &
            +(vertex_r(3)-vertex_r(1))*eta-r)
        if (result%verdict /= VERDICT_TRUE) error stop 'affine R inverse unproved'
        result = engine%zero_test(vertex_z(1)+(vertex_z(2)-vertex_z(1))*xi &
            +(vertex_z(3)-vertex_z(1))*eta-z)
        if (result%verdict /= VERDICT_TRUE) error stop 'affine Z inverse unproved'
        field = num(arena,0)
        do node = 1, 5
            write(digit,'(i0)') node+1
            centered(node) = sym(arena,'d'//digit)
            differences(node) = nodes(node+1)-nodes(1)
            field = field + centered(node)*p2(node+1)
        end do
        origin = num(arena,0)
        coefficients(1) = subs_many(diff(field,x),[x,y],origin)
        coefficients(2) = subs_many(diff(field,y),[x,y],origin)
        coefficients(3) = subs_many(diff(diff(field,x),x)/2,[x,y],origin)
        coefficients(4) = subs_many(diff(diff(field,x),y),[x,y],origin)
        coefficients(5) = subs_many(diff(diff(field,y),y)/2,[x,y],origin)
        do node = 1, 5
            result = engine%simplify(coefficients(node))
            if (.not. result%ok) error stop 'P2 interval coefficient derivation failed'
            coefficients(node) = result%value
        end do
        polynomial = coefficients(1)*x+coefficients(2)*y &
            +coefficients(3)*x**2+coefficients(4)*x*y+coefficients(5)*y**2
        result = engine%zero_test(field-polynomial)
        if (result%verdict /= VERDICT_TRUE) &
            error stop 'P2 coefficient identity unproved'
        result = engine%zero_test(diff(field,x)-diff(polynomial,x))
        if (result%verdict /= VERDICT_TRUE) error stop 'P2 xi derivative unproved'
        result = engine%zero_test(diff(field,y)-diff(polynomial,y))
        if (result%verdict /= VERDICT_TRUE) error stop 'P2 eta derivative unproved'
        recovered = subs_many(polynomial,[x,y],[xi,eta])
        oracle = subs_many(field,[x,y],[xi,eta])
        result = engine%zero_test(diff(oracle,r)-diff(recovered,r))
        if (result%verdict /= VERDICT_TRUE) &
            error stop 'P2 physical R derivative unproved'
        result = engine%zero_test(diff(oracle,z)-diff(recovered,z))
        if (result%verdict /= VERDICT_TRUE) &
            error stop 'P2 physical Z derivative unproved'
        do node = 1, 5
            coefficient_roots(node) = subs_many(coefficients(node),centered,differences)
        end do
        open(newunit=unit, &
            file=generated_path('fortfem_reference_scalar_intervals.f90'), &
            status='replace',action='write')
        write(unit,'(a)') '! Generated from canonical reference P2 basis; do not edit.'
        write(unit,'(a)') 'module fortfem_reference_scalar_intervals'
        write(unit,'(a)') 'use fortnum_interval, only: interval_t'
        write(unit,'(a)') 'implicit none'
        write(unit,'(a)') 'private'
        write(unit,'(a)') 'public :: generated_affine_p2_scalar_interval'
        write(unit,'(a)') 'contains'
        interval_spec%name = str('generated_p2_centered_coefficients_interval')
        interval_spec%generator = str( &
            'gen_reference_basis_products; '//fortsym_revision())
        interval_spec%runtime = interval_runtime('fortnum_interval','interval_t')
        interval_spec%horner = .false.
        interval_spec%args = [str('u1'),str('u2'),str('u3'),str('u4'), &
            str('u5'),str('u6')]
        interval_spec%outputs = [str('lx'),str('ly'),str('a'),str('b'),str('c')]
        call emit_interval(coefficient_roots,interval_spec)
        symbols = [sym(arena,'lx'),sym(arena,'ly'),sym(arena,'a'), &
            sym(arena,'b'),sym(arena,'c')]
        mapped = nodes(1)+symbols(1)*xi+symbols(2)*eta+symbols(3)*xi**2 &
            +symbols(4)*xi*eta+symbols(5)*eta**2
        roots_interval = [mapped,diff(mapped,r),diff(mapped,z),determinant]
        interval_spec%name = str('generated_affine_p2_polynomial_interval')
        interval_spec%horner = .true.
        interval_spec%args = [str('r'),str('z'),str('r1'),str('z1'),str('r2'), &
            str('z2'),str('r3'),str('z3'),str('u1'),str('lx'),str('ly'), &
            str('a'),str('b'),str('c')]
        interval_spec%outputs = [str('value'),str('grad_r'),str('grad_z'),str('det')]
        call emit_interval(roots_interval,interval_spec)
        write(unit,'(a)') 'pure subroutine generated_affine_p2_scalar_interval( &'
        write(unit,'(a)') &
            'r,z,r1,z1,r2,z2,r3,z3,u1,u2,u3,u4,u5,u6,value,grad_r,grad_z,det)'
        write(unit,'(a)') 'type(interval_t), intent(in) :: r,z,r1,z1,r2,z2,r3,z3'
        write(unit,'(a)') 'type(interval_t), intent(in) :: u1,u2,u3,u4,u5,u6'
        write(unit,'(a)') 'type(interval_t), intent(out) :: value,grad_r,grad_z,det'
        write(unit,'(a)') 'type(interval_t) :: lx,ly,a,b,c'
        write(unit,'(a)') 'call generated_p2_centered_coefficients_interval( &'
        write(unit,'(a)') 'u1,u2,u3,u4,u5,u6,lx,ly,a,b,c)'
        write(unit,'(a)') 'call generated_affine_p2_polynomial_interval( &'
        write(unit,'(a)') &
            'r,z,r1,z1,r2,z2,r3,z3,u1,lx,ly,a,b,c,value,grad_r,grad_z,det)'
        write(unit,'(a)') 'end subroutine generated_affine_p2_scalar_interval'
        write(unit,'(a)') 'end module fortfem_reference_scalar_intervals'
        close(unit)
        print '(a)', 'PASS: 13 exact affine P2 polynomial/map/derivative ' // &
            'identities; zero probes'
    end subroutine interval_affine_p2

    subroutine emit_interval(expressions,interval_spec)
        type(expr_t), intent(in) :: expressions(:)
        type(rigorous_kernel_spec_t), intent(in) :: interval_spec
        character(:), allocatable :: interval_code, message
        logical :: ok
        interval_code = chars(emit_rigorous_kernel(expressions,interval_spec, &
            ok=ok,message=message))
        if (.not. ok) then
            print *, message
            error stop 'P2 interval generation failed'
        end if
        write(unit,'(a)') interval_code
    end subroutine emit_interval

    subroutine nodal_scalar_coefficients(family, basis)
        character(*), intent(in) :: family
        type(expr_t), intent(in) :: basis(:)
        type(expr_t) :: field, coefficients(6), variables(2), origin(2)
        type(kernel_spec_t) :: coefficient_spec
        character(:), allocatable :: coefficient_code, message
        character(8) :: node_text, count_text
        logical :: ok
        integer :: node, coefficient

        variables = [x, y]
        origin = num(arena, 0)
        field = num(arena, 0)
        ! The caller centers on node one before invoking this product.
        do node = 2, size(basis)
            write(node_text, '(i0)') node
            field = field + sym(arena, &
                'centered_nodes('//trim(node_text)//')')*basis(node)
        end do
        coefficients(1) = subs_many(field, variables, origin)
        coefficients(2) = subs_many(diff(field, x), variables, origin)
        coefficients(3) = subs_many(diff(field, y), variables, origin)
        coefficients(4) = subs_many(diff(diff(field, x), x)/2, variables, origin)
        coefficients(5) = subs_many(diff(diff(field, x), y), variables, origin)
        coefficients(6) = subs_many(diff(diff(field, y), y)/2, variables, origin)
        do coefficient = 1, 6
            result = engine%simplify(coefficients(coefficient))
            if (.not. result%ok) error stop 'native scalar coefficient simplification failed'
            coefficients(coefficient) = result%value
        end do
        coefficient_spec%name = str('generated_reference_'//family//'_scalar_coefficients')
        coefficient_spec%module_name = str( &
            'fortfem_generated_reference_'//family//'_scalar_coefficients')
        coefficient_spec%mode = KERNEL_SUBROUTINE
        coefficient_spec%generator = str('gen_reference_basis_products')
        coefficient_spec%generator_revision = str(fortsym_revision())
        coefficient_spec%regenerate_command = str('cd tools/codegen && ./generate.sh')
        coefficient_spec%pure_procedure = .true.
        write(count_text, '(i0)') size(basis)
        coefficient_spec%args = [str('centered_nodes')]
        coefficient_spec%arg_shapes = [str('('//trim(count_text)//')')]
        coefficient_spec%outputs = [str('coefficients')]
        coefficient_spec%output_shapes = [str('(6)')]
        allocate(coefficient_spec%output_references(6))
        do coefficient = 1, 6
            write(node_text, '(i0)') coefficient
            coefficient_spec%output_references(coefficient) = &
                str('coefficients('//trim(node_text)//')')
        end do
        coefficient_code = chars(emit_kernel( &
            coefficients, coefficient_spec, ok=ok, message=message))
        if (.not. ok) then
            print *, message
            error stop 'native scalar coefficient generation failed'
        end if
        write(unit, '(a)') coefficient_code(:len(coefficient_code) - 1)
    end subroutine nodal_scalar_coefficients

    subroutine affine_geometry()
        type(expr_t)::vertices(2,3),mapped(2),jacobian(2,2),products(7)
        character(4)::row_text,column_text
        integer::row,column,k
        do column=1,3
            do row=1,2
                write(row_text,'(i0)')row;write(column_text,'(i0)')column
                vertices(row,column)=sym(arena,'vertices('//trim(row_text)//','//trim(column_text)//')')
            end do
        end do
        do row=1,2
            mapped(row)=vertices(row,1)+(vertices(row,2)-vertices(row,1))*x+ &
                (vertices(row,3)-vertices(row,1))*y
            jacobian(row,1)=diff(mapped(row),x)
            jacobian(row,2)=diff(mapped(row),y)
        end do
        products(1:2)=mapped
        products(3:6)=reshape(jacobian,[4])
        products(7)=jacobian(1,1)*jacobian(2,2)-jacobian(1,2)*jacobian(2,1)
        spec=kernel_spec_t()
        spec%name=str('generated_affine_triangle_geometry')
        spec%module_name=str('fortfem_generated_affine_triangle_geometry')
        spec%mode=KERNEL_SUBROUTINE
        spec%generator=str('gen_reference_basis_products')
        spec%generator_revision=str(fortsym_revision())
        spec%regenerate_command=str('cd tools/codegen && ./generate.sh')
        spec%pure_procedure=.true.
        spec%args=[str('xi'),str('eta'),str('vertices')]
        spec%arg_shapes=[str(''),str(''),str('(2,3)')]
        spec%outputs=[str('mapped'),str('jacobian'),str('determinant')]
        spec%output_shapes=[str('(2)'),str('(2,2)'),str('')]
        allocate(spec%output_references(7))
        spec%output_references(1:2)=[str('mapped(1)'),str('mapped(2)')]
        k=2
        do column=1,2
            do row=1,2
                k=k+1
                write(row_text,'(i0)')row;write(column_text,'(i0)')column
                spec%output_references(k)=str('jacobian('//trim(row_text)//','//trim(column_text)//')')
            end do
        end do
        spec%output_references(7)=str('determinant')
        call emit(products)
    end subroutine affine_geometry
    subroutine scalar_family(family, values)
        character(*), intent(in) :: family
        type(expr_t), intent(in) :: values(:)
        do i=1,size(values)
            roots(1)=values(i)
            roots(2)=diff(values(i),x); roots(3)=diff(values(i),y)
            roots(4)=diff(roots(2),x); roots(5)=diff(roots(2),y)
            roots(6)=diff(roots(3),y)
            call initialize(family,i,6)
            spec%outputs(1)=str('value'); spec%outputs(2)=str('gx')
            spec%outputs(3)=str('gy'); spec%outputs(4)=str('hxx')
            spec%outputs(5)=str('hxy'); spec%outputs(6)=str('hyy')
            call emit(roots)
        end do
        call dispatcher(family,size(values),.false.)
    end subroutine scalar_family

    subroutine vector_family(family, values)
        character(*), intent(in) :: family
        type(expr_t), intent(in) :: values(:,:)
        do i=1,size(values,2)
            roots(1)=values(1,i); roots(2)=values(2,i)
            roots(3)=diff(values(1,i),x)+diff(values(2,i),y)
            roots(4)=diff(values(2,i),x)-diff(values(1,i),y)
            call initialize(family,i,4)
            spec%outputs(1)=str('vx'); spec%outputs(2)=str('vy')
            spec%outputs(3)=str('divergence'); spec%outputs(4)=str('curl')
            call emit(roots(:4))
        end do
        call dispatcher(family,size(values,2),.true.)
    end subroutine vector_family

    subroutine initialize(family,basis_id,noutput)
        character(*), intent(in) :: family
        integer, intent(in) :: basis_id,noutput
        write(index_string,'(i0)') basis_id
        name = 'generated_'//family//'_'//trim(index_string)//'_jet'
        spec = kernel_spec_t()
        spec%name = str(name)
        spec%module_name = str('fortfem_'//name)
        spec%mode = KERNEL_SUBROUTINE
        spec%generator = str('gen_reference_basis_products')
        spec%generator_revision = str(fortsym_revision())
        spec%regenerate_command = str('cd tools/codegen && ./generate.sh')
        spec%pure_procedure = .true.
        allocate(spec%args(2),spec%outputs(noutput))
        spec%args(1)=str('xi'); spec%args(2)=str('eta')
    end subroutine initialize

    subroutine emit(expressions)
        type(expr_t), intent(inout) :: expressions(:)
        do j=1,size(expressions)
            result = engine%simplify(expressions(j))
            if (.not. result%ok) error stop 'reference basis simplification failed'
            expressions(j) = result%value
        end do
        code = chars(emit_kernel(expressions,spec))
        write(unit,'(a)') code(:len(code)-1)
    end subroutine emit

    subroutine dispatcher(family,n,vector)
        character(*), intent(in) :: family
        integer, intent(in) :: n
        logical, intent(in) :: vector
        integer :: k
        write(unit,'(a)') 'module fortfem_generated_'//family//'_basis'
        write(unit,'(a)') '    use, intrinsic :: iso_fortran_env, only: dp=>real64'
        do k=1,n
            write(index_string,'(i0)') k
            name='generated_'//family//'_'//trim(index_string)//'_jet'
            write(unit,'(a)') '    use fortfem_'//name//', only: '//name
        end do
        write(unit,'(a)') '    implicit none'
        write(unit,'(a)') '    private'
        write(unit,'(a)') '    public :: generated_'//family//'_jet'
        write(unit,'(a)') 'contains'
        if(vector)then
            write(unit,'(a)') '    pure subroutine generated_'//family// &
                '_jet(i,xi,eta,v,divergence,curl)'
            write(unit,'(a)') '        real(dp), intent(out) :: v(2),divergence,curl'
        else
            write(unit,'(a)') '    pure subroutine generated_'//family// &
                '_jet(i,xi,eta,v,g,h)'
            write(unit,'(a)') '        real(dp), intent(out) :: v,g(2),h(2,2)'
        end if
        write(unit,'(a)') '        integer, intent(in) :: i'
        write(unit,'(a)') '        real(dp), intent(in) :: xi,eta'
        write(unit,'(a)') '        select case(i)'
        do k=1,n
            write(index_string,'(i0)') k
            name='generated_'//family//'_'//trim(index_string)//'_jet'
            write(unit,'(a)') '        case('//trim(index_string)//')'
            if(vector)then
                write(unit,'(a)') '            call '//name// &
                    '(xi,eta,v(1),v(2),divergence,curl)'
            else
                write(unit,'(a)') '            call '//name// &
                    '(xi,eta,v,g(1),g(2),h(1,1),h(1,2),h(2,2))'
                write(unit,'(a)') '            h(2,1)=h(1,2)'
            end if
        end do
        write(unit,'(a)') '        case default'
        if(vector)then
            write(unit,'(a)') '            v=0;divergence=0;curl=0'
        else
            write(unit,'(a)') '            v=0;g=0;h=0'
        end if
        write(unit,'(a)') '        end select'
        write(unit,'(a)') '    end subroutine generated_'//family//'_jet'
        write(unit,'(a)') 'end module fortfem_generated_'//family//'_basis'
    end subroutine dispatcher
end program gen_reference_basis_products
