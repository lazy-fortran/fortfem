program gen_reference_basis_products
    use fortsym_arena, only: arena_t
    use fortsym_engine, only: engine_result_t
    use fortsym_engine_native, only: make_native_engine, native_engine_t
    use fortsym_expr, only: expr_t, sym, operator(+), operator(-), operator(*), operator(/)
    use fortsym_diff, only: diff
    use fortsym_kernel, only: emit_kernel, kernel_spec_t, KERNEL_SUBROUTINE
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
    close(unit)
contains
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
