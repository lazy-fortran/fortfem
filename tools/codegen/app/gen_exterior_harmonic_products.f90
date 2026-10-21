program gen_exterior_harmonic_products
    use fortsym_arena, only: arena_t
    use fortsym_expr, only: expr_t, sym, operator(+), operator(-), &
                            operator(*), operator(/), operator(**), sqrt, log
    use fortsym_diff, only: diff
    use fortsym_kernel, only: kernel_spec_t, emit_kernel, KERNEL_SUBROUTINE
    use fortsym_products, only: jvp, vjp
    use fortsym_engine_native, only: native_engine_t, make_native_engine
    use fortsym_engine, only: engine_result_t, VERDICT_TRUE
    use fortsym_string, only: chars, str, str_t
    use fortfem_codegen_provenance, only: generated_path, fortsym_revision
    implicit none
    type(arena_t), target :: arena
    type(native_engine_t) :: engine
    type(engine_result_t) :: identity
    type(expr_t) :: v(9), r(3), d(3), avec(3), bvec(3), numerator(3), rr, f, dotdr
    type(expr_t) :: r1, r2, length, s, rho, primitive, outputs(6), pair(1)
    type(str_t) :: names(9)
    integer :: i, j, k, unit
    call arena%init()
    engine = make_native_engine(arena)
    names = [str('r1x'), str('r1y'), str('r1z'), str('r2x'), str('r2y'), &
             str('r2z'), str('nx'), str('ny'), str('nz')]
    do i = 1, 9
        v(i) = sym(arena, chars(names(i)))
    end do
    r = v(1:3); d = v(4:6)
    rr = r(1)**2 + r(2)**2 + r(3)**2
    f = 1/sqrt(rr)
    bvec = cross(d, r)
    do i = 1, 3
        avec(i) = d(i)*f
        bvec(i) = bvec(i)/(rr*sqrt(rr))
    end do
    call prove(diff(bvec(1), r(1)) + diff(bvec(2), r(2)) + &
               diff(bvec(3), r(3)), 'Biot-Savart divergence')
    numerator = cross(d, r)
    dotdr = d(1)*r(1) + d(2)*r(2) + d(3)*r(3)
    do i = 1, 3
        j = mod(i, 3) + 1; k = mod(i + 1, 3) + 1
        call prove(diff(avec(k), r(j)) - diff(avec(j), r(k)) - bvec(i), &
                   'curl vector potential')
        ! curl B is a total source derivative; a closed wire cancels endpoints.
        ! Clear the positive denominator rr**(5/2) before the zero test.
        call prove(rr*(diff(numerator(k), r(j)) - diff(numerator(j), r(k))) - &
                   3*(numerator(k)*r(j) - numerator(j)*r(k)) - &
                   (3*r(i)*dotdr - d(i)*rr), 'closed wire curl identity')
    end do
    s = sym(arena, 's'); rho = sym(arena, 'rho')
    primitive = log(s + sqrt(s*s + rho*rho))
    call prove(diff(primitive, s) - 1/sqrt(s*s + rho*rho), &
               'straight segment vector potential primitive')
    primitive = s/(rho*rho*sqrt(s*s + rho*rho))
    call prove(diff(primitive, s) - 1/(sqrt(s*s + rho*rho)**3), &
               'straight segment magnetic field primitive')
    ! Exact straight segment: r1=x-start, r2=x-end; mu0 I = 1.
    r1 = sqrt(v(1)**2 + v(2)**2 + v(3)**2)
    r2 = sqrt(v(4)**2 + v(5)**2 + v(6)**2)
    do i = 1, 3
        d(i) = v(i) - v(i + 3)
    end do
    length = sqrt(d(1)**2 + d(2)**2 + d(3)**2)
    bvec = cross(v(1:3), v(4:6))
    do i = 1, 3
        avec(i) = d(i)/length*log((r1 + r2 + length)/(r1 + r2 - length))
        bvec(i) = bvec(i)*(r1 + r2)/ &
                  (r1*r2*(r1*r2 + v(1)*v(4) + v(2)*v(5) + v(3)*v(6)))
    end do
    outputs = [avec, bvec]
    open (newunit=unit, file=generated_path('fortfem_exterior_harmonic_products.f90'), &
          status='replace', action='write')
    call emit('current_segment', names(:6), v(:6), outputs, &
              [str('ax'), str('ay'), str('az'), str('bx'), str('by'), str('bz')])
    ! Exterior Green identity: integral_ext B_i.B_j = -int_Gamma n.(A_i x B_j).
    pair = [-(v(7)*(v(2)*v(6) - v(3)*v(5)) + &
              v(8)*(v(3)*v(4) - v(1)*v(6)) + v(9)*(v(1)*v(5) - v(2)*v(4)))]
    call emit('exterior_gram_density', names, v, pair, [str('value')])
    close (unit)
    write (*, '(a)') 'Exterior harmonic symbolic identities PASS'
contains
    function cross(a, b) result(c)
        type(expr_t), intent(in) :: a(3), b(3)
        type(expr_t) :: c(3)
        c = [a(2)*b(3) - a(3)*b(2), a(3)*b(1) - a(1)*b(3), &
             a(1)*b(2) - a(2)*b(1)]
    end function
    subroutine prove(e, label)
        type(expr_t), intent(in) :: e
        character(*), intent(in) :: label
        identity = engine%zero_test(e)
        if (identity%verdict /= VERDICT_TRUE) then
            write (*, '(a)') label
            error stop 'exterior harmonic identity'
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
        spec%generator = str('gen_exterior_harmonic_products')
        spec%generator_revision = str(fortsym_revision())
        spec%regenerate_command = str('cd tools/codegen && ./generate.sh')
        spec%mode = KERNEL_SUBROUTINE
        spec%pure_procedure = .true.
        spec%args = names
        spec%outputs = outputs
        code = chars(emit_kernel(roots, spec))
        write (unit, '(a)') code(:len(code) - 1)
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
        write (unit, '(a)') code(:len(code) - 1)
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
        write (unit, '(a)') code(:len(code) - 1)
    end subroutine
end program
