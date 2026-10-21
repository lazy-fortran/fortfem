program gen_parametric_bem_products
    use fortsym_arena, only: arena_t
    use fortsym_expr, only: expr_t, sym, operator(+), operator(-), &
        operator(*), operator(/), sqrt
    use fortsym_kernel, only: kernel_spec_t, emit_kernel, KERNEL_SUBROUTINE
    use fortsym_products, only: jvp
    use fortsym_string, only: chars, str, str_t
    use fortfem_codegen_provenance, only: generated_path, fortsym_revision
    implicit none
    type(arena_t), target :: arena
    type(expr_t) :: t(6), c(3), jac, d(6), r, f(4), g(2), a(8), curl(3)
    type(str_t) :: tn(6), dn(6), an(8)
    type(expr_t) :: pv(11), products(3)
    type(str_t) :: pn(11)
    integer :: unit, i
    call arena%init()
    tn = [str('a1'), str('a2'), str('a3'), str('b1'), str('b2'), str('b3')]
    dn = [str('d1'), str('d2'), str('d3'), str('n1'), str('n2'), str('n3')]
    an(:6) = tn
    an(7:) = [str('dxi'), str('deta')]
    do i = 1, 6
        t(i) = sym(arena, chars(tn(i)))
        d(i) = sym(arena, chars(dn(i)))
    end do
    do i = 1, 8
        a(i) = sym(arena, chars(an(i)))
    end do
    c = [t(2)*t(6)-t(3)*t(5), t(3)*t(4)-t(1)*t(6), &
        t(1)*t(5)-t(2)*t(4)]
    jac = sqrt(c(1)*c(1)+c(2)*c(2)+c(3)*c(3))
    f = [jac, c(1)/jac, c(2)/jac, c(3)/jac]
    r = sqrt(d(1)*d(1)+d(2)*d(2)+d(3)*d(3))
    ! The caller supplies 1/(4*pi), a compile-time constant, to scale G and K.
    g = [1/r, (d(1)*d(4)+d(2)*d(5)+d(3)*d(6))/(r*r*r)]
    curl = [a(4)*a(7)-a(1)*a(8), a(5)*a(7)-a(2)*a(8), &
        a(6)*a(7)-a(3)*a(8)]
    open(newunit=unit, file=generated_path('fortfem_parametric_bem_products.f90'), &
        status='replace', action='write')
    call emit('bem_frame', tn, t, f, &
        [str('jacobian'), str('n1'), str('n2'), str('n3')])
    call emit('bem_kernel', dn, d, g, [str('green'), str('double')])
    call emit('bem_curl_density', an, a, curl, &
        [str('c1'), str('c2'), str('c3')])
    pn = [str('green'),str('double'),str('jx'),str('jy'),str('phi'), &
        str('tx'),str('ty'),str('tz'),str('sx'),str('sy'),str('sz')]
    do i = 1, 11
        pv(i) = sym(arena, chars(pn(i)))
    end do
    products = [pv(1)*pv(3)*pv(4), pv(2)*pv(3)*pv(4)*pv(5), &
        pv(1)*(pv(6)*pv(9)+pv(7)*pv(10)+pv(8)*pv(11))]
    call emit('bem_pair_product', pn, pv, products, &
        [str('single_value'),str('double_value'),str('hyper_value')])
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
        integer :: k
        spec%name = str(name)
        spec%module_name = str('fortfem_generated_'//name)
        spec%generator = str('gen_parametric_bem_products')
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
    end subroutine
end program
