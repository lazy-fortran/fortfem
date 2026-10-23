program gen_tensor_fem_bem_products
    !! Laplace representation kernels for evaluating a BEM potential and its
    !! gradient off the surface: G, d_nG and their gradients in the target.
    use fortsym_arena, only: arena_t
    use fortsym_expr, only: expr_t, sym, operator(+), operator(-), &
        operator(*), operator(/), sqrt
    use fortsym_diff, only: diff
    use fortsym_kernel, only: kernel_spec_t, emit_kernel, KERNEL_SUBROUTINE
    use fortsym_products, only: jvp
    use fortsym_engine_native, only: native_engine_t, make_native_engine
    use fortsym_engine, only: engine_result_t, VERDICT_TRUE
    use fortsym_string, only: chars, str, str_t
    use fortfem_codegen_provenance, only: generated_path, fortsym_revision
    implicit none
    type(arena_t), target :: arena
    type(native_engine_t) :: engine
    type(expr_t) :: d(6), r, green, double, gg(3), gd(3), roots(8)
    type(str_t) :: names(6)
    integer :: unit, i
    call arena%init()
    engine = make_native_engine(arena)
    names = [str('d1'), str('d2'), str('d3'), str('n1'), str('n2'), str('n3')]
    do i = 1, 6
        d(i) = sym(arena, chars(names(i)))
    end do
    ! d = x - y for target x and source y with source normal n; the caller
    ! supplies 1/(4*pi).
    r = sqrt(d(1)*d(1) + d(2)*d(2) + d(3)*d(3))
    green = 1/r
    double = (d(1)*d(4) + d(2)*d(5) + d(3)*d(6))/(r*r*r)
    do i = 1, 3
        gg(i) = diff(green, d(i))
        gd(i) = diff(double, d(i))
    end do
    ! d_nG is the source normal derivative -n.grad_x G for constant n, so the
    ! harmonicity of G proved here carries over to d_nG.
    call prove(double + gg(1)*d(4) + gg(2)*d(5) + gg(3)*d(6), 'double layer')
    call prove(diff(gg(1), d(1)) + diff(gg(2), d(2)) + diff(gg(3), d(3)), &
               'single layer harmonic')
    roots = [green, double, gg, gd]
    open (newunit=unit, file=generated_path('fortfem_tensor_fem_bem_products.f90'), &
          status='replace', action='write')
    call emit('bem_field_kernel', names, d, roots, [str('green'), str('double'), &
              str('g1'), str('g2'), str('g3'), str('h1'), str('h2'), str('h3')])
    close (unit)
contains
    subroutine prove(e, label)
        type(expr_t), intent(in) :: e
        character(*), intent(in) :: label
        type(engine_result_t) :: identity
        identity = engine%zero_test(e)
        if (identity%verdict /= VERDICT_TRUE) error stop label
    end subroutine

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
        spec%generator = str('gen_tensor_fem_bem_products')
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
    end subroutine
end program
