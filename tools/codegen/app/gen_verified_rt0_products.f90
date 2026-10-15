program gen_verified_rt0_products
    use fortsym_arena, only: arena_t
    use fortsym_expr, only: expr_t, sym, num, operator(+), operator(-), &
                            operator(*), operator(/), operator(**)
    use fortsym_diff, only: diff
    use fortsym_engine_native, only: make_native_engine, native_engine_t
    use fortsym_engine, only: engine_result_t, VERDICT_TRUE
    use fortsym_rigorous_emit, only: rigorous_kernel_spec_t, interval_runtime, &
                                     emit_rigorous_kernel
    use fortsym_string, only: chars, str
    use fortfem_codegen_provenance, only: generated_path, fortsym_revision
    implicit none
    type(arena_t), target :: arena
    type(native_engine_t) :: engine
    type(engine_result_t) :: result
    type(expr_t) :: x, y, x1, y1, x2, y2, x3, y3, q1, q2, q3
    type(expr_t) :: j11, j12, j21, j22, det, sx, sy, divergence, rx, ry
    type(expr_t) :: a, f, gx, gy, s1, s2, d
    type(rigorous_kernel_spec_t) :: spec
    character(:), allocatable :: code, message
    integer :: unit
    logical :: ok

    call arena%init()
    engine = make_native_engine(arena)
    x = sym(arena, 'xi'); y = sym(arena, 'eta')
    x1 = sym(arena, 'x1'); y1 = sym(arena, 'y1')
    x2 = sym(arena, 'x2'); y2 = sym(arena, 'y2')
    x3 = sym(arena, 'x3'); y3 = sym(arena, 'y3')
    q1 = sym(arena, 'q1'); q2 = sym(arena, 'q2'); q3 = sym(arena, 'q3')
    j11 = x2 - x1; j12 = x3 - x1; j21 = y2 - y1; j22 = y3 - y1
    det = j11*j22 - j12*j21
    ! Canonical RT0 reference fields: (xi,eta-1), (xi,eta), (xi-1,eta).
    ! Their outward reference-edge fluxes are respectively (1,0,0),
    ! (0,1,0), (0,0,1), and each divergence is two.
    rx = (q1 + q2 + q3)*x - q3
    ry = (q1 + q2 + q3)*y - q1
    divergence = diff(rx, x) + diff(ry, y)
    result = engine%zero_test(divergence - 2*(q1 + q2 + q3))
    if (result%verdict /= VERDICT_TRUE) error stop 'RT0 divergence identity'
    sx = (j11*rx + j12*ry)/det
    sy = (j21*rx + j22*ry)/det
    spec%name = str('verified_rt0_geometry')
    spec%generator = str('gen_verified_rt0_products; '//fortsym_revision())
    spec%runtime = interval_runtime('fortnum_interval', 'interval_t')
    spec%args = [str('xi'), str('eta'), str('x1'), str('y1'), str('x2'), str('y2'), &
                 str('x3'), str('y3'), str('q1'), str('q2'), str('q3')]
    spec%outputs = [str('x'), str('y'), str('sigma_x'), str('sigma_y'), &
                    str('divergence'), str('determinant')]
    open (newunit=unit, file=generated_path('fortfem_verified_rt0_products.f90'), &
          status='replace')
    write (unit, '(a)') 'module fortfem_verified_rt0_products'
    write (unit, '(a)') 'implicit none'
    write (unit, '(a)') 'contains'
    call emit([x1 + j11*x + j12*y, y1 + j21*x + j22*y, sx, sy, divergence/det, det])
    a = sym(arena, 'coefficient'); f = sym(arena, 'source')
    gx = sym(arena, 'gradient_x'); gy = sym(arena, 'gradient_y')
    s1 = sym(arena, 'sigma_x'); s2 = sym(arena, 'sigma_y')
    d = sym(arena, 'divergence'); det = sym(arena, 'determinant')
    spec%name = str('verified_rt0_integrands')
    spec%args = [str('coefficient'), str('source'), str('gradient_x'), &
                 str('gradient_y'), str('sigma_x'), str('sigma_y'), str('divergence'), &
                 str('determinant')]
    spec%outputs = [str('flux_squared'), str('residual_squared')]
    call emit([det*((s1 - a*gx)**2 + (s2 - a*gy)**2)/a, det*(f + d)**2])
    write (unit, '(a)') 'end module fortfem_verified_rt0_products'
    close (unit)
contains
    subroutine emit(roots)
        type(expr_t), intent(in) :: roots(:)
        code = chars(emit_rigorous_kernel(roots, spec, ok=ok, message=message))
        if (.not. ok) then
            print *, message
            error stop 'RT0 interval generation failed'
        end if
        write (unit, '(a)') code
    end subroutine emit
end program gen_verified_rt0_products
