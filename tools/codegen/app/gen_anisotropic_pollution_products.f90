program gen_anisotropic_pollution_products
    use fortsym_arena, only: arena_t
    use fortsym_expr, only: expr_t, sym, real_expr, cos, sqrt, &
        operator(+), operator(-), operator(*), operator(/), operator(**)
    use fortsym_diff, only: diff
    use fortsym_kernel, only: emit_kernel, kernel_spec_t, KERNEL_SUBROUTINE
    use fortsym_string, only: str, chars
    use fortfem_codegen_provenance, only: fortsym_revision, generated_path
    implicit none
    type(arena_t), target :: arena
    type(kernel_spec_t) :: spec
    type(expr_t) :: x, y, pi, u, source, vx, vy, r, one
    type(expr_t) :: bx, by, dx, dy, v, gx, gy, hx, hy, w, q, numerical, exact
    character(:), allocatable :: code
    integer :: unit

    call arena%init()
    spec%mode = KERNEL_SUBROUTINE
    spec%generator = str('gen_anisotropic_pollution_products')
    spec%generator_revision = str(fortsym_revision())
    spec%regenerate_command = str('tools/codegen/generate.sh')
    spec%pure_procedure = .true.
    x = sym(arena, 'x')
    y = sym(arena, 'y')
    pi = sym(arena, 'pi')
    one = real_expr(arena, 1.0d0)
    u = cos(pi*x)*cos(pi*y)
    source = -(diff(diff(u, x), x) + diff(diff(u, y), y))
    vx = -diff(u, y)
    vy = diff(u, x)
    spec%args = [str('x'), str('y'), str('pi')]
    spec%outputs = [str('value'), str('source'), str('vx'), str('vy')]
    call emit('sovinec', [u, source, vx, vy])

    r = sqrt(x*x + y*y)
    u = one - r**3
    ! The simplified Laplacian extends continuously to the magnetic axis.
    source = real_expr(arena, 9.0d0)*r
    vx = -y
    vy = x
    spec%args = [str('x'), str('y')]
    call emit('circle', [u, source, vx, vy])

    bx = sym(arena, 'bx')
    by = sym(arena, 'by')
    dx = sym(arena, 'dx')
    dy = sym(arena, 'dy')
    spec%args = [str('bx'), str('by'), str('dx'), str('dy')]
    spec%outputs = [str('value'), str('gx'), str('gy')]
    call emit('basis', [bx*by, dx*by, bx*dy])

    gx = sym(arena, 'gx')
    gy = sym(arena, 'gy')
    hx = sym(arena, 'hx')
    hy = sym(arena, 'hy')
    w = sym(arena, 'w')
    q = sym(arena, 'q')
    v = sym(arena, 'v')
    spec%args = [str('gx'), str('gy'), str('hx'), str('hy'), &
        str('bx'), str('by'), str('w'), str('q'), str('v')]
    spec%outputs = [str('slow'), str('parallel'), str('load')]
    call emit('integrand', [w*(gx*hx + gy*hy), &
        w*(bx*gx + by*gy)*(bx*hx + by*hy), w*q*v])

    numerical = sym(arena, 'numerical')
    exact = sym(arena, 'exact')
    spec%args = [str('numerical'), str('exact'), str('w')]
    spec%outputs = [str('squared_error')]
    call emit('error', [w*(numerical - exact)**2])

    open (newunit=unit, &
        file=generated_path('fortfem_anisotropic_pollution_products.f90'), &
        status='replace', action='write')
    write (unit, '(a)') code(:len(code) - 1)
    close (unit)
contains
    subroutine emit(name, roots)
        character(*), intent(in) :: name
        type(expr_t), intent(in) :: roots(:)
        spec%name = str('pollution_'//name)
        spec%module_name = str('fortfem_generated_pollution_'//name)
        if (.not. allocated(code)) code = ''
        code = code//chars(emit_kernel(roots, spec))
    end subroutine
end program
