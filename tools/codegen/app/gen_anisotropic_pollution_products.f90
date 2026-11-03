program gen_anisotropic_pollution_products
    use fortsym_arena, only: arena_t
    use fortsym_expr, only: expr_t, sym, real_expr, cos, sqrt, &
        sin, operator(+), operator(-), operator(*), operator(/), operator(**)
    use fortsym_diff, only: diff
    use fortsym_kernel, only: emit_kernel, kernel_spec_t, KERNEL_SUBROUTINE
    use fortsym_string, only: str, chars
    use fortfem_codegen_provenance, only: fortsym_revision
    implicit none
    type(arena_t), target :: arena
    type(kernel_spec_t) :: spec
    type(expr_t) :: x, y, pi, u, source, vx, vy, r, one
    type(expr_t) :: bx, by, dx, dy, v, gx, gy, hx, hy, w, q, numerical, exact
    type(expr_t) :: psi, field_square, perturbation, along, chi, flux_x, flux_y
    character(:), allocatable :: code, directory
    character(4096) :: buffer
    integer :: unit, length, status

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

    psi = u
    vx = cos(pi*x)*sin(pi*y)
    vy = -sin(pi*x)*cos(pi*y)
    field_square = vx*vx + vy*vy
    perturbation = x*psi
    u = psi + real_expr(arena, 4.0d0)*perturbation*field_square**2
    ! Cancel the field norm before differentiation to extend across its zeros.
    along = real_expr(arena, 4.0d0)*(field_square* &
        (vx*diff(perturbation, x) + vy*diff(perturbation, y)) + &
        real_expr(arena, 2.0d0)*perturbation* &
        (vx*diff(field_square, x) + vy*diff(field_square, y)))
    chi = sym(arena, 'chi')
    flux_x = diff(u, x) + (chi - one)*vx*along
    flux_y = diff(u, y) + (chi - one)*vy*along
    source = -(diff(flux_x, x) + diff(flux_y, y))
    spec%args = [str('x'), str('y'), str('pi'), str('chi')]
    call emit('manufactured', [u, source, vx, vy])

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

    spec%args = [str('bx'), str('by'), str('gx'), str('gy')]
    spec%outputs = [str('value')]
    call emit('directional', [bx*gx + by*gy])

    numerical = sym(arena, 'numerical')
    exact = sym(arena, 'exact')
    spec%args = [str('numerical'), str('exact'), str('w')]
    spec%outputs = [str('squared_error')]
    call emit('error', [w*(numerical - exact)**2])

    call get_environment_variable('FORTFEM_TEST_CODEGEN_OUTPUT_DIR', buffer, &
        length=length, status=status)
    if (status /= 0 .or. length == 0) then
        directory = '../../test/generated/'
    else
        directory = buffer(:length)//'/'
    end if
    open (newunit=unit, &
        file=directory//'fortfem_anisotropic_pollution_products.f90', &
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
