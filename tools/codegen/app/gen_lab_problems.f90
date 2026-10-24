program gen_lab_problems
    use fortsym_arena, only: arena_t
    use fortsym_diff, only: diff
    use fortsym_expr, only: expr_t, sym, sin, cos, sqrt, &
                            operator(+), operator(-), operator(*), operator(/)
    use fortsym_kernel, only: emit_kernel, kernel_spec_t, KERNEL_SUBROUTINE
    use fortsym_string, only: str, chars
    use fortfem_codegen_provenance, only: generated_path, fortsym_revision
    implicit none
    type(arena_t), target :: arena
    type(expr_t) :: x, y, z, r, phi, u, v, a, rc, zc, s, temp
    type(expr_t) :: roots(12), tr, tp, tz, norm, iota, pi
    call arena%init()
    x = sym(arena, 'x')
    y = sym(arena, 'y')
    z = sym(arena, 'z')
    temp = sin(x) + cos(y) + sin(z) + (x + y + z)/10
    roots(1) = temp
    roots(2) = diff(temp, x)
    roots(3) = diff(temp, y)
    roots(4) = diff(temp, z)
    roots(5) = -(diff(roots(2), x) + diff(roots(3), y) + diff(roots(4), z))
    roots(6) = sin(y)
    roots(7) = sin(z)
    roots(8) = sin(x)
    call emit('scalar', ['x  ', 'y  ', 'z  '], roots(:8))
    roots(1) = sin(y)
    roots(2) = sin(z)
    roots(3) = sin(x)
    roots(4) = -cos(z)
    roots(5) = -cos(x)
    roots(6) = -cos(y)
    roots(7) = 2*sin(y)
    roots(8) = 2*sin(z)
    roots(9) = 2*sin(x)
    call emit('vector', ['x  ', 'y  ', 'z  '], roots(:9))
    u = sym(arena, 'u')
    v = sym(arena, 'v')
    phi = sym(arena, 'phi')
    roots(1) = u*sqrt(1 - v*v/2)
    roots(2) = v*sqrt(1 - u*u/2)
    call emit('square_disk', ['u  ', 'v  '], roots(:2))
    a = 1 + cos(5*phi)/10
    rc = 3 + cos(5*phi)/10
    zc = sin(5*phi)/10
    roots(1) = rc + a*u
    roots(2) = zc + 13*v/10
    roots(3) = diff(roots(1), u)
    roots(4) = diff(roots(1), v)
    roots(5) = diff(roots(1), phi)
    roots(6) = diff(roots(2), u)
    roots(7) = diff(roots(2), v)
    roots(8) = diff(roots(2), phi)
    roots(9) = u*u + v*v
    roots(10) = (2 + roots(9))/5
    roots(11) = roots(1)*a*13/10
    roots(12) = roots(1)
    call emit('geometry', ['u  ', 'v  ', 'phi'], roots)
    r = sym(arena, 'r')
    u = (r - rc)/a
    v = (z - zc)*10/13
    s = u*u + v*v
    temp = 1 - s
    roots(1) = temp
    roots(2) = diff(temp, r)
    roots(3) = diff(temp, phi)/r
    roots(4) = diff(temp, z)
    roots(5) = -(diff(roots(2), r) + roots(2)/r + &
                 diff(roots(3), phi)/r + diff(roots(4), z))
    iota = (2 + s)/5
    tr = -sin(5*phi)*(1 + u)/2 - a*iota*v
    tp = r
    tz = cos(5*phi)/2 + 13*iota*u/10
    norm = sqrt(tr*tr + tp*tp + tz*tz)
    roots(6) = tr/norm
    roots(7) = tp/norm
    roots(8) = tz/norm
    roots(9) = s
    call emit('transport', ['r  ', 'phi', 'z  '], roots(:9))
    ! pi is an input so the mathematical constant is caller controlled.
    pi = sym(arena, 'pi')
    temp = cos(pi*x)*cos(pi*y)
    roots(1) = temp
    roots(2) = diff(temp, x)
    roots(3) = diff(temp, y)
    roots(4) = -(diff(roots(2), x) + diff(roots(3), y))
    norm = sqrt(roots(2)*roots(2) + roots(3)*roots(3))
    ! Return unnormalised direction. Wrapper handles magnetic nulls.
    roots(5) = -roots(3)
    roots(6) = roots(2)
    roots(7) = x*x + y*y
    roots(8) = temp
    call emit('guenter_kernel', ['x  ', 'y  ', 'pi '], roots(:8))
contains
    subroutine emit(name, args, values)
        character(*), intent(in) :: name, args(:)
        type(expr_t), intent(in) :: values(:)
        type(kernel_spec_t) :: spec
        character(:), allocatable :: code
        character(32) :: index_text, extent
        integer :: i, unit
        spec%name = str('lab_'//name)
        spec%module_name = str('fortfem_generated_lab_'//name)
        spec%mode = KERNEL_SUBROUTINE
        spec%generator = str('gen_lab_problems')
        spec%generator_revision = str(fortsym_revision())
        spec%regenerate_command = str('cd tools/codegen && ./generate.sh')
        spec%pure_procedure = .true.
        allocate (spec%args(size(args)), spec%outputs(1), spec%output_shapes(1))
        allocate (spec%output_references(size(values)))
        do i = 1, size(args)
            spec%args(i) = str(trim(args(i)))
        end do
        spec%outputs = [str('values')]
        write (extent, '(a,i0,a)') '(', size(values), ')'
        spec%output_shapes = [str(trim(extent))]
        do i = 1, size(values)
            write (index_text, '(a,i0,a)') 'values(', i, ')'
            spec%output_references(i) = str(trim(index_text))
        end do
        code = chars(emit_kernel(values, spec))
        open (newunit=unit, file=generated_path('fortfem_lab_'//name//'.f90'), &
              status='replace', action='write')
        write (unit, '(a)') code(:len(code) - 1)
        close (unit)
    end subroutine
end program
