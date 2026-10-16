program gen_meridian_fourier_equilibration_products
    ! Fourier transport residual leaves; reconstruction is not differentiated.
    use fortsym_arena, only: arena_t
    use fortsym_engine, only: engine_result_t
    use fortsym_engine_native, only: make_native_engine, native_engine_t
    use fortsym_expr, only: expr_t, operator(*), operator(+), operator(-), &
        operator(/), sym
    use fortsym_kernel, only: emit_kernel, kernel_spec_t, KERNEL_SUBROUTINE
    use fortsym_string, only: chars, str, str_t
    use fortfem_codegen_provenance, only: fortsym_revision, generated_path
    implicit none
    type(arena_t), target :: arena
    type(native_engine_t) :: engine
    type(expr_t) :: r, n, uc, us, jc, js, fc, fs
    type(expr_t) :: a, b, c, x, y, z, u, f
    integer :: unit

    call arena%init()
    engine = make_native_engine(arena)
    open (newunit=unit, file=generated_path( &
        'fortfem_meridian_fourier_equilibration_products.f90'), &
        status='replace', action='write')
    r = sym(arena, 'r'); n = sym(arena, 'n')
    uc = sym(arena, 'uc'); us = sym(arena, 'us')
    call emit('fourier_equil_phi_gradient', &
        [str('r'), str('n'), str('uc'), str('us')], &
        [str('gc'), str('gs')], [n*us/r, -n*uc/r])

    a = sym(arena, 'a'); b = sym(arena, 'b'); c = sym(arena, 'c')
    x = sym(arena, 'x'); y = sym(arena, 'y'); z = sym(arena, 'z')
    u = sym(arena, 'u')
    ! One physical row of j = K grad p_h - f.
    call emit('fourier_equil_flux_component', &
        [str('a'), str('b'), str('c'), str('x'), str('y'), str('z'), str('u')], &
        [str('j')], [a*x + b*y + c*z - u])
    f = sym(arena, 'f')
    ! Metric-weighted scalar coefficient, meridional flux and source.
    call emit('fourier_equil_metric', &
        [str('r'), str('a'), str('x'), str('z'), str('f')], &
        [str('c'), str('px'), str('pz'), str('s')], [r*a, r*x, r*z, r*f])
    jc = sym(arena, 'jc'); js = sym(arena, 'js')
    fc = sym(arena, 'fc'); fs = sym(arena, 'fs')
    ! -div_2(R j_pol)_b = (R f0 + d_phi j_phi)_b.
    call emit('fourier_equil_source', &
        [str('r'), str('n'), str('fc'), str('fs'), str('jc'), str('js')], &
        [str('sc'), str('ss')], [r*fc + n*js, r*fs - n*jc])
    close (unit)
contains
    subroutine emit(name, arguments, outputs, expressions)
        character(*), intent(in) :: name
        type(str_t), intent(in) :: arguments(:), outputs(:)
        type(expr_t), intent(in) :: expressions(:)
        type(expr_t) :: roots(size(expressions))
        type(engine_result_t) :: result
        type(kernel_spec_t) :: spec
        character(:), allocatable :: code
        integer :: root

        do root = 1, size(expressions)
            result = engine%simplify(expressions(root))
            if (.not. result%ok) error stop "equilibrated products: simplify failed"
            roots(root) = result%value
        end do
        spec%name = str("generated_"//name)
        spec%module_name = str("fortfem_generated_"//name)
        spec%mode = KERNEL_SUBROUTINE
        spec%generator = str("gen_meridian_fourier_equilibration_products")
        spec%generator_revision = str(fortsym_revision())
        spec%regenerate_command = str("cd tools/codegen && ./generate.sh")
        spec%pure_procedure = .true.
        spec%args = arguments
        spec%outputs = outputs
        code = chars(emit_kernel(roots, spec))
        write (unit, "(a)") code(:len(code) - 1)
    end subroutine emit

end program gen_meridian_fourier_equilibration_products
