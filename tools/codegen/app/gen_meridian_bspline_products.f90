program gen_meridian_bspline_products
    !! Pointwise kernels of the per-mode B-spline curl complex
    !! (src/assembly/assembly_meridian_bspline_n.f90) and its manufactured
    !! test fields.
    !!
    !! Physical basis.  A V1 field is A = (A_R, a / R, A_Z) exp(i n phi)
    !! with the spline components ax = A_R, ay = A_Z and a = R A_phi.  Its
    !! physical curl is
    !!
    !!   B_R = i n A_Z / R - d_Z a / R,   B_phi = d_Z A_R - d_R A_Z,
    !!   B_Z = d_R a / R - i n A_R / R.
    !!
    !! Manufactured fields (test only).  With the bump
    !! b = (1 - rho^2)^6, rho^2 = ((R - R_c)^2 + (Z - Z_c)^2) / w^2,
    !! the compactly supported potential A = b (1 + zeta, i (1 + xi),
    !! 1 - xi zeta) (xi, zeta the scaled offsets), its field B = curl_n A,
    !! the current curl_n B, the gradient grad_n chi of a smooth chi, and
    !! the n = 0 toroidal field T = (0, b, 0) with its curl.  Complex fields
    !! are split into real and imaginary parts.
    use fortsym_engine, only: engine_result_t
    use fortsym_engine_native, only: make_native_engine, native_engine_t
    use fortsym_expr, only: expr_t, sym, num, real_expr, operator(+), &
        operator(-), operator(*), operator(/), operator(**)
    use fortsym_diff, only: diff
    use fortsym_kernel, only: emit_kernel, kernel_spec_t, KERNEL_SUBROUTINE
    use fortsym_string, only: chars, str, str_t
    use fortfem_codegen_provenance, only: fortsym_revision, generated_path
    implicit none

    type(native_engine_t) :: engine
    type(expr_t) :: r, z, n, rc, zc, w, one, xi, zeta, b, chi
    type(expr_t) :: are(3), aim(3), bre(3), bim(3), cre(3), cim(3)
    type(expr_t) :: gre(3), gim(3), t(3), tc(3)
    type(expr_t) :: ax, axz, ay, ayr, c, cr, cz
    character(:), allocatable :: code
    character(4096) :: buffer
    character(:), allocatable :: test_dir
    integer :: unit, ios, length, status
    type(expr_t) :: axr, axrz, axzz, ayz, ayrr, ayrz, crr, crz, czz
    type(expr_t) :: jets(7), jets_r(7), jets_z(7), bphys(5), grad(10)
    integer :: k, m

    engine = make_native_engine()
    r = sym('r')
    z = sym('z')
    n = sym('n')
    one = real_expr(1.0d0)

    ! Physical value and curl_n of a V1 basis jet.
    ax = sym('ax')
    axz = sym('ax_z')
    ay = sym('ay')
    ayr = sym('ay_r')
    c = sym('c')
    cr = sym('c_r')
    cz = sym('c_z')
    code = ''
    call emit('generated_meridian_bspline_physical', &
        'fortfem_generated_meridian_bspline_physical', &
        [str('r'), str('n'), str('ax'), str('ax_z'), str('ay'), str('ay_r'), &
        str('c'), str('c_r'), str('c_z')], &
        [str('a_phi'), str('b_r_re'), str('b_r_im'), str('b_phi'), &
        str('b_z_re'), str('b_z_im')], &
        [c/r, -cz/r, n*ay/r, axz - ayr, cr/r, -(n*ax)/r])

    ! d/dR and d/dZ of the physical curl, by the chain rule through the
    ! basis jets; the second-derivative jets are the new inputs.
    axr = sym('ax_r')
    axrz = sym('ax_rz')
    axzz = sym('ax_zz')
    ayz = sym('ay_z')
    ayrr = sym('ay_rr')
    ayrz = sym('ay_rz')
    crr = sym('c_rr')
    crz = sym('c_rz')
    czz = sym('c_zz')
    jets = [ax, axz, ay, ayr, c, cr, cz]
    jets_r = [axr, axrz, ayr, ayrr, cr, crr, crz]
    jets_z = [axz, axzz, ayz, ayrz, cz, crz, czz]
    bphys = [-cz/r, n*ay/r, axz - ayr, cr/r, -(n*ax)/r]
    do m = 1, 5
        grad(m) = diff(bphys(m), r)
        grad(5 + m) = num(0)
        do k = 1, 7
            grad(m) = grad(m) + diff(bphys(m), jets(k))*jets_r(k)
            grad(5 + m) = grad(5 + m) + diff(bphys(m), jets(k))*jets_z(k)
        end do
    end do
    call emit('generated_meridian_bspline_physical_gradient', &
        'fortfem_generated_meridian_bspline_physical_gradient', &
        [str('r'), str('n'), str('ax'), str('ax_r'), str('ax_z'), str('ax_rz'), &
        str('ax_zz'), str('ay'), str('ay_r'), str('ay_z'), str('ay_rr'), &
        str('ay_rz'), str('c'), str('c_r'), str('c_z'), str('c_rr'), &
        str('c_rz'), str('c_zz')], &
        [str('b_r_re_r'), str('b_r_im_r'), str('b_phi_r'), str('b_z_re_r'), &
        str('b_z_im_r'), str('b_r_re_z'), str('b_r_im_z'), str('b_phi_z'), &
        str('b_z_re_z'), str('b_z_im_z')], grad)
    open (newunit=unit, &
        file=generated_path('fortfem_meridian_bspline_products.f90'), &
        status='replace', action='write', iostat=ios)
    if (ios /= 0) error stop 'cannot write meridian B-spline products'
    write (unit, '(a)') code(:len(code) - 1)
    close (unit)

    ! Manufactured test fields.
    rc = sym('rc')
    zc = sym('zc')
    w = sym('w')
    xi = (r - rc)/w
    zeta = (z - zc)/w
    b = (one - xi*xi - zeta*zeta)**6
    are = [b*(one + zeta), real_expr(0.0d0), b*(one - xi*zeta)]
    aim = [real_expr(0.0d0), b*(one + xi), real_expr(0.0d0)]
    call curl_n(are, aim, bre, bim)
    call curl_n(bre, bim, cre, cim)
    chi = real_expr(0.3d0)*r*r + real_expr(0.2d0)*r*z + real_expr(0.1d0)
    gre = [diff(chi, r), real_expr(0.0d0), diff(chi, z)]
    gim = [real_expr(0.0d0), n*chi/r, real_expr(0.0d0)]
    t = [real_expr(0.0d0), b, real_expr(0.0d0)]
    tc = [-diff(b, z), real_expr(0.0d0), diff(r*b, r)/r]
    code = ''
    call emit('generated_meridian_bspline_potential', &
        'fortfem_generated_meridian_bspline_potential', &
        [str('r'), str('z'), str('n'), str('rc'), str('zc'), str('w')], &
        [str('a_r_re'), str('a_phi_im'), str('a_z_re'), &
        str('b_r_re'), str('b_r_im'), str('b_phi_re'), str('b_phi_im'), &
        str('b_z_re'), str('b_z_im'), &
        str('c_r_re'), str('c_r_im'), str('c_phi_re'), str('c_phi_im'), &
        str('c_z_re'), str('c_z_im')], &
        [are(1), aim(2), are(3), bre(1), bim(1), bre(2), bim(2), bre(3), &
        bim(3), cre(1), cim(1), cre(2), cim(2), cre(3), cim(3)])
    call emit('generated_meridian_bspline_gauge', &
        'fortfem_generated_meridian_bspline_gauge', &
        [str('r'), str('z'), str('n')], &
        [str('g_r'), str('g_phi_im'), str('g_z')], [gre(1), gim(2), gre(3)])
    call emit('generated_meridian_bspline_toroidal', &
        'fortfem_generated_meridian_bspline_toroidal', &
        [str('r'), str('z'), str('rc'), str('zc'), str('w')], &
        [str('t_phi'), str('c_r'), str('c_z')], [t(2), tc(1), tc(3)])
    call get_environment_variable('FORTFEM_TEST_CODEGEN_OUTPUT_DIR', buffer, &
        length=length, status=status)
    if (status /= 0 .or. length == 0) then
        test_dir = '../../test/generated/'
    else
        test_dir = buffer(:length)//'/'
    end if
    open (newunit=unit, &
        file=test_dir//'fortfem_meridian_bspline_manufactured.f90', &
        status='replace', action='write', iostat=ios)
    if (ios /= 0) error stop 'cannot write meridian B-spline test fields'
    write (unit, '(a)') code(:len(code) - 1)
    close (unit)

contains

    subroutine curl_n(ure, uim, vre, vim)
        !! Physical curl of u exp(i n phi), u = ure + i uim.
        type(expr_t), intent(in) :: ure(3), uim(3)
        type(expr_t), intent(out) :: vre(3), vim(3)

        vre(1) = -(n*uim(3))/r - diff(ure(2), z)
        vim(1) = (n*ure(3))/r - diff(uim(2), z)
        vre(2) = diff(ure(1), z) - diff(ure(3), r)
        vim(2) = diff(uim(1), z) - diff(uim(3), r)
        vre(3) = diff(r*ure(2), r)/r + (n*uim(1))/r
        vim(3) = diff(r*uim(2), r)/r - (n*ure(1))/r
    end subroutine curl_n

    subroutine emit(name, module_name, args, outputs, roots)
        character(*), intent(in) :: name, module_name
        type(str_t), intent(in) :: args(:), outputs(:)
        type(expr_t), intent(in) :: roots(:)

        type(kernel_spec_t) :: spec
        type(engine_result_t) :: result
        type(expr_t) :: simple(size(roots))
        integer :: i

        do i = 1, size(roots)
            simple(i) = roots(i)
            result = engine%simplify(roots(i))
            if (result%ok) simple(i) = result%value
        end do
        spec%name = str(name)
        spec%module_name = str(module_name)
        spec%mode = KERNEL_SUBROUTINE
        spec%generator = str('gen_meridian_bspline_products')
        spec%generator_revision = str(fortsym_revision())
        spec%regenerate_command = str('cd tools/codegen && ./generate.sh')
        spec%pure_procedure = .true.
        spec%args = args
        spec%outputs = outputs
        code = code//chars(emit_kernel(simple, spec))
    end subroutine emit

end program gen_meridian_bspline_products
