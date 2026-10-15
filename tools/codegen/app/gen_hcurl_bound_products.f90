program gen_hcurl_bound_products
    !! Pointwise integrands of the planar H(curl) reaction-type a posteriori
    !! energy bound (curl(nu curl u) + kappa u = f, u x n = 0): the rot leaf
    !! of scalar patch reconstructions, the patch target, the local matrix and
    !! right-hand-side contributions, and the two global defect integrands.
    use fortsym_arena, only: arena_t
    use fortsym_diff, only: diff
    use fortsym_engine, only: engine_result_t
    use fortsym_engine_native, only: make_native_engine, native_engine_t
    use fortsym_expr, only: expr_t, operator(+), operator(-), operator(*), &
        operator(/), sym, num
    use fortsym_kernel, only: emit_kernel, kernel_spec_t, KERNEL_SUBROUTINE
    use fortsym_string, only: chars, str, str_t
    use fortfem_codegen_provenance, only: fortsym_revision, generated_path
    implicit none

    type(arena_t), target :: arena
    type(native_engine_t) :: engine
    type(engine_result_t) :: result
    type(kernel_spec_t) :: spec
    type(expr_t) :: j11, j21, j12, j22, r1, r2, det
    type(expr_t) :: psi, s, g1, g2, rpsi1, rpsi2, phi_i, phi_j
    type(expr_t) :: ri1, ri2, rj1, rj2, nu, kappa, w, h, rh1, rh2
    type(expr_t) :: rot1, rot2, t1, t2, roots(2), one_root(1)
    type(expr_t) :: l(3), basis(6), pv(6), gx(6), gy(6), xr, er, outv(18)
    type(str_t) :: names(18)
    integer :: unit, ios, i
    character(2) :: label

    call arena%init()
    engine = make_native_engine(arena)
    j11 = sym(arena, "j11")
    j21 = sym(arena, "j21")
    j12 = sym(arena, "j12")
    j22 = sym(arena, "j22")
    r1 = sym(arena, "r1")
    r2 = sym(arena, "r2")
    det = j11*j22 - j12*j21
    psi = sym(arena, "psi")
    s = sym(arena, "s")
    g1 = sym(arena, "g1")
    g2 = sym(arena, "g2")
    rpsi1 = sym(arena, "rpsi1")
    rpsi2 = sym(arena, "rpsi2")
    phi_i = sym(arena, "phi_i")
    phi_j = sym(arena, "phi_j")
    ri1 = sym(arena, "ri1")
    ri2 = sym(arena, "ri2")
    rj1 = sym(arena, "rj1")
    rj2 = sym(arena, "rj2")
    nu = sym(arena, "nu")
    kappa = sym(arena, "kappa")
    w = sym(arena, "w")
    h = sym(arena, "h")
    rh1 = sym(arena, "rh1")
    rh2 = sym(arena, "rh2")

    open (newunit=unit, file=generated_path("fortfem_hcurl_bound_products.f90"), &
        status="replace", action="write", iostat=ios)
    if (ios /= 0) error stop "cannot write hcurl bound products"

    ! Reference P2 Lagrange basis (3 vertices, then the midpoints of local
    ! edges 1-2, 2-3, 3-1) and its reference gradient at (xi, eta).
    xr = sym(arena, "xi")
    er = sym(arena, "eta")
    l = [num(arena, 1) - xr - er, xr, er]
    do i = 1, 3
        basis(i) = l(i)*(num(arena, 2)*l(i) - num(arena, 1))
        basis(3 + i) = num(arena, 4)*l(i)*l(mod(i, 3) + 1)
    end do
    do i = 1, 6
        write (label, "(i0)") i
        names(i) = str("p"//trim(label))
        names(6 + i) = str("gx"//trim(label))
        names(12 + i) = str("gy"//trim(label))
        pv(i) = basis(i)
        gx(i) = diff(basis(i), xr)
        gy(i) = diff(basis(i), er)
    end do
    outv = simplified([pv, gx, gy])
    call emit_one("hcurl_p2_basis", outv, [str("xi"), str("eta")], names)

    ! rot of the physical gradient: (d_y phi, -d_x phi), grad = J^{-T} grad_ref.
    rot1 = (-j12*r1 + j11*r2)/det
    rot2 = -(j22*r1 - j21*r2)/det
    roots = simplified([rot1, rot2])
    call emit_one("hcurl_rot_leaf", roots, &
        [str("j11"), str("j21"), str("j12"), str("j22"), str("r1"), str("r2")], &
        [str("rot1"), str("rot2")])

    ! Local SPD matrix entry: w (phi_i phi_j / nu + rot_i . rot_j / kappa).
    one_root(1) = w*(phi_i*phi_j/nu + (ri1*rj1 + ri2*rj2)/kappa)
    one_root = simplified(one_root)
    call emit_one("hcurl_patch_matrix", one_root, &
        [str("w"), str("nu"), str("kappa"), str("phi_i"), str("phi_j"), &
        str("ri1"), str("ri2"), str("rj1"), str("rj2")], [str("contribution")])

    ! Local right-hand side entry: w (psi s phi_j / nu + T . rot_j / kappa).
    ! Patch target T = psi g + s rot psi (the rot-equation of the local problem).
    t1 = psi*g1 + s*rpsi1
    t2 = psi*g2 + s*rpsi2
    one_root(1) = w*(psi*s*phi_j/nu + (t1*rj1 + t2*rj2)/kappa)
    one_root = simplified(one_root)
    call emit_one("hcurl_patch_rhs", one_root, &
        [str("w"), str("nu"), str("kappa"), str("psi"), str("s"), str("g1"), &
        str("g2"), str("rpsi1"), str("rpsi2"), str("phi_j"), str("rj1"), &
        str("rj2")], [str("contribution")])

    ! Global defect integrands: flux term w (h - s)^2 / nu and field term
    ! w |g - rot h|^2 / kappa.
    roots = simplified([w*(h - s)*(h - s)/nu, &
        w*((g1 - rh1)*(g1 - rh1) + (g2 - rh2)*(g2 - rh2))/kappa])
    call emit_one("hcurl_defect", roots, &
        [str("w"), str("nu"), str("kappa"), str("h"), str("s"), str("g1"), &
        str("g2"), str("rh1"), str("rh2")], [str("flux_sq"), str("field_sq")])
    close (unit)

contains

    function simplified(expressions) result(out)
        type(expr_t), intent(in) :: expressions(:)
        type(expr_t) :: out(size(expressions))
        integer :: i
        type(engine_result_t) :: simple

        out = expressions
        do i = 1, size(expressions)
            simple = engine%simplify(expressions(i))
            if (simple%ok) out(i) = simple%value
        end do
    end function simplified

    subroutine emit_one(name, roots_in, args, outputs)
        character(*), intent(in) :: name
        type(expr_t), intent(in) :: roots_in(:)
        type(str_t), intent(in) :: args(:), outputs(:)
        character(:), allocatable :: code
        integer :: nargs

        nargs = size(args)
        spec%name = str("generated_"//name)
        spec%module_name = str("fortfem_generated_"//name)
        spec%mode = KERNEL_SUBROUTINE
        spec%generator = str("gen_hcurl_bound_products")
        spec%generator_revision = str(fortsym_revision())
        spec%regenerate_command = str("cd tools/codegen && ./generate.sh")
        spec%pure_procedure = .true.
        allocate (spec%args(nargs), spec%outputs(size(outputs)))
        spec%args = args
        spec%outputs = outputs
        code = chars(emit_kernel(roots_in, spec))
        write (unit, "(a)") code(:len(code) - 1)
        deallocate (spec%args, spec%outputs)
    end subroutine emit_one

end program gen_hcurl_bound_products
