program gen_equilibrated_nonconforming_products
    !! Pointwise leaves of the equilibrated-residual bounds for mixed RT-DG and
    !! SIPG diffusion on affine triangles (estimators/
    !! equilibrated_nonconforming_2d.f90).  Only primal leaves are emitted: the
    !! estimators are not differentiated.  A = scalar diffusion coefficient.
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
    type(expr_t) :: x1, y1, x2, y2, x3, y3, j11, j12, j21, j22, det
    type(expr_t) :: jit11, jit12, jit21, jit22, rx, ry
    type(expr_t) :: w, a, gix, giy, gjx, gjy, sx, sy, gx, gy, lx, ly
    type(expr_t) :: jump, nx, ny, phi, pl, pm, sigma, div, proj
    integer :: unit

    call arena%init()
    engine = make_native_engine(arena)
    open (newunit=unit, file=generated_path( &
        "fortfem_equilibrated_nonconforming_products.f90"), &
        status="replace", action="write")

    ! Affine geometry: J = [x2-x1 x3-x1; y2-y1 y3-y1], det J, J^(-T) and the
    ! squared edge lengths (the cell diameter is their maximum root).
    x1 = sym(arena, "x1"); y1 = sym(arena, "y1")
    x2 = sym(arena, "x2"); y2 = sym(arena, "y2")
    x3 = sym(arena, "x3"); y3 = sym(arena, "y3")
    j11 = x2 - x1; j21 = y2 - y1
    j12 = x3 - x1; j22 = y3 - y1
    det = j11*j22 - j12*j21
    call emit("eq_cell_geometry", [str("x1"), str("y1"), str("x2"), str("y2"), &
        str("x3"), str("y3")], &
        [str("det"), str("jit11"), str("jit12"), str("jit21"), str("jit22"), &
        str("len12_sq"), str("len23_sq"), str("len31_sq")], &
        [det, j22/det, -j21/det, -j12/det, j11/det, &
        j11*j11 + j21*j21, (x3 - x2)*(x3 - x2) + (y3 - y2)*(y3 - y2), &
        j12*j12 + j22*j22])

    ! Physical gradient J^(-T) * reference gradient.
    jit11 = sym(arena, "jit11"); jit12 = sym(arena, "jit12")
    jit21 = sym(arena, "jit21"); jit22 = sym(arena, "jit22")
    rx = sym(arena, "rx"); ry = sym(arena, "ry")
    call emit("eq_physical_gradient", [str("jit11"), str("jit12"), str("jit21"), &
        str("jit22"), str("rx"), str("ry")], [str("gx"), str("gy")], &
        [jit11*rx + jit12*ry, jit21*rx + jit22*ry])

    w = sym(arena, "w"); a = sym(arena, "a")
    gix = sym(arena, "gix"); giy = sym(arena, "giy")
    gjx = sym(arena, "gjx"); gjy = sym(arena, "gjy")
    sx = sym(arena, "sx"); sy = sym(arena, "sy")
    gx = sym(arena, "gx"); gy = sym(arena, "gy")
    lx = sym(arena, "lx"); ly = sym(arena, "ly")

    ! Local potential postprocessing: minimise int A^-1 |sigma + A grad s|^2,
    ! i.e. K s = -c with K_ij = int A grad phi_i . grad phi_j and
    ! c_i = int sigma . grad phi_i (w = quadrature weight times det J).
    call emit("eq_stiffness_term", [str("w"), str("a"), str("gix"), str("giy"), &
        str("gjx"), str("gjy")], [str("term")], [w*a*(gix*gjx + giy*gjy)])
    call emit("eq_flux_load_term", [str("w"), str("sx"), str("sy"), str("gix"), &
        str("giy")], [str("term")], [w*(sx*gix + sy*giy)])

    ! Weighted L2 mass entry of a scalar space.
    pl = sym(arena, "pl"); pm = sym(arena, "pm")
    call emit("eq_mass_term", [str("w"), str("a"), str("pl"), str("pm")], &
        [str("term")], [w*a*pl*pm])

    ! A^-1 |sigma + A grad s|^2, the squared flux defect of the mixed method.
    call emit("eq_flux_defect_density", [str("w"), str("a"), str("sx"), &
        str("sy"), str("gx"), str("gy")], [str("density")], &
        [w*((sx + a*gx)*(sx + a*gx) + (sy + a*gy)*(sy + a*gy))/a])

    ! A |grad v|^2, the squared energy density.
    call emit("eq_energy_density", [str("w"), str("a"), str("gx"), str("gy")], &
        [str("density")], [w*a*(gx*gx + gy*gy)])

    ! Squared divergence defect (div sigma - Pi_k f)^2.
    div = sym(arena, "div"); proj = sym(arena, "proj")
    call emit("eq_divergence_defect_density", [str("w"), str("div"), &
        str("proj")], [str("density")], [w*(div - proj)*(div - proj)])

    ! SIPG: broken flux A grad u_h, and the load flux A (grad u_h - lift),
    ! lift = DG jump lifting.
    call emit("eq_dg_fluxes", [str("a"), str("gx"), str("gy"), str("lx"), &
        str("ly")], [str("flux_x"), str("flux_y"), str("load_x"), str("load_y")], &
        [a*gx, a*gy, a*(gx - lx), a*(gy - ly)])

    ! Face term of the lifting: int_F wght A [u] (tau . n), tau = phi e_c.
    jump = sym(arena, "jump"); nx = sym(arena, "nx"); ny = sym(arena, "ny")
    phi = sym(arena, "phi")
    call emit("eq_lifting_load", [str("w"), str("a"), str("jump"), str("nx"), &
        str("ny"), str("phi")], [str("load_x"), str("load_y")], &
        [w*a*jump*nx*phi, w*a*jump*ny*phi])

    ! Squared penalty jump seminorm density sigma_F [u]^2.
    sigma = sym(arena, "sigma")
    call emit("eq_jump_density", [str("w"), str("sigma"), str("jump")], &
        [str("density")], [w*sigma*jump*jump])

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
        spec%generator = str("gen_equilibrated_nonconforming_products")
        spec%generator_revision = str(fortsym_revision())
        spec%regenerate_command = str("cd tools/codegen && ./generate.sh")
        spec%pure_procedure = .true.
        spec%args = arguments
        spec%outputs = outputs
        code = chars(emit_kernel(roots, spec))
        write (unit, "(a)") code(:len(code) - 1)
    end subroutine emit

end program gen_equilibrated_nonconforming_products
