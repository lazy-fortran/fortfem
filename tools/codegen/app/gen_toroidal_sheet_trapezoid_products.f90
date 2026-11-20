program gen_toroidal_sheet_trapezoid_products
    ! Kernels of the corrected trapezoidal rule for the energy of a surface
    ! current on a smooth toroidal surface sampled on a uniform (u, v) grid:
    !
    !   E = w^2/(8 pi) [ sum_i sum_{j /= i} k_i . k_j / |X_i - X_j| - sum_i Z_i |k_i|^2 ],
    !
    ! k = Phi_u X_v - Phi_v X_u, w = h_u h_v, and Z_i the Epstein zeta
    ! function at s = 1/2 of the lattice form Q(m, n) = |m h_u X_u + n h_v X_v|^2
    ! of point i, which corrects the punctured sum of the 1/r singularity.
    !
    ! pair_hvp: the directional derivative of the x gradient of the pair
    !   integrand (generated_toroidal_sheet_pair_* inputs) along (dx, dy).
    ! zeta_term, zeta_base: for Q with matrix (a, b; b, c), l = sqrt(ac - b^2),
    !   Z = zeta_base + sum_{(m, n) /= 0} zeta_term(m, n), the Ewald split
    !   of the form normalised to unit determinant:
    !   zeta_term = l^(-1/2) [erfc(sqrt(pi q/l))/sqrt(q/l) + erfc(sqrt(pi q*))/sqrt(q*)],
    !   q = a m^2 + 2 b m n + c n^2, q* = (c m^2 - 2 b m n + a n^2)/l,
    !   zeta_base = -4 l^(-1/2); value, gradient and Hessian in (a, b, c).
    ! diag: |k|^2 Z(a, b, c) of one point, a = h_u^2 |X_u|^2, b = h_u h_v
    !   X_u . X_v, c = h_v^2 |X_v|^2, with Z by its Taylor data in (a, b, c);
    !   value, gradient in the point inputs and its directional derivative.
    use fortsym_engine, only: engine_result_t
    use fortsym_engine_native, only: make_native_engine, native_engine_t
    use fortsym_expr, only: expr_t, operator(+), operator(-), operator(*), &
        operator(/), sqrt, sym, erfc, operator(**)
    use fortsym_diff, only: diff
    use fortsym_subs, only: subs_many
    use fortsym_products, only: jvp
    use fortsym_kernel, only: emit_kernel, kernel_spec_t, KERNEL_SUBROUTINE
    use fortsym_string, only: str_t, chars, str
    use fortfem_codegen_provenance, only: fortsym_revision, generated_path
    implicit none

    character(len=4), parameter :: px(8) = [character(len=4) :: &
        "rx", "zx", "rux", "zux", "rvx", "zvx", "pux", "pvx"]
    character(len=4), parameter :: py(8) = [character(len=4) :: &
        "ry", "zy", "ruy", "zuy", "rvy", "zvy", "puy", "pvy"]
    character(len=4), parameter :: params(7) = [character(len=4) :: &
        "t", "c1", "s1", "c2", "s2", "dc", "ds"]
    real(8), parameter :: pi = 3.14159265358979323846264338327950288_8
    character(len=4), parameter :: zj(10) = [character(len=4) :: &
        "z", "za", "zb", "zc", "zaa", "zab", "zac", "zbb", "zbc", "zcc"]
    type(native_engine_t) :: engine
    type(expr_t) :: x(8), y(8), p(7), dx(8), dy(8), k1(3), k2(3), q(3), f, g(8), &
        roots(10), a, b, c, m, n, l, qq, qs, zt, abc(3), d(8), zjet(10), zt3, h(3)
    type(str_t), allocatable :: args(:)
    character(len=8) :: names(10)
    character(:), allocatable :: code(:)
    integer :: i, unit, ios

    engine = make_native_engine()
    allocate (character(len=200000) :: code(4))

    ! pair HVP
    do i = 1, 8
        x(i) = sym(trim(px(i))); y(i) = sym(trim(py(i)))
        dx(i) = sym("d"//trim(px(i))); dy(i) = sym("d"//trim(py(i)))
    end do
    do i = 1, 7
        p(i) = sym(trim(params(i)))
    end do
    call sheet([x(1), x(3:8)], p(2), p(3), k1)
    call sheet([y(1), y(3:8)], p(4), p(5), k2)
    q(1) = (y(1) - x(1))/p(1)*p(4) + x(1)*p(6)
    q(2) = (y(1) - x(1))/p(1)*p(5) + x(1)*p(7)
    q(3) = (y(2) - x(2))/p(1)
    f = (k1(1)*k2(1) + k1(2)*k2(2) + k1(3)*k2(3))/sqrt(q(1)*q(1) + q(2)*q(2) + q(3)*q(3))
    do i = 1, 8
        g(i) = diff(f, x(i))
    end do
    roots(1:8) = jvp(g, [x, y], [dx, dy])
    args = [(str(trim(px(i))), i=1, 8), (str(trim(py(i))), i=1, 8), &
        (str(trim(params(i))), i=1, 7), (str("d"//trim(px(i))), i=1, 8), &
        (str("d"//trim(py(i))), i=1, 8)]
    names(1:8) = [character(len=8) :: "h_r", "h_z", "h_ru", "h_zu", "h_rv", "h_zv", &
        "h_pu", "h_pv"]
    code(1) = kernel("generated_ts_trap_pair_hvp", roots(1:8), names(1:8), args)

    ! zeta term and base
    a = sym("a"); b = sym("b"); c = sym("c"); m = sym("m"); n = sym("n")
    abc = [a, b, c]
    l = sqrt(a*c - b*b)
    qq = (a*m*m + 2*b*m*n + c*n*n)/l
    qs = (c*m*m - 2*b*m*n + a*n*n)/l
    zt = (erfc(sqrt(pi*qq))/sqrt(qq) &
        + erfc(sqrt(pi*qs))/sqrt(qs))/sqrt(l)
    call jets3(zt, roots)
    names = [character(len=8) :: "z", "za", "zb", "zc", "zaa", "zab", "zac", "zbb", &
        "zbc", "zcc"]
    code(2) = kernel("generated_ts_trap_zeta_term", roots, names, &
        [str("a"), str("b"), str("c"), str("m"), str("n")])
    call jets3(-4/sqrt(l), roots)
    code(3) = kernel("generated_ts_trap_zeta_base", roots, names, &
        [str("a"), str("b"), str("c")])

    ! diag: |k|^2 Z with Z as its Taylor polynomial about (a0, b0, c0)
    do i = 1, 10
        zjet(i) = sym(trim(zj(i)))
    end do
    k1(1) = x(7)*x(5) - x(8)*x(3)
    k1(2) = x(7)*x(1)
    k1(3) = x(7)*x(6) - x(8)*x(4)
    h(1) = sym("hu")*sym("hu")*(x(3)*x(3) + x(4)*x(4)) - sym("a0")
    h(2) = sym("hu")*sym("hv")*(x(3)*x(5) + x(4)*x(6)) - sym("b0")
    h(3) = sym("hv")*sym("hv")*(x(5)*x(5) + x(1)*x(1) + x(6)*x(6)) - sym("c0")
    zt3 = zjet(1) + zjet(2)*h(1) + zjet(3)*h(2) + zjet(4)*h(3) &
        + zjet(5)*h(1)*h(1)/2 + zjet(6)*h(1)*h(2) + zjet(7)*h(1)*h(3) &
        + zjet(8)*h(2)*h(2)/2 + zjet(9)*h(2)*h(3) + zjet(10)*h(3)*h(3)/2
    f = (k1(1)*k1(1) + k1(2)*k1(2) + k1(3)*k1(3))*zt3
    do i = 1, 8
        d(i) = diff(f, x(i))
    end do
    ! the HVP is formed before the Taylor point is set to the point itself
    block
        type(expr_t) :: hv(8), all(17)
        hv = jvp(d, x, dx)
        all = [f, d, hv]
        do i = 1, 17
            all(i) = subs_many(all(i), [sym("a0"), sym("b0"), sym("c0")], &
                [sym("hu")*sym("hu")*(x(3)*x(3) + x(4)*x(4)), &
                sym("hu")*sym("hv")*(x(3)*x(5) + x(4)*x(6)), &
                sym("hv")*sym("hv")*(x(5)*x(5) + x(1)*x(1) + x(6)*x(6))])
        end do
        args = [(str(trim(px(i))), i=1, 8), str("hu"), str("hv"), &
            (str(trim(zj(i))), i=1, 10), (str("d"//trim(px(i))), i=1, 8)]
        code(4) = kernel17("generated_ts_trap_diag", all, args)
    end block

    open (newunit=unit, file=generated_path("fortfem_toroidal_sheet_trapezoid_products.f90"), &
        status="replace", action="write", iostat=ios)
    if (ios /= 0) error stop "cannot write toroidal sheet trapezoid products"
    do i = 1, 4
        write (unit, "(a)") trim(code(i))
    end do
    close (unit)

contains

    subroutine sheet(v, cv, sv, k)
        !! K dA = p_u X_v - p_v X_u, v = (r, r_u, z_u, r_v, z_v, p_u, p_v).
        type(expr_t), intent(in) :: v(7), cv, sv
        type(expr_t), intent(out) :: k(3)
        type(expr_t) :: xu(3), xv(3)
        integer :: i
        xu = [v(2)*cv, v(2)*sv, v(3)]
        xv = [v(4)*cv - v(1)*sv, v(4)*sv + v(1)*cv, v(5)]
        do i = 1, 3
            k(i) = v(6)*xv(i) - v(7)*xu(i)
        end do
    end subroutine sheet

    subroutine jets3(e, r)
        !! Value, gradient and upper Hessian (row-wise) in (a, b, c).
        type(expr_t), intent(in) :: e
        type(expr_t), intent(out) :: r(10)
        type(expr_t) :: gg(3)
        integer :: i, j, k
        r(1) = e
        do i = 1, 3
            gg(i) = diff(e, abc(i))
            r(1 + i) = gg(i)
        end do
        k = 4
        do i = 1, 3
            do j = i, 3
                k = k + 1
                r(k) = diff(gg(i), abc(j))
            end do
        end do
    end subroutine jets3

    function kernel(name, rts, outs, kargs) result(text)
        character(*), intent(in) :: name
        type(expr_t), intent(in) :: rts(:)
        character(*), intent(in) :: outs(:)
        type(str_t), intent(in) :: kargs(:)
        character(:), allocatable :: text
        type(expr_t) :: r(size(rts))
        type(kernel_spec_t) :: spec
        type(engine_result_t) :: result
        integer :: i
        do i = 1, size(rts)
            r(i) = rts(i)
            result = engine%simplify(r(i))
            if (result%ok) r(i) = result%value
        end do
        spec%name = str(name)
        spec%module_name = str("fortfem_"//name)
        spec%mode = KERNEL_SUBROUTINE
        spec%generator = str("gen_toroidal_sheet_trapezoid_products")
        spec%generator_revision = str(fortsym_revision())
        spec%regenerate_command = str("cd tools/codegen && ./generate.sh")
        spec%pure_procedure = .true.
        spec%args = kargs
        allocate (spec%outputs(size(rts)))
        do i = 1, size(rts)
            spec%outputs(i) = str(trim(outs(i)))
        end do
        text = chars(emit_kernel(r, spec))
    end function kernel

    function kernel17(name, rts, kargs) result(text)
        character(*), intent(in) :: name
        type(expr_t), intent(in) :: rts(17)
        type(str_t), intent(in) :: kargs(:)
        character(:), allocatable :: text
        character(len=8) :: outs(17)
        integer :: i
        outs(1) = "value"
        do i = 1, 8
            outs(1 + i) = "g_"//trim(px(i))
            outs(9 + i) = "h_"//trim(px(i))
        end do
        text = kernel(name, rts, outs, kargs)
    end function kernel17

end program gen_toroidal_sheet_trapezoid_products
