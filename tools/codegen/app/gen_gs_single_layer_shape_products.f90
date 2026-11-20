program gen_gs_single_layer_shape_products
    ! Integrand jets for shape derivatives of the Galerkin single-layer
    ! operator of the axisymmetric Grad-Shafranov kernel
    !
    !     G = s P(m)/(4 pi),  P = (2 - m) K(m) - 2 E(m),
    !     m = 4 R R'/s^2,  s^2 = (R + R')^2 + (Z - Z')^2,  q = 1 - m.
    !
    ! The integrand of one quadrature pair is f = G |t_x| |t_y| with t the
    ! panel tangents.  Near coincidence G = -c ln(dist) + h with
    ! c = s C/(4 pi), h = s (ln(s) C + D)/(4 pi), C = (1+q) A_K - 2 A_E,
    ! D = (1+q) B_K - 2 B_E (DLMF 19.12), and ln(dist) = ell + ln|delta|:
    ! delta is the difference of the two points formed without cancellation
    ! (divided by the parameter distance on one panel, where ell is its
    ! logarithm).  Every input is linear in the panel nodes, so the value,
    ! gradient and Hessian in the inputs below give the exact first and
    ! second shape derivatives of the discrete operator.
    !
    ! Load integrand: f = psi(x) |t| for a field psi given by its value,
    ! gradient and Hessian at the point.
    use fortsym_engine, only: engine_result_t
    use fortsym_engine_native, only: make_native_engine, native_engine_t
    use fortsym_expr, only: expr_t, operator(+), operator(-), operator(*), &
        operator(/), sqrt, log, sym, pi_expr, elliptic_k, elliptic_e
    use fortsym_elliptic_log, only: elliptic_k_log_a, elliptic_k_log_b, &
        elliptic_e_log_a, elliptic_e_log_b
    use fortsym_diff, only: diff
    use fortsym_subs, only: subs_many
    use fortsym_kernel, only: emit_kernel, kernel_spec_t, KERNEL_SUBROUTINE
    use fortsym_string, only: str_t, chars, str
    use fortfem_codegen_provenance, only: fortsym_revision, generated_path
    implicit none

    integer, parameter :: n_pair = 10, n_load = 4
    character(len=3), parameter :: pair_names(n_pair) = [character(len=3) :: &
        "r", "z", "rp", "zp", "txr", "txz", "tyr", "tyz", "dr", "dz"]
    character(len=3), parameter :: load_names(n_load) = [character(len=3) :: &
        "r", "z", "txr", "txz"]
    type(native_engine_t) :: engine
    type(expr_t) :: v(n_pair), ell, dd2, uu, q, s, m, a_k, b_k, cc, dd
    type(expr_t) :: coeff, smooth, g_near, g_far, speed_x, speed_y
    type(expr_t) :: w(n_load), psi, hr, hz, f_load
    type(expr_t) :: p(6)
    character(:), allocatable :: code_near, code_far, code_load
    character(:), allocatable :: code_near0, code_near1, code_far0, code_far1
    integer :: i, ios, unit

    engine = make_native_engine()
    do i = 1, n_pair
        v(i) = sym(trim(pair_names(i)))
    end do
    ell = sym("ell")
    dd2 = (v(1) - v(3))*(v(1) - v(3)) + (v(2) - v(4))*(v(2) - v(4))
    uu = (v(1) + v(3))*(v(1) + v(3)) + (v(2) - v(4))*(v(2) - v(4))
    q = dd2/uu
    s = sqrt(uu)
    m = 4*v(1)*v(3)/uu
    speed_x = sqrt(v(5)*v(5) + v(6)*v(6))
    speed_y = sqrt(v(7)*v(7) + v(8)*v(8))

    a_k = elliptic_k_log_a(0, q)
    b_k = elliptic_k_log_b(0, q)
    cc = (1 + q)*a_k - 2*elliptic_e_log_a(q)
    dd = (1 + q)*b_k - 2*elliptic_e_log_b(q)
    coeff = s*cc/(4*pi_expr())
    smooth = s*(log(s)*cc + dd)/(4*pi_expr())
    g_near = -coeff*(ell + log(v(9)*v(9) + v(10)*v(10))/2) + smooth
    g_far = s*((2 - m)*elliptic_k(m) - 2*elliptic_e(m))/(4*pi_expr())

    code_near = jet_kernel(g_near*speed_x*speed_y, v, pair_names, &
        "generated_gs_single_layer_near_jet", [(str(trim(pair_names(i))), &
        i=1, n_pair), str("ell")])
    code_far = jet_kernel(g_far*speed_x*speed_y, v(1:8), pair_names(1:8), &
        "generated_gs_single_layer_far_jet", [(str(trim(pair_names(i))), &
        i=1, 8)])
    ! Value only (operator assembly) and value with gradient (residuals).
    code_near0 = jet_kernel(g_near*speed_x*speed_y, v, pair_names, &
        "generated_gs_single_layer_near_value", [(str(trim(pair_names(i))), &
        i=1, n_pair), str("ell")], order=0)
    code_near1 = jet_kernel(g_near*speed_x*speed_y, v, pair_names, &
        "generated_gs_single_layer_near_gradient", [(str(trim(pair_names(i))), &
        i=1, n_pair), str("ell")], order=1)
    code_far0 = jet_kernel(g_far*speed_x*speed_y, v(1:8), pair_names(1:8), &
        "generated_gs_single_layer_far_value", [(str(trim(pair_names(i))), &
        i=1, 8)], order=0)
    code_far1 = jet_kernel(g_far*speed_x*speed_y, v(1:8), pair_names(1:8), &
        "generated_gs_single_layer_far_gradient", [(str(trim(pair_names(i))), &
        i=1, 8)], order=1)

    ! psi as its second-order Taylor polynomial about the point; the jets of
    ! psi |t| at the point are exact.
    p = [sym("psi"), sym("psi_r"), sym("psi_z"), sym("psi_rr"), &
        sym("psi_rz"), sym("psi_zz")]
    do i = 1, n_load
        w(i) = sym(trim(load_names(i)))
    end do
    hr = w(1) - sym("r0")
    hz = w(2) - sym("z0")
    psi = p(1) + p(2)*hr + p(3)*hz + p(4)*hr*hr/2 + p(5)*hr*hz + p(6)*hz*hz/2
    f_load = psi*sqrt(w(3)*w(3) + w(4)*w(4))
    code_load = jet_kernel(f_load, w, load_names, &
        "generated_gs_boundary_load_jet", &
        [str("txr"), str("txz"), str("psi"), str("psi_r"), str("psi_z"), &
        str("psi_rr"), str("psi_rz"), str("psi_zz")], w(1:2), &
        [sym("r0"), sym("z0")])

    open (newunit=unit, &
        file=generated_path("fortfem_gs_single_layer_shape_products.f90"), &
        status="replace", action="write", iostat=ios)
    if (ios /= 0) error stop "cannot write GS single-layer shape products"
    write (unit, "(a)") code_near(:len(code_near) - 1)
    write (unit, "(a)") code_far(:len(code_far) - 1)
    write (unit, "(a)") code_load(:len(code_load) - 1)
    write (unit, "(a)") code_near0(:len(code_near0) - 1)
    write (unit, "(a)") code_near1(:len(code_near1) - 1)
    write (unit, "(a)") code_far0(:len(code_far0) - 1)
    write (unit, "(a)") code_far1(:len(code_far1) - 1)
    close (unit)

contains

    function jet_kernel(f, vars, names, name, args, at_old, at_new, order) &
            result(code)
        !! Value, gradient and upper-triangle Hessian of f in vars, packed
        !! row by row (H11, H12, ..., H1n, H22, ...), as a kernel of args.
        !! at_old -> at_new is substituted after differentiation.
        type(expr_t), intent(in) :: f, vars(:)
        character(*), intent(in) :: names(:), name
        type(str_t), intent(in) :: args(:)
        type(expr_t), intent(in), optional :: at_old(:), at_new(:)
        integer, intent(in), optional :: order
        character(:), allocatable :: code
        type(expr_t), allocatable :: roots(:), grad(:)
        type(kernel_spec_t) :: spec
        integer :: n, k, i, j, highest

        n = size(vars)
        highest = 2
        if (present(order)) highest = order
        allocate (grad(n))
        if (highest == 0) then
            allocate (roots(1))
        else if (highest == 1) then
            allocate (roots(1 + n))
        else
            allocate (roots(1 + n + n*(n + 1)/2))
        end if
        roots(1) = f
        if (highest >= 1) then
            do i = 1, n
                grad(i) = diff(f, vars(i))
                roots(1 + i) = grad(i)
            end do
        end if
        k = 1 + n
        if (highest == 2) then
            do i = 1, n
                do j = i, n
                    k = k + 1
                    roots(k) = diff(grad(i), vars(j))
                end do
            end do
        end if
        if (present(at_old)) then
            do k = 1, size(roots)
                roots(k) = subs_many(roots(k), at_old, at_new)
            end do
        end if
        call simplify_all(roots)

        spec%name = str(name)
        spec%module_name = str("fortfem_"//name)
        spec%mode = KERNEL_SUBROUTINE
        spec%generator = str("gen_gs_single_layer_shape_products")
        spec%generator_revision = str(fortsym_revision())
        spec%regenerate_command = str("cd tools/codegen && ./generate.sh")
        spec%pure_procedure = .true.
        allocate (spec%outputs(size(roots)))
        spec%args = args
        spec%outputs(1) = str("value")
        if (highest == 0) then
            code = chars(emit_kernel(roots, spec))
            return
        end if
        do i = 1, n
            spec%outputs(1 + i) = str("grad_"//trim(names(i)))
        end do
        if (highest == 1) then
            code = chars(emit_kernel(roots, spec))
            return
        end if
        k = 1 + n
        do i = 1, n
            do j = i, n
                k = k + 1
                spec%outputs(k) = str("hess_"//trim(names(i))//"_"// &
                    trim(names(j)))
            end do
        end do
        code = chars(emit_kernel(roots, spec))
    end function jet_kernel

    subroutine simplify_all(expressions)
        type(expr_t), intent(inout) :: expressions(:)
        type(engine_result_t) :: result
        integer :: root

        do root = 1, size(expressions)
            result = engine%simplify(expressions(root))
            if (result%ok) expressions(root) = result%value
        end do
    end subroutine simplify_all

end program gen_gs_single_layer_shape_products
