program gen_toroidal_sheet_current_products
    ! Integrand jets for the magnetic energy of a surface current on a toroidal
    ! surface X(u, v) = (R cos v, R sin v, Z), R and Z functions of the
    ! parameters, and the current K dA = Phi_u X_v - Phi_v X_u (K = n x
    ! grad_s Phi):
    !
    !   pair:  (K_x dA_x) . (K_y dA_y) / |Q|,   X_y - X_x = t Q,
    !   load:  Phi B . (X_u x X_v),             B by its value, gradient and
    !                                           Hessian in (R, Z) at the point.
    !
    ! Pair inputs: per point x, y the values (r, z, r_u, z_u, r_v, z_v, p_u,
    ! p_v), with (R_y - R_x)/t and (Z_y - Z_x)/t formed inside; parameters t,
    ! cos v, sin v at both points and dc = (cos v_y - cos v_x)/t, ds = (sin v_y
    ! - sin v_x)/t. Every input is linear in the coefficients of R, Z and Phi
    ! at fixed parameters. The energy sums pairs symmetrically, so the kernels
    ! give the value, the gradient in the x inputs and the x-x and x-y blocks
    ! of the Hessian (exact derivatives of the discrete energy).
    use fortsym_engine, only: engine_result_t
    use fortsym_engine_native, only: make_native_engine, native_engine_t
    use fortsym_expr, only: expr_t, operator(+), operator(-), operator(*), &
        operator(/), sqrt, sym
    use fortsym_diff, only: diff
    use fortsym_subs, only: subs_many
    use fortsym_kernel, only: emit_kernel, kernel_spec_t, KERNEL_SUBROUTINE
    use fortsym_string, only: str_t, chars, str
    use fortfem_codegen_provenance, only: fortsym_revision, generated_path
    implicit none

    integer, parameter :: n_self = 16, n_load = 6
    character(len=4), parameter :: self_names(n_self) = [character(len=4) :: &
        "rx", "zx", "rux", "zux", "rvx", "zvx", "pux", "pvx", &
        "ry", "zy", "ruy", "zuy", "rvy", "zvy", "puy", "pvy"]
    character(len=4), parameter :: self_params(7) = [character(len=4) :: &
        "t", "c1", "s1", "c2", "s2", "dc", "ds"]
    character(len=4), parameter :: load_names(n_load) = [character(len=4) :: &
        "r", "z", "ru", "zu", "rv", "zv"]
    type(native_engine_t) :: engine
    type(expr_t) :: v(n_self), p(7), k1(3), k2(3), q(3), f_self, dr, dz
    type(expr_t) :: w(n_load), hr, hz, br, bp, bz, f_load, bj(18)
    type(str_t), allocatable :: args(:)
    character(:), allocatable :: code(:)
    character(len=16), parameter :: bnames(18) = [character(len=16) :: &
        "b_r", "b_r_r", "b_r_z", "b_r_rr", "b_r_rz", "b_r_zz", &
        "b_p", "b_p_r", "b_p_z", "b_p_rr", "b_p_rz", "b_p_zz", &
        "b_z", "b_z_r", "b_z_z", "b_z_rr", "b_z_rz", "b_z_zz"]
    integer :: i, ios, unit, order

    engine = make_native_engine()
    do i = 1, n_self
        v(i) = sym(trim(self_names(i)))
    end do
    do i = 1, 7
        p(i) = sym(trim(self_params(i)))
    end do
    call sheet([v(1), v(3:8)], p(2), p(3), k1)
    call sheet([v(9), v(11:16)], p(4), p(5), k2)
    ! X_y - X_x = t ((R_y - R_x) cos v_y/t + R_x dc, ..., (Z_y - Z_x)/t)
    dr = (v(9) - v(1))/p(1)
    dz = (v(10) - v(2))/p(1)
    q(1) = dr*p(4) + v(1)*p(6)
    q(2) = dr*p(5) + v(1)*p(7)
    q(3) = dz
    f_self = (k1(1)*k2(1) + k1(2)*k2(2) + k1(3)*k2(3)) &
        /sqrt(q(1)*q(1) + q(2)*q(2) + q(3)*q(3))

    ! Load: B at (r, z) as its second-order Taylor polynomial about (r0, z0).
    do i = 1, n_load
        w(i) = sym(trim(load_names(i)))
    end do
    do i = 1, 18
        bj(i) = sym(trim(bnames(i)))
    end do
    hr = w(1) - sym("r0")
    hz = w(2) - sym("z0")
    br = taylor(bj(1:6))
    bp = taylor(bj(7:12))
    bz = taylor(bj(13:18))
    ! N = X_u x X_v = (-Z_u R, Z_u R_v - R_u Z_v, R_u R) in (e_R, e_phi, e_Z)
    f_load = -br*w(4)*w(1) + bp*(w(4)*w(5) - w(3)*w(6)) + bz*w(3)*w(1)

    allocate (character(len=600000) :: code(7))
    args = [(str(trim(self_names(i))), i=1, n_self), &
        (str(trim(self_params(i))), i=1, 7)]
    do order = 0, 2
        code(order + 1) = pair_kernel(f_self, v, &
            "generated_toroidal_sheet_pair_"//suffix(order), args, order)
    end do
    args = [(str(trim(load_names(i))), i=3, n_load), &
        (str(trim(bnames(i))), i=1, 18)]
    ! the point enters through the Taylor data of B; R itself is an input
    code(4) = jet_kernel(f_load, w, load_names, &
        "generated_toroidal_sheet_load_jet", [str("r"), args], 2, w(1:2), &
        [sym("r0"), sym("z0")], keep_r=.true.)

    open (newunit=unit, file=generated_path("fortfem_toroidal_sheet_products.f90"), &
        status="replace", action="write", iostat=ios)
    if (ios /= 0) error stop "cannot write toroidal sheet products"
    do i = 1, 4
        write (unit, "(a)") trim(code(i))
    end do
    close (unit)

contains

    subroutine sheet(x, c, s, k)
        !! K dA = p_u X_v - p_v X_u with X_u = (r_u c, r_u s, z_u),
        !! X_v = (r_v c - r s, r_v s + r c, z_v); x = (r, r_u, z_u, r_v, z_v,
        !! p_u, p_v).
        type(expr_t), intent(in) :: x(7), c, s
        type(expr_t), intent(out) :: k(3)
        type(expr_t) :: xu(3), xv(3)
        integer :: i

        xu(1) = x(2)*c
        xu(2) = x(2)*s
        xu(3) = x(3)
        xv(1) = x(4)*c - x(1)*s
        xv(2) = x(4)*s + x(1)*c
        xv(3) = x(5)
        do i = 1, 3
            k(i) = x(6)*xv(i) - x(7)*xu(i)
        end do
    end subroutine sheet

    function taylor(b) result(t)
        type(expr_t), intent(in) :: b(6)
        type(expr_t) :: t

        t = b(1) + b(2)*hr + b(3)*hz + b(4)*hr*hr/2 + b(5)*hr*hz + b(6)*hz*hz/2
    end function taylor

    pure function suffix(order) result(s)
        integer, intent(in) :: order
        character(:), allocatable :: s

        select case (order)
        case (0); s = "value"
        case (1); s = "gradient"
        case default; s = "jet"
        end select
    end function suffix

    function pair_kernel(f, vars, name, kargs, order) result(text)
        !! Value, gradient in the x inputs vars(1:8), the upper triangle of
        !! the x-x Hessian block (row-wise) and the x-y block (row-wise).
        type(expr_t), intent(in) :: f, vars(16)
        character(*), intent(in) :: name
        type(str_t), intent(in) :: kargs(:)
        integer, intent(in) :: order
        character(:), allocatable :: text
        type(expr_t), allocatable :: roots(:)
        type(expr_t) :: grad(8)
        type(kernel_spec_t) :: spec
        type(engine_result_t) :: result
        integer :: k, i, j, nroot

        nroot = 1
        if (order >= 1) nroot = nroot + 8
        if (order >= 2) nroot = nroot + 36 + 64
        allocate (roots(nroot), spec%outputs(nroot))
        roots(1) = f
        spec%outputs(1) = str("value")
        k = 1
        if (order >= 1) then
            do i = 1, 8
                grad(i) = diff(f, vars(i))
                k = k + 1
                roots(k) = grad(i)
                spec%outputs(k) = str("g"//trim(self_names(i)))
            end do
        end if
        if (order >= 2) then
            do i = 1, 8
                do j = i, 8
                    k = k + 1
                    roots(k) = diff(grad(i), vars(j))
                    spec%outputs(k) = str("h"//trim(self_names(i))//"_"//trim(self_names(j)))
                end do
            end do
            do i = 1, 8
                do j = 9, 16
                    k = k + 1
                    roots(k) = diff(grad(i), vars(j))
                    spec%outputs(k) = str("h"//trim(self_names(i))//"_"//trim(self_names(j)))
                end do
            end do
        end if
        do k = 1, nroot
            result = engine%simplify(roots(k))
            if (result%ok) roots(k) = result%value
        end do
        spec%name = str(name)
        spec%module_name = str("fortfem_"//name)
        spec%mode = KERNEL_SUBROUTINE
        spec%generator = str("gen_toroidal_sheet_current_products")
        spec%generator_revision = str(fortsym_revision())
        spec%regenerate_command = str("cd tools/codegen && ./generate.sh")
        spec%pure_procedure = .true.
        spec%args = kargs
        text = chars(emit_kernel(roots, spec))
    end function pair_kernel

    function jet_kernel(f, vars, names, name, kargs, order, at_old, at_new, &
            keep_r) result(text)
        !! Value, gradient and upper-triangle Hessian (row-wise) of f in vars,
        !! up to order, as a kernel of kargs; at_old -> at_new after
        !! differentiation.
        type(expr_t), intent(in) :: f, vars(:)
        character(*), intent(in) :: names(:), name
        type(str_t), intent(in) :: kargs(:)
        integer, intent(in) :: order
        type(expr_t), intent(in), optional :: at_old(:), at_new(:)
        logical, intent(in), optional :: keep_r
        character(:), allocatable :: text
        type(expr_t), allocatable :: roots(:), grad(:)
        type(kernel_spec_t) :: spec
        type(engine_result_t) :: result
        integer :: n, k, i, j, nroot

        n = size(vars)
        nroot = 1
        if (order >= 1) nroot = nroot + n
        if (order >= 2) nroot = nroot + n*(n + 1)/2
        allocate (roots(nroot), grad(n))
        roots(1) = f
        k = 1
        if (order >= 1) then
            do i = 1, n
                grad(i) = diff(f, vars(i))
                k = k + 1
                roots(k) = grad(i)
            end do
        end if
        if (order >= 2) then
            do i = 1, n
                do j = i, n
                    k = k + 1
                    roots(k) = diff(grad(i), vars(j))
                end do
            end do
        end if
        if (present(at_old)) then
            do k = 1, nroot
                roots(k) = subs_many(roots(k), at_old, at_new)
                ! the load's R is both the Taylor point and an input
                if (present(keep_r)) roots(k) = subs_many(roots(k), &
                    [sym("r0")], [sym("r")])
            end do
        end if
        do k = 1, nroot
            result = engine%simplify(roots(k))
            if (result%ok) roots(k) = result%value
        end do
        spec%name = str(name)
        spec%module_name = str("fortfem_"//name)
        spec%mode = KERNEL_SUBROUTINE
        spec%generator = str("gen_toroidal_sheet_current_products")
        spec%generator_revision = str(fortsym_revision())
        spec%regenerate_command = str("cd tools/codegen && ./generate.sh")
        spec%pure_procedure = .true.
        spec%args = kargs
        allocate (spec%outputs(nroot))
        spec%outputs(1) = str("value")
        k = 1
        if (order >= 1) then
            do i = 1, n
                k = k + 1
                spec%outputs(k) = str("grad_"//trim(names(i)))
            end do
        end if
        if (order >= 2) then
            do i = 1, n
                do j = i, n
                    k = k + 1
                    spec%outputs(k) = str("hess_"//trim(names(i))//"_"//trim(names(j)))
                end do
            end do
        end if
        text = chars(emit_kernel(roots, spec))
    end function jet_kernel

end program gen_toroidal_sheet_current_products
