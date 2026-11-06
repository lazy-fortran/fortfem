module fortfem_meridian_bspline_n
    !! Per-toroidal-mode B-spline curl complex on a rectangular meridional
    !! patch [R0, R1] x [Z0, Z1] (R0 > 0), coupled to the exact exterior of
    !! each mode by the boundary operators of
    !! fortfem_toroidal_mode_boundary_operators_2d.
    !!
    !! Fields and spaces.  A field is u(R, Z) exp(i n phi).  The variables,
    !! operators and weights are those of fortfem_meridian_feec_n (the
    !! per-mode reduction of Albert, Lainer and Biro, Comput. Phys. Commun.
    !! 277 (2022) 108401), with tensor B-spline spaces in place of the
    !! Lagrange, Nedelec and Raviart-Thomas elements.  With S the degree-p
    !! splines on an open uniform knot vector and D the degree-(p-1)
    !! splines of its derivatives (Buffa, Rivas, Sangalli and Vazquez, SIAM
    !! J. Numer. Anal. 49 (2011)):
    !!
    !!   V0  f                    S x S
    !!   V1  A = (A_pol, a)       A_pol in (D x S, S x D), a = R A_phi in S x S
    !!   V2  B = (v, s)           v in (D x S, S x D), s = B_phi in D x D
    !!   V3  g = R div B          D x D
    !!
    !! where V2 stores v = grad a - i n A_pol, so q = R B_pol = rot90 v
    !! (rot90 v = (-v_Z, v_R)).  The operators are
    !!
    !!   grad_n f = (grad f, i n f),   curl_n A = (grad a - i n A_pol,
    !!   -rot A_pol),                  div_n B = -rot v + i n s,
    !!
    !! assembled from the commuting coefficient derivatives of
    !! fortfem_assembly_bspline_2d, so curl_n grad_n = 0 and
    !! div_n curl_n = 0 hold exactly for every n.  With degree p, A is
    !! C^{p-1} and B is C^{p-2} inside the patch.  The masses are the
    !! physical L2(R dR dZ) forms (2 pi dropped):
    !!
    !!   M0 = int R |f|^2,  M1 = int R |A_pol|^2 + |a|^2 / R,
    !!   M2 = int |v|^2 / R + R |s|^2,  M3 = int |g|^2 / R.
    !!
    !! On the Cartesian patch every weight depends on R only, so each mass
    !! block is the Kronecker product of a weighted 1D mass in R and an
    !! unweighted 1D mass in Z.
    !!
    !! Numbering.  S x S: i + (j - 1) mR; D x S: i + (j - 1)(mR - 1);
    !! S x D: ncx + i + (j - 1) mR; D x D: i + (j - 1)(mR - 1), where mR, mZ
    !! are the degree-p basis counts and ncx = (mR - 1) mZ.  V1 is
    !! [A_pol (nc), a (n0)], V2 is [v (nc), s (nd)].
    !!
    !! Exterior.  The patch boundary, counter-clockwise from (R0, Z0), is a
    !! closed curve of straight P2 panels on the knot spans.  For n /= 0 the
    !! trace map is the normal-flux pairing C = C2 curl_n,
    !! C2(i, j) = int psi_i (q_j . n) ds = -int psi_i (v_j . tau) ds; for
    !! n = 0 it is the Dirichlet-trace pairing of a (see
    !! fortfem_meridian_magnetic_coupling for both closures and for the
    !! Coulomb-gauge mixed solve).
    !!
    !! Loads.  The weak source of libneo's frozen-current problem,
    !!
    !!   f(A') = (B_s, curl_n A')_V + <n x B_s, A'>_dV,
    !!
    !! is assembled from caller-supplied points x_k = (R_k, Z_k), weights
    !! w_k and complex physical vectors (R, phi, Z):
    !!
    !!   load_curl:   f_i += sum_k w_k b_k . conj(curl_n A'_i(x_k)),
    !!   load_value:  f_i += sum_k w_k g_k . conj(A'_i(x_k)).
    !!
    !! The weights carry the measure of the per-mode form, R dR dZ (or
    !! R ds on a boundary curve).  A source sampled in three dimensions at
    !! (R_k, phi_k, Z_k) with volume weights W_k enters mode n with
    !! w_k b_k = W_k B_s(x_k) exp(-i n phi_k) / (2 pi).  The points need
    !! not follow the patch cells but must lie in the patch.
    use fortfem_kinds, only: dp
    use fortfem_bspline_feec, only: evaluate_bspline_basis
    use fortfem_assembly_bspline_2d, only: build_bspline_feec_2d_operators_csc
    use fortfem_generated_meridian_bspline_physical, only: &
        generated_meridian_bspline_physical
    use fortfem_meridian_feec_n, only: meridian_op_t, meridian_op_complex, &
        meridian_zcsc_matmul, meridian_real_to_complex
    use fortfem_meridian_magnetic_coupling, only: meridian_exterior_t, &
        meridian_exterior_form, meridian_exterior_solve
    use fortfem_axisymmetric_gs_boundary_operators_2d, only: &
        gs_bem_curve_from_nodes, p2_basis
    use fortnum_quadrature, only: gauss_legendre_ab
    use fortsparse, only: csc_t, csc_z_t, csc_from_triplet, &
        csc_conjugate_transpose, fortsparse_status_t, status_ok
    implicit none
    private

    type, public :: meridian_bspline_t
        integer :: p = 0 !! degree of V0
        integer :: mr = 0, mz = 0 !! degree-p basis counts
        real(dp) :: r0 = 0.0_dp, r1 = 0.0_dp, z0 = 0.0_dp, z1 = 0.0_dp
        real(dp), allocatable :: kr(:), kz(:) !! open knot vectors
        integer :: n0 = 0, ncx = 0, nc = 0, nd = 0, n1 = 0, n2 = 0
        type(meridian_op_t) :: grad !! V0 -> V1
        type(meridian_op_t) :: curl !! V1 -> V2
        type(meridian_op_t) :: div !! V2 -> V3
        type(csc_t) :: mass(0:3) !! physical L2(R dR dZ)
    end type meridian_bspline_t

    type :: triplets_t
        integer :: n = 0
        integer, allocatable :: r(:), c(:)
        real(dp), allocatable :: v(:)
    end type triplets_t

    public :: meridian_bspline_build
    public :: meridian_bspline_local_count
    public :: meridian_bspline_v1_point
    public :: meridian_bspline_stiffness
    public :: meridian_bspline_gauge
    public :: meridian_bspline_exterior_build
    public :: meridian_bspline_solve
    public :: meridian_bspline_field_energy
    public :: meridian_bspline_load_curl
    public :: meridian_bspline_load_value
    public :: meridian_bspline_eval

    !> Point counts below which load and evaluation loops stay serial.
    integer, parameter :: par_min = 2048

contains

    subroutine meridian_bspline_build(r0, r1, z0, z1, nr, nz, p, sp, status)
        !! Complex of degree p on [r0, r1] x [z0, z1] with nr x nz uniform
        !! knot spans.
        real(dp), intent(in) :: r0, r1, z0, z1
        integer, intent(in) :: nr, nz, p
        type(meridian_bspline_t), intent(out) :: sp
        integer, intent(out) :: status

        type(csc_t) :: g, cu, rss(0:2), rdd(0:2), zss, zdd
        type(fortsparse_status_t) :: fst
        type(triplets_t) :: t
        integer :: st, w

        status = 1
        if (p < 1 .or. nr < 1 .or. nz < 1) return
        if (.not. (r0 > 0.0_dp .and. r1 > r0 .and. z1 > z0)) return
        sp%p = p
        sp%r0 = r0
        sp%r1 = r1
        sp%z0 = z0
        sp%z1 = z1
        call open_knots(r0, r1, nr, p, sp%kr)
        call open_knots(z0, z1, nz, p, sp%kz)
        sp%mr = nr + p
        sp%mz = nz + p
        sp%n0 = sp%mr*sp%mz
        sp%ncx = (sp%mr - 1)*sp%mz
        sp%nc = sp%ncx + sp%mr*(sp%mz - 1)
        sp%nd = (sp%mr - 1)*(sp%mz - 1)
        sp%n1 = sp%nc + sp%n0
        sp%n2 = sp%nc + sp%nd
        call build_bspline_feec_2d_operators_csc(sp%kr, sp%kz, p, p, g, cu, fst)
        if (.not. status_ok(fst)) return

        associate (n0 => sp%n0, nc => sp%nc, nd => sp%nd)
            ! grad_n f = (G f, i n f)
            t%n = 0
            call put(t, g, 0, 0, 1.0_dp)
            call finish(t, nc + n0, n0, sp%grad%a, st)
            if (st /= 0) return
            call put_identity(t, n0, nc, 0, 1.0_dp)
            call finish(t, nc + n0, n0, sp%grad%b, st)
            if (st /= 0) return
            ! curl_n (A_pol, a) = (G a - i n A_pol, -Cu A_pol)
            call put(t, g, 0, nc, 1.0_dp)
            call put(t, cu, nc, 0, -1.0_dp)
            call finish(t, nc + nd, nc + n0, sp%curl%a, st)
            if (st /= 0) return
            call put_identity(t, nc, 0, 0, -1.0_dp)
            call finish(t, nc + nd, nc + n0, sp%curl%b, st)
            if (st /= 0) return
            ! div_n (v, s) = -Cu v + i n s
            call put(t, cu, 0, 0, -1.0_dp)
            call finish(t, nd, nc + nd, sp%div%a, st)
            if (st /= 0) return
            call put_identity(t, nd, 0, nc, 1.0_dp)
            call finish(t, nd, nc + nd, sp%div%b, st)
            if (st /= 0) return
        end associate

        ! 1D masses: weight 0 = 1, 1 = R, 2 = 1/R
        do w = 1, 2
            call mass_1d(sp%kr, p, 0, w, rss(w), st)
            if (st /= 0) return
            call mass_1d(sp%kr, p, 1, w, rdd(w), st)
            if (st /= 0) return
        end do
        call mass_1d(sp%kz, p, 0, 0, zss, st)
        if (st /= 0) return
        call mass_1d(sp%kz, p, 1, 0, zdd, st)
        if (st /= 0) return

        call put_kron(t, zss, rss(1), 0)
        call finish(t, sp%n0, sp%n0, sp%mass(0), st)
        if (st /= 0) return
        call put_kron(t, zss, rdd(1), 0)
        call put_kron(t, zdd, rss(1), sp%ncx)
        call put_kron(t, zss, rss(2), sp%nc)
        call finish(t, sp%n1, sp%n1, sp%mass(1), st)
        if (st /= 0) return
        call put_kron(t, zss, rdd(2), 0)
        call put_kron(t, zdd, rss(2), sp%ncx)
        call put_kron(t, zdd, rdd(1), sp%nc)
        call finish(t, sp%n2, sp%n2, sp%mass(2), st)
        if (st /= 0) return
        call put_kron(t, zdd, rdd(2), 0)
        call finish(t, sp%nd, sp%nd, sp%mass(3), st)
        if (st /= 0) return
        status = 0
    end subroutine meridian_bspline_build

    pure integer function meridian_bspline_local_count(sp) result(nl)
        !! Number of V1 basis functions that do not vanish on one cell.
        type(meridian_bspline_t), intent(in) :: sp

        nl = 2*sp%p*(sp%p + 1) + (sp%p + 1)**2
    end function meridian_bspline_local_count

    subroutine meridian_bspline_v1_point(sp, n, x, idx, aval, bval, status)
        !! Global dofs idx(nl) of the V1 basis functions of the cell that
        !! contains x = (R, Z), their physical values A (aval(3, nl)) and
        !! physical curls curl_n A (bval(3, nl)), components (R, phi, Z).
        type(meridian_bspline_t), intent(in) :: sp
        integer, intent(in) :: n
        real(dp), intent(in) :: x(2)
        integer, intent(out) :: idx(:)
        complex(dp), intent(out) :: aval(:, :), bval(:, :)
        integer, intent(out) :: status

        real(dp) :: svr(0:sp%p), sdr(0:sp%p), dvr(sp%p), ddr(sp%p)
        real(dp) :: svz(0:sp%p), sdz(0:sp%p), dvz(sp%p), ddz(sp%p)
        real(dp) :: aphi, brr, bri, bph, bzr, bzi, rn, tol
        integer :: p, sr, sz, i, j, l

        status = 1
        p = sp%p
        tol = 1.0e-12_dp*max(sp%r1 - sp%r0, sp%z1 - sp%z0)
        if (x(1) < sp%r0 - tol .or. x(1) > sp%r1 + tol) return
        if (x(2) < sp%z0 - tol .or. x(2) > sp%z1 + tol) return
        sr = find_span(sp%kr, p, sp%mr, x(1))
        sz = find_span(sp%kz, p, sp%mz, x(2))
        call local_basis(sp%kr, p, sr, min(max(x(1), sp%r0), sp%r1), svr, &
            sdr, dvr, ddr, status)
        if (status /= 0) return
        call local_basis(sp%kz, p, sz, min(max(x(2), sp%z0), sp%z1), svz, &
            sdz, dvz, ddz, status)
        if (status /= 0) return
        rn = real(n, dp)
        l = 0
        do j = 0, p
            do i = 1, p
                l = l + 1
                idx(l) = sr - p - 1 + i + (sz - p + j - 1)*(sp%mr - 1)
                call generated_meridian_bspline_physical(x(1), rn, &
                    dvr(i)*svz(j), dvr(i)*sdz(j), 0.0_dp, 0.0_dp, 0.0_dp, &
                    0.0_dp, 0.0_dp, aphi, brr, bri, bph, bzr, bzi)
                aval(:, l) = [cmplx(dvr(i)*svz(j), 0.0_dp, dp), &
                    (0.0_dp, 0.0_dp), (0.0_dp, 0.0_dp)]
                bval(:, l) = [cmplx(brr, bri, dp), cmplx(bph, 0.0_dp, dp), &
                    cmplx(bzr, bzi, dp)]
            end do
        end do
        do j = 1, p
            do i = 0, p
                l = l + 1
                idx(l) = sp%ncx + sr - p + i + (sz - p + j - 2)*sp%mr
                call generated_meridian_bspline_physical(x(1), rn, 0.0_dp, &
                    0.0_dp, svr(i)*dvz(j), sdr(i)*dvz(j), 0.0_dp, 0.0_dp, &
                    0.0_dp, aphi, brr, bri, bph, bzr, bzi)
                aval(:, l) = [(0.0_dp, 0.0_dp), (0.0_dp, 0.0_dp), &
                    cmplx(svr(i)*dvz(j), 0.0_dp, dp)]
                bval(:, l) = [cmplx(brr, bri, dp), cmplx(bph, 0.0_dp, dp), &
                    cmplx(bzr, bzi, dp)]
            end do
        end do
        do j = 0, p
            do i = 0, p
                l = l + 1
                idx(l) = sp%nc + sr - p + i + (sz - p + j - 1)*sp%mr
                call generated_meridian_bspline_physical(x(1), rn, 0.0_dp, &
                    0.0_dp, 0.0_dp, 0.0_dp, svr(i)*svz(j), sdr(i)*svz(j), &
                    svr(i)*sdz(j), aphi, brr, bri, bph, bzr, bzi)
                aval(:, l) = [(0.0_dp, 0.0_dp), cmplx(aphi, 0.0_dp, dp), &
                    (0.0_dp, 0.0_dp)]
                bval(:, l) = [cmplx(brr, bri, dp), cmplx(bph, 0.0_dp, dp), &
                    cmplx(bzr, bzi, dp)]
            end do
        end do
        status = 0
    end subroutine meridian_bspline_v1_point

    subroutine meridian_bspline_stiffness(sp, n, k, status)
        !! Curl-curl form K = curl_n^H M2 curl_n.
        type(meridian_bspline_t), intent(in) :: sp
        integer, intent(in) :: n
        type(csc_z_t), intent(out) :: k
        integer, intent(out) :: status

        type(csc_z_t) :: z, zh, m, mz
        type(fortsparse_status_t) :: fst

        call meridian_op_complex(sp%curl, n, z, status)
        if (status /= 0) return
        call meridian_real_to_complex(sp%mass(2), m)
        call meridian_zcsc_matmul(m, z, mz, status)
        if (status /= 0) return
        call csc_conjugate_transpose(z, zh, fst)
        status = 1
        if (.not. status_ok(fst)) return
        call meridian_zcsc_matmul(zh, mz, k, status)
    end subroutine meridian_bspline_stiffness

    subroutine meridian_bspline_gauge(sp, n, g, status)
        !! Gauge coupling M1 grad_n (the constraint is g^H A = 0).
        type(meridian_bspline_t), intent(in) :: sp
        integer, intent(in) :: n
        type(csc_z_t), intent(out) :: g
        integer, intent(out) :: status

        type(csc_z_t) :: z, m

        call meridian_op_complex(sp%grad, n, z, status)
        if (status /= 0) return
        call meridian_real_to_complex(sp%mass(1), m)
        call meridian_zcsc_matmul(m, z, g, status)
    end subroutine meridian_bspline_gauge

    subroutine meridian_bspline_exterior_build(sp, n, ex, status, &
            panels_per_span)
        !! Exterior closure of mode n on the patch boundary, with
        !! panels_per_span (default 1) straight P2 panels per knot span.
        type(meridian_bspline_t), intent(in) :: sp
        integer, intent(in) :: n
        type(meridian_exterior_t), intent(out) :: ex
        integer, intent(out) :: status
        integer, intent(in), optional :: panels_per_span

        type(csc_z_t) :: c2, curl
        type(triplets_t) :: t
        real(dp), allocatable :: nodes(:, :), xi(:), w(:), loc(:, :)
        complex(dp), allocatable :: aval(:, :), bval(:, :)
        integer, allocatable :: idx(:)
        real(dp) :: xa(2), xb(2), x(2), tau(2), len, psi(3), dpsi(3), val
        integer :: pps, np, nb, side, s, k, panel, q, l, a, row, nl, nq, st
        integer :: first, last, step

        status = 1
        pps = 1
        if (present(panels_per_span)) pps = panels_per_span
        if (pps < 1) return
        np = 2*pps*(sp%mr - sp%p + sp%mz - sp%p)
        nb = 2*np
        allocate (nodes(2, nb))
        panel = 0
        do side = 1, 4
            select case (side)
            case (1, 2)
                first = sp%p + 1
                step = 1
            case default
                first = merge(sp%mr, sp%mz, side == 3)
                step = -1
            end select
            last = merge(sp%mr, sp%mz, side == 1 .or. side == 3)
            if (step < 0) last = sp%p + 1
            do s = first, last, step
                do k = 0, pps - 1
                    call side_point(sp, side, s, real(k, dp)/pps, xa)
                    call side_point(sp, side, s, real(k + 1, dp)/pps, xb)
                    panel = panel + 1
                    nodes(:, 2*panel - 1) = xa
                    nodes(:, 2*panel) = 0.5_dp*(xa + xb)
                end do
            end do
        end do
        ex%n = n
        call gs_bem_curve_from_nodes(nodes, ex%curve, st)
        if (st /= 0) return

        nl = meridian_bspline_local_count(sp)
        nq = sp%p + 3
        allocate (xi(nq), w(nq), idx(nl), aval(3, nl), bval(3, nl))
        allocate (loc(3, nl))
        call gauss_legendre_ab(nq, 0.0_dp, 1.0_dp, xi, w)
        t%n = 0
        do panel = 1, np
            xa = nodes(:, 2*panel - 1)
            xb = nodes(:, mod(2*panel, nb) + 1)
            len = norm2(xb - xa)
            tau = (xb - xa)/len
            loc = 0.0_dp
            ! all quadrature points of a panel lie in one cell, so idx is
            ! the same at each of them
            do q = 1, nq
                x = xa + xi(q)*(xb - xa)
                call meridian_bspline_v1_point(sp, n, x, idx, aval, bval, st)
                if (st /= 0) return
                call p2_basis(xi(q), psi, dpsi)
                do l = 1, nl
                    if (n /= 0) then
                        ! q . n = -v . tau, v = A_pol on the V1 layout
                        val = -(real(aval(1, l), dp)*tau(1) &
                            + real(aval(3, l), dp)*tau(2))
                    else
                        val = real(aval(2, l), dp)*x(1)
                    end if
                    loc(:, l) = loc(:, l) + w(q)*len*psi*val
                end do
            end do
            do l = 1, nl
                do a = 1, 3
                    if (loc(a, l) == 0.0_dp) cycle
                    row = 2*panel - 2 + a
                    if (row > nb) row = 1
                    if (n /= 0 .and. idx(l) > sp%nc) cycle
                    call push(t, row, idx(l), loc(a, l))
                end do
            end do
        end do
        if (n /= 0) then
            call finish_complex(t, nb, sp%n2, c2, st)
            if (st /= 0) return
            call meridian_op_complex(sp%curl, n, curl, st)
            if (st /= 0) return
            call meridian_zcsc_matmul(c2, curl, ex%trace, st)
            if (st /= 0) return
            call drop_zeros(ex%trace)
        else
            call finish_complex(t, nb, sp%n1, ex%trace, st)
            if (st /= 0) return
        end if
        call meridian_exterior_form(ex%curve, n, ex%form, status)
    end subroutine meridian_bspline_exterior_build

    subroutine meridian_bspline_solve(sp, ex, f, u, p, status)
        !! Coulomb-gauge FEM-BEM solve of mode ex%n for the load f(n1); u
        !! are the V1 coefficients, p the gauge multiplier in V0.
        type(meridian_bspline_t), intent(in) :: sp
        type(meridian_exterior_t), intent(in) :: ex
        complex(dp), intent(in) :: f(:)
        complex(dp), allocatable, intent(out) :: u(:), p(:)
        integer, intent(out) :: status

        type(csc_z_t) :: k, g

        call meridian_bspline_stiffness(sp, ex%n, k, status)
        if (status /= 0) return
        call meridian_bspline_gauge(sp, ex%n, g, status)
        if (status /= 0) return
        call meridian_exterior_solve(k, g, ex, f, u, p, status)
    end subroutine meridian_bspline_solve

    subroutine meridian_bspline_field_energy(sp, n, u, energy, status)
        !! Interior field energy (curl_n u, curl_n u)_M2.
        type(meridian_bspline_t), intent(in) :: sp
        integer, intent(in) :: n
        complex(dp), intent(in) :: u(:)
        real(dp), intent(out) :: energy
        integer, intent(out) :: status

        type(csc_z_t) :: k
        complex(dp), allocatable :: ku(:)
        integer :: j, e

        energy = 0.0_dp
        call meridian_bspline_stiffness(sp, n, k, status)
        if (status /= 0) return
        allocate (ku(k%nrow))
        ku = (0.0_dp, 0.0_dp)
        do j = 1, k%ncol
            do e = k%col_ptr(j), k%col_ptr(j + 1) - 1
                ku(k%row_idx(e)) = ku(k%row_idx(e)) + k%val(e)*u(j)
            end do
        end do
        energy = real(dot_product(u, ku), dp)
    end subroutine meridian_bspline_field_energy

    subroutine meridian_bspline_load_curl(sp, n, x, w, b, f, status)
        !! f_i += sum_k w_k b(:, k) . conj(curl_n A'_i(x(:, k))).
        type(meridian_bspline_t), intent(in) :: sp
        integer, intent(in) :: n
        real(dp), intent(in) :: x(:, :), w(:)
        complex(dp), intent(in) :: b(:, :)
        complex(dp), intent(inout) :: f(:)
        integer, intent(out) :: status

        call accumulate_load(sp, n, x, w, b, f, .true., status)
    end subroutine meridian_bspline_load_curl

    subroutine meridian_bspline_load_value(sp, n, x, w, g, f, status)
        !! f_i += sum_k w_k g(:, k) . conj(A'_i(x(:, k))); with g = n x B_s
        !! on a boundary quadrature this is the surface term of the source.
        type(meridian_bspline_t), intent(in) :: sp
        integer, intent(in) :: n
        real(dp), intent(in) :: x(:, :), w(:)
        complex(dp), intent(in) :: g(:, :)
        complex(dp), intent(inout) :: f(:)
        integer, intent(out) :: status

        call accumulate_load(sp, n, x, w, g, f, .false., status)
    end subroutine meridian_bspline_load_value

    subroutine meridian_bspline_eval(sp, n, u, x, a, b, status)
        !! Physical A (a(3, :)) and B = curl_n A (b(3, :)) of the V1 field u
        !! at the points x(2, :), components (R, phi, Z).
        type(meridian_bspline_t), intent(in) :: sp
        integer, intent(in) :: n
        complex(dp), intent(in) :: u(:)
        real(dp), intent(in) :: x(:, :)
        complex(dp), intent(out) :: a(:, :), b(:, :)
        integer, intent(out) :: status

        integer :: idx(meridian_bspline_local_count(sp))
        complex(dp) :: av(3, size(idx)), bv(3, size(idx))
        integer :: k, l, st, bad

        bad = 0
        !$omp parallel do if(size(x, 2) >= par_min) default(shared) &
        !$omp private(k, l, st, idx, av, bv) reduction(max:bad)
        do k = 1, size(x, 2)
            a(:, k) = (0.0_dp, 0.0_dp)
            b(:, k) = (0.0_dp, 0.0_dp)
            call meridian_bspline_v1_point(sp, n, x(:, k), idx, av, bv, st)
            if (st /= 0) then
                bad = 1
                cycle
            end if
            do l = 1, size(idx)
                a(:, k) = a(:, k) + u(idx(l))*av(:, l)
                b(:, k) = b(:, k) + u(idx(l))*bv(:, l)
            end do
        end do
        !$omp end parallel do
        status = bad
    end subroutine meridian_bspline_eval

    ! ------------------------------------------------------------------
    ! Internal helpers
    ! ------------------------------------------------------------------

    subroutine accumulate_load(sp, n, x, w, v, f, curl, status)
        type(meridian_bspline_t), intent(in) :: sp
        integer, intent(in) :: n
        real(dp), intent(in) :: x(:, :), w(:)
        complex(dp), intent(in) :: v(:, :)
        complex(dp), intent(inout) :: f(:)
        logical, intent(in) :: curl
        integer, intent(out) :: status

        integer :: idx(meridian_bspline_local_count(sp))
        complex(dp) :: av(3, size(idx)), bv(3, size(idx))
        integer :: k, l, st, bad

        status = 1
        if (size(f) /= sp%n1 .or. size(w) /= size(x, 2)) return
        if (size(v, 2) /= size(x, 2) .or. size(v, 1) /= 3) return
        bad = 0
        !$omp parallel do if(size(x, 2) >= par_min) default(shared) &
        !$omp private(k, l, st, idx, av, bv) reduction(max:bad) &
        !$omp reduction(+:f)
        do k = 1, size(x, 2)
            call meridian_bspline_v1_point(sp, n, x(:, k), idx, av, bv, st)
            if (st /= 0) then
                bad = 1
                cycle
            end if
            if (curl) av = bv
            do l = 1, size(idx)
                f(idx(l)) = f(idx(l)) + w(k)*sum(v(:, k)*conjg(av(:, l)))
            end do
        end do
        !$omp end parallel do
        status = bad
    end subroutine accumulate_load

    subroutine open_knots(a, b, n, p, k)
        real(dp), intent(in) :: a, b
        integer, intent(in) :: n, p
        real(dp), allocatable, intent(out) :: k(:)
        integer :: i

        allocate (k(n + 2*p + 1))
        k(:p + 1) = a
        do i = 1, n - 1
            k(p + 1 + i) = a + (b - a)*real(i, dp)/real(n, dp)
        end do
        k(n + p + 1:) = b
    end subroutine open_knots

    pure integer function find_span(knots, p, m, x) result(s)
        !! Span s in p+1 .. m with knots(s) <= x < knots(s+1), clamped.
        real(dp), intent(in) :: knots(:), x
        integer, intent(in) :: p, m
        integer :: lo, hi, mid

        lo = p + 1
        hi = m
        do while (hi > lo)
            mid = (lo + hi + 1)/2
            if (knots(mid) <= x) then
                lo = mid
            else
                hi = mid - 1
            end if
        end do
        s = lo
    end function find_span

    subroutine local_basis(knots, p, s, x, sv, sd, dv, dd, status)
        !! Degree-p values and derivatives sv, sd (functions s-p .. s) and
        !! degree-(p-1) ones dv, dd (reduced functions s-p .. s-1) on span s.
        real(dp), intent(in) :: knots(:), x
        integer, intent(in) :: p, s
        real(dp), intent(out) :: sv(0:p), sd(0:p), dv(p), dd(p)
        integer, intent(out) :: status
        real(dp), allocatable :: v(:), d(:)

        call evaluate_bspline_basis(knots(s - p:s + p + 1), p, x, v, d, status)
        if (status /= 0) return
        sv = v
        sd = d
        call evaluate_bspline_basis(knots(s - p + 1:s + p), p - 1, x, v, d, &
            status)
        if (status /= 0) return
        dv = v
        dd = d
    end subroutine local_basis

    subroutine mass_1d(knots, p, kind, weight, m, status)
        !! 1D mass of the degree-p (kind 0) or degree-(p-1) (kind 1) splines
        !! with weight 1 (0), x (1) or 1/x (2).
        real(dp), intent(in) :: knots(:)
        integer, intent(in) :: p, kind, weight
        type(csc_t), intent(out) :: m
        integer, intent(out) :: status

        type(triplets_t) :: t
        real(dp), allocatable :: xq(:), wq(:)
        real(dp) :: sv(0:p), sd(0:p), dv(p), dd(p), wt
        real(dp), allocatable :: f(:)
        integer :: nb, s, q, i, j, nq, nf, off

        nb = size(knots) - p - 1
        nq = p + 4
        allocate (xq(nq), wq(nq))
        nf = p + 1 - kind
        allocate (f(nf))
        t%n = 0
        do s = p + 1, nb
            if (knots(s + 1) <= knots(s)) cycle
            call gauss_legendre_ab(nq, knots(s), knots(s + 1), xq, wq)
            off = s - p - 1
            do q = 1, nq
                call local_basis(knots, p, s, xq(q), sv, sd, dv, dd, status)
                if (status /= 0) return
                if (kind == 0) then
                    f = sv
                else
                    f = dv
                end if
                select case (weight)
                case (1)
                    wt = wq(q)*xq(q)
                case (2)
                    wt = wq(q)/xq(q)
                case default
                    wt = wq(q)
                end select
                do j = 1, nf
                    do i = 1, nf
                        call push(t, off + i, off + j, wt*f(i)*f(j))
                    end do
                end do
            end do
        end do
        call finish(t, nb - kind, nb - kind, m, status)
    end subroutine mass_1d

    subroutine side_point(sp, side, s, frac, x)
        !! Point at fraction frac of knot span s along side 1 (bottom,
        !! increasing R), 2 (right, increasing Z), 3 (top, decreasing R) or
        !! 4 (left, decreasing Z).
        type(meridian_bspline_t), intent(in) :: sp
        integer, intent(in) :: side, s
        real(dp), intent(in) :: frac
        real(dp), intent(out) :: x(2)

        select case (side)
        case (1)
            x = [sp%kr(s) + frac*(sp%kr(s + 1) - sp%kr(s)), sp%z0]
        case (2)
            x = [sp%r1, sp%kz(s) + frac*(sp%kz(s + 1) - sp%kz(s))]
        case (3)
            x = [sp%kr(s + 1) - frac*(sp%kr(s + 1) - sp%kr(s)), sp%z1]
        case default
            x = [sp%r0, sp%kz(s + 1) - frac*(sp%kz(s + 1) - sp%kz(s))]
        end select
    end subroutine side_point

    subroutine push(t, r, c, v)
        type(triplets_t), intent(inout) :: t
        integer, intent(in) :: r, c
        real(dp), intent(in) :: v
        integer, allocatable :: ri(:), ci(:)
        real(dp), allocatable :: vi(:)

        if (.not. allocated(t%r)) then
            allocate (t%r(1024), t%c(1024), t%v(1024))
        else if (t%n == size(t%r)) then
            allocate (ri(2*t%n), ci(2*t%n), vi(2*t%n))
            ri(:t%n) = t%r
            ci(:t%n) = t%c
            vi(:t%n) = t%v
            call move_alloc(ri, t%r)
            call move_alloc(ci, t%c)
            call move_alloc(vi, t%v)
        end if
        t%n = t%n + 1
        t%r(t%n) = r
        t%c(t%n) = c
        t%v(t%n) = v
    end subroutine push

    subroutine put(t, a, roff, coff, scale)
        type(triplets_t), intent(inout) :: t
        type(csc_t), intent(in) :: a
        integer, intent(in) :: roff, coff
        real(dp), intent(in) :: scale
        integer :: j, e

        do j = 1, a%ncol
            do e = a%col_ptr(j), a%col_ptr(j + 1) - 1
                call push(t, roff + a%row_idx(e), coff + j, scale*a%val(e))
            end do
        end do
    end subroutine put

    subroutine put_identity(t, n, roff, coff, scale)
        type(triplets_t), intent(inout) :: t
        integer, intent(in) :: n, roff, coff
        real(dp), intent(in) :: scale
        integer :: i

        do i = 1, n
            call push(t, roff + i, coff + i, scale)
        end do
    end subroutine put_identity

    subroutine put_kron(t, az, ar, off)
        !! Block az (x) ar at diagonal offset off, R index fastest.
        type(triplets_t), intent(inout) :: t
        type(csc_t), intent(in) :: az, ar
        integer, intent(in) :: off
        integer :: jz, ez, jr, er

        do jz = 1, az%ncol
            do ez = az%col_ptr(jz), az%col_ptr(jz + 1) - 1
                do jr = 1, ar%ncol
                    do er = ar%col_ptr(jr), ar%col_ptr(jr + 1) - 1
                        call push(t, off + ar%row_idx(er) &
                            + (az%row_idx(ez) - 1)*ar%nrow, &
                            off + jr + (jz - 1)*ar%ncol, &
                            az%val(ez)*ar%val(er))
                    end do
                end do
            end do
        end do
    end subroutine put_kron

    subroutine finish(t, nrow, ncol, a, status)
        type(triplets_t), intent(inout) :: t
        integer, intent(in) :: nrow, ncol
        type(csc_t), intent(out) :: a
        integer, intent(out) :: status
        type(fortsparse_status_t) :: fst

        if (.not. allocated(t%r)) allocate (t%r(0), t%c(0), t%v(0))
        call csc_from_triplet(nrow, ncol, t%r(:t%n), t%c(:t%n), t%v(:t%n), &
            a, fst)
        t%n = 0
        status = 0
        if (.not. status_ok(fst)) status = 1
    end subroutine finish

    subroutine finish_complex(t, nrow, ncol, a, status)
        type(triplets_t), intent(inout) :: t
        integer, intent(in) :: nrow, ncol
        type(csc_z_t), intent(out) :: a
        integer, intent(out) :: status
        type(csc_t) :: ar

        call finish(t, nrow, ncol, ar, status)
        if (status /= 0) return
        call meridian_real_to_complex(ar, a)
    end subroutine finish_complex

    subroutine drop_zeros(a)
        !! Remove exact zeros left by cancellation in a sparse product.
        type(csc_z_t), intent(inout) :: a
        integer :: j, e, pos, start

        pos = 0
        do j = 1, a%ncol
            start = a%col_ptr(j)
            a%col_ptr(j) = pos + 1
            do e = start, a%col_ptr(j + 1) - 1
                if (a%val(e) == (0.0_dp, 0.0_dp)) cycle
                pos = pos + 1
                a%row_idx(pos) = a%row_idx(e)
                a%val(pos) = a%val(e)
            end do
        end do
        a%col_ptr(a%ncol + 1) = pos + 1
        a%nnz = pos
    end subroutine drop_zeros

end module fortfem_meridian_bspline_n
