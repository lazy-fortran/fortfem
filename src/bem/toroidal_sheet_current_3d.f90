module fortfem_toroidal_sheet_current_3d
    !! Magnetic energy of a divergence-free surface current on a toroidal
    !! surface, with exact first and second derivatives in the surface and
    !! current coefficients.
    !!
    !! Surface: X(u, v) = (R cos v, R sin v, Z), (R, Z)(u, v) = sum_i sum_m
    !! rz(:, i, m) N_i(u) T_m(v), N_i the P_k Lagrange basis of a closed chain
    !! of panels in u (panel p has nodes k(p-1)+1, ..., k(p-1)+k+1, the last
    !! wrapping to node 1), T_m the real Fourier basis in v with field-period
    !! symmetry: T_1 = 1, T_2n = cos(n nfp v), T_2n+1 = sin(n nfp v).
    !! Current: K = n x grad_s Phi, Phi = phi + i_t u_c + i_p v/(2 pi), phi =
    !! sum phi(i, m) N_i T_m (own mode count), u_c = (p - 1 + xi)/n_panels the
    !! chain parameter, i_t and i_p the net currents (fixed, in the units of
    !! phi). In the panel parameters K dA = Phi_xi X_v - Phi_v X_xi dxi dv.
    !!
    !!   E = (1/(8 pi)) int int (K dA)_x . (K dA)_y / |X_x - X_y|.
    !!
    !! Quadrature, fixed in the parameters (so the derivatives are those of the
    !! discrete energy), over patch pairs (panel x toroidal interval), targets
    !! on one field period and E = nfp times their sum. Separated pairs:
    !! tensor Gauss at both points, the order chosen from the distance.
    !! Touching or equal pairs: relative coordinates z = y - x, the pair
    !! domain split into rectangles in z, polar about z = 0 for those with a
    !! corner there (the radius cancels 1/|z|) and Gauss otherwise, with x
    !! integrated over the overlap of the two patches. Every integrand input is
    !! linear in the coefficients at fixed parameters, and the integrand jets
    !! are generated (tools/codegen gen_toroidal_sheet_current_products).
    !!
    !! Coefficient vector c: geometry (R then Z, node fastest, then mode),
    !! then phi (node fastest, then mode).
    use fortfem_kinds, only: dp
    use fortfem_axisymmetric_gs_boundary_operators_2d, only: gauss_unit
    use fortfem_generated_toroidal_sheet_self_value, only: &
        generated_toroidal_sheet_self_value
    use fortfem_generated_toroidal_sheet_self_gradient, only: &
        generated_toroidal_sheet_self_gradient
    use fortfem_generated_toroidal_sheet_self_jet, only: &
        generated_toroidal_sheet_self_jet
    use fortfem_generated_toroidal_sheet_load_jet, only: &
        generated_toroidal_sheet_load_jet
    implicit none
    private

    real(dp), parameter :: pi = 3.14159265358979323846264338327950288_dp
    integer, parameter :: n_self = 16

    type, public :: ts_surface_t
        integer :: degree = 0, n_panels = 0, n_nodes = 0, nfp = 1
        integer :: ntor = 0, nmode = 1, ntor_phi = 0, nmode_phi = 1
        real(dp), allocatable :: lagrange(:, :)
        real(dp), allocatable :: rz(:, :, :) !! (2, n_nodes, nmode)
        real(dp), allocatable :: phi(:, :) !! (n_nodes, nmode_phi)
        real(dp) :: i_t = 0, i_p = 0
    end type

    type, public :: ts_plan_t
        !> nv toroidal intervals per period; ng_u Gauss points per direction
        !> of the near rules (x over the overlap, the polar rule's radius and
        !> angle); separated patches take n points per direction, n from
        !> 2 to ng_far chosen by their distance for a relative error near
        !> far_tolerance.
        integer :: nv = 4, ng_u = 4, ng_v = 4, ng_far = 8, ng_t = 6, ng_s = 6
        real(dp) :: far_tolerance = 1.0e-11_dp
    end type

    !> Separated-patch quadrature points of every patch for each order.
    type :: far_points_t
        integer :: n = 0
        real(dp), allocatable :: xi(:, :), v(:, :), w(:, :) !! (n*n, patch)
        real(dp), allocatable :: u(:, :, :) !! (8, n*n, patch)
        !> basis values at the points: N, dN/dxi (k+1), T, dT/dv (geometry
        !> modes), T, dT/dv (current modes)
        real(dp), allocatable :: nb(:, :, :), dnb(:, :, :), tt(:, :, :), dtt(:, :, :)
        real(dp), allocatable :: tp(:, :, :), dtp(:, :, :)
    end type

    !> Local coefficients of each panel: values and global indices.
    type :: frame_t
        integer :: nlg = 0, nlp = 0, nl = 0
        real(dp), allocatable :: cl(:, :) !! (nl, panel)
        integer, allocatable :: gidx(:, :) !! (nl, panel)
    end type

    !> Linear functional of the coefficients: sum coef(j) c(idx(j)) + shift.
    type :: functional_t
        integer :: n = 0
        integer, allocatable :: idx(:)
        real(dp), allocatable :: coef(:)
        real(dp) :: shift = 0
    end type

    abstract interface
        subroutine ts_field(r, z, v, b, db, d2b)
            !! (B_R, B_phi, B_Z) at (R, Z, phi = v), derivatives in (R, Z)
            !! and second derivatives (RR, RZ, ZZ).
            import :: dp
            real(dp), intent(in) :: r, z, v
            real(dp), intent(out) :: b(3), db(3, 2), d2b(3, 3)
        end subroutine ts_field
    end interface

    public :: ts_field
    public :: ts_surface_create, ts_size, ts_coefficients, ts_set_coefficients
    public :: ts_energy, ts_load, ts_point

contains

    subroutine ts_surface_create(degree, n_panels, nfp, ntor, ntor_phi, s, status)
        integer, intent(in) :: degree, n_panels, nfp, ntor, ntor_phi
        type(ts_surface_t), intent(out) :: s
        integer, intent(out) :: status
        real(dp) :: vander(degree + 1, degree + 1)
        integer :: a, j

        status = 1
        if (degree < 1 .or. degree > 4 .or. n_panels < 3 .or. nfp < 1) return
        if (ntor < 0 .or. ntor_phi < 0) return
        s%degree = degree; s%n_panels = n_panels; s%n_nodes = degree*n_panels
        s%nfp = nfp; s%ntor = ntor; s%nmode = 2*ntor + 1
        s%ntor_phi = ntor_phi; s%nmode_phi = 2*ntor_phi + 1
        do a = 1, degree + 1
            do j = 1, degree + 1
                vander(a, j) = (real(a - 1, dp)/degree)**(j - 1)
            end do
        end do
        call invert_small(vander, status)
        if (status /= 0) return
        s%lagrange = transpose(vander)
        allocate (s%rz(2, s%n_nodes, s%nmode), s%phi(s%n_nodes, s%nmode_phi))
        s%rz = 0; s%phi = 0
        status = 0
    end subroutine ts_surface_create

    pure integer function ts_size(s) result(n)
        type(ts_surface_t), intent(in) :: s
        n = 2*s%n_nodes*s%nmode + s%n_nodes*s%nmode_phi
    end function ts_size

    function ts_coefficients(s) result(c)
        type(ts_surface_t), intent(in) :: s
        real(dp) :: c(ts_size(s))
        integer :: ng
        ng = 2*s%n_nodes*s%nmode
        c(1:ng/2) = reshape(s%rz(1, :, :), [ng/2])
        c(ng/2 + 1:ng) = reshape(s%rz(2, :, :), [ng/2])
        c(ng + 1:) = reshape(s%phi, [s%n_nodes*s%nmode_phi])
    end function ts_coefficients

    subroutine ts_set_coefficients(s, c)
        type(ts_surface_t), intent(inout) :: s
        real(dp), intent(in) :: c(:)
        integer :: ng
        ng = 2*s%n_nodes*s%nmode
        s%rz(1, :, :) = reshape(c(1:ng/2), [s%n_nodes, s%nmode])
        s%rz(2, :, :) = reshape(c(ng/2 + 1:ng), [s%n_nodes, s%nmode])
        s%phi = reshape(c(ng + 1:), [s%n_nodes, s%nmode_phi])
    end subroutine ts_set_coefficients

    ! ------------------------------------------------------------------
    ! Point functionals
    ! ------------------------------------------------------------------

    pure integer function node_of(s, p, a)
        type(ts_surface_t), intent(in) :: s
        integer, intent(in) :: p, a
        node_of = (p - 1)*s%degree + a
        if (node_of > s%n_nodes) node_of = node_of - s%n_nodes
    end function node_of

    pure subroutine basis(s, xi, n, dn)
        type(ts_surface_t), intent(in) :: s
        real(dp), intent(in) :: xi
        real(dp), intent(out) :: n(:), dn(:)
        integer :: a, j
        n = 0; dn = 0
        do a = 1, s%degree + 1
            do j = s%degree + 1, 1, -1
                n(a) = n(a)*xi + s%lagrange(a, j)
                if (j > 1) dn(a) = dn(a)*xi + (j - 1)*s%lagrange(a, j)
            end do
        end do
    end subroutine basis

    pure subroutine trig(nmode, nfp, v, t, dt)
        integer, intent(in) :: nmode, nfp
        real(dp), intent(in) :: v
        real(dp), intent(out) :: t(nmode), dt(nmode)
        integer :: n
        t(1) = 1; dt(1) = 0
        do n = 1, (nmode - 1)/2
            t(2*n) = cos(n*nfp*v); dt(2*n) = -n*nfp*sin(n*nfp*v)
            t(2*n + 1) = sin(n*nfp*v); dt(2*n + 1) = n*nfp*cos(n*nfp*v)
        end do
    end subroutine trig

    subroutine point_functionals(s, p, xi, v, f)
        !! f(1:7) = (R, R_xi, Z_xi, R_v, Z_v, Phi_xi, Phi_v) at panel p,
        !! parameter xi (any real: the panel's polynomial continued) and v.
        type(ts_surface_t), intent(in) :: s
        integer, intent(in) :: p
        real(dp), intent(in) :: xi, v
        type(functional_t), intent(inout) :: f(7)
        real(dp) :: n(s%degree + 1), dn(s%degree + 1)
        real(dp) :: t(s%nmode), dt(s%nmode), tp(s%nmode_phi), dtp(s%nmode_phi)
        integer :: a, m, k, ng, node, ng1, lg, lp

        call basis(s, xi, n, dn)
        call trig(s%nmode, s%nfp, v, t, dt)
        call trig(s%nmode_phi, s%nfp, v, tp, dtp)
        ng1 = s%n_nodes*s%nmode
        ng = 2*ng1
        ! storage is kept across calls (allocated once per work array)
        lg = (s%degree + 1)*s%nmode
        lp = (s%degree + 1)*s%nmode_phi
        do k = 1, 7
            if (.not. allocated(f(k)%idx)) allocate (f(k)%idx(max(lg, lp)), &
                f(k)%coef(max(lg, lp)))
            f(k)%n = merge(lg, lp, k <= 5)
            f(k)%shift = 0
        end do
        do m = 1, s%nmode
            do a = 1, s%degree + 1
                node = node_of(s, p, a)
                k = a + (s%degree + 1)*(m - 1)
                f(1)%idx(k) = node + s%n_nodes*(m - 1); f(1)%coef(k) = n(a)*t(m)
                f(2)%idx(k) = f(1)%idx(k); f(2)%coef(k) = dn(a)*t(m)
                f(3)%idx(k) = ng1 + f(1)%idx(k); f(3)%coef(k) = dn(a)*t(m)
                f(4)%idx(k) = f(1)%idx(k); f(4)%coef(k) = n(a)*dt(m)
                f(5)%idx(k) = ng1 + f(1)%idx(k); f(5)%coef(k) = n(a)*dt(m)
            end do
        end do
        do m = 1, s%nmode_phi
            do a = 1, s%degree + 1
                node = node_of(s, p, a)
                k = a + (s%degree + 1)*(m - 1)
                f(6)%idx(k) = ng + node + s%n_nodes*(m - 1); f(6)%coef(k) = dn(a)*tp(m)
                f(7)%idx(k) = f(6)%idx(k); f(7)%coef(k) = n(a)*dtp(m)
            end do
        end do
        f(6)%shift = s%i_t/s%n_panels
        f(7)%shift = s%i_p/(2*pi)
    end subroutine point_functionals

    pure real(dp) function apply(f, c) result(x)
        type(functional_t), intent(in) :: f
        real(dp), intent(in) :: c(:)
        x = dot_product(f%coef(1:f%n), c(f%idx(1:f%n))) + f%shift
    end function apply

    subroutine ts_point(s, p, xi, v, x, n)
        !! Position and the unnormalised normal X_xi x X_v (Cartesian).
        type(ts_surface_t), intent(in) :: s
        integer, intent(in) :: p
        real(dp), intent(in) :: xi, v
        real(dp), intent(out) :: x(3), n(3)
        type(functional_t) :: f(7)
        real(dp) :: c(ts_size(s)), r, z, ru, zu, rv, zv, xu(3), xv(3)
        type(functional_t) :: fz
        c = ts_coefficients(s)
        call point_functionals(s, p, xi, v, f)
        r = apply(f(1), c); ru = apply(f(2), c); zu = apply(f(3), c)
        rv = apply(f(4), c); zv = apply(f(5), c)
        fz = f(1)
        fz%idx = fz%idx + s%n_nodes*s%nmode
        z = apply(fz, c)
        x = [r*cos(v), r*sin(v), z]
        xu = [ru*cos(v), ru*sin(v), zu]
        xv = [rv*cos(v) - r*sin(v), rv*sin(v) + r*cos(v), zv]
        n = [xu(2)*xv(3) - xu(3)*xv(2), xu(3)*xv(1) - xu(1)*xv(3), &
            xu(1)*xv(2) - xu(2)*xv(1)]
    end subroutine ts_point

    ! ------------------------------------------------------------------
    ! Self energy
    ! ------------------------------------------------------------------

    subroutine ts_energy(s, plan, energy, gradient, hessian)
        !! E, and when present its gradient and Hessian in ts_coefficients(s).
        !! The double sum is invariant under exchanging the two points with a
        !! shift by whole periods, and both the separated and the near rules
        !! are mirror symmetric, so the source side of the gradient and of
        !! the own Hessian blocks equals the target side: these are taken
        !! twice, a cross block C enters as 2 C or 2 C^T and the Hessian is
        !! symmetrised at the end.
        type(ts_surface_t), intent(in) :: s
        type(ts_plan_t), intent(in) :: plan
        real(dp), intent(out) :: energy
        real(dp), intent(out), optional :: gradient(:), hessian(:, :)
        type(far_points_t), allocatable :: far(:)
        type(frame_t) :: fr
        real(dp), allocatable :: c(:), lg(:), lh(:, :), centre(:, :), radius(:)
        real(dp), allocatable :: gl(:), lxx(:, :)
        real(dp) :: le
        integer :: order, nc, ng, nh, p, b, q, bb, nv_all, target, dq, dbb

        c = ts_coefficients(s)
        nc = size(c)
        order = 0
        if (present(gradient)) order = 1
        if (present(hessian)) order = 2
        ng = 0; nh = 0
        if (order >= 1) ng = nc
        if (order == 2) nh = nc
        allocate (lg(ng), lh(nh, nh), source=0.0_dp)
        nv_all = plan%nv*s%nfp
        call make_frame(s, c, fr)
        call far_tables(s, plan, c, far, centre, radius)
        le = 0
        !$omp parallel do schedule(dynamic, 1) default(shared) &
        !$omp& private(target, p, b, q, bb, dq, dbb, gl, lxx) reduction(+:le, lg, lh)
        do target = 1, plan%nv*s%n_panels
            p = 1 + mod(target - 1, s%n_panels)
            b = 1 + (target - 1)/s%n_panels
            if (.not. allocated(gl)) allocate (gl(fr%nl), lxx(fr%nl, fr%nl))
            gl = 0
            if (order == 2) lxx = 0
            call far_target(s, plan, fr, far, centre, radius, p, b, order, le, gl, &
                lxx, lh)
            do dbb = -1, 1
                bb = modulo(b - 1 + dbb, nv_all) + 1
                if (centred(bb - b, nv_all) /= dbb) cycle
                do dq = -1, 1
                    q = modulo(p - 1 + dq, s%n_panels) + 1
                    if (centred(q - p, s%n_panels) /= dq) cycle
                    call patch_pair(s, plan, fr, p, b, q, bb, order, le, gl, lxx, lh)
                end do
            end do
            if (order >= 1) lg(fr%gidx(:, p)) = lg(fr%gidx(:, p)) + 2*gl
            if (order == 2) lh(fr%gidx(:, p), fr%gidx(:, p)) = &
                lh(fr%gidx(:, p), fr%gidx(:, p)) + 2*lxx
        end do
        !$omp end parallel do
        energy = s%nfp*le/(8*pi)
        if (present(gradient)) gradient = s%nfp*lg/(8*pi)
        if (present(hessian)) hessian = s%nfp*(lh + transpose(lh))/(16*pi)
    end subroutine ts_energy

    pure integer function centred(d, n)
        !! d modulo n in [-n/2, n/2).
        integer, intent(in) :: d, n
        centred = modulo(d + n/2, n) - n/2
    end function centred

    subroutine make_frame(s, c, fr)
        !! Local coefficients of each panel: R and Z (node a fastest, then
        !! mode), then Phi; their global indices and values.
        type(ts_surface_t), intent(in) :: s
        real(dp), intent(in) :: c(:)
        type(frame_t), intent(out) :: fr
        integer :: q, a, m, k1, l, node, zoff
        k1 = s%degree + 1
        fr%nlg = k1*s%nmode; fr%nlp = k1*s%nmode_phi
        fr%nl = 2*fr%nlg + fr%nlp
        zoff = s%n_nodes*s%nmode
        allocate (fr%gidx(fr%nl, s%n_panels), fr%cl(fr%nl, s%n_panels))
        do q = 1, s%n_panels
            do a = 1, k1
                node = node_of(s, q, a)
                do m = 1, s%nmode
                    l = a + k1*(m - 1)
                    fr%gidx(l, q) = node + s%n_nodes*(m - 1)
                    fr%gidx(fr%nlg + l, q) = zoff + fr%gidx(l, q)
                end do
                do m = 1, s%nmode_phi
                    fr%gidx(2*fr%nlg + a + k1*(m - 1), q) = 2*zoff + node + s%n_nodes*(m - 1)
                end do
            end do
            fr%cl(:, q) = c(fr%gidx(:, q))
        end do
    end subroutine make_frame

    pure subroutine point_vectors(s, n, dn, t, dt, tp, dtp, gv, pv)
        !! Basis vectors of the point inputs in the local coefficients:
        !! gv(:, 1:3) = N T, dN T, N dT (R and Z), pv(:, 1:2) = dN T, N dT (Phi).
        type(ts_surface_t), intent(in) :: s
        real(dp), intent(in) :: n(:), dn(:), t(:), dt(:), tp(:), dtp(:)
        real(dp), intent(out) :: gv(:, :), pv(:, :)
        integer :: a, m, l, k1
        k1 = s%degree + 1
        do m = 1, s%nmode
            do a = 1, k1
                l = a + k1*(m - 1)
                gv(l, 1) = n(a)*t(m); gv(l, 2) = dn(a)*t(m); gv(l, 3) = n(a)*dt(m)
            end do
        end do
        do m = 1, s%nmode_phi
            do a = 1, k1
                l = a + k1*(m - 1)
                pv(l, 1) = dn(a)*tp(m); pv(l, 2) = n(a)*dtp(m)
            end do
        end do
    end subroutine point_vectors

    pure subroutine point_at(s, xi, v, gv, pv)
        type(ts_surface_t), intent(in) :: s
        real(dp), intent(in) :: xi, v
        real(dp), intent(out) :: gv(:, :), pv(:, :)
        real(dp) :: n(s%degree + 1), dn(s%degree + 1)
        real(dp) :: t(s%nmode), dt(s%nmode), tp(s%nmode_phi), dtp(s%nmode_phi)
        call basis(s, xi, n, dn)
        call trig(s%nmode, s%nfp, v, t, dt)
        call trig(s%nmode_phi, s%nfp, v, tp, dtp)
        call point_vectors(s, n, dn, t, dt, tp, dtp, gv, pv)
    end subroutine point_at

    pure function point_inputs(s, fr, q, gv, pv) result(u)
        !! (R, Z, R_xi, Z_xi, R_v, Z_v, Phi_xi, Phi_v) on panel q, with the
        !! net-current ramps in Phi.
        type(ts_surface_t), intent(in) :: s
        type(frame_t), intent(in) :: fr
        integer, intent(in) :: q
        real(dp), intent(in) :: gv(:, :), pv(:, :)
        real(dp) :: u(8)
        integer :: g
        g = fr%nlg
        u(1) = dot_product(gv(:, 1), fr%cl(1:g, q))
        u(2) = dot_product(gv(:, 1), fr%cl(g + 1:2*g, q))
        u(3) = dot_product(gv(:, 2), fr%cl(1:g, q))
        u(4) = dot_product(gv(:, 2), fr%cl(g + 1:2*g, q))
        u(5) = dot_product(gv(:, 3), fr%cl(1:g, q))
        u(6) = dot_product(gv(:, 3), fr%cl(g + 1:2*g, q))
        u(7) = dot_product(pv(:, 1), fr%cl(2*g + 1:, q)) + s%i_t/s%n_panels
        u(8) = dot_product(pv(:, 2), fr%cl(2*g + 1:, q)) + s%i_p/(2*pi)
    end function point_inputs

    pure subroutine add_gradient(nlg, gv, pv, gp, w, gl)
        !! gl += w B^T gp for one point.
        integer, intent(in) :: nlg
        real(dp), intent(in) :: gv(:, :), pv(:, :), gp(8), w
        real(dp), intent(inout) :: gl(:)
        gl(1:nlg) = gl(1:nlg) + w*(gp(1)*gv(:, 1) + gp(3)*gv(:, 2) + gp(5)*gv(:, 3))
        gl(nlg + 1:2*nlg) = gl(nlg + 1:2*nlg) &
            + w*(gp(2)*gv(:, 1) + gp(4)*gv(:, 2) + gp(6)*gv(:, 3))
        gl(2*nlg + 1:) = gl(2*nlg + 1:) + w*(gp(7)*pv(:, 1) + gp(8)*pv(:, 2))
    end subroutine add_gradient

    pure subroutine add_block(nlg, ga, pa, gb, pb, h, w, l)
        !! l += w B_a^T h B_b, B the point inputs in local coefficients: the
        !! inputs of R are 1, 3, 5 (vectors ga(:, 1:3)), of Z 2, 4, 6, of Phi
        !! 7, 8 (pa(:, 1:2)).
        integer, intent(in) :: nlg
        real(dp), intent(in) :: ga(:, :), pa(:, :), gb(:, :), pb(:, :), h(8, 8), w
        real(dp), intent(inout) :: l(:, :)
        integer, parameter :: rows(3, 3) = reshape([1, 3, 5, 2, 4, 6, 7, 8, 0], [3, 3])
        integer, parameter :: nk(3) = [3, 3, 2]
        real(dp) :: colv(3), wk(3)
        integer :: ca, cb, ka, jj, lb, ob, nlp
        nlp = size(pa, 1)
        do cb = 1, 3
            lb = merge(nlg, nlp, cb < 3)
            ob = (cb - 1)*nlg
            do jj = 1, lb
                if (cb < 3) then
                    colv = gb(jj, 1:3)
                else
                    colv(1:2) = pb(jj, 1:2)
                end if
                do ca = 1, 3
                    do ka = 1, nk(ca)
                        wk(ka) = w*dot_product(h(rows(ka, ca), rows(1:nk(cb), cb)), &
                            colv(1:nk(cb)))
                    end do
                    if (ca < 3) then
                        l((ca - 1)*nlg + 1:ca*nlg, ob + jj) = l((ca - 1)*nlg + 1:ca*nlg, &
                            ob + jj) + wk(1)*ga(:, 1) + wk(2)*ga(:, 2) + wk(3)*ga(:, 3)
                    else
                        l(2*nlg + 1:, ob + jj) = l(2*nlg + 1:, ob + jj) + wk(1)*pa(:, 1) &
                            + wk(2)*pa(:, 2)
                    end if
                end do
            end do
        end do
    end subroutine add_block

    subroutine far_tables(s, plan, c, far, centre, radius)
        !! Quadrature points of every patch (panel q, interval bb of the torus,
        !! patch index q + n_panels (bb - 1)) for the orders 2..ng_far with
        !! their inputs u = (R, Z, R_xi, Z_xi, R_v, Z_v, Phi_xi, Phi_v), and
        !! each patch's centre and radius.
        type(ts_surface_t), intent(in) :: s
        type(ts_plan_t), intent(in) :: plan
        real(dp), intent(in) :: c(:)
        type(far_points_t), allocatable, intent(out) :: far(:)
        real(dp), allocatable, intent(out) :: centre(:, :), radius(:)
        type(functional_t) :: f(7)
        real(dp) :: hv, g(plan%ng_far), gw(plan%ng_far), xx(3), nn(3)
        integer :: n, q, bb, i, j, k, patch, npatch, zoff
        hv = 2*pi/(s%nfp*plan%nv)
        npatch = plan%nv*s%nfp*s%n_panels
        zoff = s%n_nodes*s%nmode
        allocate (far(2:plan%ng_far), centre(3, npatch), radius(npatch))
        do n = 2, plan%ng_far
            call gauss_unit(n, g(:n), gw(:n))
            far(n)%n = n
            allocate (far(n)%xi(n*n, npatch), far(n)%v(n*n, npatch), &
                far(n)%w(n*n, npatch), far(n)%u(8, n*n, npatch))
            allocate (far(n)%nb(s%degree + 1, n*n, npatch), &
                far(n)%dnb(s%degree + 1, n*n, npatch), far(n)%tt(s%nmode, n*n, npatch), &
                far(n)%dtt(s%nmode, n*n, npatch), far(n)%tp(s%nmode_phi, n*n, npatch), &
                far(n)%dtp(s%nmode_phi, n*n, npatch))
            !$omp parallel do default(shared) private(patch, q, bb, i, j, k, f)
            do patch = 1, npatch
                q = 1 + mod(patch - 1, s%n_panels)
                bb = 1 + (patch - 1)/s%n_panels
                do j = 1, n
                    do i = 1, n
                        k = i + n*(j - 1)
                        far(n)%xi(k, patch) = g(i)
                        far(n)%v(k, patch) = (bb - 1 + g(j))*hv
                        far(n)%w(k, patch) = gw(i)*gw(j)*hv
                        call point_functionals(s, q, g(i), far(n)%v(k, patch), f)
                        call basis(s, g(i), far(n)%nb(:, k, patch), far(n)%dnb(:, k, patch))
                        call trig(s%nmode, s%nfp, far(n)%v(k, patch), far(n)%tt(:, k, patch), &
                            far(n)%dtt(:, k, patch))
                        call trig(s%nmode_phi, s%nfp, far(n)%v(k, patch), &
                            far(n)%tp(:, k, patch), far(n)%dtp(:, k, patch))
                        far(n)%u(:, k, patch) = [apply(f(1), c), &
                            apply_offset(f(1), c, zoff), apply(f(2), c), &
                            apply(f(3), c), apply(f(4), c), apply(f(5), c), &
                            apply(f(6), c), apply(f(7), c)]
                    end do
                end do
            end do
            !$omp end parallel do
        end do
        do patch = 1, npatch
            q = 1 + mod(patch - 1, s%n_panels)
            bb = 1 + (patch - 1)/s%n_panels
            call ts_point(s, q, 0.5_dp, (bb - 0.5_dp)*hv, centre(:, patch), nn)
            radius(patch) = 0
            do j = 0, 2
                do i = 0, 2
                    call ts_point(s, q, 0.5_dp*i, (bb - 1 + 0.5_dp*j)*hv, xx, nn)
                    radius(patch) = max(radius(patch), norm2(xx - centre(:, patch)))
                end do
            end do
        end do
    end subroutine far_tables

    pure integer function far_order(plan, distance, size) result(n)
        !! Gauss points per direction for separated patches: the tensor rule's
        !! error behaves like (size/(2 distance))^(2n).
        type(ts_plan_t), intent(in) :: plan
        real(dp), intent(in) :: distance, size
        real(dp) :: ratio
        ratio = distance/max(size, tiny(1.0_dp))
        if (ratio <= 1) then
            n = plan%ng_far
        else
            n = ceiling(-log10(plan%far_tolerance)/(2*log10(2*ratio)))
            n = min(plan%ng_far, max(2, n))
        end if
    end function far_order

    subroutine far_target(s, plan, fr, far, centre, radius, p, b, order, le, gl, &
            lxx, lh)
        !! All separated source patches of the target patch (p, b). For each
        !! order n the target points of that order meet every source patch
        !! that takes it, with the inputs from the tables. Per target point the
        !! gradient and own block go to the local gl, lxx, the cross block
        !! through a row block in the global coefficients into lh as 2 C^T.
        type(ts_surface_t), intent(in) :: s
        type(ts_plan_t), intent(in) :: plan
        type(frame_t), intent(in) :: fr
        real(dp), intent(in) :: centre(:, :), radius(:)
        type(far_points_t), intent(in) :: far(2:)
        integer, intent(in) :: p, b, order
        real(dp), intent(inout) :: le, gl(:), lxx(:, :), lh(:, :)
        real(dp), allocatable :: cross(:, :), col(:)
        logical, allocatable :: use(:, :)
        real(dp) :: out(1 + n_self + n_self*(n_self + 1)/2), ux(8), uy(8), wx, wy
        real(dp) :: gp(16), hxx(8, 8), hxy(8, 8), sxx(8, 8), gx(8)
        real(dp) :: gvx(fr%nlg, 3), pvx(fr%nlp, 2)
        integer :: nv_all, npatch, tp, patch, n, k, l, q, bb, nc, lx

        nv_all = plan%nv*s%nfp
        npatch = nv_all*s%n_panels
        nc = size(lh, 1)
        tp = p + s%n_panels*(b - 1)
        ! which separated patches take which order
        allocate (use(npatch, 2:plan%ng_far))
        use = .false.
        do patch = 1, npatch
            q = 1 + mod(patch - 1, s%n_panels)
            bb = 1 + (patch - 1)/s%n_panels
            if (abs(centred(q - p, s%n_panels)) <= 1 .and. &
                abs(centred(bb - b, nv_all)) <= 1) cycle
            n = far_order(plan, norm2(centre(:, patch) - centre(:, tp)), &
                radius(patch) + radius(tp))
            use(patch, n) = .true.
        end do
        if (order == 2) allocate (cross(8, nc), col(nc))
        do n = 2, plan%ng_far
            if (.not. any(use(:, n))) cycle
            do k = 1, n*n
                ux = far(n)%u(:, k, tp)
                wx = far(n)%w(k, tp)
                gx = 0
                if (order == 2) then
                    sxx = 0; cross = 0
                end if
                do patch = 1, npatch
                    if (.not. use(patch, n)) cycle
                    do l = 1, n*n
                        uy = far(n)%u(:, l, patch)
                        wy = far(n)%w(l, patch)
                        call far_leaf(ux, uy, far(n)%v(k, tp), far(n)%v(l, patch), &
                            order, out)
                        le = le + wx*wy*out(1)
                        if (order == 0) cycle
                        call leaf_to_points(out, 1.0_dp, order, gp, hxx, hxy)
                        gx = gx + wy*gp(1:8)
                        if (order < 2) cycle
                        sxx = sxx + wy*hxx
                        call add_cross(s, far(n), l, patch, wy*hxy, cross)
                    end do
                end do
                if (order == 0) cycle
                call point_vectors(s, far(n)%nb(:, k, tp), far(n)%dnb(:, k, tp), &
                    far(n)%tt(:, k, tp), far(n)%dtt(:, k, tp), far(n)%tp(:, k, tp), &
                    far(n)%dtp(:, k, tp), gvx, pvx)
                call add_gradient(fr%nlg, gvx, pvx, gx, wx, gl)
                if (order < 2) cycle
                call add_block(fr%nlg, gvx, pvx, gvx, pvx, sxx, wx, lxx)
                do lx = 1, fr%nlg
                    col = gvx(lx, 1)*cross(1, :) + gvx(lx, 2)*cross(3, :) &
                        + gvx(lx, 3)*cross(5, :)
                    lh(:, fr%gidx(lx, p)) = lh(:, fr%gidx(lx, p)) + 2*wx*col
                    col = gvx(lx, 1)*cross(2, :) + gvx(lx, 2)*cross(4, :) &
                        + gvx(lx, 3)*cross(6, :)
                    lh(:, fr%gidx(fr%nlg + lx, p)) = lh(:, fr%gidx(fr%nlg + lx, p)) &
                        + 2*wx*col
                end do
                do lx = 1, fr%nlp
                    col = pvx(lx, 1)*cross(7, :) + pvx(lx, 2)*cross(8, :)
                    lh(:, fr%gidx(2*fr%nlg + lx, p)) = lh(:, fr%gidx(2*fr%nlg + lx, p)) &
                        + 2*wx*col
                end do
            end do
        end do
    end subroutine far_target

    subroutine far_leaf(ux, uy, vx, vy, order, out)
        !! The self integrand of a separated pair (t = 1) from point inputs.
        real(dp), intent(in) :: ux(8), uy(8), vx, vy
        integer, intent(in) :: order
        real(dp), intent(out) :: out(:)
        real(dp) :: u(n_self), c1, s1, c2, s2
        u = [ux(1), ux(3:8), uy(1), uy(3:8), uy(1) - ux(1), uy(2) - ux(2)]
        c1 = cos(vx); s1 = sin(vx); c2 = cos(vy); s2 = sin(vy)
        call leaf(u, c1, s1, c2, s2, c2 - c1, s2 - s1, order, out)
    end subroutine far_leaf

    subroutine leaf(u, c1, s1, c2, s2, dc, ds, order, out)
        real(dp), intent(in) :: u(n_self), c1, s1, c2, s2, dc, ds
        integer, intent(in) :: order
        real(dp), intent(out) :: out(:)
        select case (order)
        case (0)
            call generated_toroidal_sheet_self_value(u(1), u(2), u(3), u(4), u(5), &
                u(6), u(7), u(8), u(9), u(10), u(11), u(12), u(13), u(14), u(15), &
                u(16), c1, s1, c2, s2, dc, ds, out(1))
        case (1)
            call generated_toroidal_sheet_self_gradient(u(1), u(2), u(3), u(4), &
                u(5), u(6), u(7), u(8), u(9), u(10), u(11), u(12), u(13), u(14), &
                u(15), u(16), c1, s1, c2, s2, dc, ds, out(1), out(2), out(3), &
                out(4), out(5), out(6), out(7), out(8), out(9), out(10), out(11), &
                out(12), out(13), out(14), out(15), out(16), out(17))
        case default
            call self_jet(u, c1, s1, c2, s2, dc, ds, out)
        end select
    end subroutine leaf

    pure subroutine leaf_to_points(out, t, order, gp, hxx, hxy)
        !! Gradient of the leaf (inputs x 7, y 7, (R_y - R_x)/t, (Z_y - Z_x)/t)
        !! in the point inputs (x: R, Z, R_xi, Z_xi, R_v, Z_v, Phi_xi, Phi_v;
        !! then y), and the x-x and x-y Hessian blocks.
        real(dp), intent(in) :: out(:), t
        integer, intent(in) :: order
        real(dp), intent(out) :: gp(16), hxx(8, 8), hxy(8, 8)
        real(dp) :: he(n_self, n_self), cm(2, 16)
        integer :: lm(2, 16), nm(16), i, j, k, a, c
        nm = 1; lm = 0; cm = 0
        lm(1, 1) = 1; cm(1, 1) = 1; lm(2, 1) = 15; cm(2, 1) = -1/t; nm(1) = 2
        lm(1, 2) = 16; cm(1, 2) = -1/t
        lm(1, 9) = 8; cm(1, 9) = 1; lm(2, 9) = 15; cm(2, 9) = 1/t; nm(9) = 2
        lm(1, 10) = 16; cm(1, 10) = 1/t
        do i = 3, 8
            lm(1, i) = i - 1; cm(1, i) = 1
            lm(1, i + 8) = i + 6; cm(1, i + 8) = 1
        end do
        do i = 1, 16
            gp(i) = 0
            do a = 1, nm(i)
                gp(i) = gp(i) + cm(a, i)*out(1 + lm(a, i))
            end do
        end do
        hxx = 0; hxy = 0
        if (order < 2) return
        k = 1 + n_self
        do i = 1, n_self
            do j = i, n_self
                k = k + 1
                he(i, j) = out(k); he(j, i) = out(k)
            end do
        end do
        do j = 1, 16
            do i = 1, 8
                do c = 1, nm(j)
                    do a = 1, nm(i)
                        if (j <= 8) then
                            hxx(i, j) = hxx(i, j) + cm(a, i)*cm(c, j)*he(lm(a, i), lm(c, j))
                        else
                            hxy(i, j - 8) = hxy(i, j - 8) &
                                + cm(a, i)*cm(c, j)*he(lm(a, i), lm(c, j))
                        end if
                    end do
                end do
            end do
        end do
    end subroutine leaf_to_points

    subroutine add_cross(s, far, l, patch, hxy, cross)
        !! cross(i, :) += sum_j hxy(i, j) J_y(j, :) for source point l of the
        !! patch, J_y from the tabulated basis values.
        type(ts_surface_t), intent(in) :: s
        type(far_points_t), intent(in) :: far
        integer, intent(in) :: l, patch
        real(dp), intent(in) :: hxy(8, 8)
        real(dp), intent(inout) :: cross(:, :)
        integer :: q, a, m, node, idx, zoff, poff
        real(dp) :: n, dn, t, dt
        q = 1 + mod(patch - 1, s%n_panels)
        zoff = s%n_nodes*s%nmode
        poff = 2*zoff
        do a = 1, s%degree + 1
            node = node_of(s, q, a)
            n = far%nb(a, l, patch); dn = far%dnb(a, l, patch)
            do m = 1, s%nmode
                t = far%tt(m, l, patch); dt = far%dtt(m, l, patch)
                idx = node + s%n_nodes*(m - 1)
                cross(:, idx) = cross(:, idx) + hxy(:, 1)*n*t + hxy(:, 3)*dn*t &
                    + hxy(:, 5)*n*dt
                cross(:, zoff + idx) = cross(:, zoff + idx) + hxy(:, 2)*n*t &
                    + hxy(:, 4)*dn*t + hxy(:, 6)*n*dt
            end do
            do m = 1, s%nmode_phi
                idx = poff + node + s%n_nodes*(m - 1)
                cross(:, idx) = cross(:, idx) + hxy(:, 7)*dn*far%tp(m, l, patch) &
                    + hxy(:, 8)*n*far%dtp(m, l, patch)
            end do
        end do
    end subroutine add_cross

    subroutine patch_pair(s, plan, fr, p, b, q, bb, order, le, gl, lxx, lh)
        !! Touching or equal patches: target (panel p, toroidal interval b of
        !! the first period) and source (q, bb over the torus); in the
        !! target's continued parameters the source occupies [o1, o1 + 1] x
        !! [o2, o2 + h]. Relative coordinates z = y - x; per direction the
        !! overlap length of x and y = x + z is the tent w - |z_i - o_i|, so z
        !! is split at o_i into rectangles, those with a corner at z = 0 take a
        !! polar rule about it (the radius cancels 1/|X_y - X_x|), the others
        !! tensor Gauss, and x runs over the overlap with tensor Gauss. The
        !! cross block goes into lh as 2 C.
        type(ts_surface_t), intent(in) :: s
        type(ts_plan_t), intent(in) :: plan
        type(frame_t), intent(in) :: fr
        integer, intent(in) :: p, b, q, bb, order
        real(dp), intent(inout) :: le, gl(:), lxx(:, :), lh(:, :)
        real(dp) :: hv, o(2), w(2), lo(2), hi(2), corner(2, 5), d(2), z(2)
        real(dp) :: gu(plan%ng_u), gw(plan%ng_u)
        real(dp) :: tq(plan%ng_t), tw(plan%ng_t), sq(plan%ng_s), sw(plan%ng_s)
        real(dp) :: det, rho, wz
        real(dp) :: gvx(fr%nlg, 3), pvx(fr%nlp, 2), gvy(fr%nlg, 3), pvy(fr%nlp, 2)
        real(dp), allocatable :: lxy(:, :)
        integer :: nv_all, dp_, db, k, l, m, ri, rj

        hv = 2*pi/(s%nfp*plan%nv)
        nv_all = plan%nv*s%nfp
        dp_ = centred(q - p, s%n_panels)
        db = centred(bb - b, nv_all)
        call gauss_unit(plan%ng_u, gu, gw)
        call gauss_unit(plan%ng_t, tq, tw)
        call gauss_unit(plan%ng_s, sq, sw)
        if (order == 2) allocate (lxy(fr%nl, fr%nl), source=0.0_dp)
        o = [real(dp_, dp), db*hv]
        w = [1.0_dp, hv]
        ! z rectangles [o_i - w_i, o_i] and [o_i, o_i + w_i]
        do rj = 0, 1
            do ri = 0, 1
                lo = [o(1) - w(1) + ri*w(1), o(2) - w(2) + rj*w(2)]
                hi = lo + w
                if (abs(lo(1)) < 1.0e-12_dp .or. abs(hi(1)) < 1.0e-12_dp) then
                    if (abs(lo(2)) < 1.0e-12_dp*hv .or. abs(hi(2)) < 1.0e-12_dp*hv) then
                        ! polar about the corner z = 0
                        corner(:, 1) = lo
                        corner(:, 2) = [hi(1), lo(2)]
                        corner(:, 3) = hi
                        corner(:, 4) = [lo(1), hi(2)]
                        corner(:, 5) = corner(:, 1)
                        do k = 1, 4
                            d = corner(:, k + 1) - corner(:, k)
                            det = corner(1, k)*d(2) - corner(2, k)*d(1)
                            if (abs(det) <= 1.0e-12_dp*hv) cycle
                            do l = 1, plan%ng_s
                                do m = 1, plan%ng_t
                                    rho = tq(m)
                                    z = rho*(corner(:, k) + sq(l)*d)
                                    wz = tw(m)*sw(l)*det
                                    call overlap(z, rho, wz)
                                end do
                            end do
                        end do
                        cycle
                    end if
                end if
                ! z = 0 is not on this rectangle: tensor Gauss
                do m = 1, plan%ng_t
                    do l = 1, plan%ng_t
                        z = lo + [tq(l), tq(m)]*w
                        call overlap(z, 1.0_dp, tw(l)*tw(m)*w(1)*w(2))
                    end do
                end do
            end do
        end do
        if (order == 2) lh(fr%gidx(:, p), fr%gidx(:, q)) = &
            lh(fr%gidx(:, p), fr%gidx(:, q)) + 2*lxy
    contains
        subroutine overlap(z, t, wz)
            !! x over the overlap of the target patch and the source patch
            !! shifted by -z; X_y - X_x = t Q (t the polar radius, or 1).
            real(dp), intent(in) :: z(2), t, wz
            real(dp) :: xl(2), xh(2), xx(2), ux(8), uy(8), u(n_self), vy, wt
            real(dp) :: out(1 + n_self + n_self*(n_self + 1)/2), gp(16), hxx(8, 8), hxy(8, 8)
            integer :: i1, i2
            xl = [max(0.0_dp, o(1) - z(1)), max((b - 1)*hv, (b - 1)*hv + o(2) - z(2))]
            xh = [min(1.0_dp, o(1) + 1 - z(1)), min(b*hv, (b - 1)*hv + o(2) + hv - z(2))]
            if (any(xh <= xl)) return
            do i2 = 1, plan%ng_u
                do i1 = 1, plan%ng_u
                    xx = xl + [gu(i1), gu(i2)]*(xh - xl)
                    wt = wz*gw(i1)*gw(i2)*(xh(1) - xl(1))*(xh(2) - xl(2))
                    vy = xx(2) + z(2)
                    call point_at(s, xx(1), xx(2), gvx, pvx)
                    call point_at(s, xx(1) + z(1) - o(1), vy, gvy, pvy)
                    ux = point_inputs(s, fr, p, gvx, pvx)
                    uy = point_inputs(s, fr, q, gvy, pvy)
                    u = [ux(1), ux(3:8), uy(1), uy(3:8), (uy(1) - ux(1))/t, &
                        (uy(2) - ux(2))/t]
                    ! (cos v_y - cos v_x)/t and (sin v_y - sin v_x)/t without
                    ! cancellation
                    call leaf(u, cos(xx(2)), sin(xx(2)), cos(vy), sin(vy), &
                        -2*sin(0.5_dp*z(2))*sin(0.5_dp*(vy + xx(2)))/t, &
                        2*sin(0.5_dp*z(2))*cos(0.5_dp*(vy + xx(2)))/t, order, out)
                    le = le + wt*out(1)
                    if (order == 0) cycle
                    call leaf_to_points(out, t, order, gp, hxx, hxy)
                    call add_gradient(fr%nlg, gvx, pvx, gp(1:8), wt, gl)
                    if (order < 2) cycle
                    call add_block(fr%nlg, gvx, pvx, gvx, pvx, hxx, wt, lxx)
                    call add_block(fr%nlg, gvx, pvx, gvy, pvy, hxy, wt, lxy)
                end do
            end do
        end subroutine overlap
    end subroutine patch_pair

    pure real(dp) function apply_offset(f, c, offset) result(x)
        type(functional_t), intent(in) :: f
        real(dp), intent(in) :: c(:)
        integer, intent(in) :: offset
        x = dot_product(f%coef(1:f%n), c(f%idx(1:f%n) + offset))
    end function apply_offset


    subroutine self_jet(u, c1, s1, c2, s2, dc, ds, out)
        real(dp), intent(in) :: u(n_self), c1, s1, c2, s2, dc, ds
        real(dp), intent(out) :: out(:)
        call generated_toroidal_sheet_self_jet(u(1), u(2), u(3), u(4), u(5), u(6), &
            u(7), u(8), u(9), u(10), u(11), u(12), u(13), u(14), u(15), u(16), c1, &
            s1, c2, s2, dc, ds, out(1), out(2), out(3), out(4), out(5), out(6), &
            out(7), out(8), out(9), out(10), out(11), out(12), out(13), out(14), &
            out(15), out(16), out(17), out(18), out(19), out(20), out(21), out(22), &
            out(23), out(24), out(25), out(26), out(27), out(28), out(29), out(30), &
            out(31), out(32), out(33), out(34), out(35), out(36), out(37), out(38), &
            out(39), out(40), out(41), out(42), out(43), out(44), out(45), out(46), &
            out(47), out(48), out(49), out(50), out(51), out(52), out(53), out(54), &
            out(55), out(56), out(57), out(58), out(59), out(60), out(61), out(62), &
            out(63), out(64), out(65), out(66), out(67), out(68), out(69), out(70), &
            out(71), out(72), out(73), out(74), out(75), out(76), out(77), out(78), &
            out(79), out(80), out(81), out(82), out(83), out(84), out(85), out(86), &
            out(87), out(88), out(89), out(90), out(91), out(92), out(93), out(94), &
            out(95), out(96), out(97), out(98), out(99), out(100), out(101), &
            out(102), out(103), out(104), out(105), out(106), out(107), out(108), &
            out(109), out(110), out(111), out(112), out(113), out(114), out(115), &
            out(116), out(117), out(118), out(119), out(120), out(121), out(122), &
            out(123), out(124), out(125), out(126), out(127), out(128), out(129), &
            out(130), out(131), out(132), out(133), out(134), out(135), out(136), &
            out(137), out(138), out(139), out(140), out(141), out(142), out(143), &
            out(144), out(145), out(146), out(147), out(148), out(149), out(150), &
            out(151), out(152), out(153))
    end subroutine self_jet

    ! ------------------------------------------------------------------
    ! External field load
    ! ------------------------------------------------------------------

    subroutine ts_load(s, plan, field, value, gradient, hessian)
        !! int Phi B . (X_xi x X_v) dxi dv over the torus cut at u_c = 0 and
        !! v = 0, Phi = phi + i_t u_c + i_p v/(2 pi) with the ramps of the net
        !! currents, and when present its gradient and Hessian in the
        !! coefficients. Gauss per patch.
        type(ts_surface_t), intent(in) :: s
        type(ts_plan_t), intent(in) :: plan
        procedure(ts_field) :: field
        real(dp), intent(out) :: value
        real(dp), intent(out), optional :: gradient(:), hessian(:, :)
        real(dp), allocatable :: c(:)
        real(dp) :: gu(plan%ng_u), gw(plan%ng_u), gv(plan%ng_v), gvw(plan%ng_v)
        real(dp) :: hv, xi, v, w, phi, out(1 + 6 + 21), b(3), db(3, 2), d2b(3, 3)
        real(dp) :: u(5), gr(6), he(6, 6), r, z
        type(functional_t) :: f(7), fphi, fz, lin(6)
        integer :: p, bb, i, j, k, a, l, m, ng1, kk

        c = ts_coefficients(s)
        call gauss_unit(plan%ng_u, gu, gw)
        call gauss_unit(plan%ng_v, gv, gvw)
        hv = 2*pi/(s%nfp*plan%nv)
        ng1 = s%n_nodes*s%nmode
        value = 0
        if (present(gradient)) gradient = 0
        if (present(hessian)) hessian = 0
        do bb = 1, plan%nv*s%nfp
            do p = 1, s%n_panels
                do j = 1, plan%ng_v
                    do i = 1, plan%ng_u
                        xi = gu(i); v = (bb - 1 + gv(j))*hv
                        w = gw(i)*gvw(j)*hv
                        call point_functionals(s, p, xi, v, f)
                        fz = f(1)
                        fz%idx = fz%idx + ng1
                        ! phi itself: the value functional of the current
                        fphi = phi_value(s, p, xi, v)
                        r = apply(f(1), c); z = apply(fz, c)
                        call field(r, z, v, b, db, d2b)
                        ! leaf variables (r, z, r_xi, z_xi, r_v, z_v)
                        lin = [f(1), fz, f(2), f(3), f(4), f(5)]
                        u = [r, apply(f(2), c), apply(f(3), c), apply(f(4), c), &
                            apply(f(5), c)]
                        call generated_toroidal_sheet_load_jet(u(1), u(2), u(3), &
                            u(4), u(5), b(1), db(1, 1), db(1, 2), d2b(1, 1), &
                            d2b(1, 2), d2b(1, 3), b(2), db(2, 1), db(2, 2), &
                            d2b(2, 1), d2b(2, 2), d2b(2, 3), b(3), db(3, 1), &
                            db(3, 2), d2b(3, 1), d2b(3, 2), d2b(3, 3), out(1), &
                            out(2), out(3), out(4), out(5), out(6), out(7), out(8), &
                            out(9), out(10), out(11), out(12), out(13), out(14), &
                            out(15), out(16), out(17), out(18), out(19), out(20), &
                            out(21), out(22), out(23), out(24), out(25), out(26), &
                            out(27), out(28))
                        ! out: B.N, its gradient in (r, z, r_xi, z_xi, r_v, z_v)
                        ! (z through B's Taylor data) and the Hessian
                        phi = apply(fphi, c) + s%i_t*(p - 1 + xi)/s%n_panels &
                            + s%i_p*v/(2*pi)
                        value = value + w*phi*out(1)
                        if (.not. present(gradient)) cycle
                        call unpack_load(out, gr, he)
                        ! geometry inputs: r, ru, zu, rv, zv and z through B
                        do kk = 1, 6
                            do k = 1, lin(kk)%n
                                gradient(lin(kk)%idx(k)) = gradient(lin(kk)%idx(k)) &
                                    + w*phi*gr(kk)*lin(kk)%coef(k)
                            end do
                        end do
                        do k = 1, fphi%n
                            gradient(fphi%idx(k)) = gradient(fphi%idx(k)) &
                                + w*out(1)*fphi%coef(k)
                        end do
                        if (.not. present(hessian)) cycle
                        do kk = 1, 6
                            do l = 1, 6
                                do k = 1, lin(kk)%n
                                    do a = 1, lin(l)%n
                                        hessian(lin(kk)%idx(k), lin(l)%idx(a)) = &
                                            hessian(lin(kk)%idx(k), lin(l)%idx(a)) &
                                            + w*phi*he(kk, l)*lin(kk)%coef(k) &
                                            *lin(l)%coef(a)
                                    end do
                                end do
                            end do
                            do k = 1, lin(kk)%n
                                do m = 1, fphi%n
                                    hessian(lin(kk)%idx(k), fphi%idx(m)) = &
                                        hessian(lin(kk)%idx(k), fphi%idx(m)) &
                                        + w*gr(kk)*lin(kk)%coef(k)*fphi%coef(m)
                                    hessian(fphi%idx(m), lin(kk)%idx(k)) = &
                                        hessian(fphi%idx(m), lin(kk)%idx(k)) &
                                        + w*gr(kk)*lin(kk)%coef(k)*fphi%coef(m)
                                end do
                            end do
                        end do
                    end do
                end do
            end do
        end do
    end subroutine ts_load

    subroutine unpack_load(out, gr, he)
        real(dp), intent(in) :: out(:)
        real(dp), intent(out) :: gr(6), he(6, 6)
        integer :: i, j, k
        gr = out(2:7)
        k = 7
        do i = 1, 6
            do j = i, 6
                k = k + 1
                he(i, j) = out(k); he(j, i) = out(k)
            end do
        end do
    end subroutine unpack_load

    function phi_value(s, p, xi, v) result(f)
        type(ts_surface_t), intent(in) :: s
        integer, intent(in) :: p
        real(dp), intent(in) :: xi, v
        type(functional_t) :: f
        real(dp) :: n(s%degree + 1), dn(s%degree + 1), tp(s%nmode_phi), dtp(s%nmode_phi)
        integer :: a, m, k
        call basis(s, xi, n, dn)
        call trig(s%nmode_phi, s%nfp, v, tp, dtp)
        f%n = (s%degree + 1)*s%nmode_phi
        allocate (f%idx(f%n), f%coef(f%n))
        do m = 1, s%nmode_phi
            do a = 1, s%degree + 1
                k = a + (s%degree + 1)*(m - 1)
                f%idx(k) = 2*s%n_nodes*s%nmode + node_of(s, p, a) + s%n_nodes*(m - 1)
                f%coef(k) = n(a)*tp(m)
            end do
        end do
    end function phi_value

    pure subroutine invert_small(a, status)
        real(dp), intent(inout) :: a(:, :)
        integer, intent(out) :: status
        real(dp) :: aug(size(a, 1), 2*size(a, 1)), row(2*size(a, 1)), piv
        integer :: n, i, k, p
        n = size(a, 1)
        aug = 0
        aug(:, 1:n) = a
        do i = 1, n
            aug(i, n + i) = 1
        end do
        status = 1
        do k = 1, n
            p = k - 1 + maxloc(abs(aug(k:n, k)), dim=1)
            if (aug(p, k) == 0) return
            row = aug(k, :); aug(k, :) = aug(p, :); aug(p, :) = row
            piv = aug(k, k)
            aug(k, :) = aug(k, :)/piv
            do i = 1, n
                if (i /= k) aug(i, :) = aug(i, :) - aug(i, k)*aug(k, :)
            end do
        end do
        a = aug(:, n + 1:2*n)
        status = 0
    end subroutine invert_small

end module fortfem_toroidal_sheet_current_3d
