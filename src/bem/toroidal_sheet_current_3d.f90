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
        type(ts_surface_t), intent(in) :: s
        type(ts_plan_t), intent(in) :: plan
        real(dp), intent(out) :: energy
        real(dp), intent(out), optional :: gradient(:), hessian(:, :)
        type(far_points_t), allocatable :: far(:)
        real(dp), allocatable :: c(:), lg(:), lh(:, :), centre(:, :), radius(:)
        real(dp) :: le
        integer :: order, nc, ng, nh, p, b, q, bb, nv_all, npatch, target

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
        npatch = nv_all*s%n_panels
        call far_tables(s, plan, c, far, centre, radius)
        le = 0
        !$omp parallel do schedule(dynamic, 1) default(shared) &
        !$omp& private(target, p, b, q, bb) reduction(+:le, lg, lh)
        do target = 1, plan%nv*s%n_panels
            p = 1 + mod(target - 1, s%n_panels)
            b = 1 + (target - 1)/s%n_panels
            call far_target(s, plan, c, far, centre, radius, p, b, order, le, lg, lh)
            do bb = 1, nv_all
                do q = 1, s%n_panels
                    call patch_pair(s, plan, c, p, b, q, bb, order, le, lg, lh)
                end do
            end do
        end do
        !$omp end parallel do
        energy = s%nfp*le/(8*pi)
        if (present(gradient)) gradient = s%nfp*lg/(8*pi)
        if (present(hessian)) hessian = s%nfp*lh/(8*pi)
    end subroutine ts_energy

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

    subroutine far_target(s, plan, c, far, centre, radius, p, b, order, le, lg, lh)
        !! All separated source patches of the target patch (p, b). For each
        !! order n the target points of that order meet every source patch
        !! that takes it; the inputs come from the tables, the gradient and
        !! the point's own Hessian blocks are accumulated per point in input
        !! space and scattered once, the cross blocks through a dense row
        !! block per target point.
        type(ts_surface_t), intent(in) :: s
        type(ts_plan_t), intent(in) :: plan
        real(dp), intent(in) :: c(:), centre(:, :), radius(:)
        type(far_points_t), intent(in) :: far(2:)
        integer, intent(in) :: p, b, order
        real(dp), intent(inout) :: le, lg(:), lh(:, :)
        real(dp), allocatable :: gx(:), gy(:, :, :), hxx(:, :), hyy(:, :, :, :)
        real(dp), allocatable :: cross(:, :)
        logical, allocatable :: use(:, :)
        real(dp) :: out(1 + n_self + n_self*(n_self + 1)/2), hv, ux(8), uy(8), lin(16, 16)
        real(dp) :: hp(16, 16), gp(16), wx, wy
        integer :: nv_all, npatch, tp, patch, n, k, l, dp_, db, q, bb, nc
        type(functional_t) :: f(7)

        hv = 2*pi/(s%nfp*plan%nv)
        nv_all = plan%nv*s%nfp
        npatch = nv_all*s%n_panels
        nc = size(c)
        tp = p + s%n_panels*(b - 1)
        ! which separated patches take which order
        allocate (use(npatch, 2:plan%ng_far))
        use = .false.
        do patch = 1, npatch
            q = 1 + mod(patch - 1, s%n_panels)
            bb = 1 + (patch - 1)/s%n_panels
            dp_ = modulo(q - p + s%n_panels/2, s%n_panels) - s%n_panels/2
            db = modulo(bb - b + nv_all/2, nv_all) - nv_all/2
            if (abs(dp_) <= 1 .and. abs(db) <= 1) cycle
            n = far_order(plan, norm2(centre(:, patch) - centre(:, tp)), &
                radius(patch) + radius(tp))
            use(patch, n) = .true.
        end do
        ! input-space accumulators of the source points (gradient and own
        ! Hessian block), per order
        if (order >= 1) then
            allocate (gy(8, plan%ng_far**2, npatch), gx(8))
            if (order == 2) then
                allocate (hyy(8, 8, plan%ng_far**2, npatch), hxx(8, 8), cross(8, nc))
                hyy = 0
            end if
            gy = 0
        end if
        lin = 0
        do n = 2, plan%ng_far
            if (.not. any(use(:, n))) cycle
            do k = 1, n*n
                ux = far(n)%u(:, k, tp)
                wx = far(n)%w(k, tp)
                if (order >= 1) gx = 0
                if (order == 2) then
                    hxx = 0; cross = 0
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
                        ! leaf inputs (x 7, y 7, dr, dz) in point inputs (x 8, y 8)
                        call leaf_to_points(out, gp, hp, order)
                        gx = gx + wx*wy*gp(1:8)
                        gy(:, l, patch) = gy(:, l, patch) + wx*wy*gp(9:16)
                        if (order < 2) cycle
                        hxx = hxx + wx*wy*hp(1:8, 1:8)
                        hyy(:, :, l, patch) = hyy(:, :, l, patch) + wx*wy*hp(9:16, 9:16)
                        call add_cross(s, far(n), l, patch, wx*wy*hp(1:8, 9:16), cross)
                    end do
                end do
                if (order == 0) cycle
                ! scatter the target point
                call point_functionals(s, p, far(n)%xi(k, tp), far(n)%v(k, tp), f)
                call scatter_point(s, f, gx, lg)
                if (order == 2) then
                    call scatter_point_hessian(s, f, hxx, lh)
                    call scatter_cross(s, f, cross, lh)
                end if
            end do
        end do
        if (order == 0) return
        ! scatter the source points
        do patch = 1, npatch
            q = 1 + mod(patch - 1, s%n_panels)
            do n = 2, plan%ng_far
                if (.not. use(patch, n)) cycle
                do l = 1, n*n
                    call point_functionals(s, q, far(n)%xi(l, patch), far(n)%v(l, patch), f)
                    call scatter_point(s, f, gy(:, l, patch), lg)
                    if (order == 2) call scatter_point_hessian(s, f, hyy(:, :, l, patch), lh)
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
        select case (order)
        case (0)
            call generated_toroidal_sheet_self_value(u(1), u(2), u(3), u(4), u(5), &
                u(6), u(7), u(8), u(9), u(10), u(11), u(12), u(13), u(14), u(15), &
                u(16), c1, s1, c2, s2, c2 - c1, s2 - s1, out(1))
        case (1)
            call generated_toroidal_sheet_self_gradient(u(1), u(2), u(3), u(4), &
                u(5), u(6), u(7), u(8), u(9), u(10), u(11), u(12), u(13), u(14), &
                u(15), u(16), c1, s1, c2, s2, c2 - c1, s2 - s1, out(1), out(2), &
                out(3), out(4), out(5), out(6), out(7), out(8), out(9), out(10), &
                out(11), out(12), out(13), out(14), out(15), out(16), out(17))
        case default
            call self_jet(u, c1, s1, c2, s2, c2 - c1, s2 - s1, out)
        end select
    end subroutine far_leaf

    subroutine leaf_to_points(out, gp, hp, order)
        !! Gradient and Hessian of the leaf (inputs x 7, y 7, dr, dz) in the
        !! point inputs (x: R, Z, R_xi, Z_xi, R_v, Z_v, Phi_xi, Phi_v; then y).
        real(dp), intent(in) :: out(:)
        integer, intent(in) :: order
        real(dp), intent(out) :: gp(16), hp(16, 16)
        real(dp) :: m(16, 16), he(16, 16)
        integer :: i, j, k
        m = 0
        ! leaf input i <- point variables
        m(1, 1) = 1
        do i = 2, 7
            m(i, i + 1) = 1
            m(i + 7, i + 9) = 1
        end do
        m(8, 9) = 1
        m(15, 9) = 1; m(15, 1) = -1
        m(16, 10) = 1; m(16, 2) = -1
        gp = matmul(out(2:17), m)
        hp = 0
        if (order < 2) return
        k = 17
        do i = 1, 16
            do j = i, 16
                k = k + 1
                he(i, j) = out(k); he(j, i) = out(k)
            end do
        end do
        hp = matmul(transpose(m), matmul(he, m))
    end subroutine leaf_to_points

    subroutine point_rows(s, f, rows)
        !! The 8 point inputs as functionals: rows(1:8) = (R, Z, R_xi, Z_xi,
        !! R_v, Z_v, Phi_xi, Phi_v) from the 7 point functionals.
        type(ts_surface_t), intent(in) :: s
        type(functional_t), intent(in) :: f(7)
        type(functional_t), intent(inout) :: rows(8)
        integer :: k, map(8)
        map = [1, 1, 2, 3, 4, 5, 6, 7]
        do k = 1, 8
            rows(k)%n = f(map(k))%n
            rows(k)%idx = f(map(k))%idx(1:f(map(k))%n)
            rows(k)%coef = f(map(k))%coef(1:f(map(k))%n)
        end do
        rows(2)%idx = rows(2)%idx + s%n_nodes*s%nmode
    end subroutine point_rows

    subroutine scatter_point(s, f, g8, lg)
        type(ts_surface_t), intent(in) :: s
        type(functional_t), intent(in) :: f(7)
        real(dp), intent(in) :: g8(8)
        real(dp), intent(inout) :: lg(:)
        type(functional_t) :: rows(8)
        integer :: i, a
        call point_rows(s, f, rows)
        do i = 1, 8
            do a = 1, rows(i)%n
                lg(rows(i)%idx(a)) = lg(rows(i)%idx(a)) + g8(i)*rows(i)%coef(a)
            end do
        end do
    end subroutine scatter_point

    subroutine scatter_point_hessian(s, f, h8, lh)
        type(ts_surface_t), intent(in) :: s
        type(functional_t), intent(in) :: f(7)
        real(dp), intent(in) :: h8(8, 8)
        real(dp), intent(inout) :: lh(:, :)
        type(functional_t) :: rows(8)
        integer :: i, j, a, b
        call point_rows(s, f, rows)
        do j = 1, 8
            do i = 1, 8
                if (h8(i, j) == 0) cycle
                do b = 1, rows(j)%n
                    do a = 1, rows(i)%n
                        lh(rows(i)%idx(a), rows(j)%idx(b)) = &
                            lh(rows(i)%idx(a), rows(j)%idx(b)) &
                            + h8(i, j)*rows(i)%coef(a)*rows(j)%coef(b)
                    end do
                end do
            end do
        end do
    end subroutine scatter_point_hessian

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

    subroutine scatter_cross(s, f, cross, lh)
        !! lh += J_x^T cross + (J_x^T cross)^T.
        type(ts_surface_t), intent(in) :: s
        type(functional_t), intent(in) :: f(7)
        real(dp), intent(in) :: cross(:, :)
        real(dp), intent(inout) :: lh(:, :)
        type(functional_t) :: rows(8)
        integer :: i, a
        call point_rows(s, f, rows)
        do i = 1, 8
            do a = 1, rows(i)%n
                lh(rows(i)%idx(a), :) = lh(rows(i)%idx(a), :) + rows(i)%coef(a)*cross(i, :)
                lh(:, rows(i)%idx(a)) = lh(:, rows(i)%idx(a)) + rows(i)%coef(a)*cross(i, :)
            end do
        end do
    end subroutine scatter_cross

    subroutine patch_pair(s, plan, c, p, b, q, bb, order, le, lg, lh)
        !! Target patch (panel p, toroidal interval b of the first period) and
        !! source patch (q, bb over the torus). In the target's continued
        !! parameters the source occupies [o1, o1 + 1] x [o2, o2 + h].
        !! Separated patches: tensor Gauss in both. Touching or equal patches:
        !! relative coordinates z = y - x; per direction the overlap length of
        !! x and y = x + z is the tent w - |z_i - o_i|, so z is split at o_i
        !! into rectangles, those with a corner at z = 0 take a polar rule
        !! about it (the radius cancels 1/|X_y - X_x|), the others tensor
        !! Gauss, and x runs over the overlap with tensor Gauss.
        type(ts_surface_t), intent(in) :: s
        type(ts_plan_t), intent(in) :: plan
        real(dp), intent(in) :: c(:)
        integer, intent(in) :: p, b, q, bb, order
        real(dp), intent(inout) :: le, lg(:), lh(:, :)
        real(dp) :: hv, o(2), w(2), lo(2), hi(2), corner(2, 5), d(2), z(2)
        real(dp) :: gu(plan%ng_u), gw(plan%ng_u)
        real(dp) :: tq(plan%ng_t), tw(plan%ng_t), sq(plan%ng_s), sw(plan%ng_s)
        real(dp) :: det, rho, wz
        integer :: nv_all, dp_, db, k, l, m, ri, rj
        type(functional_t) :: fx(7), fy(7)

        hv = 2*pi/(s%nfp*plan%nv)
        nv_all = plan%nv*s%nfp
        dp_ = modulo(q - p + s%n_panels/2, s%n_panels) - s%n_panels/2
        db = modulo(bb - b + nv_all/2, nv_all) - nv_all/2
        call gauss_unit(plan%ng_u, gu, gw)
        if (abs(dp_) > 1 .or. abs(db) > 1) return
        call gauss_unit(plan%ng_t, tq, tw)
        call gauss_unit(plan%ng_s, sq, sw)
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
    contains
        subroutine overlap(z, t, wz)
            !! x over the overlap of the target patch and the source patch
            !! shifted by -z; X_y - X_x = t Q (t the polar radius, or 1).
            real(dp), intent(in) :: z(2), t, wz
            real(dp) :: xl(2), xh(2), xx(2)
            integer :: i1, i2
            xl = [max(0.0_dp, o(1) - z(1)), max((b - 1)*hv, (b - 1)*hv + o(2) - z(2))]
            xh = [min(1.0_dp, o(1) + 1 - z(1)), min(b*hv, (b - 1)*hv + o(2) + hv - z(2))]
            if (any(xh <= xl)) return
            do i2 = 1, plan%ng_u
                do i1 = 1, plan%ng_u
                    xx = xl + [gu(i1), gu(i2)]*(xh - xl)
                    call point_functionals(s, p, xx(1), xx(2), fx)
                    call source_point(s, c, fx, xx, q, xx(1) + z(1) - o(1), &
                        xx(2) + z(2), t, wz*gw(i1)*gw(i2)*(xh(1) - xl(1)) &
                        *(xh(2) - xl(2)), order, le, lg, lh, fy)
                end do
            end do
        end subroutine overlap
    end subroutine patch_pair

    subroutine source_point(s, c, fx, x0, q, xi_y, v_y, t, w, order, le, lg, lh, fy)
        !! One source at (xi_y, v_y) of panel q for the target fx; X_y - X_x
        !! = t Q (t = 1 for far pairs). fy is work storage.
        type(ts_surface_t), intent(in) :: s
        real(dp), intent(in) :: c(:), x0(2), xi_y, v_y, t, w
        type(functional_t), intent(in) :: fx(7)
        integer, intent(in) :: q, order
        real(dp), intent(inout) :: le, lg(:), lh(:, :)
        type(functional_t), intent(inout) :: fy(7)
        real(dp) :: u(n_self), out(1 + n_self + n_self*(n_self + 1)/2), f
        real(dp) :: gr(n_self), he(n_self, n_self), c1, s1, c2, s2, dc, ds
        integer :: i, j, k, bidx, zoff

        call point_functionals(s, q, xi_y, v_y, fy)
        zoff = s%n_nodes*s%nmode
        do i = 1, 7
            u(i) = apply(fx(i), c)
            u(7 + i) = apply(fy(i), c)
        end do
        u(15) = (u(8) - u(1))/t
        u(16) = (apply_offset(fy(1), c, zoff) - apply_offset(fx(1), c, zoff))/t
        c1 = cos(x0(2)); s1 = sin(x0(2)); c2 = cos(v_y); s2 = sin(v_y)
        ! (cos v_y - cos v_x)/t without cancellation
        dc = -2*sin(0.5_dp*(v_y - x0(2)))*sin(0.5_dp*(v_y + x0(2)))/t
        ds = 2*sin(0.5_dp*(v_y - x0(2)))*cos(0.5_dp*(v_y + x0(2)))/t
        select case (order)
        case (0)
            call generated_toroidal_sheet_self_value(u(1), u(2), u(3), u(4), u(5), &
                u(6), u(7), u(8), u(9), u(10), u(11), u(12), u(13), u(14), u(15), &
                u(16), c1, s1, c2, s2, dc, ds, f)
            le = le + w*f
            return
        case (1)
            call generated_toroidal_sheet_self_gradient(u(1), u(2), u(3), u(4), &
                u(5), u(6), u(7), u(8), u(9), u(10), u(11), u(12), u(13), u(14), &
                u(15), u(16), c1, s1, c2, s2, dc, ds, out(1), out(2), out(3), &
                out(4), out(5), out(6), out(7), out(8), out(9), out(10), out(11), &
                out(12), out(13), out(14), out(15), out(16), out(17))
        case default
            call self_jet(u, c1, s1, c2, s2, dc, ds, out)
        end select
        le = le + w*out(1)
        gr = out(2:1 + n_self)
        do i = 1, n_self
            call scatter_gradient(i, w*gr(i))
        end do
        if (order < 2) return
        bidx = 1 + n_self
        do i = 1, n_self
            do j = i, n_self
                bidx = bidx + 1
                he(i, j) = out(bidx); he(j, i) = out(bidx)
            end do
        end do
        do j = 1, n_self
            do i = 1, n_self
                if (he(i, j) /= 0) call scatter_hessian(i, j, w*he(i, j))
            end do
        end do
    contains
        subroutine input(i, part, idx, coef, n)
            !! Coefficients of input i: part 1 (and part 2 for the
            !! differences dr, dz: x enters with -1/t).
            integer, intent(in) :: i, part
            integer, intent(out) :: n
            integer, intent(out) :: idx(:)
            real(dp), intent(out) :: coef(:)
            select case (i)
            case (1:7)
                n = merge(fx(i)%n, 0, part == 1)
                if (n > 0) then
                    idx(:n) = fx(i)%idx(:n); coef(:n) = fx(i)%coef(:n)
                end if
            case (8:14)
                n = merge(fy(i - 7)%n, 0, part == 1)
                if (n > 0) then
                    idx(:n) = fy(i - 7)%idx(:n); coef(:n) = fy(i - 7)%coef(:n)
                end if
            case default
                k = merge(0, zoff, i == 15)
                if (part == 1) then
                    n = fy(1)%n
                    idx(:n) = fy(1)%idx(:n) + k; coef(:n) = fy(1)%coef(:n)/t
                else
                    n = fx(1)%n
                    idx(:n) = fx(1)%idx(:n) + k; coef(:n) = -fx(1)%coef(:n)/t
                end if
            end select
        end subroutine input

        subroutine scatter_gradient(i, gi)
            integer, intent(in) :: i
            real(dp), intent(in) :: gi
            integer :: idx(size(fx(1)%idx)), n, a, part
            real(dp) :: coef(size(fx(1)%idx))
            do part = 1, 2
                call input(i, part, idx, coef, n)
                do a = 1, n
                    lg(idx(a)) = lg(idx(a)) + gi*coef(a)
                end do
            end do
        end subroutine scatter_gradient

        subroutine scatter_hessian(i, j, hij)
            integer, intent(in) :: i, j
            real(dp), intent(in) :: hij
            integer :: ii(size(fx(1)%idx)), jj(size(fx(1)%idx)), ni, nj, a, b2, pi_, pj
            real(dp) :: ci(size(fx(1)%idx)), cj(size(fx(1)%idx))
            do pj = 1, 2
                call input(j, pj, jj, cj, nj)
                if (nj == 0) cycle
                do pi_ = 1, 2
                    call input(i, pi_, ii, ci, ni)
                    do b2 = 1, nj
                        do a = 1, ni
                            lh(ii(a), jj(b2)) = lh(ii(a), jj(b2)) + hij*ci(a)*cj(b2)
                        end do
                    end do
                end do
            end do
        end subroutine scatter_hessian
    end subroutine source_point

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
        !! int phi B . (X_xi x X_v) dxi dv over the torus (phi the single-valued
        !! part), and when present its gradient and Hessian in the
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
                        phi = apply(fphi, c)
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
