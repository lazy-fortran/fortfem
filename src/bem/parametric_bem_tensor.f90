module fortfem_parametric_bem_tensor
    !! Tensor-product trace spaces for the scalar Laplace Galerkin BEM of
    !! fortfem_parametric_laplace_bem_3d on quadrilateral parametric panels.
    !! The surface contract, kernels, singular treatment and DtN builder are
    !! the existing ones: `bem_surface_t%evaluate` is called on the reference
    !! square [0,1]^2, coincident and close pairs use a Duffy fan about the
    !! target point or about its closest point on the source panel, far pairs
    !! use the tensor Gauss product rule, the double layer integrates
    !! d_nG (u_y-u_x) and restores -u_x/2, and the hypersingular form uses
    !! surface curls. Potentials are continuous Q_p (Lagrange at caller nodes,
    !! normally the GLL nodes of a volume trace); fluxes are discontinuous
    !! Q_pf at Gauss nodes. The result fills `bem_operators_t`, so
    !! `build_exterior_dtn` applies unchanged.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortfem_kinds, only: dp
    use fortfem_parametric_laplace_bem_3d, only: bem_surface_t, bem_operators_t
    use fortfem_generated_bem_frame, only: bem_frame
    use fortfem_generated_bem_kernel, only: bem_kernel
    use fortfem_generated_bem_field_kernel, only: bem_field_kernel
    use fortfem_gll_line, only: gauss_rule, lagrange_jet
    implicit none
    private
    public :: bem_tensor_options_t, bem_tensor_space_t, build_bem_tensor_space
    public :: assemble_parametric_bem_tensor, bem_tensor_exterior_field
    public :: bem_tensor_flux_l2_error, bem_tensor_boundary_load
    public :: bem_tensor_trace_l2_error, boundary_function

    type :: bem_tensor_options_t
        !! Gauss points per direction for target and far source integration
        !! (0: trace degree + 3) and per fan coordinate (0: trace degree + 5).
        integer :: target_points = 0, fan_points = 0
        real(dp) :: near_factor = 1.5_dp
    end type

    type :: bem_tensor_space_t
        integer :: trace_degree = 0, flux_degree = 0
        integer :: npanel = 0, ntrace = 0, nflux = 0
        real(dp), allocatable :: trace_nodes(:), flux_nodes(:)
        !! trace_dofs(i + (p+1)*(j-1), panel): global trace index at node (i,j).
        integer, allocatable :: trace_dofs(:, :)
    end type

    type :: rule_t
        integer :: n = 0
        real(dp), allocatable :: uv(:, :), w(:), phi(:, :), du(:, :), dv(:, :)
        real(dp), allocatable :: psi(:, :)
    end type

    type :: panel_points_t
        real(dp), allocatable :: x(:, :, :), a(:, :, :), b(:, :, :)
        real(dp), allocatable :: normal(:, :, :), jacobian(:, :)
    end type

    abstract interface
        function boundary_function(x, normal) result(value)
            import :: dp
            real(dp), intent(in) :: x(3), normal(3)
            real(dp) :: value
        end function
    end interface

    real(dp), parameter :: scale = 1.0_dp/(4.0_dp*acos(-1.0_dp))
    real(dp), parameter :: square(2, 4) = reshape([0.0_dp, 0.0_dp, 1.0_dp, &
        0.0_dp, 1.0_dp, 1.0_dp, 0.0_dp, 1.0_dp], [2, 4])
contains
    subroutine build_bem_tensor_space(trace_nodes, flux_degree, trace_dofs, &
                                      space, status)
        !! trace_nodes on [0,1], increasing, including both end points so that
        !! neighbouring panels share edge nodes; flux_degree >= 0.
        real(dp), intent(in) :: trace_nodes(:)
        integer, intent(in) :: flux_degree, trace_dofs(:, :)
        type(bem_tensor_space_t), intent(out) :: space
        integer, intent(out) :: status
        real(dp), allocatable :: x(:), w(:)
        integer :: n
        status = 1
        n = size(trace_nodes)
        if (n < 2 .or. flux_degree < 0) return
        if (abs(trace_nodes(1)) > 1.e-14_dp) return
        if (abs(trace_nodes(n) - 1) > 1.e-14_dp) return
        if (any(trace_nodes(2:) <= trace_nodes(:n - 1))) return
        if (size(trace_dofs, 1) /= n*n .or. size(trace_dofs, 2) < 4) return
        if (any(trace_dofs < 1)) return
        call gauss_rule(flux_degree + 1, x, w, status)
        if (status /= 0) return
        space%trace_degree = n - 1
        space%flux_degree = flux_degree
        space%npanel = size(trace_dofs, 2)
        space%ntrace = maxval(trace_dofs)
        space%nflux = (flux_degree + 1)**2*space%npanel
        space%trace_nodes = trace_nodes
        space%flux_nodes = (x + 1)/2
        space%trace_dofs = trace_dofs
        status = 0
    end subroutine

    pure subroutine basis_at(space, uv, phi, du, dv, psi)
        type(bem_tensor_space_t), intent(in) :: space
        real(dp), intent(in) :: uv(2)
        real(dp), intent(out) :: phi(:), du(:), dv(:), psi(:)
        real(dp) :: bu(space%trace_degree + 1), bdu(space%trace_degree + 1)
        real(dp) :: bv(space%trace_degree + 1), bdv(space%trace_degree + 1)
        real(dp) :: fu(space%flux_degree + 1), fdu(space%flux_degree + 1)
        real(dp) :: fv(space%flux_degree + 1), fdv(space%flux_degree + 1)
        integer :: i, j, n, m
        n = space%trace_degree + 1
        m = space%flux_degree + 1
        call lagrange_jet(space%trace_nodes, uv(1), bu, bdu)
        call lagrange_jet(space%trace_nodes, uv(2), bv, bdv)
        call lagrange_jet(space%flux_nodes, uv(1), fu, fdu)
        call lagrange_jet(space%flux_nodes, uv(2), fv, fdv)
        do j = 1, n
            do i = 1, n
                phi(i + n*(j - 1)) = bu(i)*bv(j)
                du(i + n*(j - 1)) = bdu(i)*bv(j)
                dv(i + n*(j - 1)) = bu(i)*bdv(j)
            end do
        end do
        do j = 1, m
            do i = 1, m
                psi(i + m*(j - 1)) = fu(i)*fv(j)
            end do
        end do
    end subroutine

    subroutine build_rule(space, n, rule, status)
        type(bem_tensor_space_t), intent(in) :: space
        integer, intent(in) :: n
        type(rule_t), intent(out) :: rule
        integer, intent(out) :: status
        real(dp), allocatable :: x(:), w(:)
        integer :: i, j, t, nt, nf
        call gauss_rule(n, x, w, status)
        if (status /= 0) return
        nt = (space%trace_degree + 1)**2
        nf = (space%flux_degree + 1)**2
        rule%n = n*n
        allocate (rule%uv(2, n*n), rule%w(n*n), rule%phi(nt, n*n), &
                  rule%du(nt, n*n), rule%dv(nt, n*n), rule%psi(nf, n*n))
        do j = 1, n
            do i = 1, n
                t = i + n*(j - 1)
                rule%uv(:, t) = [(x(i) + 1)/2, (x(j) + 1)/2]
                rule%w(t) = w(i)*w(j)/4
                call basis_at(space, rule%uv(:, t), rule%phi(:, t), rule%du(:, t), &
                              rule%dv(:, t), rule%psi(:, t))
            end do
        end do
    end subroutine

    subroutine sample_panels(surface, rule, npanel, points, status)
        class(bem_surface_t), intent(in) :: surface
        type(rule_t), intent(in) :: rule
        integer, intent(in) :: npanel
        type(panel_points_t), intent(out) :: points
        integer, intent(out) :: status
        integer :: k, t, errors, local
        allocate (points%x(3, rule%n, npanel), points%a(3, rule%n, npanel), &
                  points%b(3, rule%n, npanel), points%normal(3, rule%n, npanel), &
                  points%jacobian(rule%n, npanel))
        errors = 0
        !$omp parallel do default(shared) private(k,t,local) reduction(+:errors)
        do k = 1, npanel
            do t = 1, rule%n
                call sample_point(surface, k, rule%uv(:, t), points%x(:, t, k), &
                                  points%a(:, t, k), points%b(:, t, k), &
                                  points%normal(:, t, k), points%jacobian(t, k), local)
                if (local /= 0) errors = errors + 1
            end do
        end do
        !$omp end parallel do
        status = merge(0, 2, errors == 0)
    end subroutine

    pure subroutine sample_point(surface, panel, uv, x, a, b, normal, jacobian, status)
        class(bem_surface_t), intent(in) :: surface
        integer, intent(in) :: panel
        real(dp), intent(in) :: uv(2)
        real(dp), intent(out) :: x(3), a(3), b(3), normal(3), jacobian
        integer, intent(out) :: status
        real(dp) :: norm_a, norm_b
        normal = 0
        jacobian = 0
        call surface%evaluate(panel, uv, x, a, b, status)
        if (status /= 0) return
        status = 2
        if (.not. all(ieee_is_finite(x))) return
        if (.not. all(ieee_is_finite(a))) return
        if (.not. all(ieee_is_finite(b))) return
        norm_a = norm2(a)
        norm_b = norm2(b)
        if (min(norm_a, norm_b) <= tiny(1.0_dp)) return
        if (abs(dot_product(a, b)) >= &
            (1.0_dp - 64.0_dp*epsilon(1.0_dp))*norm_a*norm_b) return
        call bem_frame(a(1), a(2), a(3), b(1), b(2), b(3), jacobian, &
                       normal(1), normal(2), normal(3))
        status = 0
    end subroutine

    subroutine assemble_parametric_bem_tensor(surface, space, options, op, status)
        !! Dense V (flux x flux), K and mass (flux x trace), W (trace x trace).
        !! Status 1 rejects inputs, 2 geometry or nonfinite values. `area` is
        !! the panel area. `genus` is not set here; the caller owns topology.
        class(bem_surface_t), intent(in) :: surface
        type(bem_tensor_space_t), intent(in) :: space
        type(bem_tensor_options_t), intent(in) :: options
        type(bem_operators_t), intent(out) :: op
        integer, intent(out) :: status
        type(rule_t) :: rule
        type(panel_points_t) :: pts
        real(dp), allocatable :: centre(:, :), diameter(:), fan_x(:), fan_w(:)
        integer :: nq, nl, np, nt, nf, i, k, errors
        status = 1
        if (space%npanel < 4 .or. space%ntrace < 1) return
        if (options%near_factor < 0) return
        nq = options%target_points
        if (nq == 0) nq = space%trace_degree + 3
        nl = options%fan_points
        if (nl == 0) nl = space%trace_degree + 5
        if (nq < 1 .or. nl < 1) return
        np = space%npanel
        nt = (space%trace_degree + 1)**2
        nf = (space%flux_degree + 1)**2
        call build_rule(space, nq, rule, status)
        if (status /= 0) return
        call gauss_rule(nl, fan_x, fan_w, status)
        if (status /= 0) return
        fan_x = (fan_x + 1)/2
        fan_w = fan_w/2
        call sample_panels(surface, rule, np, pts, status)
        if (status /= 0) return
        status = 2
        allocate (centre(3, np), diameter(np))
        allocate (op%v(space%nflux, space%nflux), op%k(space%nflux, space%ntrace), &
                  op%mass(space%nflux, space%ntrace), &
                  op%w(space%ntrace, space%ntrace), &
                  op%area(np))
        op%v = 0
        op%k = 0
        op%mass = 0
        op%w = 0
        do k = 1, np
            centre(:, k) = matmul(pts%x(:, :, k), rule%w)/sum(rule%w)
            diameter(k) = 0
            do i = 1, rule%n
                diameter(k) = max(diameter(k), norm2(pts%x(:, i, k) - centre(:, k)))
            end do
            op%area(k) = dot_product(rule%w, pts%jacobian(:, k))
            do i = 1, rule%n
                op%mass(nf*(k - 1) + 1:nf*k, space%trace_dofs(:, k)) = &
                    op%mass(nf*(k - 1) + 1:nf*k, space%trace_dofs(:, k)) + &
                    rule%w(i)*pts%jacobian(i, k)* &
                    spread(rule%psi(:, i), 2, nt)*spread(rule%phi(:, i), 1, nf)
            end do
        end do
        errors = 0
        !$omp parallel default(shared) reduction(+:errors)
        block
            real(dp), allocatable :: vloc(:, :), kloc(:, :), kself(:, :), wloc(:, :)
            real(dp), allocatable :: vs(:), ks(:), ws(:, :), ct(:, :)
            real(dp) :: k0, wt
            integer, allocatable :: rows(:), cols(:)
            integer :: j, t, a, b, local
            logical :: near
            allocate (vloc(nf, nf), kloc(nf, nt), kself(nf, nt), wloc(nt, nt), &
                      vs(nf), ks(nt), ws(3, nt), ct(3, nt), rows(nf), cols(nf))
            !$omp do schedule(dynamic, 1)
            do i = 1, np
                rows = [(nf*(i - 1) + a, a=1, nf)]
                do j = 1, np
                    near = norm2(centre(:, i) - centre(:, j)) < &
                           options%near_factor*(diameter(i) + diameter(j))
                    vloc = 0
                    kloc = 0
                    kself = 0
                    wloc = 0
                    do t = 1, rule%n
                        if (near) then
                            call source_fan(surface, space, j, pts%x(:, t, i), &
                                            rule%uv(:, t), i == j, fan_x, fan_w, &
                                            vs, ks, k0, ws, local)
                            if (local /= 0) errors = errors + 1
                        else
                            call source_regular(rule, pts, j, pts%x(:, t, i), &
                                                vs, ks, k0, ws)
                        end if
                        wt = rule%w(t)*pts%jacobian(t, i)
                        do a = 1, nt
                            ct(:, a) = pts%b(:, t, i)*rule%du(a, t) - &
                                       pts%a(:, t, i)*rule%dv(a, t)
                        end do
                        do b = 1, nf
                            vloc(:, b) = vloc(:, b) + wt*rule%psi(:, t)*vs(b)
                        end do
                        do b = 1, nt
                            kloc(:, b) = kloc(:, b) + wt*rule%psi(:, t)*ks(b)
                            kself(:, b) = kself(:, b) + &
                                          wt*k0*rule%psi(:, t)*rule%phi(b, t)
                            do a = 1, nt
                                wloc(a, b) = wloc(a, b) + rule%w(t)* &
                                             dot_product(ct(:, a), ws(:, b))
                            end do
                        end do
                    end do
                    cols = [(nf*(j - 1) + a, a=1, nf)]
                    op%v(rows, cols) = op%v(rows, cols) + vloc
                    op%k(rows, space%trace_dofs(:, j)) = &
                        op%k(rows, space%trace_dofs(:, j)) + kloc
                    op%k(rows, space%trace_dofs(:, i)) = &
                        op%k(rows, space%trace_dofs(:, i)) - kself
                    do b = 1, nt
                        do a = 1, nt
                            !$omp atomic update
                            op%w(space%trace_dofs(a, i), space%trace_dofs(b, j)) = &
                                op%w(space%trace_dofs(a, i), &
                                     space%trace_dofs(b, j)) + wloc(a, b)
                        end do
                    end do
                end do
            end do
            !$omp end do
        end block
        !$omp end parallel
        if (errors /= 0) return
        op%k = op%k - 0.5_dp*op%mass
        op%v = 0.5_dp*(op%v + transpose(op%v))
        op%w = 0.5_dp*(op%w + transpose(op%w))
        op%adjoint = transpose(op%k)
        if (.not. all(ieee_is_finite(op%v))) return
        if (.not. all(ieee_is_finite(op%k))) return
        if (.not. all(ieee_is_finite(op%w))) return
        status = 0
    end subroutine

    pure subroutine accumulate(xt, x, a, b, normal, jacobian, weight, phi, du, dv, &
                               psi, vs, ks, k0, ws)
        !! weight is the reference weight; curl densities carry the Jacobian.
        real(dp), intent(in) :: xt(3), x(3), a(3), b(3), normal(3), jacobian, weight
        real(dp), intent(in) :: phi(:), du(:), dv(:), psi(:)
        real(dp), intent(inout) :: vs(:), ks(:), k0, ws(:, :)
        real(dp) :: d(3), green, double, g, h
        integer :: m
        d = xt - x
        call bem_kernel(d(1), d(2), d(3), normal(1), normal(2), normal(3), &
                        green, double)
        g = scale*green*weight
        h = scale*double*weight*jacobian
        vs = vs + g*jacobian*psi
        ks = ks + h*phi
        k0 = k0 + h
        do m = 1, size(phi)
            ws(:, m) = ws(:, m) + g*(b*du(m) - a*dv(m))
        end do
    end subroutine

    pure subroutine source_regular(rule, pts, panel, xt, vs, ks, k0, ws)
        type(rule_t), intent(in) :: rule
        type(panel_points_t), intent(in) :: pts
        integer, intent(in) :: panel
        real(dp), intent(in) :: xt(3)
        real(dp), intent(out) :: vs(:), ks(:), k0, ws(:, :)
        integer :: s
        vs = 0
        ks = 0
        k0 = 0
        ws = 0
        do s = 1, rule%n
            call accumulate(xt, pts%x(:, s, panel), pts%a(:, s, panel), &
                            pts%b(:, s, panel), pts%normal(:, s, panel), &
                            pts%jacobian(s, panel), rule%w(s), rule%phi(:, s), &
                            rule%du(:, s), rule%dv(:, s), rule%psi(:, s), vs, ks, &
                            k0, ws)
        end do
    end subroutine

    pure subroutine source_fan(surface, space, panel, xt, uvt, same, nodes, weights, &
                               vs, ks, k0, ws, status)
        !! Four Duffy triangles from the apex to the square edges; the radial
        !! coordinate r = s^2 clusters nodes at the apex.
        class(bem_surface_t), intent(in) :: surface
        type(bem_tensor_space_t), intent(in) :: space
        integer, intent(in) :: panel
        real(dp), intent(in) :: xt(3), uvt(2), nodes(:), weights(:)
        logical, intent(in) :: same
        real(dp), intent(out) :: vs(:), ks(:), k0, ws(:, :)
        integer, intent(out) :: status
        real(dp) :: apex(2), e1(2), e2(2), det, uv(2), r, factor
        real(dp) :: x(3), a(3), b(3), normal(3), jacobian
        real(dp) :: phi(size(ks)), du(size(ks)), dv(size(ks)), psi(size(vs))
        integer :: edge, i, j
        vs = 0
        ks = 0
        k0 = 0
        ws = 0
        status = 0
        apex = uvt
        if (.not. same) then
            call closest_square_point(surface, panel, xt, apex, status)
            if (status /= 0) return
        end if
        do edge = 1, 4
            e1 = square(:, edge) - apex
            e2 = square(:, modulo(edge, 4) + 1) - apex
            det = abs(e1(1)*e2(2) - e1(2)*e2(1))
            if (det < 64.0_dp*epsilon(1.0_dp)) cycle
            do i = 1, size(nodes)
                r = nodes(i)**2
                do j = 1, size(nodes)
                    uv = apex + r*((1.0_dp - nodes(j))*e1 + nodes(j)*e2)
                    call sample_point(surface, panel, uv, x, a, b, normal, jacobian, &
                                      status)
                    if (status /= 0) return
                    call basis_at(space, uv, phi, du, dv, psi)
                    factor = weights(i)*weights(j)*2.0_dp*nodes(i)**3*det
                    call accumulate(xt, x, a, b, normal, jacobian, factor, phi, du, &
                                    dv, psi, vs, ks, k0, ws)
                end do
            end do
        end do
    end subroutine

    pure subroutine closest_square_point(surface, panel, x, uv, status)
        class(bem_surface_t), intent(in) :: surface
        integer, intent(in) :: panel
        real(dp), intent(in) :: x(3)
        real(dp), intent(out) :: uv(2)
        integer, intent(out) :: status
        real(dp) :: p(3), a(3), b(3), d(3), aa, ab, bb, det, step(2)
        integer :: iteration
        uv = [0.5_dp, 0.5_dp]
        do iteration = 1, 12
            call surface%evaluate(panel, uv, p, a, b, status)
            if (status /= 0) return
            status = 2
            if (.not. all(ieee_is_finite(p))) return
            d = p - x
            aa = dot_product(a, a)
            ab = dot_product(a, b)
            bb = dot_product(b, b)
            det = aa*bb - ab*ab
            if (det <= tiny(1.0_dp)) return
            step = [(bb*dot_product(a, d) - ab*dot_product(b, d))/det, &
                    (aa*dot_product(b, d) - ab*dot_product(a, d))/det]
            uv = min(1.0_dp, max(0.0_dp, uv - step))
        end do
        status = 0
    end subroutine

    subroutine bem_tensor_exterior_field(surface, space, trace, flux, points, &
                                         potential, gradient, status, quadrature)
        !! Green representation u(x) = int u d_nG - G q for x off the surface.
        !! For exterior x, trace and flux are the exterior Dirichlet data and
        !! their derivative along the normal out of the bounded region; for
        !! exact exterior data the integral vanishes inside that region.
        !! Regular tensor Gauss only: points must be a few panels away.
        class(bem_surface_t), intent(in) :: surface
        type(bem_tensor_space_t), intent(in) :: space
        real(dp), intent(in) :: trace(:), flux(:), points(:, :)
        real(dp), intent(out) :: potential(:), gradient(:, :)
        integer, intent(out) :: status
        integer, intent(in), optional :: quadrature
        type(rule_t) :: rule
        type(panel_points_t) :: pts
        real(dp), allocatable :: ut(:, :), qt(:, :)
        integer :: n, k, t, nf
        status = 1
        if (size(trace) /= space%ntrace .or. size(flux) /= space%nflux) return
        if (size(points, 1) /= 3 .or. size(potential) /= size(points, 2)) return
        if (any(shape(gradient) /= [3, size(points, 2)])) return
        n = 2*space%trace_degree + 4
        if (present(quadrature)) n = quadrature
        call build_rule(space, n, rule, status)
        if (status /= 0) return
        call sample_panels(surface, rule, space%npanel, pts, status)
        if (status /= 0) return
        nf = (space%flux_degree + 1)**2
        allocate (ut(rule%n, space%npanel), qt(rule%n, space%npanel))
        do k = 1, space%npanel
            do t = 1, rule%n
                ut(t, k) = dot_product(rule%phi(:, t), trace(space%trace_dofs(:, k)))
                qt(t, k) = dot_product(rule%psi(:, t), flux(nf*(k - 1) + 1:nf*k))
            end do
        end do
        !$omp parallel do default(shared) schedule(static)
        do k = 1, size(points, 2)
            call representation(points(:, k), potential(k), gradient(:, k))
        end do
        !$omp end parallel do
        status = 0
    contains
        subroutine representation(x, u, g)
            real(dp), intent(in) :: x(3)
            real(dp), intent(out) :: u, g(3)
            real(dp) :: d(3), green, double, gg(3), gd(3), w
            integer :: panel, s
            u = 0
            g = 0
            do panel = 1, space%npanel
                do s = 1, rule%n
                    d = x - pts%x(:, s, panel)
                    call bem_field_kernel(d(1), d(2), d(3), pts%normal(1, s, panel), &
                                          pts%normal(2, s, panel), &
                                          pts%normal(3, s, panel), green, double, &
                                          gg(1), gg(2), gg(3), gd(1), gd(2), gd(3))
                    w = scale*rule%w(s)*pts%jacobian(s, panel)
                    u = u + w*(ut(s, panel)*double - qt(s, panel)*green)
                    g = g + w*(ut(s, panel)*gd - qt(s, panel)*gg)
                end do
            end do
        end subroutine
    end subroutine

    subroutine bem_tensor_flux_l2_error(surface, space, flux, exact, error, &
                                        status, quadrature)
        !! L2(surface) distance between the discrete flux and exact(x, normal).
        class(bem_surface_t), intent(in) :: surface
        type(bem_tensor_space_t), intent(in) :: space
        real(dp), intent(in) :: flux(:)
        procedure(boundary_function) :: exact
        real(dp), intent(out) :: error
        integer, intent(out) :: status
        integer, intent(in), optional :: quadrature
        type(rule_t) :: rule
        type(panel_points_t) :: pts
        integer :: n, k, t, nf
        error = huge(1.0_dp)
        status = 1
        if (size(flux) /= space%nflux) return
        n = 2*space%trace_degree + 4
        if (present(quadrature)) n = quadrature
        call build_rule(space, n, rule, status)
        if (status /= 0) return
        call sample_panels(surface, rule, space%npanel, pts, status)
        if (status /= 0) return
        nf = (space%flux_degree + 1)**2
        error = 0
        do k = 1, space%npanel
            do t = 1, rule%n
                error = error + rule%w(t)*pts%jacobian(t, k)*(dot_product( &
                        rule%psi(:, t), flux(nf*(k - 1) + 1:nf*k)) - &
                        exact(pts%x(:, t, k), pts%normal(:, t, k)))**2
            end do
        end do
        error = sqrt(error)
    end subroutine

    subroutine bem_tensor_trace_l2_error(surface, space, trace, exact, error, &
                                         status, quadrature)
        !! L2(surface) distance between the continuous trace and exact.
        class(bem_surface_t), intent(in) :: surface
        type(bem_tensor_space_t), intent(in) :: space
        real(dp), intent(in) :: trace(:)
        procedure(boundary_function) :: exact
        real(dp), intent(out) :: error
        integer, intent(out) :: status
        integer, intent(in), optional :: quadrature
        type(rule_t) :: rule
        type(panel_points_t) :: pts
        integer :: n, k, t
        error = huge(1.0_dp)
        status = 1
        if (size(trace) /= space%ntrace) return
        n = 2*space%trace_degree + 4
        if (present(quadrature)) n = quadrature
        call build_rule(space, n, rule, status)
        if (status /= 0) return
        call sample_panels(surface, rule, space%npanel, pts, status)
        if (status /= 0) return
        error = 0
        do k = 1, space%npanel
            do t = 1, rule%n
                error = error + rule%w(t)*pts%jacobian(t, k)*(dot_product( &
                        rule%phi(:, t), trace(space%trace_dofs(:, k))) - &
                        exact(pts%x(:, t, k), pts%normal(:, t, k)))**2
            end do
        end do
        error = sqrt(error)
    end subroutine

    subroutine bem_tensor_boundary_load(surface, space, g, load, status, quadrature)
        !! load(i) = int g(x, n) phi_i over the surface, in trace numbering.
        class(bem_surface_t), intent(in) :: surface
        type(bem_tensor_space_t), intent(in) :: space
        procedure(boundary_function) :: g
        real(dp), allocatable, intent(out) :: load(:)
        integer, intent(out) :: status
        integer, intent(in), optional :: quadrature
        type(rule_t) :: rule
        type(panel_points_t) :: pts
        integer :: n, k, t
        n = space%trace_degree + 3
        if (present(quadrature)) n = quadrature
        call build_rule(space, n, rule, status)
        if (status /= 0) return
        call sample_panels(surface, rule, space%npanel, pts, status)
        if (status /= 0) return
        allocate (load(space%ntrace))
        load = 0
        do k = 1, space%npanel
            do t = 1, rule%n
                load(space%trace_dofs(:, k)) = load(space%trace_dofs(:, k)) + &
                    rule%w(t)*pts%jacobian(t, k)*rule%phi(:, t)* &
                    g(pts%x(:, t, k), pts%normal(:, t, k))
            end do
        end do
    end subroutine
end module
