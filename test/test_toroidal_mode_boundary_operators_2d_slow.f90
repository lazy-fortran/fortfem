program test_toroidal_mode_boundary_operators_2d_slow
    !! Independent oracles for the toroidal-mode boundary operators on curved
    !! P2 panels (n = 0, 1, 3, 6):
    !!  - structure: K' = K^T, symmetry and definiteness of V, W and S_n,
    !!    K 1 = -1/2 and W 1 = 0 for n = 0,
    !!  - Calderon identities on the Cauchy data of an interior-source field
    !!    (V t = (K - 1/2) u and W u = -(1/2 + K') t) and the Galerkin
    !!    identity V W = 1/4 - K^2,
    !!  - the exterior Dirichlet problem for that field: Neumann data, energy
    !!    and observed convergence orders,
    !!  - n = 1 against the existing Grad-Shafranov operators through
    !!    psi = R A (energies of the same physical field),
    !!  - the Johnson-Nedelec form and the Neumann-to-Dirichlet map.
    use check, only: check_condition, check_summary
    use fortfem_kinds, only: dp
    use fortfem_toroidal_mode_kernel, only: tm_green, tm_green_conormal
    use fortfem_axisymmetric_gs_kernel, only: gs_green_gradient
    use fortfem_axisymmetric_gs_boundary_operators_2d, only: &
        gs_bem_curve_t, gs_bem_operators_t, gs_bem_curve_from_nodes, &
        gs_bem_assemble, gs_bem_steklov_symmetric, gs_bem_dual_norm_squared, &
        gs_bem_panel_point, gs_bem_flux_value, gs_bem_trace_node_index
    use fortfem_toroidal_mode_boundary_operators_2d, only: tm_bem_assemble, &
        tm_bem_steklov_symmetric, tm_bem_steklov_jn, tm_bem_exterior_neumann, &
        tm_bem_neumann_to_dirichlet
    implicit none

    real(dp), parameter :: pi = acos(-1.0_dp)
    real(dp), parameter :: circle_center(2) = [3.0_dp, 0.5_dp]
    real(dp), parameter :: circle_source(2) = [3.25_dp, 0.62_dp]
    real(dp), parameter :: d_source(2) = [3.3_dp, 0.3_dp]
    integer, parameter :: orders(4) = [0, 1, 3, 6]
    logical :: all_passed

    interface
        subroutine dsyev(jobz, uplo, n, a, lda, w, work, lwork, info)
            import :: dp
            character(len=1), intent(in) :: jobz, uplo
            integer, intent(in) :: n, lda, lwork
            real(dp), intent(inout) :: a(lda, *)
            real(dp), intent(out) :: w(*), work(*)
            integer, intent(out) :: info
        end subroutine dsyev
    end interface

    all_passed = .true.

    call test_structure(2, 0)
    call test_structure(2, 1)
    call test_structure(2, 3)
    call test_structure(2, 6)
    call test_cauchy_data(1)
    call test_cauchy_data(3)
    call test_cauchy_data(6)
    call test_calderon(2, 0)
    call test_calderon(2, 1)
    call test_calderon(2, 3)
    call test_calderon(2, 6)
    call test_exterior_dirichlet(1, 0)
    call test_exterior_dirichlet(1, 1)
    call test_exterior_dirichlet(1, 3)
    call test_exterior_dirichlet(2, 6)
    call test_exterior_dirichlet(1, 6)
    call test_matches_grad_shafranov()
    call test_jn_and_ntd(3)

    call check_summary("Toroidal-mode boundary operators")
    if (.not. all_passed) error stop 1

contains

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

    pure function curve_point(shape, theta) result(x)
        !! 1: circle of radius 1 about (3, 0.5); 2: D shape
        !! R = 3 + cos(theta + 0.4 sin theta), Z = 0.5 + 1.4 sin(theta).
        integer, intent(in) :: shape
        real(dp), intent(in) :: theta
        real(dp) :: x(2)

        if (shape == 1) then
            x = circle_center + [cos(theta), sin(theta)]
        else
            x = [3.0_dp + cos(theta + 0.4_dp*sin(theta)), &
                0.5_dp + 1.4_dp*sin(theta)]
        end if
    end function curve_point

    subroutine make_curve(shape, np, curve)
        integer, intent(in) :: shape, np
        type(gs_bem_curve_t), intent(out) :: curve
        real(dp), allocatable :: nodes(:, :)
        integer :: j, status

        allocate(nodes(2, 2*np))
        do j = 1, 2*np
            nodes(:, j) = curve_point(shape, pi*real(j - 1, dp)/np)
        end do
        call gs_bem_curve_from_nodes(nodes, curve, status)
        if (status /= 0) error stop "curve construction failed"
    end subroutine make_curve

    function source_point(shape) result(y)
        integer, intent(in) :: shape
        real(dp) :: y(2)

        if (shape == 1) then
            y = circle_source
        else
            y = d_source
        end if
    end function source_point

    subroutine gauss_unit(n, x, w)
        use fortnum_quadrature, only: gauss_legendre
        integer, intent(in) :: n
        real(dp), intent(out) :: x(n), w(n)

        call gauss_legendre(n, x, w)
        x = 0.5_dp*(x + 1.0_dp)
        w = 0.5_dp*w
    end subroutine gauss_unit

    ! ------------------------------------------------------------------
    ! Exact exterior data of the n-mode field of an interior source
    ! ------------------------------------------------------------------

    subroutine field_data(n, src, x, tau, u, t)
        !! Trace u = G_n(x, src) and t = R n . grad G_n, n = (tau_Z, -tau_R).
        integer, intent(in) :: n
        real(dp), intent(in) :: src(2), x(2), tau(2)
        real(dp), intent(out) :: u, t
        real(dp) :: nrm(2), dsrc, dboth

        nrm = [tau(2), -tau(1)]
        u = tm_green(n, x(1), x(2), src(1), src(2))
        call tm_green_conormal(n, x(1), x(2), src(1), src(2), nrm, nrm, t, &
            dsrc, dboth)
    end subroutine field_data

    subroutine nodal_trace(n, curve, src, u)
        integer, intent(in) :: n
        type(gs_bem_curve_t), intent(in) :: curve
        real(dp), intent(in) :: src(2)
        real(dp), allocatable, intent(out) :: u(:)
        real(dp) :: x(2), tau(2), sp, t, xi
        integer :: p, k

        allocate(u(2*curve%n_panels))
        do p = 1, curve%n_panels
            do k = 1, 3
                xi = 0.5_dp*(k - 1)
                call gs_bem_panel_point(curve, p, xi, x, tau, sp)
                call field_data(n, src, x, tau, u(gs_bem_trace_node_index( &
                    curve, p, k)), t)
            end do
        end do
    end subroutine nodal_trace

    subroutine nodal_flux(n, curve, deg, src, t)
        !! Flux dofs: values of the exact conormal flux at the equispaced
        !! nodes of each panel.
        integer, intent(in) :: n, deg
        type(gs_bem_curve_t), intent(in) :: curve
        real(dp), intent(in) :: src(2)
        real(dp), allocatable, intent(out) :: t(:)
        real(dp) :: x(2), tau(2), sp, u, xi
        integer :: p, k

        allocate(t((deg + 1)*curve%n_panels))
        do p = 1, curve%n_panels
            do k = 1, deg + 1
                xi = real(k - 1, dp)/max(deg, 1)
                call gs_bem_panel_point(curve, p, xi, x, tau, sp)
                call field_data(n, src, x, tau, u, t((p - 1)*(deg + 1) + k))
            end do
        end do
    end subroutine nodal_flux

    function eigenvalues(a) result(w)
        real(dp), intent(in) :: a(:, :)
        real(dp), allocatable :: w(:)
        real(dp), allocatable :: b(:, :), work(:)
        integer :: n, info

        n = size(a, 1)
        b = a
        allocate(w(n), work(max(1, 3*n)))
        call dsyev('N', 'U', n, b, n, w, work, size(work), info)
        if (info /= 0) w = -huge(1.0_dp)
    end function eigenvalues

    ! ------------------------------------------------------------------
    ! Tests
    ! ------------------------------------------------------------------

    subroutine test_structure(shape, n)
        integer, intent(in) :: shape, n
        type(gs_bem_curve_t) :: curve
        type(gs_bem_operators_t) :: ops
        real(dp), allocatable :: ones_x(:), r(:), s(:, :), ev(:)
        real(dp) :: wnorm, scale
        integer :: status
        character(len=24) :: tag

        write (tag, '(a,i0)') "n=", n
        call make_curve(shape, 32, curve)
        call tm_bem_assemble(curve, n, 1, ops, status)
        call record_condition(status == 0, "assembly status, "//trim(tag))
        allocate(ones_x(ops%nx))
        ones_x = 1.0_dp
        if (n == 0) then
            r = matmul(ops%k, ones_x) + 0.5_dp*matmul(ops%m, ones_x)
            call record_condition(maxval(abs(r)) < 1.0e-8_dp, &
                "K 1 = -1/2 (interior jump), "//trim(tag))
            call record_condition(maxval(abs(matmul(ops%w, ones_x))) &
                < 1.0e-10_dp*maxval(abs(ops%w)), "W 1 = 0, "//trim(tag))
        end if
        call record_condition(maxval(abs(ops%kt - transpose(ops%k))) &
            < 1.0e-12_dp*maxval(abs(ops%k)), &
            "K' assembled from d_nx G equals K^T, "//trim(tag))
        scale = maxval(abs(ops%v))
        call record_condition(maxval(abs(ops%v - transpose(ops%v))) &
            < 1.0e-12_dp*scale, "V symmetric, "//trim(tag))
        wnorm = maxval(abs(ops%w))
        call record_condition(maxval(abs(ops%w - transpose(ops%w))) &
            < 1.0e-6_dp*wnorm, "W symmetric (Maue-type form), "//trim(tag))
        ev = eigenvalues(ops%v)
        call record_condition(ev(1) > 0.0_dp, "V positive definite, "//trim(tag))
        ev = eigenvalues(0.5_dp*(ops%w + transpose(ops%w)))
        if (n == 0) then
            call record_condition(ev(1) > -1.0e-7_dp*ev(size(ev)), &
                "W positive semidefinite, "//trim(tag))
        else
            call record_condition(ev(1) > 0.0_dp, &
                "W positive definite, "//trim(tag))
        end if
        call tm_bem_steklov_symmetric(ops, s, status)
        ev = eigenvalues(s)
        call record_condition(status == 0 .and. ev(1) > 0.0_dp, &
            "S_n symmetric positive definite, "//trim(tag))
    end subroutine test_structure

    subroutine test_cauchy_data(n)
        !! Calderon identities on the exact Cauchy data (u, t) of an exterior
        !! field: V t = (K - 1/2) u and W u = -(1/2 + K') t, as relative
        !! residuals of the Galerkin vectors on three meshes.
        integer, intent(in) :: n
        integer, parameter :: levels(3) = [12, 24, 48]
        type(gs_bem_curve_t) :: curve
        type(gs_bem_operators_t) :: ops
        real(dp), allocatable :: u(:), t(:), r1(:), r2(:), a1(:), a2(:)
        real(dp) :: e1(3), e2(3), src(2)
        integer :: i, status
        character(len=16) :: tag

        write (tag, '(a,i0)') "n=", n
        src = source_point(2)
        do i = 1, 3
            call make_curve(2, levels(i), curve)
            call tm_bem_assemble(curve, n, 1, ops, status)
            call nodal_trace(n, curve, src, u)
            call nodal_flux(n, curve, 1, src, t)
            a1 = matmul(ops%v, t)
            r1 = a1 - matmul(ops%k - 0.5_dp*ops%m, u)
            a2 = matmul(ops%w, u)
            r2 = a2 + matmul(0.5_dp*transpose(ops%m) + ops%kt, t)
            e1(i) = norm2(r1)/norm2(a1)
            e2(i) = norm2(r2)/norm2(a2)
        end do
        print '(a,a,3es11.3,a,f6.2,a,3es11.3,a,f6.2)', "  Cauchy ", tag, &
            e1, " rate", log(e1(2)/e1(3))/log(2.0_dp), "  hyper", e2, &
            " rate", log(e2(2)/e2(3))/log(2.0_dp)
        ! The flux is interpolated linearly, so the residuals converge with
        ! order 2.
        call record_condition(e1(3) < e1(2) .and. e1(2) < e1(1) .and. &
            log(e1(2)/e1(3))/log(2.0_dp) > 1.8_dp .and. e1(3) < 1.0e-2_dp, &
            "V t = (K - 1/2) u on Cauchy data, "//trim(tag))
        call record_condition(e2(3) < e2(2) .and. e2(2) < e2(1) .and. &
            log(e2(2)/e2(3))/log(2.0_dp) > 1.8_dp .and. e2(3) < 5.0e-2_dp, &
            "W u = -(1/2 + K') t on Cauchy data, "//trim(tag))
    end subroutine test_cauchy_data

    subroutine test_calderon(shape, n)
        !! <W u, v> against <M^T V^{-1}(1/2 M - K) P (1/2 + K) u, v>,
        !! P the L2 projection onto the trace space (V W = 1/4 - K^2).
        integer, intent(in) :: shape, n
        integer, parameter :: levels(3) = [12, 24, 48]
        type(gs_bem_curve_t) :: curve
        type(gs_bem_operators_t) :: ops
        real(dp), allocatable :: u(:), v(:), g(:), c(:, :), mx(:, :)
        real(dp) :: err(3), lhs, rhs, src(2), srcv(2)
        integer :: i, status, info, ipiv(512)
        character(len=16) :: tag

        write (tag, '(a,i0)') "n=", n
        src = source_point(shape)
        srcv = src + [-0.15_dp, 0.1_dp]
        do i = 1, 3
            call make_curve(shape, levels(i), curve)
            call tm_bem_assemble(curve, n, 2, ops, status)
            call nodal_trace(n, curve, src, u)
            call nodal_trace(n, curve, srcv, v)
            mx = ops%mxx
            g = matmul(0.5_dp*ops%mxx + ops%kxx, u)
            call solve_dense(mx, g, info, ipiv)
            c = reshape(matmul(0.5_dp*ops%m - ops%k, g), [ops%ny, 1])
            call solve_spd(ops%v, c)
            rhs = dot_product(matmul(transpose(ops%m), c(:, 1)), v)
            lhs = dot_product(matmul(ops%w, u), v)
            err(i) = abs(lhs - rhs)/abs(lhs)
        end do
        print '(a,a,3es11.3,a,f6.2)', "  Calderon ", tag, err, " rate", &
            log(err(2)/err(3))/log(2.0_dp)
        call record_condition(err(2) < 0.1_dp*err(1), &
            "Galerkin Calderon V W = 1/4 - K^2 converges, "//trim(tag))
        call record_condition(err(3) < 1.0e-5_dp, &
            "Calderon identity small on the fine mesh, "//trim(tag))
    end subroutine test_calderon

    subroutine solve_dense(a, b, info, ipiv)
        real(dp), intent(inout) :: a(:, :), b(:)
        integer, intent(out) :: info, ipiv(:)
        interface
            subroutine dgesv(n, nrhs, a, lda, ipiv, b, ldb, info)
                import :: dp
                integer, intent(in) :: n, nrhs, lda, ldb
                real(dp), intent(inout) :: a(lda, *), b(ldb, *)
                integer, intent(out) :: ipiv(*), info
            end subroutine dgesv
        end interface

        call dgesv(size(a, 1), 1, a, size(a, 1), ipiv, b, size(a, 1), info)
    end subroutine solve_dense

    subroutine solve_spd(a, b)
        real(dp), intent(in) :: a(:, :)
        real(dp), intent(inout) :: b(:, :)
        real(dp), allocatable :: f(:, :)
        integer :: info
        interface
            subroutine dposv(uplo, n, nrhs, a, lda, b, ldb, info)
                import :: dp
                character(len=1), intent(in) :: uplo
                integer, intent(in) :: n, nrhs, lda, ldb
                real(dp), intent(inout) :: a(lda, *), b(ldb, *)
                integer, intent(out) :: info
            end subroutine dposv
        end interface

        f = a
        call dposv('L', size(a, 1), size(b, 2), f, size(a, 1), b, size(a, 1), &
            info)
    end subroutine solve_spd

    subroutine neumann_l2_error(n, curve, ops, src, tcoef, err)
        !! L2(Gamma) error of the computed flux against the exact one.
        integer, intent(in) :: n
        type(gs_bem_curve_t), intent(in) :: curve
        type(gs_bem_operators_t), intent(in) :: ops
        real(dp), intent(in) :: src(2), tcoef(:)
        real(dp), intent(out) :: err
        integer, parameter :: nq = 8
        real(dp) :: gx(nq), gw(nq), x(2), tau(2), sp, u, t, th
        integer :: p, k

        call gauss_unit(nq, gx, gw)
        err = 0.0_dp
        do p = 1, curve%n_panels
            do k = 1, nq
                call gs_bem_panel_point(curve, p, gx(k), x, tau, sp)
                call field_data(n, src, x, tau, u, t)
                th = gs_bem_flux_value(ops%flux_degree, p, gx(k), tcoef)
                err = err + gw(k)*sp*(th - t)**2
            end do
        end do
        err = sqrt(err)
    end subroutine neumann_l2_error

    subroutine exact_energy(n, curve, src, energy)
        !! Exterior energy of the field = -int t u ds.
        integer, intent(in) :: n
        type(gs_bem_curve_t), intent(in) :: curve
        real(dp), intent(in) :: src(2)
        real(dp), intent(out) :: energy
        integer, parameter :: nq = 10
        real(dp) :: gx(nq), gw(nq), x(2), tau(2), sp, u, t
        integer :: p, k

        call gauss_unit(nq, gx, gw)
        energy = 0.0_dp
        do p = 1, curve%n_panels
            do k = 1, nq
                call gs_bem_panel_point(curve, p, gx(k), x, tau, sp)
                call field_data(n, src, x, tau, u, t)
                energy = energy - gw(k)*sp*t*u
            end do
        end do
    end subroutine exact_energy

    subroutine test_exterior_dirichlet(shape, n)
        integer, intent(in) :: shape, n
        integer, parameter :: levels(3) = [12, 24, 48]
        type(gs_bem_curve_t) :: curve
        type(gs_bem_operators_t) :: ops
        real(dp), allocatable :: u(:), t(:), s(:, :)
        real(dp) :: err(3), eerr(3), energy, src(2), rate_t, rate_e
        integer :: i, status
        character(len=40) :: tag

        write (tag, '(a,i0,a,i0)') "shape ", shape, ", n=", n
        src = source_point(shape)
        do i = 1, 3
            call make_curve(shape, levels(i), curve)
            call tm_bem_assemble(curve, n, 1, ops, status)
            call nodal_trace(n, curve, src, u)
            call tm_bem_exterior_neumann(ops, u, t, status)
            call neumann_l2_error(n, curve, ops, src, t, err(i))
            call tm_bem_steklov_symmetric(ops, s, status)
            call exact_energy(n, curve, src, energy)
            eerr(i) = abs(dot_product(u, matmul(s, u)) - energy)/energy
        end do
        rate_t = log(err(2)/err(3))/log(2.0_dp)
        rate_e = log(eerr(2)/eerr(3))/log(2.0_dp)
        print '(a,a,3es11.3,a,f6.2,a,3es11.3,a,f6.2)', "  Neumann ", trim(tag), &
            err, " rate", rate_t, "  energy", eerr, " rate", rate_e
        call record_condition(err(3) < err(2) .and. err(2) < err(1), &
            "Neumann data converge monotonically, "//trim(tag))
        call record_condition(rate_t > 1.7_dp, &
            "Neumann L2 order 2 (flux degree 1), "//trim(tag))
        call record_condition(eerr(3) < eerr(2) .and. rate_e > 2.0_dp .and. &
            eerr(3) < 1.0e-4_dp, &
            "exterior energy <S u, u> converges, "//trim(tag))
    end subroutine test_exterior_dirichlet

    subroutine test_matches_grad_shafranov()
        !! The GS operators act on psi = R A_phi with weight 1/R.  For the
        !! field A = G_1(., y0) of a ring at y0, psi = R0 R A solves the GS
        !! equation (G_GS = R R0 G_1), and with Gamma-integrals of exterior
        !! energies E_GS = E_1 - int_Gamma A^2 n_R ds (the term from
        !! int d_R(A^2)).  Both sides come from independently assembled
        !! operators (different kernels, different hypersingular forms).
        integer, parameter :: levels(3) = [12, 24, 48]
        type(gs_bem_curve_t) :: curve
        type(gs_bem_operators_t) :: ops_t, ops_g
        real(dp), allocatable :: u(:), psi(:), s_t(:, :), s_g(:, :)
        real(dp) :: diff(3), e_t, e_g, bdry, src(2)
        integer :: i, status

        src = source_point(2)
        do i = 1, 3
            call make_curve(2, levels(i), curve)
            call tm_bem_assemble(curve, 1, 1, ops_t, status)
            call gs_bem_assemble(curve, 1, ops_g, status)
            call nodal_trace(1, curve, src, u)
            call psi_trace(curve, src, psi)
            call tm_bem_steklov_symmetric(ops_t, s_t, status)
            call gs_bem_steklov_symmetric(ops_g, s_g, status)
            e_t = dot_product(u, matmul(s_t, u))
            e_g = dot_product(psi, matmul(s_g, psi))/src(1)**2
            call boundary_term(curve, src, bdry)
            diff(i) = abs(e_g - (e_t - bdry))/abs(e_t)
        end do
        print '(a,3es11.3,a,f6.2)', "  n=1 vs Grad-Shafranov energy ", diff, &
            " rate", log(diff(2)/diff(3))/log(2.0_dp)
        call record_condition(diff(3) < diff(2) .and. diff(2) < diff(1) .and. &
            diff(3) < 1.0e-5_dp, &
            "n = 1 energy equals the Grad-Shafranov energy via psi = R A")
    end subroutine test_matches_grad_shafranov

    subroutine psi_trace(curve, src, psi)
        type(gs_bem_curve_t), intent(in) :: curve
        real(dp), intent(in) :: src(2)
        real(dp), allocatable, intent(out) :: psi(:)
        real(dp) :: x(2), tau(2), sp, grad(4)
        integer :: p, k

        allocate(psi(2*curve%n_panels))
        do p = 1, curve%n_panels
            do k = 1, 3
                call gs_bem_panel_point(curve, p, 0.5_dp*(k - 1), x, tau, sp)
                call gs_green_gradient(x(1), x(2), src(1), src(2), &
                    psi(gs_bem_trace_node_index(curve, p, k)), grad)
            end do
        end do
    end subroutine psi_trace

    subroutine boundary_term(curve, src, val)
        !! int_Gamma A^2 n_R ds with A = G_1(., src), n_R = tau_Z.
        type(gs_bem_curve_t), intent(in) :: curve
        real(dp), intent(in) :: src(2)
        real(dp), intent(out) :: val
        integer, parameter :: nq = 10
        real(dp) :: gx(nq), gw(nq), x(2), tau(2), sp, a
        integer :: p, k

        call gauss_unit(nq, gx, gw)
        val = 0.0_dp
        do p = 1, curve%n_panels
            do k = 1, nq
                call gs_bem_panel_point(curve, p, gx(k), x, tau, sp)
                a = tm_green(1, x(1), x(2), src(1), src(2))
                val = val + gw(k)*sp*a*a*tau(2)
            end do
        end do
    end subroutine boundary_term

    subroutine test_jn_and_ntd(n)
        integer, intent(in) :: n
        integer, parameter :: levels(3) = [12, 24, 48]
        type(gs_bem_curve_t) :: curve
        type(gs_bem_operators_t) :: ops
        real(dp), allocatable :: u(:), s(:, :), sj(:, :), nd(:, :), x(:)
        real(dp) :: diff(3), esym, ejn, res, dual
        integer :: i, status

        do i = 1, 3
            call make_curve(2, levels(i), curve)
            call tm_bem_assemble(curve, n, 1, ops, status)
            call nodal_trace(n, curve, d_source, u)
            call tm_bem_steklov_symmetric(ops, s, status)
            call tm_bem_steklov_jn(ops, sj, status)
            esym = dot_product(u, matmul(s, u))
            ejn = dot_product(u, matmul(sj, u))
            diff(i) = abs(esym - ejn)/esym
        end do
        print '(a,3es11.3)', "  JN vs symmetric energy difference ", diff
        call record_condition(diff(3) < diff(2) .and. diff(2) < diff(1), &
            "Johnson-Nedelec energy converges to the symmetric form")

        call tm_bem_neumann_to_dirichlet(s, nd, status)
        call record_condition(status == 0, "NtD map from Cholesky inverse")
        x = matmul(nd, matmul(s, u))
        res = maxval(abs(x - u))/maxval(abs(u))
        call record_condition(res < 1.0e-8_dp, "N_h S_h u = u")
        dual = gs_bem_dual_norm_squared(nd, matmul(s, u))
        call record_condition(abs(dual - esym)/esym < 1.0e-8_dp, &
            "dual norm of S u equals the energy <S u, u>")
    end subroutine test_jn_and_ntd

end program test_toroidal_mode_boundary_operators_2d_slow
