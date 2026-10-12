program test_meridian_feec_gs_n0
    !! n = 0 reduces to the axisymmetric Grad-Shafranov operator
    !! L u = -div((1/R) grad u):
    !!  - the n = 0 curl-curl form of A = (0, R A_phi = a) equals the
    !!    directly assembled GS stiffness int (1/R) grad a . grad v,
    !!  - for the filament field u = G(x, y0) with y0 in the hole of a curved
    !!    annulus, the FEM energy u^T K u converges to the annulus energy
    !!    E_ext(inner circle) - E_ext(outer circle) of the existing GS
    !!    boundary-operator module (exact flux quadrature) and agrees with
    !!    the Steklov-Poincare energies u^T S u of the BEM operators.
    !$ use omp_lib, only: omp_get_max_threads, omp_set_num_threads
    use check, only: check_condition, check_summary
    use fortfem_kinds, only: dp
    use fortfem_axisymmetric_gs_kernel, only: gs_green, gs_green_gradient
    use fortfem_axisymmetric_gs_boundary_operators_2d, only: gs_bem_curve_t, &
        gs_bem_operators_t, gs_bem_curve_from_nodes, gs_bem_assemble, &
        gs_bem_steklov_symmetric, gs_bem_panel_point, gs_bem_trace_node_index
    use fortfem_meridian_curved_mesh, only: meridian_mesh_t
    use fortfem_meridian_feec_n, only: meridian_feec_t, meridian_feec_build, &
        meridian_stiffness, meridian_gs_stiffness, meridian_lagrange_nodes
    use fortsparse, only: csc_t, csc_z_t
    use meridian_test_meshes, only: make_annulus_mesh
    implicit none

    real(dp), parameter :: pi = acos(-1.0_dp)
    real(dp), parameter :: ctr(2) = [3.0_dp, 0.0_dp]
    real(dp), parameter :: src(2) = [3.1_dp, 0.1_dp]
    real(dp), parameter :: r1 = 0.4_dp, r2 = 1.0_dp
    logical :: all_passed

    all_passed = .true.
    !$ call omp_set_num_threads(min(4, omp_get_max_threads()))
    call test_block_identity()
    call test_energy_convergence()
    call check_summary("Meridian per-n FEEC: n = 0 Grad-Shafranov agreement")
    if (.not. all_passed) error stop 1

contains

    subroutine record(c, d)
        logical, intent(in) :: c
        character(len=*), intent(in) :: d
        all_passed = all_passed .and. c
        call check_condition(c, d)
    end subroutine record

    function quad_form(a, u) result(e)
        !! u^T A u for a real CSC matrix.
        type(csc_t), intent(in) :: a
        real(dp), intent(in) :: u(:)
        real(dp) :: e
        integer :: j, p

        e = 0.0_dp
        do j = 1, a%ncol
            do p = a%col_ptr(j), a%col_ptr(j + 1) - 1
                e = e + u(a%row_idx(p))*a%val(p)*u(j)
            end do
        end do
    end function quad_form

    subroutine test_block_identity()
        type(meridian_mesh_t) :: mm
        type(meridian_feec_t) :: fs
        type(csc_t) :: kgs
        type(csc_z_t) :: kc
        real(dp), allocatable :: xy(:, :), u(:)
        complex(dp), allocatable :: e(:)
        real(dp) :: egs, ecc
        integer :: k, status, i, j, p
        character(len=40) :: tag

        do k = 1, 3
            call make_annulus_mesh(3, 12, r1, r2, ctr(1), k, mm)
            call meridian_feec_build(mm, k, fs, status)
            call meridian_gs_stiffness(fs, kgs, status)
            call meridian_stiffness(fs, 1, 0, kc, status)
            call meridian_lagrange_nodes(fs, xy)
            allocate (u(fs%n0), e(fs%sp(1)%ndof))
            do i = 1, fs%n0
                u(i) = gs_green(xy(1, i), xy(2, i), src(1), src(2))
            end do
            e = (0.0_dp, 0.0_dp)
            e(fs%nn + 1:) = u
            egs = quad_form(kgs, u)
            ecc = 0.0_dp
            do j = 1, kc%ncol
                do p = kc%col_ptr(j), kc%col_ptr(j + 1) - 1
                    ecc = ecc + real(conjg(e(kc%row_idx(p)))*kc%val(p)*e(j), dp)
                end do
            end do
            write (tag, '(a,i0)') "order ", k
            print '(a,a,2es16.8)', "  GS energy vs n=0 curl-curl ", trim(tag), &
                egs, ecc
            call record(abs(egs - ecc) < 1.0e-12_dp*abs(egs), &
                "n = 0 curl-curl on (0, a) equals GS stiffness, "//trim(tag))
            deallocate (u, e)
        end do
    end subroutine test_block_identity

    subroutine circle_curve(radius, np, curve)
        real(dp), intent(in) :: radius
        integer, intent(in) :: np
        type(gs_bem_curve_t), intent(out) :: curve
        real(dp) :: nodes(2, 2*np)
        integer :: j, status

        do j = 1, 2*np
            nodes(:, j) = ctr + radius*[cos(pi*(j - 1)/np), sin(pi*(j - 1)/np)]
        end do
        call gs_bem_curve_from_nodes(nodes, curve, status)
        if (status /= 0) error stop "curve failed"
    end subroutine circle_curve

    subroutine exterior_energies(radius, e_exact, e_bem)
        !! Exterior energy of u = G(., src) outside the circle: -int t u ds
        !! with exact flux, and u^T S u from the BEM operators.
        real(dp), intent(in) :: radius
        real(dp), intent(out) :: e_exact, e_bem
        integer, parameter :: np = 48, nq = 10
        type(gs_bem_curve_t) :: curve
        type(gs_bem_operators_t) :: ops
        real(dp), allocatable :: u(:), s(:, :)
        real(dp) :: gx(nq), gw(nq), x(2), tau(2), sp, g, t, grad(4), xi
        integer :: p, q, status, kk

        call circle_curve(radius, np, curve)
        call gauss_nodes(nq, gx, gw)
        e_exact = 0.0_dp
        do p = 1, np
            do q = 1, nq
                call gs_bem_panel_point(curve, p, gx(q), x, tau, sp)
                call gs_green_gradient(x(1), x(2), src(1), src(2), g, grad)
                t = (tau(2)*grad(1) - tau(1)*grad(2))/x(1)
                e_exact = e_exact - gw(q)*sp*t*g
            end do
        end do
        call gs_bem_assemble(curve, 1, ops, status)
        allocate (u(ops%nx))
        do p = 1, np
            do kk = 1, 3
                xi = 0.5_dp*(kk - 1)
                call gs_bem_panel_point(curve, p, xi, x, tau, sp)
                u(gs_bem_trace_node_index(curve, p, kk)) = &
                    gs_green(x(1), x(2), src(1), src(2))
            end do
        end do
        call gs_bem_steklov_symmetric(ops, s, status)
        e_bem = dot_product(u, matmul(s, u))
    end subroutine exterior_energies

    subroutine gauss_nodes(n, x, w)
        use fortnum_quadrature, only: gauss_legendre_ab
        integer, intent(in) :: n
        real(dp), intent(out) :: x(n), w(n)
        call gauss_legendre_ab(n, 0.0_dp, 1.0_dp, x, w)
    end subroutine gauss_nodes

    subroutine test_energy_convergence()
        type(meridian_mesh_t) :: mm
        type(meridian_feec_t) :: fs
        type(csc_t) :: kgs
        real(dp), allocatable :: xy(:, :), u(:)
        real(dp) :: ein_x, ein_b, eout_x, eout_b, e_ref, e_bem, err(3)
        real(dp) :: rate, efem
        integer :: k, lev, status, i
        character(len=40) :: tag

        call exterior_energies(r1, ein_x, ein_b)
        call exterior_energies(r2, eout_x, eout_b)
        e_ref = ein_x - eout_x
        e_bem = ein_b - eout_b
        print '(a,2es18.10)', "  annulus energy: exact flux quadrature, BEM S ", &
            e_ref, e_bem
        call record(abs(e_bem - e_ref) < 2.0e-3_dp*abs(e_ref), &
            "BEM Steklov energies reproduce the annulus energy")
        do k = 1, 3
            do lev = 1, 3
                call make_annulus_mesh(2**(lev - 1)*2, 2**(lev - 1)*8, r1, r2, &
                    ctr(1), k, mm)
                call meridian_feec_build(mm, k, fs, status)
                call meridian_gs_stiffness(fs, kgs, status)
                call meridian_lagrange_nodes(fs, xy)
                allocate (u(fs%n0))
                do i = 1, fs%n0
                    u(i) = gs_green(xy(1, i), xy(2, i), src(1), src(2))
                end do
                efem = quad_form(kgs, u)
                err(lev) = abs(efem - e_ref)/abs(e_ref)
                deallocate (u)
            end do
            rate = log(err(2)/err(3))/log(2.0_dp)
            write (tag, '(a,i0)') "order ", k
            print '(a,a,3es11.3,a,f6.2)', "  FEM energy error ", trim(tag), err, &
                " rate", rate
            call record(err(3) < err(2) .and. err(2) < err(1), &
                "FEM GS energy converges monotonically, "//trim(tag))
            call record(rate > min(real(2*k, dp) - 0.8_dp, 4.4_dp), &
                "FEM GS energy converges at order 2k, "//trim(tag))
            if (k == 3) call record(abs(efem - e_bem) < 1.0e-5_dp*e_ref, &
                "FEM GS energy agrees with the BEM operator energy, "//trim(tag))
        end do
    end subroutine test_energy_convergence

end program test_meridian_feec_gs_n0
