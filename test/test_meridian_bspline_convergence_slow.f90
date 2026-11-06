program test_meridian_bspline_convergence_slow
    !! Convergence of the per-mode B-spline FEM-BEM solve.
    !!  1. Manufactured compactly supported B (n = 0 with toroidal flux,
    !!     n = 2): L2 order of B at least p - 1/4, and the field outside the
    !!     support converges to zero.
    !!  2. Circular filament, n = 0: the potential converges to the elliptic
    !!     integral potential.
    !!  3. Spline (degree 2) and Nedelec (order 2) solutions of the same
    !!     per-mode problem, a compact current whose field reaches the
    !!     boundary, agree under refinement: field at interior points, the
    !!     boundary flux moments and the energies.
    !!  4. The exterior energy of the patch [2, 3] x [-0.5, 0.5] equals the
    !!     field energy between it and a larger patch plus the exterior
    !!     energy of the larger patch.
    use check, only: check_condition, check_summary
    use fortfem_kinds, only: dp
    use fortfem_meridian_bspline_n, only: meridian_bspline_t, &
        meridian_bspline_build, meridian_bspline_exterior_build, &
        meridian_bspline_solve, meridian_bspline_load_value, &
        meridian_bspline_eval, meridian_bspline_field_energy
    use fortfem_meridian_curved_mesh, only: meridian_mesh_t
    use fortfem_meridian_feec_n, only: meridian_feec_t, meridian_feec_build, &
        meridian_stiffness, meridian_op_complex, meridian_zcsc_matmul, &
        meridian_real_to_complex, meridian_load, meridian_point_basis_t, &
        meridian_point_basis_init, meridian_point_basis
    use fortfem_meridian_magnetic_coupling, only: meridian_boundary_t, &
        meridian_boundary_build, meridian_exterior_t, meridian_exterior_build, &
        meridian_exterior_solve, meridian_exterior_energy
    use fortsparse, only: csc_z_t
    use meridian_bspline_fields, only: bump_fields, manufactured_solve, &
        loop_solve, patch_quadrature
    use meridian_test_meshes, only: make_bubble_rectangle_mesh
    implicit none

    logical :: all_passed
    integer :: mode

    all_passed = .true.
    call test_manufactured_order()
    call test_filament()
    call test_nedelec_agreement()
    call test_boundary_invariance()
    call check_summary("Meridian per-n B-spline convergence")
    if (.not. all_passed) error stop 1

contains

    subroutine record(c, d)
        logical, intent(in) :: c
        character(len=*), intent(in) :: d
        all_passed = all_passed .and. c
        call check_condition(c, d)
    end subroutine record

    subroutine test_manufactured_order()
        real(dp) :: eb(3), eo(3), ext, pmax
        integer :: p, n, l, status

        do p = 2, 3
            do n = 0, 2, 2
                do l = 1, 3
                    call manufactured_solve(4*2**l, p, n, .false., eb(l), &
                        eo(l), ext, pmax, status)
                    call record(status == 0, "manufactured solve")
                end do
                print '(a,i0,a,i0,a,3es10.2,a,3es10.2)', "  p = ", p, &
                    ", n = ", n, ": |B_h - B| ", eb, "; outside ", eo
                print '(a,2f6.2)', "    orders", log(eb(1:2)/eb(2:3))/log(2.0_dp)
                call record(log(eb(2)/eb(3))/log(2.0_dp) > p - 0.25_dp, &
                    "order of B")
                call record(eo(3) < 0.1_dp*eo(2) .and. eo(3) < 1.0e-3_dp*eb(1), &
                    "zero outside the support")
            end do
        end do
    end subroutine test_manufactured_order

    subroutine test_filament()
        real(dp) :: err(3)
        integer :: l, status

        do l = 1, 3
            call loop_solve(4*2**l, 2, 2.5_dp, err(l), status)
            call record(status == 0, "filament solve")
        end do
        print '(a,3es10.2)', "  filament, p = 2: max error of A_phi ", err
        print '(a,2f6.2)', "    orders", log(err(1:2)/err(2:3))/log(2.0_dp)
        call record(err(3) < 0.1_dp*err(2) .and. err(3) < 2.0e-4_dp, &
            "filament potential converges")
    end subroutine test_filament

    subroutine bump_current(x, v)
        !! Current density mu0 J = curl_n A of the manufactured potential.
        real(dp), intent(in) :: x(2)
        complex(dp), intent(out) :: v(3)
        complex(dp) :: a(3), c(3)

        call bump_fields(x, mode, a, v, c)
    end subroutine bump_current

    subroutine spline_problem(r0, r1, z0, z1, nr, nz, p, n, sp, ex, u, wint, &
            wext, status)
        real(dp), intent(in) :: r0, r1, z0, z1
        integer, intent(in) :: nr, nz, p, n
        type(meridian_bspline_t), intent(out) :: sp
        type(meridian_exterior_t), intent(out) :: ex
        complex(dp), allocatable, intent(out) :: u(:)
        real(dp), intent(out) :: wint, wext
        integer, intent(out) :: status

        real(dp), allocatable :: x(:, :), w(:)
        complex(dp), allocatable :: f(:), g(:, :), pm(:)
        integer :: k

        call meridian_bspline_build(r0, r1, z0, z1, nr, nz, p, sp, status)
        if (status /= 0) return
        call meridian_bspline_exterior_build(sp, n, ex, status)
        if (status /= 0) return
        call patch_quadrature(sp, p + 4, x, w)
        allocate (g(3, size(w)), f(sp%n1))
        mode = n
        do k = 1, size(w)
            call bump_current(x(:, k), g(:, k))
        end do
        f = (0.0_dp, 0.0_dp)
        call meridian_bspline_load_value(sp, n, x, w, g, f, status)
        if (status /= 0) return
        call meridian_bspline_solve(sp, ex, f, u, pm, status)
        if (status /= 0) return
        call meridian_bspline_field_energy(sp, n, u, wint, status)
        wext = meridian_exterior_energy(ex, u)
    end subroutine spline_problem

    subroutine test_nedelec_agreement()
        type(meridian_bspline_t) :: sp
        type(meridian_exterior_t) :: exs, exn
        type(meridian_mesh_t) :: mm
        type(meridian_feec_t) :: fs
        type(meridian_boundary_t) :: bd
        type(csc_z_t) :: k, g, gr, m1
        complex(dp), allocatable :: us(:), un(:), pn(:), f(:), as(:, :), &
            bs(:, :), bn(:, :), ts(:), tn(:)
        real(dp) :: xs(2, 9), ws, we, wn, wne, db(3), dt(3), dw(3)
        integer :: l, nx, n, status, i, j

        do j = 1, 3
            do i = 1, 3
                xs(:, i + 3*(j - 1)) = [2.0_dp + 0.25_dp*i, &
                    -0.5_dp + 0.25_dp*j]
            end do
        end do
        do n = 0, 2, 2
            do l = 1, 3
                nx = 4*2**l
                call spline_problem(2.0_dp, 3.0_dp, -0.5_dp, 0.5_dp, nx, nx, 2, &
                    n, sp, exs, us, ws, we, status)
                call record(status == 0, "spline solve")
                allocate (as(3, 9), bs(3, 9), bn(3, 9))
                call meridian_bspline_eval(sp, n, us, xs, as, bs, status)

                call make_bubble_rectangle_mesh(nx, 1, mm)
                call meridian_feec_build(mm, 2, fs, status)
                call meridian_boundary_build(fs, bd, status)
                call meridian_exterior_build(fs, bd, n, exn, status)
                call record(status == 0, "Nedelec exterior")
                call meridian_stiffness(fs, 1, n, k, status)
                call meridian_op_complex(fs%grad, n, gr, status)
                call meridian_real_to_complex(fs%mass(1), m1)
                call meridian_zcsc_matmul(m1, gr, g, status)
                mode = n
                call meridian_load(fs, 1, bump_current, f, status)
                call meridian_exterior_solve(k, g, exn, f, un, pn, status)
                call record(status == 0, "Nedelec solve")
                wne = meridian_exterior_energy(exn, un)
                wn = real(dot_product(f, un), dp) - wne
                call nedelec_field(fs, n, un, xs, bn)
                db(l) = maxval(abs(bs - bn))/maxval(abs(bn))
                call trace_moments(exs, us, ts)
                call trace_moments(exn, un, tn)
                dt(l) = maxval(abs(ts - tn))/maxval(abs(tn))
                dw(l) = abs(ws + we - wn - wne)/(wn + wne)
                deallocate (as, bs, bn)
            end do
            print '(a,i0,a,3es10.2)', "  n = ", n, &
                ": spline - Nedelec, B at points ", db
            print '(a,3es10.2)', "    boundary flux moments ", dt
            print '(a,3es10.2)', "    total energy ", dw
            call record(db(3) < 0.4_dp*db(2) .and. db(2) < 0.4_dp*db(1), &
                "spline and Nedelec fields agree under refinement")
            call record(dt(3) < 0.3_dp*dt(2), &
                "spline and Nedelec boundary traces agree under refinement")
            call record(dw(3) < 0.3_dp*dw(2) .and. dw(3) < 1.0e-3_dp, &
                "spline and Nedelec energies agree under refinement")
        end do
    end subroutine test_nedelec_agreement

    subroutine trace_moments(ex, u, t)
        !! T u, the flux (n /= 0) or Dirichlet (n = 0) moments on the curve.
        type(meridian_exterior_t), intent(in) :: ex
        complex(dp), intent(in) :: u(:)
        complex(dp), allocatable, intent(out) :: t(:)
        integer :: j, e

        allocate (t(ex%trace%nrow))
        t = (0.0_dp, 0.0_dp)
        do j = 1, ex%trace%ncol
            do e = ex%trace%col_ptr(j), ex%trace%col_ptr(j + 1) - 1
                t(ex%trace%row_idx(e)) = t(ex%trace%row_idx(e)) &
                    + ex%trace%val(e)*u(j)
            end do
        end do
    end subroutine trace_moments

    subroutine nedelec_field(fs, n, u, x, b)
        !! B = curl_n u of the Nedelec field at points of an affine mesh.
        type(meridian_feec_t), intent(in) :: fs
        integer, intent(in) :: n
        complex(dp), intent(in) :: u(:)
        real(dp), intent(in) :: x(:, :)
        complex(dp), intent(out) :: b(:, :)

        type(meridian_point_basis_t) :: pb
        type(csc_z_t) :: c
        complex(dp), allocatable :: v(:)
        real(dp), allocatable :: phi(:, :)
        real(dp) :: e1(2), e2(2), d(2), det, xi, eta, xx(2), x0(2)
        integer :: k, t, j, st, e

        call meridian_point_basis_init(fs, pb, st)
        call meridian_op_complex(fs%curl, n, c, st)
        allocate (v(c%nrow), phi(3, fs%sp(2)%nloc))
        v = (0.0_dp, 0.0_dp)
        do j = 1, c%ncol
            do e = c%col_ptr(j), c%col_ptr(j + 1) - 1
                v(c%row_idx(e)) = v(c%row_idx(e)) + c%val(e)*u(j)
            end do
        end do
        do k = 1, size(x, 2)
            b(:, k) = (0.0_dp, 0.0_dp)
            do t = 1, fs%mm%topology%n_triangles
                ! vertices of the affine element from its own geometry map
                call meridian_point_basis(fs, pb, 2, t, 0.0_dp, 0.0_dp, x0, &
                    phi, st)
                call meridian_point_basis(fs, pb, 2, t, 1.0_dp, 0.0_dp, xx, &
                    phi, st)
                e1 = xx - x0
                call meridian_point_basis(fs, pb, 2, t, 0.0_dp, 1.0_dp, xx, &
                    phi, st)
                e2 = xx - x0
                d = x(:, k) - x0
                det = e1(1)*e2(2) - e1(2)*e2(1)
                xi = (d(1)*e2(2) - d(2)*e2(1))/det
                eta = (e1(1)*d(2) - e1(2)*d(1))/det
                if (min(xi, eta, 1.0_dp - xi - eta) < -1.0e-12_dp) cycle
                call meridian_point_basis(fs, pb, 2, t, xi, eta, xx, phi, st)
                do j = 1, fs%sp(2)%nloc
                    b(:, k) = b(:, k) + v(fs%sp(2)%dof(j, t)) &
                        *real(fs%sp(2)%sgn(j, t), dp)*phi(:, j)
                end do
                exit
            end do
        end do
    end subroutine nedelec_field

    subroutine test_boundary_invariance()
        type(meridian_bspline_t) :: s1, s2
        type(meridian_exterior_t) :: e1, e2
        complex(dp), allocatable :: u1(:), u2(:), a(:, :), b(:, :)
        real(dp), allocatable :: x(:, :), w(:)
        real(dp) :: wi1, we1, wi2, we2, wring, rel(2)
        integer :: n, l, status, k

        do n = 0, 2, 2
            do l = 1, 2
                call spline_problem(2.0_dp, 3.0_dp, -0.5_dp, 0.5_dp, 8*l, 8*l, &
                    3, n, s1, e1, u1, wi1, we1, status)
                call spline_problem(1.75_dp, 3.25_dp, -0.75_dp, 0.75_dp, 12*l, &
                    12*l, 3, n, s2, e2, u2, wi2, we2, status)
                call record(status == 0, "two patches solved")
                call patch_quadrature(s2, 6, x, w)
                allocate (a(3, size(w)), b(3, size(w)))
                call meridian_bspline_eval(s2, n, u2, x, a, b, status)
                wring = 0.0_dp
                do k = 1, size(w)
                    if (x(1, k) > 2.0_dp .and. x(1, k) < 3.0_dp .and. &
                        abs(x(2, k)) < 0.5_dp) cycle
                    wring = wring + w(k)*sum(abs(b(:, k))**2)
                end do
                rel(l) = abs(we1 - wring - we2)/we1
                deallocate (a, b)
            end do
            print '(a,i0,a,2es10.2,a,es10.2)', "  n = ", n, &
                ": |W_ext(inner) - W(ring) - W_ext(outer)| / W_ext(inner)", &
                rel, ", W_ext / W_int ", we1/wi1
            call record(rel(2) < 1.0e-3_dp .and. rel(2) < rel(1), &
                "exterior energy invariant under moving the boundary")
        end do
    end subroutine test_boundary_invariance

end program test_meridian_bspline_convergence_slow
