program test_equilibrated_hcurl_meridian_n
    !! Patchwise equilibrated bound for the gauge-fixed per-n meridian
    !! curl-curl problem (fortfem_equilibrated_hcurl_meridian_n).
    !!   1. Exact recovery: A_h a degree-2 polynomial field in the discrete
    !!      space on a straight mesh with continuous rot A_h, J = curl_n^*
    !!      curl_n A_h in closed form: the bound vanishes to rounding.
    !!   2. Manufactured A* (A_phi = 0, zero tangential trace) on a curved
    !!      annulus, J = curl curl A* by finite differences of the analytic
    !!      field, Galerkin A_h in V1 with a = 0, for several n: the bound is
    !!      >= the exact energy error ||curl_n (A* - A_h)|| on every mesh, the
    !!      efficiency index stays below a measured cap and the bound
    !!      converges at the energy rate.
    !!   3. A non-Galerkin A_h (L2 projection of A*): bound >= error as well.
    !$ use omp_lib, only: omp_get_max_threads, omp_set_num_threads
    use check, only: check_condition, check_summary
    use fortfem_equilibrated_hcurl_meridian_n, only: HCURL_GLOBAL, HCURL_OK, &
        HCURL_NOT_GAUGE_FIXED, HCURL_PATCHWISE, meridian_hcurl_equilibrated_bound, &
        meridian_hcurl_estimate_t
    use fortfem_kinds, only: dp
    use fortfem_meridian_curved_mesh, only: meridian_mesh_t
    use fortfem_meridian_feec_n, only: meridian_feec_build, meridian_feec_t, &
        meridian_error_norm, meridian_load, meridian_op_complex, &
        meridian_real_to_complex, meridian_stiffness
    use fortsparse, only: csc_matvec, csc_z_t
    use meridian_test_fields, only: ann_c, ann_r1, ann_r2, fun_a0star, &
        fun_b0star, fun_j0star, mode_n, project, rect_chi, solve_z, zero_rows
    use meridian_test_meshes, only: make_annulus_mesh, make_bubble_rectangle_mesh
    implicit none

    complex(dp), parameter :: ca1 = (1.0_dp, 0.5_dp), ca2 = (0.8_dp, -0.3_dp)
    real(dp) :: cyl_r0 = 0.0_dp !! > 0: cylinder chart (metric radius R0)
    logical :: all_passed

    all_passed = .true.
    !$ call omp_set_num_threads(min(4, omp_get_max_threads()))
    ann_c = 3.0_dp
    ann_r1 = 0.4_dp
    ann_r2 = 1.0_dp
    call exact_recovery(3, 3, 2, 0.0_dp)
    call exact_recovery(3, 4, 1, 0.0_dp)
    call exact_recovery(3, 3, 2, 1.7_dp)
    rect_chi = .true.
    call run_manufactured(1, 2, [1, 3], 5, .true., HCURL_GLOBAL, 2.6_dp, .true.)
    call run_manufactured(1, 3, [2], 5, .true., HCURL_GLOBAL, 1.05_dp, .true.)
    call run_manufactured(2, 3, [1, 3], 5, .true., HCURL_GLOBAL, 1.8_dp, .true.)
    call run_manufactured(2, 4, [2], 4, .true., HCURL_GLOBAL, 1.05_dp, .true.)
    call run_manufactured(3, 4, [2], 3, .true., HCURL_GLOBAL, 1.3_dp, .true.)
    call run_manufactured(2, 3, [3], 4, .false., HCURL_GLOBAL, 0.0_dp, .false.)
    call run_manufactured(2, 3, [1, 3, 10], 4, .true., HCURL_PATCHWISE, 0.0_dp, &
        .false.)
    call check_summary("Meridian per-n H(curl) equilibrated bound")
    if (.not. all_passed) error stop 1

contains

    subroutine record(c, d)
        logical, intent(in) :: c
        character(len=*), intent(in) :: d
        all_passed = all_passed .and. c
        call check_condition(c, d)
    end subroutine record

    ! ---- exact recovery -------------------------------------------------

    subroutine poly_a(x, v)
        !! Degree-2 poloidal field with linear rot A (physical components).
        real(dp), intent(in) :: x(2)
        complex(dp), intent(out) :: v(3)
        real(dp) :: r, z

        r = x(1)
        z = x(2)
        v(1) = ca1*(1.0_dp + 0.5_dp*z + 0.3_dp*r + 0.2_dp*r*z - 0.4_dp*z*z)
        v(2) = (0.0_dp, 0.0_dp)
        v(3) = ca2*(-0.7_dp + 0.4_dp*r + 0.6_dp*z + 0.3_dp*r*r - 0.5_dp*r*z)
    end subroutine poly_a

    subroutine poly_j(x, v)
        !! mu0 J_pol = (1/R) (-rot90 grad(R w) + n^2 A / R), w = rot A, the
        !! strong curl curl of poly_a for gauge a = 0 (poloidal part).
        real(dp), intent(in) :: x(2)
        complex(dp), intent(out) :: v(3)
        complex(dp) :: w, wr, wz, a(3)
        real(dp) :: r, z, nn

        r = x(1)
        z = x(2)
        nn = real(mode_n, dp)
        call poly_a(x, a)
        w = ca2*(0.4_dp + 0.6_dp*r - 0.5_dp*z) &
            - ca1*(0.5_dp + 0.2_dp*r - 0.8_dp*z)
        ! w = rot A = d_R A_Z - d_Z A_R and its constant gradient
        wr = ca2*(0.6_dp) - ca1*(0.2_dp)
        wz = ca2*(-0.5_dp) + ca1*(0.8_dp)
        v(2) = (0.0_dp, 0.0_dp)
        if (cyl_r0 > 0.0_dp) then
            ! constant metric radius: J = -rot90 grad w + n^2 A / R0^2
            v(1) = wz + nn*nn*a(1)/cyl_r0**2
            v(3) = -wr + nn*nn*a(3)/cyl_r0**2
        else
            v(1) = (r*wz + nn*nn*a(1)/r)/r
            v(3) = (-(w + r*wr) + nn*nn*a(3)/r)/r
        end if
    end subroutine poly_j

    subroutine exact_recovery(k, m, n, r0)
        integer, intent(in) :: k, m, n
        real(dp), intent(in) :: r0
        type(meridian_mesh_t) :: mm
        type(meridian_feec_t) :: fs, fsd
        type(meridian_hcurl_estimate_t) :: est
        complex(dp), allocatable :: u(:)
        character(len=80) :: tag
        integer :: status, method
        logical :: ok

        mode_n = n
        call make_annulus_mesh(2, 10, ann_r1, ann_r2, ann_c, 1, mm)
        cyl_r0 = r0
        if (r0 > 0.0_dp) then
            call meridian_feec_build(mm, k, fs, status, quad_degree=2*m + 8, &
                cylinder_r0=r0)
            if (status /= 0) error stop "build failed"
            call meridian_feec_build(mm, m, fsd, status, quad_degree=2*m + 8, &
                cylinder_r0=r0)
        else
            call meridian_feec_build(mm, k, fs, status, quad_degree=2*m + 8)
            if (status /= 0) error stop "build failed"
            call meridian_feec_build(mm, m, fsd, status, quad_degree=2*m + 8)
        end if
        if (status /= 0) error stop "build failed"
        call project(fs, 1, poly_a, u, ok)
        if (.not. ok) error stop "projection failed"
        do method = HCURL_GLOBAL, HCURL_PATCHWISE
            call meridian_hcurl_equilibrated_bound(fs, fsd, n, u, poly_j, est, &
                status, method)
            write (tag, '(a,i0,a,i0,a,i0,a,i0,a,f4.1)') "exact recovery k=", k, &
                " m=", m, " n=", n, " method=", method, " R0=", r0
            print '(a,a,2es11.3)', "  ", trim(tag), est%bound, est%field_norm
            call record(status == HCURL_OK, "status ok, "//trim(tag))
            call record(est%bound < 1.0e-9_dp*est%field_norm, &
                "bound vanishes for exactly equilibrated A_h, "//trim(tag))
            call record(est%field_norm > 1.0e-3_dp, "nontrivial field, "//trim(tag))
        end do
        ! a field with a nonzero a component is rejected
        u(fs%nn + 1) = (1.0_dp, 0.0_dp)
        call meridian_hcurl_equilibrated_bound(fs, fsd, n, u, poly_j, est, status)
        call record(status == HCURL_NOT_GAUGE_FIXED, "gauge a = 0 enforced")
    end subroutine exact_recovery

    ! ---- manufactured problem -------------------------------------------

    subroutine galerkin(fs, n, u)
        type(meridian_feec_t), intent(in) :: fs
        integer, intent(in) :: n
        complex(dp), allocatable, intent(out) :: u(:)
        type(csc_z_t) :: kz
        complex(dp), allocatable :: rhs(:)
        logical, allocatable :: mask(:)
        logical :: ok
        integer :: st

        mode_n = n
        call meridian_stiffness(fs, 1, n, kz, st)
        call meridian_load(fs, 1, fun_j0star, rhs, st)
        allocate (mask(fs%sp(1)%ndof))
        mask = fs%sp(1)%bnd
        mask(fs%nn + 1:) = .true.
        where (mask) rhs = (0.0_dp, 0.0_dp)
        call zero_rows(kz, mask)
        allocate (u(size(rhs)))
        call solve_z(kz, rhs, u, ok)
        if (.not. ok) error stop "solve failed"
    end subroutine galerkin

    subroutine energy_error(fs, n, u, err)
        type(meridian_feec_t), intent(in) :: fs
        integer, intent(in) :: n
        complex(dp), intent(in) :: u(:)
        real(dp), intent(out) :: err
        type(csc_z_t) :: op
        complex(dp), allocatable :: w(:)
        integer :: st

        call meridian_op_complex(fs%curl, n, op, st)
        w = csc_matvec(op, u)
        call meridian_error_norm(fs, 2, w, fun_b0star, err, st)
    end subroutine energy_error

    subroutine run_manufactured(k, m, ns, nlev, galerkin_ah, method, cap, strict)
        !! strict: assert the efficiency cap and the energy rate (otherwise
        !! only the guaranteed bound >= error).
        integer, intent(in) :: k, m, ns(:), nlev, method
        logical, intent(in) :: galerkin_ah, strict
        real(dp), intent(in) :: cap
        type(meridian_mesh_t) :: mm
        type(meridian_feec_t) :: fs, fsd
        type(meridian_hcurl_estimate_t) :: est
        complex(dp), allocatable :: u(:)
        real(dp) :: err(8), bnd(8), eff(8), rate
        integer :: i, lev, status, qd
        logical :: ok, bound_ok
        character(len=90) :: tag

        qd = 2*m + 6
        do i = 1, size(ns)
            bound_ok = .true.
            do lev = 1, nlev
                call make_bubble_rectangle_mesh(4*2**(lev - 1), k, mm)
                call meridian_feec_build(mm, k, fs, status, quad_degree=qd)
                if (status /= 0) error stop "build failed"
                call meridian_feec_build(mm, m, fsd, status, quad_degree=qd)
                if (status /= 0) error stop "build failed"
                mode_n = ns(i)
                if (galerkin_ah) then
                    call galerkin(fs, ns(i), u)
                else
                    call project(fs, 1, fun_a0star, u, ok)
                    if (.not. ok) error stop "projection failed"
                    u(fs%nn + 1:) = (0.0_dp, 0.0_dp)
                    where (fs%sp(1)%bnd) u = (0.0_dp, 0.0_dp)
                end if
                call energy_error(fs, ns(i), u, err(lev))
                call meridian_hcurl_equilibrated_bound(fs, fsd, ns(i), u, &
                    fun_j0star, est, status, method)
                if (status /= HCURL_OK) error stop "estimator failed"
                bnd(lev) = est%bound
                eff(lev) = bnd(lev)/err(lev)
                bound_ok = bound_ok .and. bnd(lev) >= err(lev)
            end do
            rate = log(bnd(nlev - 1)/bnd(nlev))/log(2.0_dp)
            write (tag, '(a,i0,a,i0,a,i0,a,l1,a,i0)') "k=", k, " m=", m, &
                " n=", ns(i), " galerkin=", galerkin_ah, " method=", method
            print '(a,a,a,f6.2)', "  ", trim(tag), "  bound rate", rate
            print '(a,5f8.3)', "    efficiency index", eff(1:nlev)
            print '(a,5es10.2)', "    exact energy error", err(1:nlev)
            call record(bound_ok, "bound >= exact energy error on every mesh, " &
                //trim(tag))
            if (.not. strict) cycle
            call record(all(eff(1:nlev) < cap), &
                "efficiency index below cap, "//trim(tag))
            call record(rate > real(k, dp) - 0.35_dp, &
                "bound converges at the energy rate, "//trim(tag))
        end do
    end subroutine run_manufactured

end program test_equilibrated_hcurl_meridian_n
