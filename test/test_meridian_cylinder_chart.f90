program test_meridian_cylinder_chart
    !! Periodic-cylinder chart of the per-n FEEC complex (cylinder_r0 > 0):
    !! a straight cylinder of axial period 2 pi R0 with the cross-section
    !! (x, y) in the (R, Z) slots.
    !!  - Exact-sequence identities curl_n grad_n = 0 and div_n curl_n = 0 on
    !!    a disc that contains the axis, orders 1-3, geometry orders 1-3,
    !!    n = 0, 1, 3 (the identities are metric free, so they must survive
    !!    negative x1 and the axis).
    !!  - The chart scaling of the masses: M_0 ~ R0, M_3 ~ 1/R0, and the
    !!    volume R0 * pi a^2 of the cylinder.
    !!  - An analytic vacuum field in an annulus a < r < b,
    !!      psi = [I_m(k r) + c K_m(k r)] exp(i m theta) exp(i k z),  k = n / R0
    !!    (r^m + c r^-m for k = 0),
    !!    solved from its Dirichlet data with the per-n stiffness
    !!    grad_n^H M_1 grad_n: the potential converges at order k + 1 and
    !!    B = grad psi = (psi_x, i k psi, psi_y) from grad_n at order k, with
    !!    the axial component i n psi / R0.
    !$ use omp_lib, only: omp_get_max_threads, omp_set_num_threads
    use check, only: check_condition, check_summary
    use fortfem_kinds, only: dp
    use fortfem_mesh_2d, only: mesh_2d_t
    use fortfem_meridian_curved_mesh, only: meridian_mesh_t, &
        meridian_mesh_from_map
    use fortfem_triangle_lagrange_arbitrary_order, only: &
        initialize_triangle_lagrange_basis, triangle_lagrange_basis_t, &
        triangle_lagrange_nodes
    use fortfem_meridian_feec_n, only: meridian_feec_t, meridian_feec_build, &
        meridian_op_complex, meridian_zcsc_matmul, meridian_error_norm, &
        meridian_stiffness, meridian_lagrange_nodes
    use fortnum_special_bessel, only: bessel_in, bessel_kn
    use fortsparse, only: csc_matvec, csc_z_t
    use meridian_test_fields, only: solve_z, zero_rows
    use meridian_test_meshes, only: make_annulus_mesh
    implicit none

    real(dp), parameter :: pi = acos(-1.0_dp)
    real(dp), parameter :: ra = 0.5_dp, rb = 1.0_dp
    ! Mode of the analytic vacuum field (module-level for the callbacks).
    integer :: vm = 1
    real(dp) :: vk = 1.0_dp, vc = 0.0_dp
    logical :: all_passed

    all_passed = .true.
    !$ call omp_set_num_threads(min(4, omp_get_max_threads()))
    call identities_on_disc()
    call chart_scaling()
    call run_vacuum(1, 1, 2, 0.0_dp, 1.7_dp, 3)
    call run_vacuum(2, 1, 2, 0.0_dp, 1.7_dp, 3)
    call run_vacuum(3, 1, 2, 0.0_dp, 1.7_dp, 2)
    call run_vacuum(1, 3, 2, 0.3_dp, 0.8_dp, 3)
    call run_vacuum(2, 0, 1, 0.3_dp, 1.0_dp, 3)
    call check_summary("Meridian per-n FEEC: periodic-cylinder chart")
    if (.not. all_passed) error stop 1

contains

    subroutine record(c, d)
        logical, intent(in) :: c
        character(len=*), intent(in) :: d
        all_passed = all_passed .and. c
        call check_condition(c, d)
    end subroutine record

    ! ------------------------------------------------------------------
    ! Disc mesh containing the axis
    ! ------------------------------------------------------------------

    subroutine disc_map(p, x)
        real(dp), intent(in) :: p(2)
        real(dp), intent(out) :: x(2)

        x = [p(1)*cos(p(2)), p(1)*sin(p(2))]
    end subroutine disc_map

    subroutine make_disc_mesh(nr, nth, radius, geo, mm)
        !! Disc of the given radius: a centre vertex, nr rings, nth sectors.
        integer, intent(in) :: nr, nth, geo
        real(dp), intent(in) :: radius
        type(meridian_mesh_t), intent(out) :: mm
        type(mesh_2d_t) :: topo
        real(dp), allocatable :: pv(:, :, :)
        real(dp) :: dr, dth, x(2)
        integer :: i, j, jn, k, status

        dr = radius/real(nr, dp)
        dth = 2.0_dp*pi/real(nth, dp)
        topo%n_vertices = 1 + nr*nth
        topo%n_triangles = nth + 2*(nr - 1)*nth
        topo%has_triangles = .true.
        allocate (topo%vertices(2, topo%n_vertices))
        allocate (topo%triangles(3, topo%n_triangles))
        allocate (pv(2, 3, topo%n_triangles))
        topo%vertices(:, 1) = 0.0_dp
        do i = 1, nr
            do j = 0, nth - 1
                call disc_map([i*dr, j*dth], x)
                topo%vertices(:, 1 + (i - 1)*nth + j + 1) = x
            end do
        end do
        k = 0
        do j = 0, nth - 1
            jn = mod(j + 1, nth)
            k = k + 1
            topo%triangles(:, k) = [1, 1 + j + 1, 1 + jn + 1]
            pv(:, :, k) = reshape([0.0_dp, j*dth, dr, j*dth, dr, (j + 1)*dth], &
                [2, 3])
        end do
        do i = 1, nr - 1
            do j = 0, nth - 1
                jn = mod(j + 1, nth)
                k = k + 1
                topo%triangles(:, k) = [ring_id(i, j, nth), &
                    ring_id(i + 1, j, nth), ring_id(i + 1, jn, nth)]
                pv(:, :, k) = reshape([i*dr, j*dth, (i + 1)*dr, j*dth, &
                    (i + 1)*dr, (j + 1)*dth], [2, 3])
                k = k + 1
                topo%triangles(:, k) = [ring_id(i, j, nth), &
                    ring_id(i + 1, jn, nth), ring_id(i, jn, nth)]
                pv(:, :, k) = reshape([i*dr, j*dth, (i + 1)*dr, (j + 1)*dth, &
                    i*dr, (j + 1)*dth], [2, 3])
            end do
        end do
        call meridian_mesh_from_map(topo, pv, geo, mm, status, disc_map)
        if (status /= 0) error stop "disc mesh failed"
        call curve_centre_triangles(mm, nth, dr, dth)
    end subroutine make_disc_mesh

    subroutine curve_centre_triangles(mm, nth, dr, dth)
        !! The first nth triangles touch the centre.  A polar parameter map
        !! would turn their centre-to-ring edge into a spiral, so their
        !! geometry nodes are set from the lattice (a, b) of the triangle
        !! (centre, ring j, ring j + 1) as r = (a + b) dr, theta = theta_j +
        !! b / (a + b) dth: radial edges and the circular arc, conforming with
        !! the neighbouring triangles.
        type(meridian_mesh_t), intent(inout) :: mm
        integer, intent(in) :: nth
        real(dp), intent(in) :: dr, dth
        type(triangle_lagrange_basis_t) :: basis
        real(dp), allocatable :: lat(:, :)
        real(dp) :: a, b, s, rho
        integer :: t, i, status

        call initialize_triangle_lagrange_basis(mm%geo_order, basis, status)
        call triangle_lagrange_nodes(basis, lat)
        do t = 1, nth
            do i = 1, size(lat, 2)
                a = lat(1, i)
                b = lat(2, i)
                if (a + b < 1.0e-14_dp) then
                    mm%nodes(:, i, t) = 0.0_dp
                    cycle
                end if
                s = b/(a + b)
                rho = (a + b)*dr
                mm%nodes(:, i, t) = rho*[cos((t - 1)*dth + s*dth), &
                    sin((t - 1)*dth + s*dth)]
            end do
        end do
    end subroutine curve_centre_triangles

    ! ------------------------------------------------------------------
    ! Identities and chart scaling
    ! ------------------------------------------------------------------

    function max_abs(z) result(m)
        type(csc_z_t), intent(in) :: z
        real(dp) :: m
        m = 0.0_dp
        if (z%nnz > 0) m = maxval(abs(z%val))
    end function max_abs

    subroutine identities_on_disc()
        type(meridian_mesh_t) :: mm
        type(meridian_feec_t) :: fs
        type(csc_z_t) :: g, c, d, cg, dc
        integer, parameter :: ns(3) = [0, 1, 3]
        character(len=60) :: tag
        real(dp) :: e1, e2, s1, s2
        integer :: geo, k, i, status

        do geo = 1, 3
            call make_disc_mesh(2, 6, 1.0_dp, geo, mm)
            do k = 1, 3
                call meridian_feec_build(mm, k, fs, status, cylinder_r0=2.5_dp)
                write (tag, '(a,i0,a,i0)') "disc geo ", geo, " order ", k
                call record(status == 0, "build "//trim(tag))
                if (status /= 0) cycle
                e1 = 0.0_dp
                e2 = 0.0_dp
                s1 = 0.0_dp
                s2 = 0.0_dp
                do i = 1, size(ns)
                    call meridian_op_complex(fs%grad, ns(i), g, status)
                    call meridian_op_complex(fs%curl, ns(i), c, status)
                    call meridian_op_complex(fs%div, ns(i), d, status)
                    call meridian_zcsc_matmul(c, g, cg, status)
                    call meridian_zcsc_matmul(d, c, dc, status)
                    e1 = max(e1, max_abs(cg))
                    e2 = max(e2, max_abs(dc))
                    s1 = max(s1, max_abs(g))
                    s2 = max(s2, max_abs(c))
                end do
                call record(e1 < 1.0e-12_dp*max(1.0_dp, s1*s2), &
                    "curl_n grad_n = 0 on the cylinder chart, "//trim(tag))
                call record(e2 < 1.0e-12_dp*max(1.0_dp, s2), &
                    "div_n curl_n = 0 on the cylinder chart, "//trim(tag))
                call record(s1 > 0.5_dp .and. s2 > 0.5_dp, &
                    "operators are nontrivial, "//trim(tag))
            end do
        end do
        ! The default chart still rejects a mesh through the axis x1 <= 0.
        call make_disc_mesh(2, 6, 1.0_dp, 1, mm)
        call meridian_feec_build(mm, 1, fs, status)
        call record(status /= 0, "axisymmetric chart rejects R <= 0 meshes")
        call meridian_feec_build(mm, 1, fs, status, cylinder_r0=-1.0_dp)
        call record(status /= 0, "cylinder_r0 <= 0 is rejected")
    end subroutine identities_on_disc

    pure function ring_id(i, j, nth) result(id)
        integer, intent(in) :: i, j, nth
        integer :: id
        id = 1 + (i - 1)*nth + j + 1
    end function ring_id

    subroutine chart_scaling()
        type(meridian_mesh_t) :: mm
        type(meridian_feec_t) :: f1, f2
        real(dp) :: v1, v2, vol
        integer :: status, geo

        geo = 3
        call make_disc_mesh(3, 12, 1.0_dp, geo, mm)
        call meridian_feec_build(mm, 2, f1, status, cylinder_r0=1.0_dp)
        call record(status == 0, "build R0 = 1")
        call meridian_feec_build(mm, 2, f2, status, cylinder_r0=2.5_dp)
        call record(status == 0, "build R0 = 2.5")
        v1 = sum(f1%mass(0)%val)
        v2 = sum(f2%mass(0)%val)
        vol = pi
        print '(a,2es14.6)', "  mass(0) volume / (R0 pi a^2):", v1/vol, v2/(2.5_dp*vol)
        call record(abs(v1/vol - 1.0_dp) < 2.0e-4_dp, &
            "mass(0) sums to R0 pi a^2 (R0 = 1)")
        call record(abs(v2/(2.5_dp*vol) - 1.0_dp) < 2.0e-4_dp, &
            "mass(0) sums to R0 pi a^2 (R0 = 2.5)")
        call record(abs(v2/v1 - 2.5_dp) < 1.0e-12_dp, "mass(0) scales as R0")
        v1 = sum(f1%mass(3)%val)
        v2 = sum(f2%mass(3)%val)
        call record(abs(v1/v2 - 2.5_dp) < 1.0e-10_dp, "mass(3) scales as 1/R0")
    end subroutine chart_scaling

    ! ------------------------------------------------------------------
    ! Analytic vacuum field in an annulus
    ! ------------------------------------------------------------------

    subroutine radial_profile(r, f, df)
        !! F(r) = I_m(k r) + c K_m(k r) and F'(r).
        real(dp), intent(in) :: r
        real(dp), intent(out) :: f, df
        real(dp) :: x

        x = vk*r
        if (vk == 0.0_dp) then
            ! k = 0: the 2D harmonic functions r^m and r^-m (m >= 1).
            f = r**vm + vc*r**(-vm)
            df = real(vm, dp)*(r**(vm - 1) - vc*r**(-vm - 1))
            return
        end if
        f = bessel_in(vm, x) + vc*bessel_kn(vm, x)
        if (vm == 0) then
            df = vk*(bessel_in(1, x) - vc*bessel_kn(1, x))
        else
            df = 0.5_dp*vk*(bessel_in(vm - 1, x) + bessel_in(vm + 1, x) &
                - vc*(bessel_kn(vm - 1, x) + bessel_kn(vm + 1, x)))
        end if
    end subroutine radial_profile

    subroutine psi_value(x, psi)
        real(dp), intent(in) :: x(2)
        complex(dp), intent(out) :: psi
        real(dp) :: r, th, f, df

        r = hypot(x(1), x(2))
        th = atan2(x(2), x(1))
        call radial_profile(r, f, df)
        psi = f*exp(cmplx(0.0_dp, real(vm, dp)*th, dp))
    end subroutine psi_value

    subroutine fun_psi(x, v)
        real(dp), intent(in) :: x(2)
        complex(dp), intent(out) :: v(3)
        v = (0.0_dp, 0.0_dp)
        call psi_value(x, v(1))
    end subroutine fun_psi

    subroutine fun_bvac(x, v)
        !! B = grad psi in the slots (x, phi, y): (psi_x, i k psi, psi_y).
        real(dp), intent(in) :: x(2)
        complex(dp), intent(out) :: v(3)
        real(dp) :: r, th, f, df
        complex(dp) :: e, im

        r = hypot(x(1), x(2))
        th = atan2(x(2), x(1))
        call radial_profile(r, f, df)
        e = exp(cmplx(0.0_dp, real(vm, dp)*th, dp))
        im = cmplx(0.0_dp, real(vm, dp), dp)
        v(1) = e*(df*cos(th) - im*f*sin(th)/r)
        v(3) = e*(df*sin(th) + im*f*cos(th)/r)
        v(2) = cmplx(0.0_dp, vk, dp)*f*e
    end subroutine fun_bvac

    subroutine solve_vacuum(fs, n, e_pot, e_b)
        !! Dirichlet problem grad_n^H M_1 grad_n u = 0, u = psi on the walls.
        type(meridian_feec_t), intent(in) :: fs
        integer, intent(in) :: n
        real(dp), intent(out) :: e_pot, e_b
        type(csc_z_t) :: kz, op
        real(dp), allocatable :: xy(:, :)
        complex(dp), allocatable :: u(:), ub(:), rhs(:), w(:)
        complex(dp) :: psi
        logical :: ok
        integer :: i, st

        call meridian_stiffness(fs, 0, n, kz, st)
        call meridian_lagrange_nodes(fs, xy)
        allocate (ub(fs%n0), u(fs%n0))
        ub = (0.0_dp, 0.0_dp)
        do i = 1, fs%n0
            if (fs%sp(0)%bnd(i)) then
                call psi_value(xy(:, i), psi)
                ub(i) = psi
            end if
        end do
        rhs = -csc_matvec(kz, ub)
        where (fs%sp(0)%bnd) rhs = ub
        call zero_rows(kz, fs%sp(0)%bnd)
        call solve_z(kz, rhs, u, ok)
        if (.not. ok) error stop "vacuum solve failed"
        call meridian_error_norm(fs, 0, u, fun_psi, e_pot, st)
        call meridian_op_complex(fs%grad, n, op, st)
        w = csc_matvec(op, u)
        call meridian_error_norm(fs, 1, w, fun_bvac, e_b, st)
    end subroutine solve_vacuum

    subroutine run_vacuum(m, n, k, cmix, r0, nlev)
        integer, intent(in) :: m, n, k, nlev
        real(dp), intent(in) :: cmix, r0
        type(meridian_mesh_t) :: mm
        type(meridian_feec_t) :: fs
        real(dp) :: ep(3), eb(3), rp, rb_rate
        integer :: lev, nr, nth, status
        character(len=70) :: tag

        vm = m
        vc = cmix
        vk = real(n, dp)/r0
        do lev = 1, nlev
            nr = 2**(lev - 1)*2
            nth = 2**(lev - 1)*12
            call make_annulus_mesh(nr, nth, ra, rb, 0.0_dp, k, mm)
            call meridian_feec_build(mm, k, fs, status, cylinder_r0=r0)
            if (status /= 0) error stop "build failed"
            call solve_vacuum(fs, n, ep(lev), eb(lev))
        end do
        rp = log(ep(nlev - 1)/ep(nlev))/log(2.0_dp)
        rb_rate = log(eb(nlev - 1)/eb(nlev))/log(2.0_dp)
        write (tag, '(a,i0,a,i0,a,i0,a,f4.1)') "vacuum m ", m, " n ", n, &
            " order ", k, " R0 ", r0
        print '(a,a,a,2f7.2,a,2es10.2)', "  ", trim(tag), " rates psi/B", &
            rp, rb_rate, "  fine errors", ep(nlev), eb(nlev)
        call record(rp > real(k + 1, dp) - 0.4_dp, &
            "potential converges at order k + 1, "//trim(tag))
        call record(rb_rate > real(k, dp) - 0.4_dp, &
            "B = grad_n psi (incl. i n psi / R0) converges at order k, "// &
            trim(tag))
    end subroutine run_vacuum

end program test_meridian_cylinder_chart
