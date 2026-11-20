program test_toroidal_sheet_current_3d
    !! Energy of a surface current on a toroidal surface and its derivatives.
    !! Oracles:
    !!  - a net poloidal current i_p on a circular torus (R0, a) is a
    !!    toroidal-field coil: B = i_p/(2 pi R) inside and 0 outside, so
    !!    E = (1/2) int B^2 dV = (i_p^2/2) (R0 - sqrt(R0^2 - a^2));
    !!  - the load of the net-current ramps alone (phi = 0) in the field
    !!    B = (1 + cos v, 0, b0): int (i_t u + i_p v/2 pi) B.N = i_t 2 pi
    !!    (b0 U_r - U_z) - i_p pi A, U_r = int u R R_u, U_z = int u R Z_u,
    !!    A = oint R Z_u along the chain;
    !!  - the gradient and Hessian of the self energy and of the field load
    !!    against central differences on a non-axisymmetric surface (nfp = 2).
    use check, only: check_condition, check_summary
    use fortfem_kinds, only: dp
    use fortfem_toroidal_sheet_current_3d, only: ts_surface_t, ts_plan_t, &
        ts_surface_create, ts_size, ts_coefficients, ts_set_coefficients, &
        ts_energy, ts_load
    implicit none

    real(dp), parameter :: pi = acos(-1.0_dp)
    logical :: all_passed

    all_passed = .true.
    call test_toroidal_field_coil()
    call test_ramp_load()
    call test_derivatives()
    call check_summary("Toroidal sheet current energy")
    if (.not. all_passed) error stop 1

contains

    subroutine record(condition, description, err)
        logical, intent(in) :: condition
        character(*), intent(in) :: description
        real(dp), intent(in) :: err
        all_passed = all_passed .and. condition
        call check_condition(condition, description)
        write (*, "(a,es10.3)") "   relative error ", err
    end subroutine record

    subroutine circle(s, r0, a, shift)
        !! R = r0 + a cos(theta), Z = a sin(theta) at the chain nodes, theta =
        !! 2 pi u; shift adds a small non-axisymmetric part.
        type(ts_surface_t), intent(inout) :: s
        real(dp), intent(in) :: r0, a, shift
        integer :: i
        real(dp) :: th
        s%rz = 0
        do i = 1, s%n_nodes
            th = 2*pi*(i - 1)/s%n_nodes
            s%rz(1, i, 1) = r0 + a*cos(th)
            s%rz(2, i, 1) = a*sin(th)
            if (s%nmode >= 3) then
                s%rz(1, i, 2) = shift*cos(th)
                s%rz(2, i, 3) = shift*sin(th)
            end if
        end do
    end subroutine circle

    subroutine test_toroidal_field_coil()
        type(ts_surface_t) :: s
        type(ts_plan_t) :: plan
        real(dp) :: e, exact, err
        integer :: st, level
        real(dp), parameter :: r0 = 3.0_dp, a = 1.0_dp

        do level = 1, 2
            call ts_surface_create(3, 8*level, 1, 0, 0, s, st)
            call circle(s, r0, a, 0.0_dp)
            s%i_p = 1
            plan%nv = 8*level
            call ts_energy(s, plan, e)
            ! E = (i_p^2/(4 pi)) int int dR dZ/R over the P3 section,
            ! int int dR dZ/R = oint ln R dZ (Gauss along the panels)
            exact = section_integral(s)/(4*pi)
            err = abs(e - exact)/exact
            write (*, "(a,i0,a,2es20.12)") "   level ", level, ": E, exact ", e, exact
        end do
        call record(err < 1.0e-6_dp, "poloidal current energy is the toroidal-field coil's", &
            err)
        write (*, "(a,es10.3)") "   section vs circle ", &
            abs(exact - 0.5_dp*(r0 - sqrt(r0**2 - a**2)))/exact
    end subroutine test_toroidal_field_coil

    function section_integral(s) result(total)
        !! oint ln R dZ over the P_k chain (counter-clockwise), 20 Gauss points
        !! per panel.
        type(ts_surface_t), intent(in) :: s
        real(dp) :: total, x(20), w(20), xi, r, dz, n(s%degree + 1), dn(s%degree + 1)
        integer :: p, q, i, j, node
        call gauss01(x, w)
        total = 0
        do p = 1, s%n_panels
            do q = 1, 20
                xi = x(q)
                n = 0; dn = 0
                do i = 1, s%degree + 1
                    do j = s%degree + 1, 1, -1
                        n(i) = n(i)*xi + s%lagrange(i, j)
                        if (j > 1) dn(i) = dn(i)*xi + (j - 1)*s%lagrange(i, j)
                    end do
                end do
                r = 0; dz = 0
                do i = 1, s%degree + 1
                    node = (p - 1)*s%degree + i
                    if (node > s%n_nodes) node = node - s%n_nodes
                    r = r + n(i)*s%rz(1, node, 1)
                    dz = dz + dn(i)*s%rz(2, node, 1)
                end do
                total = total + w(q)*log(r)*dz
            end do
        end do
    end function section_integral

    subroutine gauss01(x, w)
        use fortnum_quadrature, only: gauss_legendre
        real(dp), intent(out) :: x(:), w(:)
        call gauss_legendre(size(x), x, w)
        x = 0.5_dp*(x + 1); w = 0.5_dp*w
    end subroutine gauss01

    subroutine test_ramp_load()
        type(ts_surface_t) :: s
        type(ts_plan_t) :: plan
        real(dp) :: e, exact, ur, uz, area, x(20), w(20), xi, uc, r, dr, dz
        real(dp) :: n(4), dn(4)
        integer :: st, p, q, i, j, node
        call ts_surface_create(3, 8, 1, 0, 0, s, st)
        call circle(s, 3.0_dp, 1.0_dp, 0.0_dp)
        s%i_t = 0.7_dp; s%i_p = 1.3_dp
        plan%nv = 8; plan%ng_v = 6
        call ts_load(s, plan, ramp_field, e)
        call gauss01(x, w)
        ur = 0; uz = 0; area = 0
        do p = 1, s%n_panels
            do q = 1, 20
                xi = x(q); uc = (p - 1 + xi)/s%n_panels
                n = 0; dn = 0
                do i = 1, 4
                    do j = 4, 1, -1
                        n(i) = n(i)*xi + s%lagrange(i, j)
                        if (j > 1) dn(i) = dn(i)*xi + (j - 1)*s%lagrange(i, j)
                    end do
                end do
                r = 0; dr = 0; dz = 0
                do i = 1, 4
                    node = (p - 1)*3 + i
                    if (node > s%n_nodes) node = node - s%n_nodes
                    r = r + n(i)*s%rz(1, node, 1)
                    dr = dr + dn(i)*s%rz(1, node, 1)
                    dz = dz + dn(i)*s%rz(2, node, 1)
                end do
                ur = ur + w(q)*uc*r*dr
                uz = uz + w(q)*uc*r*dz
                area = area + w(q)*r*dz
            end do
        end do
        exact = s%i_t*2*pi*(0.4_dp*ur - uz) - s%i_p*pi*area
        call record(abs(e - exact) < 1.0e-12_dp*abs(exact), &
            "load of the net-current ramps", abs(e - exact)/abs(exact))
    end subroutine test_ramp_load

    subroutine ramp_field(r, z, v, b, db, d2b)
        real(dp), intent(in) :: r, z, v
        real(dp), intent(out) :: b(3), db(3, 2), d2b(3, 3)
        b = [1 + cos(v), 0.0_dp, 0.4_dp]
        db = 0; d2b = 0
    end subroutine ramp_field

    subroutine field(r, z, v, b, db, d2b)
        !! A smooth test field with all jets: B = (z r, 1/r + sin v z, r^2 - z).
        real(dp), intent(in) :: r, z, v
        real(dp), intent(out) :: b(3), db(3, 2), d2b(3, 3)
        b = [z*r, 1/r + sin(v)*z, r*r - z]
        db(:, 1) = [z, -1/r**2, 2*r]
        db(:, 2) = [r, sin(v), -1.0_dp]
        d2b(:, 1) = [0.0_dp, 2/r**3, 2.0_dp]
        d2b(:, 2) = [1.0_dp, 0.0_dp, 0.0_dp]
        d2b(:, 3) = [0.0_dp, 0.0_dp, 0.0_dp]
    end subroutine field

    subroutine test_derivatives()
        type(ts_surface_t) :: s, t
        type(ts_plan_t) :: plan
        real(dp), allocatable :: c(:), dir(:), g(:), h(:, :), gp(:), gm(:), hp(:, :)
        real(dp) :: e, ep, em, step, err, fd
        integer :: st, n, i, kind

        call ts_surface_create(2, 6, 2, 1, 1, s, st)
        call circle(s, 3.0_dp, 0.8_dp, 0.1_dp)
        s%i_t = 0.7_dp
        s%i_p = 1.3_dp
        do i = 1, s%n_nodes
            s%phi(i, :) = [0.2_dp*sin(1.1_dp*i), 0.1_dp*cos(0.7_dp*i), 0.05_dp*sin(0.3_dp*i)]
        end do
        plan%nv = 4; plan%ng_u = 3; plan%ng_v = 3; plan%ng_far = 3
        plan%ng_t = 6; plan%ng_s = 6
        n = ts_size(s)
        allocate (dir(n), g(n), h(n, n), gp(n), gm(n), hp(n, n))
        do i = 1, n
            dir(i) = 0.01_dp*cos(1.9_dp*i + 0.3_dp)
        end do
        c = ts_coefficients(s)
        step = 1.0e-5_dp
        do kind = 1, 2
            if (kind == 1) then
                call ts_energy(s, plan, e, g, h)
            else
                call ts_load(s, plan, field, e, g, h)
            end if
            t = s
            call ts_set_coefficients(t, c + step*dir)
            if (kind == 1) then
                call ts_energy(t, plan, ep, gp, hp)
            else
                call ts_load(t, plan, field, ep, gp, hp)
            end if
            call ts_set_coefficients(t, c - step*dir)
            if (kind == 1) then
                call ts_energy(t, plan, em, gm, hp)
            else
                call ts_load(t, plan, field, em, gm, hp)
            end if
            fd = (ep - em)/(2*step)
            err = abs(fd - dot_product(g, dir))/sum(abs(g*dir))
            call record(err < 1.0e-7_dp, merge("self energy gradient", "load gradient       ", &
                kind == 1), err)
            err = maxval(abs((gp - gm)/(2*step) - matmul(h, dir)))/maxval(abs(matmul(h, dir)))
            call record(err < 1.0e-6_dp, merge("self energy Hessian ", "load Hessian        ", &
                kind == 1), err)
            err = maxval(abs(h - transpose(h)))/maxval(abs(h))
            call record(err < 1.0e-12_dp, merge("self Hessian symmetric", "load Hessian symmetric", &
                kind == 1), err)
        end do
    end subroutine test_derivatives

end program test_toroidal_sheet_current_3d
