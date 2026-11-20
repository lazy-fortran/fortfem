program test_toroidal_sheet_trapezoid
    !! Corrected trapezoidal rule for the energy of a surface current.
    !! Oracles:
    !!  - a net poloidal current i_p on a circular torus (R0, a) is a
    !!    toroidal-field coil: E = (i_p^2/2) (R0 - sqrt(R0^2 - a^2)); the rule
    !!    converges at third order (nfp = 3, so the period symmetry is used);
    !!  - the gradient and the Hessian-vector product against central
    !!    differences on a non-axisymmetric surface (nfp = 2), and the
    !!    current block and the dense Hessian against the product, with one
    !!    level and with the two-level partition.
    use check, only: check_condition, check_summary
    use fortfem_kinds, only: dp
    use fortfem_toroidal_sheet_trapezoid, only: tst_grid_t, tst_energy, tst_hvp, &
        tst_current_modes, tst_hessian
    implicit none

    real(dp), parameter :: pi = acos(-1.0_dp)
    logical :: all_passed

    all_passed = .true.
    call test_toroidal_field_coil()
    call test_derivatives(0.0_dp)
    call test_two_level()
    call test_derivatives(2.5_dp)
    call check_summary("Toroidal sheet trapezoidal rule")
    if (.not. all_passed) error stop 1

contains

    subroutine record(condition, description, err)
        logical, intent(in) :: condition
        character(*), intent(in) :: description
        real(dp), intent(in) :: err
        all_passed = all_passed .and. condition
        call check_condition(condition, description)
        write (*, "(a,es10.3)") "   error ", err
    end subroutine record

    subroutine torus(g, u, shift)
        type(tst_grid_t), intent(in) :: g
        real(dp), intent(out) :: u(:, :, :)
        real(dp), intent(in) :: shift
        real(dp) :: th, v
        integer :: i, j
        do j = 1, g%nv
            v = 2*pi*(j - 1)/(g%nfp*g%nv)
            do i = 1, g%nu
                th = 2*pi*(i - 1)/g%nu
                u(1, i, j) = 3 + cos(th) + shift*cos(th - g%nfp*v)
                u(2, i, j) = sin(th) + shift*sin(th - g%nfp*v)
                u(3, i, j) = -2*pi*(sin(th) + shift*sin(th - g%nfp*v))
                u(4, i, j) = 2*pi*(cos(th) + shift*cos(th - g%nfp*v))
                u(5, i, j) = g%nfp*shift*sin(th - g%nfp*v)
                u(6, i, j) = -g%nfp*shift*cos(th - g%nfp*v)
                u(7, i, j) = 0
                u(8, i, j) = 1/(2*pi)
            end do
        end do
    end subroutine torus

    subroutine test_toroidal_field_coil()
        type(tst_grid_t) :: g
        real(dp), allocatable :: u(:, :, :)
        real(dp) :: e, exact, err(2)
        integer :: level
        exact = 0.5_dp*(3 - sqrt(8.0_dp))
        do level = 1, 2
            g = tst_grid_t(nu=16*level, nv=16*level, nfp=3)
            allocate (u(8, g%nu, g%nv))
            call torus(g, u, 0.0_dp)
            call tst_energy(g, u, e)
            err(level) = abs(e - exact)/exact
            write (*, "(a,i0,a,2es20.12)") "   nu ", g%nu, ": E, exact ", e, exact
            deallocate (u)
        end do
        call record(err(2) < 1.0e-4_dp .and. err(1)/err(2) > 7, &
            "toroidal-field coil energy at third order", err(2))
    end subroutine test_toroidal_field_coil

    subroutine test_two_level()
        !! The toroidal-field coil with the two-level partition (r0 = 3 at 64 x
        !! 64 points, nfp = 3): the error stays at the one-level value (it
        !! grows to 1.3e-4 at r0 = 1.2, where the step spans 3 coarse cells).
        type(tst_grid_t) :: g
        real(dp), allocatable :: u(:, :, :)
        real(dp) :: e1, e2, exact
        exact = 0.5_dp*(3 - sqrt(8.0_dp))
        g = tst_grid_t(nu=64, nv=64, nfp=3)
        allocate (u(8, 64, 64))
        call torus(g, u, 0.0_dp)
        call tst_energy(g, u, e1)
        ! r0/2 is about 7 coarse spacings: the coarse sum resolves the step
        g%r0 = 3.0_dp
        call tst_energy(g, u, e2)
        write (*, "(a,2es12.4)") "   one-level, two-level error ", abs(e1 - exact)/exact, &
            abs(e2 - exact)/exact
        call record(abs(e2 - exact)/exact < 2*abs(e1 - exact)/exact, &
            "two-level partition keeps the accuracy", abs(e2 - e1)/exact)
    end subroutine test_two_level

    subroutine test_derivatives(r0)
        real(dp), intent(in) :: r0
        type(tst_grid_t) :: g
        real(dp), allocatable :: u(:, :, :), du(:, :, :), gr(:, :, :), gp(:, :, :), &
            gm(:, :, :), hv(:, :, :)
        real(dp) :: e, ep, em, h, err
        integer :: i, j, c
        g = tst_grid_t(nu=12, nv=8, nfp=2, r0=r0)
        write (*, "(a,f5.2)") "   r0 ", r0
        allocate (u(8, 12, 8), du(8, 12, 8), gr(8, 12, 8), gp(8, 12, 8), gm(8, 12, 8), &
            hv(8, 12, 8))
        call torus(g, u, 0.1_dp)
        do j = 1, 8
            do i = 1, 12
                u(7, i, j) = 0.3_dp*sin(1.3_dp*i + 0.4_dp*j)
                u(8, i, j) = u(8, i, j) + 0.2_dp*cos(0.7_dp*i - 1.1_dp*j)
                do c = 1, 8
                    du(c, i, j) = 0.01_dp*cos(1.9_dp*c + 2.3_dp*i + 0.7_dp*j)
                end do
            end do
        end do
        call tst_energy(g, u, e, gr)
        call tst_hvp(g, u, du, hv)
        h = 1.0e-5_dp
        call tst_energy(g, u + h*du, ep, gp)
        call tst_energy(g, u - h*du, em, gm)
        err = abs((ep - em)/(2*h) - sum(gr*du))/sum(abs(gr*du))
        call record(err < 1.0e-8_dp, "gradient", err)
        err = maxval(abs((gp - gm)/(2*h) - hv))/maxval(abs(hv))
        call record(err < 1.0e-7_dp, "Hessian-vector product", err)
        call test_modes(g, u)
        block
            real(dp), allocatable :: hd(:, :, :, :)
            allocate (hd(8, 96, 8, 96))
            call tst_hessian(g, u, hd)
            err = maxval(abs(reshape(matmul(reshape(hd, [768, 768]), reshape(du, [768])), &
                [8, 12, 8]) - hv))/maxval(abs(hv))
            call record(err < 1.0e-12_dp, "dense Hessian", err)
        end block
    end subroutine test_derivatives

    subroutine test_modes(g, u)
        !! The current block in the basis phi(i, n) T_n(v) (2 modes) against
        !! the product along a basis direction.
        type(tst_grid_t), intent(in) :: g
        real(dp), intent(in) :: u(:, :, :)
        real(dp) :: du(12, 12), t(3, 8), dt(3, 8), m(36, 36), phi(12, 3), dq(8, 12, 8), &
            hp(8, 12, 8), v, y(36), err
        integer :: i, k, j
        du = 0
        do k = 1, 12
            do i = 1, 12
                if (i /= k) du(i, k) = 2*pi*0.5_dp*(-1)**mod(abs(i - k), 2) &
                    /tan((i - k)*pi/12)
            end do
        end do
        do j = 1, 8
            v = 2*pi*(j - 1)/(2*8)
            t(:, j) = [1.0_dp, cos(2*v), sin(2*v)]
            dt(:, j) = [0.0_dp, -2*sin(2*v), 2*cos(2*v)]
        end do
        call tst_current_modes(g, u, du, t, dt, m)
        do k = 1, 3
            do i = 1, 12
                phi(i, k) = cos(0.7_dp*i + 1.3_dp*k)
            end do
        end do
        dq = 0
        dq(7, :, :) = matmul(du, matmul(phi, t))
        dq(8, :, :) = matmul(phi, dt)
        call tst_hvp(g, u, dq, hp)
        y = reshape(matmul(matmul(transpose(du), hp(7, :, :)), transpose(t)) &
            + matmul(hp(8, :, :), transpose(dt)), [36])
        err = maxval(abs(matmul(m, reshape(phi, [36])) - y))/maxval(abs(y))
        call record(err < 1.0e-12_dp, "current block in the mode basis", err)
    end subroutine test_modes

end program test_toroidal_sheet_trapezoid
