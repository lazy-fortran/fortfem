program test_toroidal_sheet_trapezoid
    !! Corrected trapezoidal rule for the energy of a surface current.
    !! Oracles:
    !!  - a net poloidal current i_p on a circular torus (R0, a) is a
    !!    toroidal-field coil: E = (i_p^2/2) (R0 - sqrt(R0^2 - a^2)); the rule
    !!    converges at third order (nfp = 3, so the period symmetry is used);
    !!  - the gradient and the Hessian-vector product against central
    !!    differences on a non-axisymmetric surface (nfp = 2), and the
    !!    current block of the Hessian against the product.
    use check, only: check_condition, check_summary
    use fortfem_kinds, only: dp
    use fortfem_toroidal_sheet_trapezoid, only: tst_grid_t, tst_energy, tst_hvp, &
        tst_current_form
    implicit none

    real(dp), parameter :: pi = acos(-1.0_dp)
    logical :: all_passed

    all_passed = .true.
    call test_toroidal_field_coil()
    call test_derivatives()
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

    subroutine test_derivatives()
        type(tst_grid_t) :: g
        real(dp), allocatable :: u(:, :, :), du(:, :, :), gr(:, :, :), gp(:, :, :), &
            gm(:, :, :), hv(:, :, :)
        real(dp) :: e, ep, em, h, err
        integer :: i, j, c
        g = tst_grid_t(nu=12, nv=8, nfp=2)
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
        block
            real(dp) :: form(2, 96, 2, 96), dq(8, 12, 8), hp(8, 12, 8), y(2, 96)
            dq = 0
            dq(7:8, :, :) = du(7:8, :, :)
            call tst_hvp(g, u, dq, hp)
            call tst_current_form(g, u, form)
            y = reshape(matmul(reshape(form, [192, 192]), reshape(dq(7:8, :, :), [192])), &
                [2, 96])
            err = maxval(abs(y - reshape(hp(7:8, :, :), [2, 96])))/maxval(abs(hp(7:8, :, :)))
            call record(err < 1.0e-12_dp, "current block of the Hessian", err)
        end block
    end subroutine test_derivatives

end program test_toroidal_sheet_trapezoid
