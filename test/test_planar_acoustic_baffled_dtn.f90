program test_planar_acoustic_baffled_dtn
    !! Baffled planar DtN form (period_sample_count = lfac*m) against the
    !! Rayleigh integral of a finite interface in an infinite hard wall.
    !!
    !! Oracle (time convention exp(-i omega t), lossy fluid k = omega/c with
    !! Im(k) > 0): a normal displacement u on the line y = 0, zero on the
    !! baffle, radiates into y > 0 the wall pressure
    !!   p(x) = -(i/2) rho omega**2 int H0(k|x - x'|) u(x') dx',
    !! the Neumann Green's function -2 (i/4) H0 with d p/d y = rho omega**2 u;
    !! its Fourier symbol is the DtN multiplier -i rho omega**2/k_y. The
    !! Galerkin entries for P1 hats f, g on the sample grid are then
    !!   F(f, g) = -(i/2) rho omega**2 int H0(k|s|) C(s) ds,
    !!   C(s) = int f(x) g(x - s) dx,
    !! evaluated with H0(z) = -(2i/pi) K0(-i z), C exact by Gauss on its linear
    !! pieces and s = h t**3 near the logarithmic singularity at s = 0.
    !! Compared: the piston force sum_ab F_ab (all hats), the load on the end
    !! sample sum_b F_1b, and the corner entry F_1m. The periodic form
    !! (lfac = 1) makes the two end samples neighbours (F_1m = F_12); the
    !! baffled form removes that wrap-around and converges to the oracle as
    !! lfac grows, up to the sampling error of the interface, which falls when
    !! h is halved. With m = 48 samples over two wavelengths and 5 % loss the
    !! end-load error is 0.94, 0.11, 1.4e-2, 8.9e-4, 1.1e-3 for lfac = 1, 2, 4,
    !! 8, 16 (force 4.3e-2 down to 4.0e-5), and 5.0e-4 at lfac = 16, m = 96.
    !! The images of the periodic DtN decay only through the fluid loss, so a
    !! lossless fluid needs a larger lfac.
    use check, only: check_condition, check_summary
    use fortfem_boundary, only: assemble_planar_acoustic_displacement_dtn_form
    use fortfem_kinds, only: dp
    use fortnum_special_complex_bessel, only: bessel_k_complex
    use fortnum_status, only: fortnum_status_t
    implicit none

    real(dp), parameter :: pi = acos(-1.0_dp)
    real(dp), parameter :: omega = 2.0_dp*pi, rho = 1.3_dp
    real(dp), parameter :: width = 2.0_dp
    complex(dp), parameter :: c = (1.0_dp, -0.05_dp)
    integer, parameter :: factors(5) = [1, 2, 4, 8, 16]
    complex(dp) :: k
    logical :: all_passed

    all_passed = .true.
    k = omega/c
    call compare(48, .true.)
    call compare(96, .false.)
    call check_rejections()
    call check_summary("Baffled planar acoustic DtN against the Rayleigh integral")
    if (.not. all_passed) error stop 1

contains

    subroutine compare(m, check_lfac)
        integer, intent(in) :: m
        logical, intent(in) :: check_lfac

        complex(dp), allocatable :: form(:, :)
        real(dp), allocatable :: piston(:), first(:), last(:)
        complex(dp) :: ref_force, ref_end, ref_corner
        real(dp) :: h, force_error(5), end_error(5), corner_error(5), scale
        integer :: i, status, total
        real(dp), save :: coarse_end_error = -1.0_dp

        h = width/real(m, dp)
        allocate(form(m, m), piston(m), first(m), last(m))
        piston = 1.0_dp
        first = 0.0_dp
        first(1) = 1.0_dp
        last = 0.0_dp
        last(m) = 1.0_dp
        ref_force = rayleigh(piston, piston, h)
        ref_end = rayleigh(first, piston, h)
        ref_corner = rayleigh(first, last, h)
        do i = 1, size(factors)
            total = factors(i)*m
            call assemble_planar_acoustic_displacement_dtn_form(m, omega, c, &
                rho, h*real(total, dp), total/2, form, status, &
                period_sample_count=total)
            call record(status == 0, "baffled form assembles")
            scale = abs(form(1, 2))
            force_error(i) = abs(sum(form) - ref_force)/abs(ref_force)
            end_error(i) = abs(sum(form(1, :)) - ref_end)/abs(ref_end)
            corner_error(i) = abs(form(1, m) - ref_corner)/scale
            write (*, "(a,i4,a,i3,a,3es10.3)") "m ", m, " lfac ", factors(i), &
                ": force, end load, corner/|F12| errors ", force_error(i), &
                end_error(i), corner_error(i)
        end do
        if (check_lfac) then
            call assemble_planar_acoustic_displacement_dtn_form(m, omega, c, &
                rho, width, m/2, form, status)
            call record(status == 0 .and. abs(form(1, m) - form(1, 2)) <= &
                1.0e-12_dp*scale, "periodic form couples the end samples as neighbours")
            call record(corner_error(1) > 0.5_dp, &
                "periodic corner entry is an end-sample artefact")
            call record(corner_error(4) < 2.0e-3_dp .and. &
                corner_error(5) < 1.0e-4_dp, &
                "baffled corner entry matches the Rayleigh integral")
            call record(all(force_error(2:4) < 0.2_dp*force_error(1:3)) .and. &
                all(end_error(2:4) < 0.2_dp*end_error(1:3)), &
                "force and end load converge geometrically up to lfac = 8")
            call record(force_error(5) < 1.0e-4_dp .and. &
                end_error(5) < 2.0e-3_dp, &
                "lfac = 16 matches the Rayleigh force and end load")
            coarse_end_error = end_error(5)
        else
            call record(end_error(5) < 0.6_dp*coarse_end_error, &
                "the remaining end-load error is sampling error that falls with h")
        end if
    end subroutine compare

    subroutine check_rejections()
        complex(dp) :: form(8, 8)
        integer :: status

        call assemble_planar_acoustic_displacement_dtn_form(8, omega, c, rho, &
            1.0_dp, 3, form, status, period_sample_count=7)
        call record(status /= 0, "period shorter than the interface is rejected")
        call assemble_planar_acoustic_displacement_dtn_form(8, omega, c, rho, &
            1.0_dp, 9, form, status, period_sample_count=16)
        call record(status /= 0, "modes beyond the extended Nyquist are rejected")
    end subroutine check_rejections

    complex(dp) function rayleigh(f, g, h) result(value)
        !! -(i/2) rho omega**2 int H0(k|s|) C(s) ds for nodal hats f, g on
        !! samples 0..m-1 of spacing h.
        real(dp), intent(in) :: f(:), g(:), h

        integer, parameter :: order = 10
        real(dp) :: nodes(order), weights(order), s, t
        integer :: j, q, m

        m = size(f)
        call gauss_legendre(order, nodes, weights)
        value = (0.0_dp, 0.0_dp)
        do j = -m - 1, m
            do q = 1, order
                t = 0.5_dp*(nodes(q) + 1.0_dp)
                if (j == 0) then
                    s = h*t**3
                    value = value + 0.5_dp*weights(q)*3.0_dp*h*t**2* &
                        hankel0(s)*correlation(f, g, h, s)
                else if (j == -1) then
                    s = -h*t**3
                    value = value + 0.5_dp*weights(q)*3.0_dp*h*t**2* &
                        hankel0(-s)*correlation(f, g, h, s)
                else
                    s = h*(real(j, dp) + t)
                    value = value + 0.5_dp*weights(q)*h*hankel0(abs(s))* &
                        correlation(f, g, h, s)
                end if
            end do
        end do
        value = -(0.0_dp, 0.5_dp)*rho*omega**2*value
    end function rayleigh

    complex(dp) function hankel0(r)
        real(dp), intent(in) :: r

        complex(dp) :: k0
        type(fortnum_status_t) :: status

        call bessel_k_complex(0, -(0.0_dp, 1.0_dp)*k*r, .false., k0, status)
        if (status%code /= 0) error stop "K0 failed"
        hankel0 = -(0.0_dp, 2.0_dp)/pi*k0
    end function hankel0

    real(dp) function correlation(f, g, h, s)
        !! int f(x) g(x - s) dx, exact: two-point Gauss on each linear piece.
        real(dp), intent(in) :: f(:), g(:), h, s

        real(dp) :: shift, pieces(3), a, b, mid, half, x
        integer :: cell, piece, q

        correlation = 0.0_dp
        shift = s/h - floor(s/h)
        do cell = -1, size(f)
            pieces = [real(cell, dp), real(cell, dp) + shift, &
                real(cell + 1, dp)]
            do piece = 1, 2
                a = pieces(piece)
                b = pieces(piece + 1)
                if (b - a <= 0.0_dp) cycle
                mid = 0.5_dp*(a + b)
                half = 0.5_dp*(b - a)
                do q = -1, 1, 2
                    x = mid + real(q, dp)*half/sqrt(3.0_dp)
                    correlation = correlation + half*h*nodal(f, x)* &
                        nodal(g, x - s/h)
                end do
            end do
        end do
    end function correlation

    pure real(dp) function nodal(v, x)
        !! P1 interpolant of v(1:m) at samples 0..m-1, zero beyond.
        real(dp), intent(in) :: v(:), x

        integer :: j
        real(dp) :: w, left, right

        j = floor(x)
        w = x - real(j, dp)
        left = 0.0_dp
        right = 0.0_dp
        if (j >= 0 .and. j < size(v)) left = v(j + 1)
        if (j + 1 >= 0 .and. j + 1 < size(v)) right = v(j + 2)
        nodal = (1.0_dp - w)*left + w*right
    end function nodal

    subroutine gauss_legendre(n, x, w)
        integer, intent(in) :: n
        real(dp), intent(out) :: x(n), w(n)

        real(dp) :: z, p0, p1, p2, dp1
        integer :: i, j, it

        do i = 1, n
            z = cos(pi*(real(i, dp) - 0.25_dp)/(real(n, dp) + 0.5_dp))
            do it = 1, 100
                p0 = 1.0_dp
                p1 = z
                do j = 2, n
                    p2 = (real(2*j - 1, dp)*z*p1 - real(j - 1, dp)*p0)/real(j, dp)
                    p0 = p1
                    p1 = p2
                end do
                dp1 = real(n, dp)*(z*p1 - p0)/(z*z - 1.0_dp)
                z = z - p1/dp1
                if (abs(p1/dp1) < 1.0e-15_dp) exit
            end do
            x(i) = z
            w(i) = 2.0_dp/((1.0_dp - z*z)*dp1*dp1)
        end do
    end subroutine gauss_legendre

    subroutine record(condition, description)
        logical, intent(in) :: condition
        character(*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record

end program test_planar_acoustic_baffled_dtn
