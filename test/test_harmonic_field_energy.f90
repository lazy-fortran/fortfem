program test_harmonic_field_energy
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
    use fortfem_kinds, only: dp
    use check, only: check_condition, check_summary
    use fortnum_quadrature, only: gauss_legendre_ab
    use fortfem_harmonic_field_energy, only: assemble_field_gram, &
        couple_harmonic_sector
    use fortfem_generated_field_inner_product, only: field_inner_product
    use fortfem_generated_field_inner_product_jvp, only: field_inner_product_jvp
    use fortfem_generated_field_inner_product_vjp, only: field_inner_product_vjp
    implicit none
    integer, parameter :: nr = 12, nt = 32, np = 32, nq = nr*nt*np
    real(dp), parameter :: pi = acos(-1.0_dp), radius = 3.0_dp
    real(dp) :: fields(3, nq, 2), weights(nq), gram(2, 2), nodes(nr), rw(nr)
    real(dp) :: scalar(1, 1), cross(1, 1), harmonic(1, 1), coordinates(1, 1)
    real(dp) :: target(1), u(1), coefficient(1), load(1), energy, expected
    real(dp) :: r, theta, phi, cylindrical_r, volume, henergy, flux
    integer :: i, j, k, q, status
    logical :: passed

    passed = .true.
    call gauss_legendre_ab(nr, 0.0_dp, 1.0_dp, nodes, rw)
    q = 0
    do k = 1, np
        phi = 2*pi*real(k, dp)/np
        do j = 1, nt
            theta = 2*pi*real(j, dp)/nt
            do i = 1, nr
                q = q+1
                r = nodes(i)
                cylindrical_r = radius+r*cos(theta)
                weights(q) = rw(i)*r*cylindrical_r*(2*pi/nt)*(2*pi/np)
                fields(:, q, 1) = [1.0_dp, 0.0_dp, 0.0_dp]
                ! H = grad(phi)+grad(x), an intentionally non-tangent
                ! representative of the toroidal cohomology class.
                fields(:, q, 2) = [-sin(phi)/cylindrical_r+1, &
                    cos(phi)/cylindrical_r, 0.0_dp]
            end do
        end do
    end do
    volume = 2*pi*pi*radius
    flux = 2*pi/(radius+sqrt(radius*radius-1))
    henergy = 2*pi*flux
    call assemble_field_gram(weights, fields, gram, status)
    call record(status == 0, 'physical field Gram assembly succeeds')
    call record(abs(gram(1, 1)/volume-1) < 1.e-12_dp, 'solid torus volume')
    call record(abs(gram(1, 2)/volume-1) < 1.e-12_dp, &
        'nonzero scalar/harmonic cross energy equals volume')
    call record(abs(gram(2, 2)/(volume+henergy)-1) < 1.e-12_dp, &
        'harmonic energy equals period times meridional flux plus volume')
    scalar = gram(1:1, 1:1)
    cross = gram(1:1, 2:2)
    harmonic = gram(2:2, 2:2)
    coordinates = 2*pi
    target = 4*pi
    u = -2.0_dp
    call couple()
    call record(status == 0 .and. abs(coefficient(1)-2) < 1.e-13_dp, &
        'prescribed toroidal circulation fixes physical amplitude')
    call record(abs(energy/(4*henergy)-1) < 1.e-11_dp, &
        'scalar correction removes exact part with full cross energy')
    call record(abs(load(1)) < 1.e-11_dp*volume, &
        'corrected tangent harmonic field is orthogonal to scalar x')
    coordinates = flux
    target = -3*flux
    u = 0.75_dp
    call couple()
    expected = (u(1)-3)**2*volume+9*henergy
    call record(status == 0 .and. abs(coefficient(1)+3) < 1.e-13_dp, &
        'prescribed meridional cut flux uses its own physical units')
    call record(abs(energy/expected-1) < 1.e-12_dp, &
        'flux prescribed sector retains signed cross energy')
    coordinates = 0.0_dp
    call couple()
    call record(status /= 0, 'singular coordinate map is rejected')
    coordinates = flux
    harmonic = -1.0_dp
    call couple()
    call record(status /= 0, 'negative harmonic metric is rejected')
    harmonic = gram(2:2, 2:2)
    weights(1) = -1.0_dp
    call assemble_field_gram(weights, fields, gram, status)
    call record(status /= 0, 'negative volume weight is rejected')
    weights(1) = ieee_value(0.0_dp, ieee_quiet_nan)
    call assemble_field_gram(weights, fields, gram, status)
    call record(status /= 0, 'nonfinite quadrature weight is rejected')
    call multiple_sectors()
    call derivatives()
    call check_summary('Physical harmonic field energy on a solid torus')
    if (.not. passed) error stop 1
contains
    subroutine multiple_sectors()
        real(dp) :: s(1, 1), c(1, 2), g(2, 2), p(2, 2), t(2), a(2), f(1), e
        s = 2.0_dp
        c = reshape([0.5_dp, -0.5_dp], [1, 2])
        g = reshape([3.0_dp, 1.0_dp, 1.0_dp, 2.0_dp], [2, 2])
        p = reshape([1.0_dp, 2.0_dp, 3.0_dp, 1.0_dp], [2, 2])
        t = [-1.0_dp, 3.0_dp]
        call couple_harmonic_sector(s, c, g, p, t, [2.0_dp], &
            a, f, e, status)
        call record(status == 0 .and. maxval(abs(a-[2.0_dp, -1.0_dp])) &
            < 1.e-13_dp, 'coupled coordinate map handles multiple sectors')
        call record(abs(e-24) < 1.e-12_dp .and. abs(f(1)-5.5_dp) < 1.e-12_dp, &
            'multiple sectors include off-diagonal and scalar cross energies')
        g(1, 2) = 2.0_dp
        call couple_harmonic_sector(s, c, g, p, t, [2.0_dp], &
            a, f, e, status)
        call record(status /= 0, 'nonsymmetric harmonic metric is rejected')
        g(1, 2) = 1.0_dp
        c = 100.0_dp
        call couple_harmonic_sector(s, c, g, p, t, [-2.0_dp], &
            a, f, e, status)
        call record(status /= 0 .and. e == 0, &
            'inconsistent blocks giving negative total energy are rejected')
    end subroutine
    subroutine derivatives()
        real(dp) :: v(7), d(7), bars(7), jd, plus, minus, step
        v = [2.0_dp, 1.0_dp, -3.0_dp, 2.0_dp, 4.0_dp, -2.0_dp, 1.0_dp]
        d = [0.1_dp, -0.2_dp, 0.3_dp, 0.4_dp, -0.5_dp, 0.6_dp, -0.7_dp]
        call field_inner_product_jvp(v(1), v(2), v(3), v(4), v(5), v(6), &
            v(7), d(1), d(2), d(3), d(4), d(5), d(6), d(7), jd)
        call field_inner_product_vjp(v(1), v(2), v(3), v(4), v(5), v(6), &
            v(7), 1.0_dp, bars(1), bars(2), bars(3), bars(4), bars(5), &
            bars(6), bars(7))
        call record(abs(jd-dot_product(bars, d)) < 1.e-12_dp, &
            'generated field contraction JVP/VJP satisfy duality')
        step = 1.e-5_dp
        v = v+step*d
        call field_inner_product(v(1), v(2), v(3), v(4), v(5), v(6), v(7), plus)
        v = v-2*step*d
        call field_inner_product(v(1), v(2), v(3), v(4), v(5), v(6), v(7), minus)
        call record(abs((plus-minus)/(2*step)-jd) < 1.e-8_dp, &
            'generated derivative agrees with independent central difference')
    end subroutine
    subroutine couple()
        call couple_harmonic_sector(scalar, cross, harmonic, coordinates, &
            target, u, coefficient, load, energy, status)
    end subroutine
    subroutine record(ok, label)
        logical, intent(in) :: ok
        character(*), intent(in) :: label
        passed = passed .and. ok
        call check_condition(ok, label)
    end subroutine
end program
