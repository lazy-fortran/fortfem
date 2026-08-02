program test_team3_neutral_benchmark
    !! Independent analytic oracle for the TEAM-3-shaped neutral fixture.
    !! It checks identities from a separately written stream-function oracle,
    !! rather than inspecting generated gallery files.
    use check, only: check_condition, check_summary
    use fortfem_kinds, only: dp
    implicit none

    real(dp), parameter :: step = 2.0e-5_dp
    real(dp), parameter :: tolerance = 5.0e-6_dp
    real(dp), parameter :: sample_points(2, 7) = reshape([ &
        -1.05_dp, -0.72_dp, -0.67_dp, -0.15_dp, 0.0_dp, 0.22_dp, &
         0.76_dp, -0.35_dp,  0.00_dp,  0.41_dp, 0.66_dp, -0.54_dp, &
         0.83_dp, 0.23_dp], [2, 7])
    real(dp) :: divergence_error, ampere_error, energy, bx, by, jz
    real(dp) :: x, y, az, ax, ay, axx, ayy
    integer :: sample

    divergence_error = 0.0_dp
    ampere_error = 0.0_dp
    energy = 0.0_dp
    do sample = 1, size(sample_points, 2)
        x = sample_points(1, sample)
        y = sample_points(2, sample)
        call analytic_field(x, y, bx, by, jz)
        divergence_error = max(divergence_error, abs( &
            (analytic_bx(x + step, y) - analytic_bx(x - step, y))/(2.0_dp*step) + &
            (analytic_by(x, y + step) - analytic_by(x, y - step))/(2.0_dp*step)))
        call analytic_stream_derivatives(x, y, az, ax, ay, axx, ayy)
        ampere_error = max(ampere_error, abs(jz + axx + ayy))
        energy = energy + 0.5_dp*(bx**2 + by**2)
    end do

    call check_condition(energy > 0.0_dp, &
        "TEAM-3-shaped manufactured field has positive magnetic energy")
    call check_condition(divergence_error < tolerance, &
        "stream-function field is numerically divergence free")
    call check_condition(ampere_error < 2.0e-13_dp, &
        "Ampere current equals the analytical stream-function Laplacian")
    call check_condition(abs(analytic_bx(-0.15_dp, 0.41_dp)) > 1.0e-8_dp .and. &
        abs(analytic_by(-0.15_dp, 0.41_dp)) > 1.0e-8_dp, &
        "air-gap probe has both field components")
    call check_summary("TEAM-3 neutral benchmark oracle")

contains

    pure subroutine analytic_stream_derivatives(x, y, az, ax, ay, axx, ayy)
        real(dp), intent(in) :: x, y
        real(dp), intent(out) :: az, ax, ay, axx, ayy
        real(dp) :: value, dx, dy, dxx, dyy

        call component(x, y, -0.40_dp, 0.0_dp, 3.0_dp, 2.0_dp, 2.2_dp, &
            value, dx, dy, dxx, dyy)
        az = value
        ax = dx
        ay = dy
        axx = dxx
        ayy = dyy
        call component(x, y, 0.32_dp, 0.0_dp, 1.2_dp, 2.5_dp, 3.1_dp, &
            value, dx, dy, dxx, dyy)
        az = az + 0.38_dp*value
        ax = ax + 0.38_dp*dx
        ay = ay + 0.38_dp*dy
        axx = axx + 0.38_dp*dxx
        ayy = ayy + 0.38_dp*dyy
        call component(x, y, 0.0_dp, 0.0_dp, 0.35_dp, 0.45_dp, 0.0_dp, &
            value, dx, dy, dxx, dyy)
        az = az + 0.10_dp*value
        ax = ax + 0.10_dp*dx
        ay = ay + 0.10_dp*dy
        axx = axx + 0.10_dp*dxx
        ayy = ayy + 0.10_dp*dyy
    end subroutine analytic_stream_derivatives

    pure subroutine analytic_field(x, y, bx, by, jz)
        real(dp), intent(in) :: x, y
        real(dp), intent(out) :: bx, by, jz
        real(dp) :: az, ax, ay, axx, ayy

        call analytic_stream_derivatives(x, y, az, ax, ay, axx, ayy)
        bx = ay
        by = -ax
        jz = -(axx + ayy)
    end subroutine analytic_field

    pure real(dp) function analytic_bx(x, y)
        real(dp), intent(in) :: x, y
        real(dp) :: bx, by, jz

        call analytic_field(x, y, bx, by, jz)
        analytic_bx = bx
    end function analytic_bx

    pure real(dp) function analytic_by(x, y)
        real(dp), intent(in) :: x, y
        real(dp) :: bx, by, jz

        call analytic_field(x, y, bx, by, jz)
        analytic_by = by
    end function analytic_by

    pure subroutine component(x, y, center_x, center_y, alpha, beta, wave, &
            value, dx, dy, dxx, dyy)
        real(dp), intent(in) :: x, y, center_x, center_y, alpha, beta, wave
        real(dp), intent(out) :: value, dx, dy, dxx, dyy
        real(dp) :: xx, yy, envelope, cosine, sine

        xx = x - center_x
        yy = y - center_y
        envelope = exp(-alpha*xx**2 - beta*yy**2)
        cosine = cos(wave*xx)
        sine = sin(wave*xx)
        value = envelope*cosine
        dx = envelope*(-2.0_dp*alpha*xx*cosine - wave*sine)
        dy = envelope*(-2.0_dp*beta*yy*cosine)
        dxx = envelope*((4.0_dp*alpha**2*xx**2 - 2.0_dp*alpha - wave**2)* &
            cosine + 4.0_dp*alpha*wave*xx*sine)
        dyy = envelope*(4.0_dp*beta**2*yy**2 - 2.0_dp*beta)*cosine
    end subroutine component

end program test_team3_neutral_benchmark
