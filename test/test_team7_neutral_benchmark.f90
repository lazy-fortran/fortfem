program test_team7_neutral_benchmark
    !! Independent analytic oracle for the TEAM-7-shaped neutral fixture.
    !!
    !! This test deliberately derives the field from a stream function and
    !! checks the two identities used by the gallery: div(B)=0 and
    !! curl(B)=-laplacian(A).  It does not inspect generated files or repeat
    !! the example's sampling loop.
    use check, only: check_condition, check_summary
    use fortfem_kinds, only: dp
    implicit none

    real(dp), parameter :: tolerance = 4.0e-6_dp
    real(dp), parameter :: sample_points(2, 5) = reshape([ &
        -0.85_dp, -0.62_dp, -0.30_dp, -0.18_dp, 0.0_dp, -0.35_dp, &
         0.25_dp, -0.12_dp, 0.72_dp,  0.54_dp], [2, 5])
    real(dp), parameter :: step = 2.0e-5_dp
    real(dp) :: divergence_error, curl_error, energy, bx, by, jz
    real(dp) :: x, y, a_re, a_im, ax_re, ay_re, axx_re, ayy_re
    real(dp) :: ax_im, ay_im, axx_im, ayy_im
    integer :: sample

    divergence_error = 0.0_dp
    curl_error = 0.0_dp
    energy = 0.0_dp
    do sample = 1, size(sample_points, 2)
        x = sample_points(1, sample)
        y = sample_points(2, sample)
        call analytic_field(x, y, bx, by, jz)
        divergence_error = max(divergence_error, abs( &
            (analytic_bx(x + step, y) - analytic_bx(x - step, y))/(2.0_dp*step) + &
            (analytic_by(x, y + step) - analytic_by(x, y - step))/(2.0_dp*step)))

        call analytic_stream_derivatives(x, y, a_re, a_im, ax_re, ay_re, &
            axx_re, ayy_re, ax_im, ay_im, axx_im, ayy_im)
        curl_error = max(curl_error, abs(jz + axx_re + ayy_re))
        energy = energy + 0.5_dp*(bx**2 + by**2)
    end do

    call check_condition(energy > 0.0_dp, &
        "TEAM-7-shaped manufactured field has positive magnetic energy")
    call check_condition(divergence_error < tolerance, &
        "stream-function field is numerically divergence free")
    call check_condition(curl_error < tolerance, &
        "Ampere current agrees with the analytical stream-function Laplacian")
    call check_condition(abs(analytic_bx(0.25_dp, -0.12_dp)) > 1.0e-8_dp .and. &
        abs(analytic_by(0.25_dp, -0.12_dp)) > 1.0e-8_dp, &
        "plate probe contains both tangential field components")
    call check_summary("TEAM-7 neutral benchmark oracle")

contains

    pure subroutine analytic_stream_derivatives(x, y, a_re, a_im, ax_re, &
            ay_re, axx_re, ayy_re, ax_im, ay_im, axx_im, ayy_im)
        real(dp), intent(in) :: x, y
        real(dp), intent(out) :: a_re, a_im, ax_re, ay_re, axx_re, ayy_re
        real(dp), intent(out) :: ax_im, ay_im, axx_im, ayy_im
        real(dp) :: value, dx, dy, dxx, dyy

        call component(x, y, -0.42_dp, 0.56_dp, 18.0_dp, 18.0_dp, 0.0_dp, &
            value, dx, dy, dxx, dyy)
        a_re = value
        ax_re = dx
        ay_re = dy
        axx_re = dxx
        ayy_re = dyy
        call component(x, y, 0.0_dp, -0.34_dp, 1.3_dp, 5.0_dp, &
            2.0_dp*acos(-1.0_dp)/1.4_dp, value, dx, dy, dxx, dyy)
        a_re = a_re + 0.42_dp*value
        ax_re = ax_re + 0.42_dp*dx
        ay_re = ay_re + 0.42_dp*dy
        axx_re = axx_re + 0.42_dp*dxx
        ayy_re = ayy_re + 0.42_dp*dyy
        a_im = 0.27_dp*value
        ax_im = 0.27_dp*dx
        ay_im = 0.27_dp*dy
        axx_im = 0.27_dp*dxx
        ayy_im = 0.27_dp*dyy
    end subroutine analytic_stream_derivatives

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

    pure subroutine analytic_field(x, y, bx, by, jz)
        real(dp), intent(in) :: x, y
        real(dp), intent(out) :: bx, by, jz
        real(dp) :: a_re, a_im, ax_re, ay_re, axx_re, ayy_re
        real(dp) :: ax_im, ay_im, axx_im, ayy_im

        call analytic_stream_derivatives(x, y, a_re, a_im, ax_re, ay_re, &
            axx_re, ayy_re, ax_im, ay_im, axx_im, ayy_im)
        bx = ay_re
        by = -ax_re
        jz = -(axx_re + ayy_re)
    end subroutine analytic_field

    pure real(dp) function analytic_bx(x, y)
        real(dp), intent(in) :: x, y
        real(dp) :: a_re, a_im, ax_re, ay_re, axx_re, ayy_re
        real(dp) :: ax_im, ay_im, axx_im, ayy_im

        call analytic_stream_derivatives(x, y, a_re, a_im, ax_re, ay_re, &
            axx_re, ayy_re, ax_im, ay_im, axx_im, ayy_im)
        analytic_bx = ay_re
    end function analytic_bx

    pure real(dp) function analytic_by(x, y)
        real(dp), intent(in) :: x, y
        real(dp) :: a_re, a_im, ax_re, ay_re, axx_re, ayy_re
        real(dp) :: ax_im, ay_im, axx_im, ayy_im

        call analytic_stream_derivatives(x, y, a_re, a_im, ax_re, ay_re, &
            axx_re, ayy_re, ax_im, ay_im, axx_im, ayy_im)
        analytic_by = -ax_re
    end function analytic_by

end program test_team7_neutral_benchmark
