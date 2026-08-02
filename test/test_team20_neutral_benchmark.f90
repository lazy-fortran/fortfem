program test_team20_neutral_benchmark
    !! Independent analytic oracle for the TEAM-20-shaped neutral fixture.
    !! The vector-potential identities are written independently of the gallery
    !! sampling loop and check a nonzero 3-D field and div(curl(A)) numerically.
    use check, only: check_condition, check_summary
    use fortfem_kinds, only: dp
    implicit none

    real(dp), parameter :: step = 2.0e-5_dp
    real(dp), parameter :: tolerance = 3.0e-5_dp
    real(dp), parameter :: sample_points(3, 6) = reshape([ &
        -0.72_dp, -0.31_dp, -0.47_dp, 0.00_dp, 0.00_dp, 0.00_dp, &
         0.34_dp,  0.28_dp,  0.33_dp, 0.61_dp, -0.42_dp, 0.52_dp, &
        -0.55_dp,  0.41_dp, -0.19_dp, 0.12_dp, -0.68_dp, 0.74_dp], [3, 6])
    real(dp) :: divergence_error, energy, x, y, z, bx, by, bz
    integer :: sample

    divergence_error = 0.0_dp
    energy = 0.0_dp
    do sample = 1, size(sample_points, 2)
        x = sample_points(1, sample)
        y = sample_points(2, sample)
        z = sample_points(3, sample)
        call analytic_field(x, y, z, bx, by, bz)
        divergence_error = max(divergence_error, abs( &
            (analytic_bx(x + step, y, z) - analytic_bx(x - step, y, z))/(2.0_dp*step) + &
            (analytic_by(x, y + step, z) - analytic_by(x, y - step, z))/(2.0_dp*step) + &
            (analytic_bz(x, y, z + step) - analytic_bz(x, y, z - step))/(2.0_dp*step)))
        energy = energy + 0.5_dp*(bx**2 + by**2 + bz**2)
    end do

    call check_condition(energy > 0.0_dp, &
        "TEAM-20-shaped manufactured 3-D field has positive energy")
    call check_condition(divergence_error < tolerance, &
        "curl of vector potential is numerically divergence free")
    call analytic_field(0.0_dp, 0.0_dp, 0.0_dp, bx, by, bz)
    call check_condition(abs(bx) + abs(by) + abs(bz) > 1.0e-8_dp, &
        "central solenoid probe has a nonzero magnetic field")
    call check_condition(abs(analytic_bz(0.34_dp, 0.28_dp, 0.33_dp)) > 1.0e-8_dp, &
        "off-axis probe contains the axial field component")
    call check_summary("TEAM-20 neutral benchmark oracle")

contains

    pure subroutine analytic_field(x, y, z, bx, by, bz)
        real(dp), intent(in) :: x, y, z
        real(dp), intent(out) :: bx, by, bz
        real(dp) :: ax, ay, az, d_ax_x, d_ax_y, d_ax_z
        real(dp) :: d_ay_x, d_ay_y, d_ay_z, d_az_x, d_az_y, d_az_z
        real(dp) :: e_x, e_y, e_z, e_p, xp, pi

        pi = acos(-1.0_dp)
        xp = x + 0.52_dp
        e_x = exp(-1.8_dp*x*x - 2.1_dp*y*y)
        e_y = exp(-2.2_dp*x*x - 1.5_dp*y*y)
        e_z = exp(-3.2_dp*x*x - 3.8_dp*y*y)
        e_p = exp(-8.0_dp*xp*xp - 11.2_dp*y*y)
        ax = 0.12_dp*e_x*sin(0.5_dp*pi*z)
        ay = 0.25_dp*e_y*cos(0.5_dp*pi*z)
        az = e_z*cos(0.5_dp*pi*z) + 0.18_dp*e_p*cos(pi*z)
        d_ax_x = -3.6_dp*x*e_x*0.12_dp*sin(0.5_dp*pi*z)
        d_ax_y = -4.2_dp*y*e_x*0.12_dp*sin(0.5_dp*pi*z)
        d_ax_z = 0.12_dp*e_x*0.5_dp*pi*cos(0.5_dp*pi*z)
        d_ay_x = -4.4_dp*x*e_y*0.25_dp*cos(0.5_dp*pi*z)
        d_ay_y = -3.0_dp*y*e_y*0.25_dp*cos(0.5_dp*pi*z)
        d_ay_z = -0.25_dp*e_y*0.5_dp*pi*sin(0.5_dp*pi*z)
        d_az_x = -6.4_dp*x*e_z*cos(0.5_dp*pi*z) - &
            16.0_dp*xp*0.18_dp*e_p*cos(pi*z)
        d_az_y = -7.6_dp*y*e_z*cos(0.5_dp*pi*z) - &
            22.4_dp*y*0.18_dp*e_p*cos(pi*z)
        d_az_z = -0.5_dp*pi*e_z*sin(0.5_dp*pi*z) - &
            0.18_dp*pi*e_p*sin(pi*z)
        bx = d_az_y - d_ay_z
        by = d_ax_z - d_az_x
        bz = d_ay_x - d_ax_y
    end subroutine analytic_field

    pure real(dp) function analytic_bx(x, y, z)
        real(dp), intent(in) :: x, y, z
        real(dp) :: bx, by, bz
        call analytic_field(x, y, z, bx, by, bz)
        analytic_bx = bx
    end function analytic_bx

    pure real(dp) function analytic_by(x, y, z)
        real(dp), intent(in) :: x, y, z
        real(dp) :: bx, by, bz
        call analytic_field(x, y, z, bx, by, bz)
        analytic_by = by
    end function analytic_by

    pure real(dp) function analytic_bz(x, y, z)
        real(dp), intent(in) :: x, y, z
        real(dp) :: bx, by, bz
        call analytic_field(x, y, z, bx, by, bz)
        analytic_bz = bz
    end function analytic_bz

end program test_team20_neutral_benchmark
