program team13_neutral_benchmark
    !! Solution-first, license-safe TEAM-13-shaped electromagnetic fixture.
    !!
    !! The geometry is a small supplied-array surrogate for the coil, channels,
    !! and plate layout.  The stream function is manufactured, so the example
    !! tests field sampling, probes, energy, and plotting without embedding a
    !! TEAM reader, nonlinear B-H curve, or restricted reference data.
    use fortfem_kinds, only: dp
    use fortplot, only: colorbar, contourf, figure, plot, quiver, savefig, &
        title, xlabel, ylabel
    implicit none

    integer, parameter :: nx = 61, ny = 45
    integer, parameter :: arrow_nx = 16, arrow_ny = 12
    integer, parameter :: arrow_count = arrow_nx*arrow_ny
    real(dp), parameter :: pi = acos(-1.0_dp)
    character(*), parameter :: output_directory = &
        "output/example/team13_neutral_benchmark"
    ! FortPlot stores raster fields as (ny,nx); keeping that convention here
    ! also makes the example work with older released FortPlot versions.
    real(dp) :: x_grid(nx), y_grid(ny), bmag(ny, nx)
    real(dp) :: arrow_x(arrow_count), arrow_y(arrow_count)
    real(dp) :: arrow_u(arrow_count), arrow_v(arrow_count)
    real(dp) :: bx, by, x, y, psi
    real(dp) :: bx_plus, bx_minus, by_plus, by_minus
    real(dp) :: bx_yplus, bx_yminus, by_yplus, by_yminus
    real(dp) :: coil, channel, plate, energy, divergence_error
    real(dp) :: probe_x(nx), probe_b(nx)
    integer :: ix, iy, arrow, unit

    call execute_command_line("mkdir -p "//output_directory)
    divergence_error = 0.0_dp
    energy = 0.0_dp
    do ix = 1, nx
        x_grid(ix) = -1.0_dp + 2.0_dp*real(ix - 1, dp)/real(nx - 1, dp)
    end do
    do iy = 1, ny
        y_grid(iy) = -0.8_dp + 1.6_dp*real(iy - 1, dp)/real(ny - 1, dp)
    end do

    do iy = 1, ny
        y = y_grid(iy)
        do ix = 1, nx
            x = x_grid(ix)
            call manufactured_field(x, y, psi, bx, by, coil, channel, plate)
            bmag(iy, ix) = sqrt(bx**2 + by**2)
            energy = energy + 0.5_dp*bmag(iy, ix)**2
            ! The channel uses |x| and |y| to suggest a benchmark slot.  Its
            ! manufactured field is smooth away from those coordinate axes,
            ! so exclude the kink lines from this finite-difference diagnostic.
            if (ix > 1 .and. ix < nx .and. iy > 1 .and. iy < ny .and. &
                abs(x) > 2.0e-5_dp .and. abs(y) > 2.0e-5_dp) then
                call manufactured_field(x + 1.0e-5_dp, y, psi, bx_plus, &
                    by_plus, coil, channel, plate)
                call manufactured_field(x - 1.0e-5_dp, y, psi, bx_minus, &
                    by_minus, coil, channel, plate)
                call manufactured_field(x, y + 1.0e-5_dp, psi, bx_yplus, &
                    by_yplus, coil, channel, plate)
                call manufactured_field(x, y - 1.0e-5_dp, psi, bx_yminus, &
                    by_yminus, coil, channel, plate)
                divergence_error = max(divergence_error, abs( &
                    (bx_plus - bx_minus)/(2.0e-5_dp) + &
                    (by_yplus - by_yminus)/(2.0e-5_dp)))
            end if
        end do
    end do
    energy = energy/real(nx*ny, dp)

    arrow = 0
    do iy = 1, arrow_ny
        y = -0.75_dp + 1.5_dp*real(iy - 1, dp)/real(arrow_ny - 1, dp)
        do ix = 1, arrow_nx
            x = -0.9_dp + 1.8_dp*real(ix - 1, dp)/real(arrow_nx - 1, dp)
            arrow = arrow + 1
            call manufactured_field(x, y, psi, bx, by, coil, channel, plate)
            arrow_x(arrow) = x
            arrow_y(arrow) = y
            arrow_u(arrow) = 0.16_dp*bx/max(1.0e-12_dp, sqrt(bx**2 + by**2))
            arrow_v(arrow) = 0.16_dp*by/max(1.0e-12_dp, sqrt(bx**2 + by**2))
        end do
    end do

    do ix = 1, nx
        x = x_grid(ix)
        probe_x(ix) = x
        call manufactured_field(x, 0.30_dp, psi, bx, by, coil, channel, plate)
        probe_b(ix) = sqrt(bx**2 + by**2)
    end do
    open (newunit=unit, file=output_directory//"/solution.csv", &
        status="replace", action="write")
    write (unit, "(a)") "x,y,Bx,By,B_magnitude"
    do iy = 1, ny
        do ix = 1, nx
            x = x_grid(ix)
            y = y_grid(iy)
            call manufactured_field(x, y, psi, bx, by, coil, channel, plate)
            write (unit, "(*(es24.16,:,','))") x, y, bx, by, sqrt(bx**2 + by**2)
        end do
    end do
    close (unit)
    open (newunit=unit, file=output_directory//"/probe.csv", &
        status="replace", action="write")
    write (unit, "(a)") "x,y,field_magnitude"
    do ix = 1, nx
        write (unit, "(*(es24.16,:,','))") probe_x(ix), 0.30_dp, probe_b(ix)
    end do
    close (unit)
    open (newunit=unit, file=output_directory//"/diagnostics.csv", &
        status="replace", action="write")
    write (unit, "(a,es24.16)") "mean_magnetic_energy,", energy
    write (unit, "(a,es24.16)") "divergence_proxy,", divergence_error
    close (unit)

    ! Physical solution first: field magnitude with vectors and benchmark shape.
    call figure(figsize=[8.0_dp, 6.4_dp])
    call contourf(x_grid, y_grid, bmag, cmap="viridis", show_colorbar=.true.)
    call colorbar(label="|B| (supplied-array units)")
    call quiver(arrow_x, arrow_y, arrow_u, arrow_v, scale=1.0_dp, &
        scale_units="xy", angles="xy", color=[0.0_dp, 0.0_dp, 0.0_dp], &
        width=0.004_dp, headwidth=4.0_dp)
    ! Keep a line-segment fallback so vector directions remain visible in
    ! minimal backends that do not render quiver heads.
    do arrow = 1, arrow_count
        call plot([arrow_x(arrow), arrow_x(arrow) + arrow_u(arrow)], &
            [arrow_y(arrow), arrow_y(arrow) + arrow_v(arrow)], &
            color=[0.0_dp, 0.0_dp, 0.0_dp], linewidth=1.0_dp)
    end do
    call plot([-0.85_dp, -0.85_dp, 0.85_dp, 0.85_dp, -0.85_dp], &
        [-0.55_dp, -0.25_dp, -0.25_dp, -0.55_dp, -0.55_dp], &
        color=[1.0_dp, 0.3_dp, 0.1_dp], linewidth=2.0_dp)
    call plot([-0.75_dp, -0.75_dp, 0.75_dp, 0.75_dp, -0.75_dp], &
        [0.22_dp, 0.42_dp, 0.42_dp, 0.22_dp, 0.22_dp], &
        color=[1.0_dp, 0.8_dp, 0.1_dp], linewidth=2.0_dp)
    call xlabel("x")
    call ylabel("y")
    call title("TEAM-13-shaped neutral magnetostatic field")
    call savefig(output_directory//"/solution.png")

    call figure(figsize=[7.0_dp, 4.5_dp])
    call plot(probe_x, probe_b, color=[0.12_dp, 0.35_dp, 0.75_dp], &
        linewidth=2.0_dp)
    call xlabel("probe x at y=0.30")
    call ylabel("|B|")
    call title("TEAM neutral probe response")
    call savefig(output_directory//"/probe.png")

contains

    subroutine manufactured_field(x, y, psi, bx, by, coil, channel, plate)
        real(dp), intent(in) :: x, y
        real(dp), intent(out) :: psi, bx, by, coil, channel, plate
        real(dp) :: coil_x, coil_y, channel_x, channel_y, plate_x, plate_y
        real(dp) :: radius, gaussian

        coil_x = x + 0.38_dp
        coil_y = y
        radius = coil_x**2 + coil_y**2
        coil = exp(-18.0_dp*radius)
        channel_x = abs(x) - 0.73_dp
        channel_y = abs(y) - 0.32_dp
        channel = exp(-55.0_dp*(channel_x**2 + channel_y**2))
        plate_x = x - 0.12_dp
        plate_y = y - 0.02_dp
        plate = exp(-18.0_dp*plate_x**2 - 5.0_dp*plate_y**2)
        gaussian = exp(-2.0_dp*(x**2 + y**2))
        psi = coil + 0.25_dp*channel + 0.18_dp*plate + 0.06_dp*gaussian
        bx = -( -36.0_dp*coil_y*coil + &
            0.25_dp*(-110.0_dp*sign(1.0_dp, y)*channel_y*channel) + &
            0.18_dp*(-10.0_dp*plate_y*plate) - 0.24_dp*y*gaussian)
        by = -36.0_dp*coil_x*coil + &
            0.25_dp*(-110.0_dp*sign(1.0_dp, x)*channel_x*channel) + &
            0.18_dp*(-36.0_dp*plate_x*plate) - 0.24_dp*x*gaussian
    end subroutine manufactured_field

end program team13_neutral_benchmark
