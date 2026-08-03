---
title: team20_neutral_benchmark Example
---

# team20_neutral_benchmark Example

# TEAM-20-shaped neutral magnetostatic benchmark

This is a small, solution-first foundation fixture shaped like TEAM Problem
20: a three-dimensional solenoid with a central pole and an enclosing yoke.
The field is generated from a smooth manufactured vector potential, so the
curl construction supplies an analytical divergence-free solution and a fast
plotting/diagnostics path.

The fixture is deliberately **not** an exact TEAM-20 reproduction. It does
not contain a TEAM reader, the nonlinear steel B--H curve, the workshop mesh,
force measurements, or redistributed reference data. Those inputs belong in
the sister benchmark-data repository after provenance and license review. The
supplied outlines are only a visual geometry surrogate. `solution.png` is a
2-D pole/yoke cut with field vectors; `solution_3d.png` is a 3-D solenoid
surface with sparse vector segments; `probe.png` is a one-dimensional cut.
CSV and diagnostic files are generated under
`output/example/team20_neutral_benchmark`; images are not committed.

The provenance targets are the [TEAM workshop report](https://www.osti.gov/servlets/purl/7179128),
the [TEAM-20 static-force description](https://www.simscale.com/docs/validation-cases/team-20-magnetostatics/),
and the [Bíró/Ostergaard TEAM-20 proceedings reference](https://ansyshelp.ansys.com/public/Views/Secured/corp/v252/en/ans_vm/Hlp_V_VM233.html).

Exact TEAM-20 mesh, B--H, force, and probe arrays remain ineligible for
redistribution in FortFEM.  A reviewed sister-repository artifact can be
validated and rendered through the optional adapter described in
[`benchmark/external_oracles`](../../benchmark/external_oracles/README.md).

Run it with:

```text
fo run --example team20_neutral_benchmark
```

## Usage

```bash
fpm run --example team20_neutral_benchmark
```

## Source Code

```fortran
program team20_neutral_benchmark
    !! Solution-first, license-safe TEAM-20-shaped magnetostatic fixture.
    !!
    !! TEAM 20 is a 3-D static-force problem with a solenoid, steel pole, and
    !! yoke. This fixture keeps that topology as supplied geometry and uses a
    !! smooth manufactured vector potential. It is deliberately not an exact
    !! TEAM reader, nonlinear B-H model, or redistribution of TEAM data.
    use fortfem_kinds, only: dp
    use fortplot, only: add_3d_plot, add_parametric_surface, colorbar, contourf, &
        figure, legend, plot, quiver, savefig, title, xlabel, ylabel
    implicit none

    integer, parameter :: nx = 81, ny = 61
    integer, parameter :: ntheta = 48, nz_surface = 25
    integer, parameter :: arrow_nx = 15, arrow_ny = 11
    integer, parameter :: arrow_count = arrow_nx*arrow_ny
    real(dp), parameter :: pi = acos(-1.0_dp)
    character(*), parameter :: output_directory = &
        "output/example/team20_neutral_benchmark"
    real(dp) :: x_grid(nx), y_grid(ny), bmag_slice(ny, nx)
    real(dp) :: bx_slice(ny, nx), by_slice(ny, nx), bz_slice(ny, nx)
    real(dp) :: arrow_x(arrow_count), arrow_y(arrow_count)
    real(dp) :: arrow_u(arrow_count), arrow_v(arrow_count)
    real(dp) :: surface_x(nz_surface, ntheta), surface_y(nz_surface, ntheta)
    real(dp) :: surface_z(nz_surface, ntheta), surface_b(nz_surface, ntheta)
    real(dp) :: x, y, z, ax, ay, az, bx, by, bz, bmag, energy
    real(dp) :: max_divergence, step, bx_plus, bx_minus, by_plus, by_minus
    real(dp) :: bz_plus, bz_minus, scale, theta, radius
    real(dp) :: pole_x(5), pole_y(5), yoke_x(97), yoke_y(97)
    real(dp) :: coil_x(97), coil_y(97), probe_x(nx), probe_b(nx)
    integer :: ix, iy, iz, itheta, arrow, unit, command_status

    call execute_command_line("mkdir -p "//output_directory, &
        exitstat=command_status)
    if (command_status /= 0) error stop "cannot create output directory"

    do ix = 1, nx
        x_grid(ix) = -1.30_dp + 2.60_dp*real(ix - 1, dp)/real(nx - 1, dp)
    end do
    do iy = 1, ny
        y_grid(iy) = -1.00_dp + 2.00_dp*real(iy - 1, dp)/real(ny - 1, dp)
    end do

    energy = 0.0_dp
    max_divergence = 0.0_dp
    step = 2.0e-5_dp
    do iy = 1, ny
        y = y_grid(iy)
        do ix = 1, nx
            x = x_grid(ix)
            call manufactured_field(x, y, 0.0_dp, ax, ay, az, bx, by, bz)
            bx_slice(iy, ix) = bx
            by_slice(iy, ix) = by
            bz_slice(iy, ix) = bz
            bmag_slice(iy, ix) = sqrt(bx**2 + by**2 + bz**2)
            energy = energy + 0.5_dp*bmag_slice(iy, ix)**2
            if (ix > 1 .and. ix < nx .and. iy > 1 .and. iy < ny) then
                max_divergence = max(max_divergence, abs( &
                    (analytic_bx(x + step, y, 0.0_dp) - &
                     analytic_bx(x - step, y, 0.0_dp))/(2.0_dp*step) + &
                    (analytic_by(x, y + step, 0.0_dp) - &
                     analytic_by(x, y - step, 0.0_dp))/(2.0_dp*step) + &
                    (analytic_bz(x, y, step) - analytic_bz(x, y, -step))/ &
                    (2.0_dp*step)))
            end if
        end do
    end do
    energy = energy/real(nx*ny, dp)

    scale = max(maxval(bmag_slice), 1.0e-12_dp)
    arrow = 0
    do iy = 1, arrow_ny
        y = -0.80_dp + 1.60_dp*real(iy - 1, dp)/real(arrow_ny - 1, dp)
        do ix = 1, arrow_nx
            x = -1.08_dp + 2.16_dp*real(ix - 1, dp)/real(arrow_nx - 1, dp)
            arrow = arrow + 1
            call manufactured_field(x, y, 0.0_dp, ax, ay, az, bx, by, bz)
            arrow_x(arrow) = x
            arrow_y(arrow) = y
            arrow_u(arrow) = 0.14_dp*bx/scale
            arrow_v(arrow) = 0.14_dp*by/scale
        end do
    end do

    do ix = 1, nx
        probe_x(ix) = x_grid(ix)
        call manufactured_field(x_grid(ix), 0.0_dp, 0.0_dp, ax, ay, az, &
            bx, by, bz)
        probe_b(ix) = sqrt(bx**2 + by**2 + bz**2)
    end do

    ! Write the solved/manufactured field before diagnostics or plotting.
    open (newunit=unit, file=output_directory//"/solution.csv", &
        status="replace", action="write")
    write (unit, "(a)") "x,y,z,Bx,By,Bz,B_magnitude"
    do iy = 1, ny
        do ix = 1, nx
            write (unit, "(*(es24.16,:,','))") x_grid(ix), y_grid(iy), 0.0_dp, &
                bx_slice(iy, ix), by_slice(iy, ix), bz_slice(iy, ix), &
                bmag_slice(iy, ix)
        end do
    end do
    close (unit)
    open (newunit=unit, file=output_directory//"/probe.csv", &
        status="replace", action="write")
    write (unit, "(a)") "x,y,z,B_magnitude"
    do ix = 1, nx
        write (unit, "(*(es24.16,:,','))") probe_x(ix), 0.0_dp, 0.0_dp, probe_b(ix)
    end do
    close (unit)
    open (newunit=unit, file=output_directory//"/diagnostics.csv", &
        status="replace", action="write")
    write (unit, "(a,es24.16)") "mean_magnetic_energy,", energy
    write (unit, "(a,es24.16)") "divergence_proxy,", max_divergence
    write (unit, "(a,es24.16)") "max_field_magnitude,", maxval(bmag_slice)
    close (unit)

    ! Supplied TEAM-20-shaped solenoid/pole/yoke outlines on the solution cut.
    pole_x = [-0.24_dp, 0.24_dp, 0.24_dp, -0.24_dp, -0.24_dp]
    pole_y = [-0.72_dp, -0.72_dp, 0.72_dp, 0.72_dp, -0.72_dp]
    do ix = 1, size(yoke_x)
        theta = 2.0_dp*pi*real(ix - 1, dp)/real(size(yoke_x) - 1, dp)
        yoke_x(ix) = 0.66_dp*cos(theta)
        yoke_y(ix) = 0.88_dp*sin(theta)
        coil_x(ix) = 0.40_dp + 0.19_dp*cos(theta)
        coil_y(ix) = 0.00_dp + 0.70_dp*sin(theta)
    end do

    ! Primary physical view: field magnitude, vectors, and benchmark layout.
    call figure(figsize=[8.6_dp, 6.4_dp])
    call contourf(x_grid, y_grid, bmag_slice, cmap="viridis", show_colorbar=.true.)
    call colorbar(label="|B| (manufactured static-force field)")
    call quiver(arrow_x, arrow_y, arrow_u, arrow_v, scale=1.0_dp, &
        scale_units="xy", angles="xy", color=[0.0_dp, 0.0_dp, 0.0_dp], &
        width=0.004_dp, headwidth=4.0_dp)
    do arrow = 1, arrow_count
        call plot([arrow_x(arrow), arrow_x(arrow) + arrow_u(arrow)], &
            [arrow_y(arrow), arrow_y(arrow) + arrow_v(arrow)], &
            color=[0.0_dp, 0.0_dp, 0.0_dp], linewidth=1.0_dp)
    end do
    call plot(pole_x, pole_y, color=[0.95_dp, 0.95_dp, 0.95_dp], linewidth=2.5_dp)
    call plot(yoke_x, yoke_y, color=[1.0_dp, 0.80_dp, 0.10_dp], linewidth=2.5_dp)
    call plot(coil_x, coil_y, color=[1.0_dp, 0.25_dp, 0.10_dp], linewidth=2.5_dp)
    call xlabel("x (pole/yoke cut)")
    call ylabel("y")
    call title("TEAM-20-shaped neutral solenoid field")
    call savefig(output_directory//"/solution.png")

    ! Three-dimensional physical surface: a cylindrical coil envelope lifted by
    ! sampled field magnitude, with sparse 3-D vector segments.
    radius = 0.74_dp
    do iz = 1, nz_surface
        z = -0.92_dp + 1.84_dp*real(iz - 1, dp)/real(nz_surface - 1, dp)
        do itheta = 1, ntheta
            theta = 2.0_dp*pi*real(itheta - 1, dp)/real(ntheta - 1, dp)
            x = radius*cos(theta)
            y = radius*sin(theta)
            call manufactured_field(x, y, z, ax, ay, az, bx, by, bz)
            surface_x(iz, itheta) = x
            surface_y(iz, itheta) = y
            surface_b(iz, itheta) = sqrt(bx**2 + by**2 + bz**2)
            surface_z(iz, itheta) = z + 0.22_dp*surface_b(iz, itheta)/scale
        end do
    end do
    call figure(figsize=[8.4_dp, 6.4_dp])
    call add_parametric_surface(surface_x, surface_y, surface_z, &
        color="royalblue", alpha=0.72_dp, linewidth=0.0_dp, filled=.true., &
        row_stride=2, column_stride=2, label="|B| height on solenoid")
    do iz = 2, nz_surface - 1, 4
        do itheta = 1, ntheta - 1, 4
            theta = 2.0_dp*pi*real(itheta - 1, dp)/real(ntheta - 1, dp)
            x = radius*cos(theta)
            y = radius*sin(theta)
            z = -0.92_dp + 1.84_dp*real(iz - 1, dp)/real(nz_surface - 1, dp)
            call manufactured_field(x, y, z, ax, ay, az, bx, by, bz)
            call add_3d_plot([x, x + 0.13_dp*bx/scale], &
                [y, y + 0.13_dp*by/scale], &
                [z + 0.22_dp*sqrt(bx**2 + by**2 + bz**2)/scale, &
                 z + 0.22_dp*sqrt(bx**2 + by**2 + bz**2)/scale + 0.13_dp*bz/scale], &
                color="black", linewidth=1.0_dp)
        end do
    end do
    call title("TEAM-20-shaped neutral 3-D field surface")
    call legend()
    call savefig(output_directory//"/solution_3d.png")

    call figure(figsize=[8.0_dp, 4.8_dp])
    call plot(probe_x, probe_b, color=[0.12_dp, 0.35_dp, 0.75_dp], linewidth=2.0_dp)
    call xlabel("x at y=z=0")
    call ylabel("|B|")
    call title("TEAM-20 neutral axial probe response")
    call savefig(output_directory//"/probe.png")

contains

    pure subroutine manufactured_field(x, y, z, ax, ay, az, bx, by, bz)
        real(dp), intent(in) :: x, y, z
        real(dp), intent(out) :: ax, ay, az, bx, by, bz
        real(dp) :: e_x, e_y, e_z, e_p, d_ax_x, d_ax_y, d_ax_z
        real(dp) :: d_ay_x, d_ay_y, d_ay_z, d_az_x, d_az_y, d_az_z
        real(dp) :: xp

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
    end subroutine manufactured_field

    pure real(dp) function analytic_bx(x, y, z)
        real(dp), intent(in) :: x, y, z
        real(dp) :: ax, ay, az, by, bz
        call manufactured_field(x, y, z, ax, ay, az, analytic_bx, by, bz)
    end function analytic_bx

    pure real(dp) function analytic_by(x, y, z)
        real(dp), intent(in) :: x, y, z
        real(dp) :: ax, ay, az, bx, bz
        call manufactured_field(x, y, z, ax, ay, az, bx, analytic_by, bz)
    end function analytic_by

    pure real(dp) function analytic_bz(x, y, z)
        real(dp), intent(in) :: x, y, z
        real(dp) :: ax, ay, az, bx, by
        call manufactured_field(x, y, z, ax, ay, az, bx, by, analytic_bz)
    end function analytic_bz

end program team20_neutral_benchmark
```

## Generated Plots

### primary.png

![primary.png](../../media/examples/team20_neutral_benchmark/primary.png)

### probe.png

![probe.png](../../media/examples/team20_neutral_benchmark/probe.png)

### solution.png

![solution.png](../../media/examples/team20_neutral_benchmark/solution.png)

### solution_3d.png

![solution_3d.png](../../media/examples/team20_neutral_benchmark/solution_3d.png)

---

[← Back to all examples](../index.html)
