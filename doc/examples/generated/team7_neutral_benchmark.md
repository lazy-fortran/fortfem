---
title: team7_neutral_benchmark Example
---

# team7_neutral_benchmark Example

# TEAM-7-shaped neutral eddy-current benchmark

This solution-first fixture is a small, license-safe 2-D analogue of the
coil/plate eddy-current layout used in the TEAM benchmark family.  A complex
manufactured stream function defines a divergence-free magnetic field,
`B = (dA/dy,-dA/dx)`, and the corresponding complex Ampere current
`J_z = -laplacian(A)`.  The first plot, `solution.png`, shows the field
magnitude, vector directions, coil outline, conducting plate, and slot.  The
additional `current.png` and `probe.png` plots expose the induced-current
response and a probe cut.  `solution.csv`, `probe.csv`, and `diagnostics.csv`
are generated beside the plots.

The geometry, source amplitudes, and material response are intentionally
manufactured.  This example does **not** embed a TEAM reader, nonlinear B-H
curve, proprietary geometry, or external reference data; an exact benchmark
comparison belongs in the separate benchmark-data repository when its source
terms and license permit redistribution.

The provenance target is the [TEAM workshop report](https://www.osti.gov/servlets/purl/7179128)
and the public [TEAM problem catalogue](https://www.osti.gov/biblio/7179128).  The
fixture is intended to exercise FortFEM's 2-D vector-field plotting, analytic
manufactured forcing, and diagnostic hooks before a specialized solver is
added.

Exact TEAM-7 arrays are intentionally out of tree.  If a separately licensed
sister-repository artifact is available, use the provenance/checksum-gated
adapter documented in [`benchmark/external_oracles`](../../benchmark/external_oracles/README.md).

Run it with:

```text
fo run --example team7_neutral_benchmark
```

## Usage

```bash
fpm run --example team7_neutral_benchmark
```

## Source Code

```fortran
program team7_neutral_benchmark
    !! Solution-first, license-safe TEAM-7-shaped eddy-current fixture.
    !!
    !! A complex stream function manufactures a divergence-free magnetic field
    !! and its Ampere current on a small 2-D plate/coil layout.  The geometry,
    !! source amplitudes, and material response are deliberately supplied
    !! arrays: this is a foundation fixture, not a TEAM reader or a claim to
    !! reproduce restricted benchmark data.
    use fortfem_kinds, only: dp
    use fortplot, only: colorbar, contourf, figure, legend, plot, savefig, &
        title, xlabel, ylabel
    implicit none

    integer, parameter :: nx = 81, ny = 61
    integer, parameter :: arrow_nx = 17, arrow_ny = 13
    integer, parameter :: arrow_count = arrow_nx*arrow_ny
    integer, parameter :: outline_count = 65
    real(dp), parameter :: pi = acos(-1.0_dp)
    character(*), parameter :: output_directory = &
        "output/example/team7_neutral_benchmark"
    real(dp) :: x_grid(nx), y_grid(ny)
    real(dp) :: bmag(ny, nx), jmag(ny, nx), bx_field(ny, nx)
    real(dp) :: by_field(ny, nx), bx_im_field(ny, nx), by_im_field(ny, nx)
    real(dp) :: jz_field(ny, nx)
    real(dp) :: arrow_x(arrow_count), arrow_y(arrow_count)
    real(dp) :: arrow_u(arrow_count), arrow_v(arrow_count)
    real(dp) :: coil_x(outline_count), coil_y(outline_count)
    real(dp) :: plate_x(5), plate_y(5), slot_x(5), slot_y(5)
    real(dp) :: line_x(2), line_y(2)
    real(dp) :: probe_x(nx), probe_b(nx), probe_j(nx)
    real(dp) :: x, y, a_re, a_im, bx_re, by_re, bx_im, by_im
    real(dp) :: jz_re, jz_im, energy, current_energy, max_divergence
    real(dp), parameter :: divergence_step = 2.0e-5_dp
    real(dp) :: bx_plus, bx_minus, by_plus, by_minus
    real(dp) :: bx_yplus, bx_yminus, by_yplus, by_yminus
    real(dp) :: scale, angle
    integer :: ix, iy, arrow, unit

    call execute_command_line("mkdir -p "//output_directory)
    do ix = 1, nx
        x_grid(ix) = -1.2_dp + 2.4_dp*real(ix - 1, dp)/real(nx - 1, dp)
    end do
    do iy = 1, ny
        y_grid(iy) = -0.95_dp + 1.9_dp*real(iy - 1, dp)/real(ny - 1, dp)
    end do

    energy = 0.0_dp
    current_energy = 0.0_dp
    max_divergence = 0.0_dp
    do iy = 1, ny
        y = y_grid(iy)
        do ix = 1, nx
            x = x_grid(ix)
            call manufactured_field(x, y, a_re, a_im, bx_re, by_re, bx_im, &
                by_im, jz_re, jz_im)
            bx_field(iy, ix) = bx_re
            by_field(iy, ix) = by_re
            bx_im_field(iy, ix) = bx_im
            by_im_field(iy, ix) = by_im
            jz_field(iy, ix) = jz_re
            bmag(iy, ix) = sqrt(bx_re**2 + by_re**2 + bx_im**2 + by_im**2)
            jmag(iy, ix) = sqrt(jz_re**2 + jz_im**2)
            energy = energy + 0.5_dp*bmag(iy, ix)**2
            current_energy = current_energy + 0.5_dp*jmag(iy, ix)**2
        end do
    end do
    energy = energy/real(nx*ny, dp)
    current_energy = current_energy/real(nx*ny, dp)

    ! A finite-difference divergence check is a cheap run-time diagnostic;
    ! the manufactured stream-function construction is exact analytically.
    do iy = 2, ny - 1
        y = y_grid(iy)
        do ix = 2, nx - 1
            x = x_grid(ix)
            call manufactured_field(x + divergence_step, y, a_re, a_im, &
                bx_plus, by_plus, bx_im, by_im, jz_re, jz_im)
            call manufactured_field(x - divergence_step, y, a_re, a_im, &
                bx_minus, by_minus, bx_im, by_im, jz_re, jz_im)
            call manufactured_field(x, y + divergence_step, a_re, a_im, &
                bx_yplus, by_yplus, bx_im, by_im, jz_re, jz_im)
            call manufactured_field(x, y - divergence_step, a_re, a_im, &
                bx_yminus, by_yminus, bx_im, by_im, jz_re, jz_im)
            max_divergence = max(max_divergence, abs( &
                (bx_plus - bx_minus)/(2.0_dp*divergence_step) + &
                (by_yplus - by_yminus)/(2.0_dp*divergence_step)))
        end do
    end do

    arrow = 0
    scale = max(maxval(sqrt(bx_field**2 + by_field**2)), 1.0e-12_dp)
    do iy = 1, arrow_ny
        y = -0.82_dp + 1.64_dp*real(iy - 1, dp)/real(arrow_ny - 1, dp)
        do ix = 1, arrow_nx
            x = -1.05_dp + 2.10_dp*real(ix - 1, dp)/real(arrow_nx - 1, dp)
            arrow = arrow + 1
            call manufactured_field(x, y, a_re, a_im, bx_re, by_re, bx_im, &
                by_im, jz_re, jz_im)
            arrow_x(arrow) = x
            arrow_y(arrow) = y
            arrow_u(arrow) = 0.12_dp*bx_re/scale
            arrow_v(arrow) = 0.12_dp*by_re/scale
        end do
    end do

    do ix = 1, nx
        x = x_grid(ix)
        probe_x(ix) = x
        call manufactured_field(x, -0.23_dp, a_re, a_im, bx_re, by_re, &
            bx_im, by_im, jz_re, jz_im)
        probe_b(ix) = sqrt(bx_re**2 + by_re**2 + bx_im**2 + by_im**2)
        probe_j(ix) = sqrt(jz_re**2 + jz_im**2)
    end do

    open (newunit=unit, file=output_directory//"/solution.csv", &
        status="replace", action="write")
    write (unit, "(a)") "x,y,Bx_re,By_re,Bx_im,By_im,B_magnitude,Jz_magnitude"
    do iy = 1, ny
        do ix = 1, nx
            write (unit, "(*(es24.16,:,','))") x_grid(ix), y_grid(iy), &
                bx_field(iy, ix), by_field(iy, ix), bx_im_field(iy, ix), &
                by_im_field(iy, ix), &
                bmag(iy, ix), jmag(iy, ix)
        end do
    end do
    close (unit)
    open (newunit=unit, file=output_directory//"/probe.csv", &
        status="replace", action="write")
    write (unit, "(a)") "x,y,B_magnitude,Jz_magnitude"
    do ix = 1, nx
        write (unit, "(*(es24.16,:,','))") probe_x(ix), -0.23_dp, &
            probe_b(ix), probe_j(ix)
    end do
    close (unit)
    open (newunit=unit, file=output_directory//"/diagnostics.csv", &
        status="replace", action="write")
    write (unit, "(a,es24.16)") "mean_magnetic_energy,", energy
    write (unit, "(a,es24.16)") "mean_current_energy,", current_energy
    write (unit, "(a,es24.16)") "divergence_proxy,", max_divergence
    close (unit)

    do ix = 1, outline_count
        angle = 2.0_dp*pi*real(ix - 1, dp)/real(outline_count - 1, dp)
        coil_x(ix) = -0.42_dp + 0.24_dp*cos(angle)
        coil_y(ix) = 0.56_dp + 0.12_dp*sin(angle)
    end do
    plate_x = [-1.0_dp, 1.0_dp, 1.0_dp, -1.0_dp, -1.0_dp]
    plate_y = [-0.66_dp, -0.66_dp, -0.04_dp, -0.04_dp, -0.66_dp]
    slot_x = [-0.20_dp, 0.20_dp, 0.20_dp, -0.20_dp, -0.20_dp]
    slot_y = [-0.66_dp, -0.66_dp, -0.36_dp, -0.36_dp, -0.66_dp]

    ! Physical solution first: field magnitude, vectors, and TEAM-like layout.
    call figure(figsize=[8.6_dp, 6.3_dp])
    call contourf(x_grid, y_grid, bmag, cmap="viridis", show_colorbar=.true.)
    call colorbar(label="|B| (manufactured phasor amplitude)")
    do arrow = 1, arrow_count
        line_x = [arrow_x(arrow), arrow_x(arrow) + arrow_u(arrow)]
        line_y = [arrow_y(arrow), arrow_y(arrow) + arrow_v(arrow)]
        call plot(line_x, line_y, &
            color=[0.0_dp, 0.0_dp, 0.0_dp], linewidth=1.0_dp)
    end do
    call plot(coil_x, coil_y, color=[1.0_dp, 0.3_dp, 0.1_dp], linewidth=2.5_dp)
    call plot(plate_x, plate_y, color=[1.0_dp, 0.85_dp, 0.1_dp], linewidth=2.5_dp)
    call plot(slot_x, slot_y, color=[0.95_dp, 0.95_dp, 0.95_dp], linewidth=2.5_dp)
    call xlabel("x")
    call ylabel("y")
    call title("TEAM-7-shaped neutral eddy-current field")
    call savefig(output_directory//"/solution.png")

    call figure(figsize=[8.0_dp, 5.8_dp])
    call contourf(x_grid, y_grid, jmag, cmap="magma", show_colorbar=.true.)
    call colorbar(label="|J_z| (manufactured Ampere current)")
    call plot(coil_x, coil_y, color=[1.0_dp, 0.9_dp, 0.2_dp], linewidth=2.0_dp)
    call plot(plate_x, plate_y, color=[1.0_dp, 0.9_dp, 0.2_dp], linewidth=2.0_dp)
    call xlabel("x")
    call ylabel("y")
    call title("TEAM-7-shaped induced-current magnitude")
    call savefig(output_directory//"/current.png")

    call figure(figsize=[8.0_dp, 4.8_dp])
    call plot(probe_x, probe_b, label="|B|", color=[0.12_dp, 0.35_dp, 0.75_dp], &
        linewidth=2.0_dp)
    call plot(probe_x, probe_j, label="|J_z|", color=[0.8_dp, 0.25_dp, 0.1_dp], &
        linewidth=2.0_dp)
    call xlabel("x at y=-0.23")
    call ylabel("manufactured amplitude")
    call title("TEAM neutral probe response")
    call legend()
    call savefig(output_directory//"/probe.png")

contains

    pure subroutine manufactured_field(x, y, a_re, a_im, bx_re, by_re, bx_im, &
            by_im, jz_re, jz_im)
        real(dp), intent(in) :: x, y
        real(dp), intent(out) :: a_re, a_im, bx_re, by_re, bx_im, by_im
        real(dp), intent(out) :: jz_re, jz_im
        real(dp) :: value, dx_value, dy_value, dxx_value, dyy_value
        real(dp) :: ax_re, ay_re, axx_re, ayy_re
        real(dp) :: ax_im, ay_im, axx_im, ayy_im

        call component(x, y, -0.42_dp, 0.56_dp, 18.0_dp, 18.0_dp, 0.0_dp, &
            a_re, ax_re, ay_re, axx_re, ayy_re)
        call component(x, y, 0.0_dp, -0.34_dp, 1.3_dp, 5.0_dp, 2.0_dp*pi/1.4_dp, &
            value, dx_value, dy_value, dxx_value, dyy_value)
        a_re = a_re + 0.42_dp*value
        ax_re = ax_re + 0.42_dp*dx_value
        ay_re = ay_re + 0.42_dp*dy_value
        axx_re = axx_re + 0.42_dp*dxx_value
        ayy_re = ayy_re + 0.42_dp*dyy_value
        a_im = 0.27_dp*value
        ax_im = 0.27_dp*dx_value
        ay_im = 0.27_dp*dy_value
        axx_im = 0.27_dp*dxx_value
        ayy_im = 0.27_dp*dyy_value
        bx_re = ay_re
        by_re = -ax_re
        bx_im = ay_im
        by_im = -ax_im
        jz_re = -(axx_re + ayy_re)
        jz_im = -(axx_im + ayy_im)
    end subroutine manufactured_field

    pure subroutine component(x, y, center_x, center_y, alpha, beta, wave, &
            value, dx_value, dy_value, dxx_value, dyy_value)
        real(dp), intent(in) :: x, y, center_x, center_y, alpha, beta, wave
        real(dp), intent(out) :: value, dx_value, dy_value, dxx_value, dyy_value
        real(dp) :: xx, yy, envelope, cosine, sine

        xx = x - center_x
        yy = y - center_y
        envelope = exp(-alpha*xx**2 - beta*yy**2)
        cosine = cos(wave*xx)
        sine = sin(wave*xx)
        value = envelope*cosine
        dx_value = envelope*(-2.0_dp*alpha*xx*cosine - wave*sine)
        dy_value = envelope*(-2.0_dp*beta*yy*cosine)
        dxx_value = envelope*((4.0_dp*alpha**2*xx**2 - 2.0_dp*alpha - wave**2)* &
            cosine + 4.0_dp*alpha*wave*xx*sine)
        dyy_value = envelope*(4.0_dp*beta**2*yy**2 - 2.0_dp*beta)*cosine
    end subroutine component

end program team7_neutral_benchmark
```

## Generated Plots

### primary.png

![primary.png](../../media/examples/team7_neutral_benchmark/primary.png)

### current.png

![current.png](../../media/examples/team7_neutral_benchmark/current.png)

### probe.png

![probe.png](../../media/examples/team7_neutral_benchmark/probe.png)

### solution.png

![solution.png](../../media/examples/team7_neutral_benchmark/solution.png)

---

[← Back to all examples](../index.html)
