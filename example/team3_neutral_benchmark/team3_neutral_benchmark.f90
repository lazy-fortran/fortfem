program team3_neutral_benchmark
    !! Solution-first, license-safe TEAM-3-shaped neutral fixture.
    !!
    !! A smooth manufactured A_z gives B=(dA/dy,-dA/dx) and
    !! J_z=-laplacian(A_z).  The arrays and geometry are deliberately supplied
    !! surrogates, not TEAM input data or a nonlinear material model.
    use fortfem_kinds, only: dp
    use fortplot, only: add_3d_plot, add_parametric_surface, colorbar, contourf, &
        figure, legend, plot, quiver, savefig, title, xlabel, ylabel
    implicit none

    integer, parameter :: nx = 81, ny = 61
    integer, parameter :: arrow_nx = 17, arrow_ny = 13
    integer, parameter :: arrow_count = arrow_nx*arrow_ny
    integer, parameter :: outline_count = 65
    integer, parameter :: surface_stride = 2
    real(dp), parameter :: pi = acos(-1.0_dp)
    character(*), parameter :: output_directory = &
        "output/example/team3_neutral_benchmark"
    real(dp) :: x_grid(nx), y_grid(ny)
    real(dp) :: az_field(ny, nx), bmag(ny, nx), bx_field(ny, nx), by_field(ny, nx)
    real(dp) :: jz_field(ny, nx), mu_r_field(ny, nx), source_field(ny, nx)
    real(dp) :: arrow_x(arrow_count), arrow_y(arrow_count)
    real(dp) :: arrow_u(arrow_count), arrow_v(arrow_count)
    real(dp) :: core_x(9), core_y(9), coil_x(outline_count), coil_y(outline_count)
    real(dp) :: gap_x(2), gap_y(2), probe_x(nx), probe_b(nx), probe_j(nx)
    real(dp) :: x, y, az, ax, ay, axx, ayy, bx, by, jz
    real(dp) :: energy, current_energy, max_divergence, source_energy
    real(dp) :: bx_plus, bx_minus, by_plus, by_minus, step, scale, angle
    real(dp) :: surface_x(ny, nx), surface_y(ny, nx), surface_z(ny, nx)
    real(dp) :: field_scale, arrow_length, arrow_end(3)
    integer :: ix, iy, arrow, unit, command_status

    call execute_command_line("mkdir -p "//output_directory, &
        exitstat=command_status)
    if (command_status /= 0) error stop "cannot create output directory"
    call initialize_gallery_sequence()

    do ix = 1, nx
        x_grid(ix) = -1.30_dp + 2.60_dp*real(ix - 1, dp)/real(nx - 1, dp)
    end do
    do iy = 1, ny
        y_grid(iy) = -1.00_dp + 2.00_dp*real(iy - 1, dp)/real(ny - 1, dp)
    end do

    energy = 0.0_dp
    current_energy = 0.0_dp
    source_energy = 0.0_dp
    max_divergence = 0.0_dp
    do iy = 1, ny
        y = y_grid(iy)
        do ix = 1, nx
            x = x_grid(ix)
            call manufactured_field(x, y, az, ax, ay, axx, ayy, bx, by, jz)
            az_field(iy, ix) = az
            bx_field(iy, ix) = bx
            by_field(iy, ix) = by
            jz_field(iy, ix) = jz
            bmag(iy, ix) = sqrt(bx**2 + by**2)
            mu_r_field(iy, ix) = material_response(x, y)
            source_field(iy, ix) = source_response(x, y)
            energy = energy + 0.5_dp*mu_r_field(iy, ix)*bmag(iy, ix)**2
            current_energy = current_energy + 0.5_dp*jz**2
            source_energy = source_energy + source_field(iy, ix)**2
            if (ix > 1 .and. ix < nx .and. iy > 1 .and. iy < ny) then
                step = 2.0e-5_dp
                call manufactured_field(x + step, y, az, ax, ay, axx, ayy, &
                    bx_plus, by_plus, jz)
                call manufactured_field(x - step, y, az, ax, ay, axx, ayy, &
                    bx_minus, by_minus, jz)
                call manufactured_field(x, y + step, az, ax, ay, axx, ayy, &
                    bx, by_plus, jz)
                call manufactured_field(x, y - step, az, ax, ay, axx, ayy, &
                    bx, by_minus, jz)
                max_divergence = max(max_divergence, abs( &
                    (bx_plus - bx_minus)/(2.0_dp*step) + &
                    (by_plus - by_minus)/(2.0_dp*step)))
            end if
        end do
    end do
    energy = energy/real(nx*ny, dp)
    current_energy = current_energy/real(nx*ny, dp)
    source_energy = source_energy/real(nx*ny, dp)

    arrow = 0
    scale = max(maxval(bmag), 1.0e-12_dp)
    do iy = 1, arrow_ny
        y = -0.82_dp + 1.64_dp*real(iy - 1, dp)/real(arrow_ny - 1, dp)
        do ix = 1, arrow_nx
            x = -1.13_dp + 2.26_dp*real(ix - 1, dp)/real(arrow_nx - 1, dp)
            arrow = arrow + 1
            call manufactured_field(x, y, az, ax, ay, axx, ayy, bx, by, jz)
            arrow_x(arrow) = x
            arrow_y(arrow) = y
            arrow_u(arrow) = 0.14_dp*bx/scale
            arrow_v(arrow) = 0.14_dp*by/scale
        end do
    end do

    do ix = 1, nx
        x = x_grid(ix)
        probe_x(ix) = x
        call manufactured_field(x, 0.0_dp, az, ax, ay, axx, ayy, bx, by, jz)
        probe_b(ix) = sqrt(bx**2 + by**2)
        probe_j(ix) = abs(jz)
    end do

    ! A solution CSV is intentionally written before any plots or diagnostics.
    open (newunit=unit, file=output_directory//"/solution.csv", &
        status="replace", action="write")
    write (unit, "(a)") "x,y,Az,Bx,By,B_magnitude,Jz"
    do iy = 1, ny
        do ix = 1, nx
            write (unit, "(*(es24.16,:,','))") x_grid(ix), y_grid(iy), &
                az_field(iy, ix), bx_field(iy, ix), by_field(iy, ix), bmag(iy, ix), &
                jz_field(iy, ix)
        end do
    end do
    close (unit)
    open (newunit=unit, file=output_directory//"/material.csv", &
        status="replace", action="write")
    write (unit, "(a)") "x,y,mu_r,source"
    do iy = 1, ny
        do ix = 1, nx
            write (unit, "(*(es24.16,:,','))") x_grid(ix), y_grid(iy), &
                mu_r_field(iy, ix), source_field(iy, ix)
        end do
    end do
    close (unit)

    core_x = [-0.90_dp, 0.55_dp, 0.55_dp, -0.45_dp, -0.45_dp, 0.55_dp, &
        0.55_dp, -0.90_dp, -0.90_dp]
    core_y = [-0.70_dp, -0.70_dp, -0.40_dp, -0.40_dp, 0.40_dp, 0.40_dp, &
        0.70_dp, 0.70_dp, -0.70_dp]
    do ix = 1, outline_count
        angle = 2.0_dp*pi*real(ix - 1, dp)/real(outline_count - 1, dp)
        coil_x(ix) = -0.62_dp + 0.15_dp*cos(angle)
        coil_y(ix) = 0.00_dp + 0.31_dp*sin(angle)
    end do
    gap_x = [0.55_dp, 0.55_dp]
    gap_y = [-0.40_dp, 0.40_dp]
    open (newunit=unit, file=output_directory//"/geometry.csv", &
        status="replace", action="write")
    write (unit, "(a)") "feature,index,x,y"
    do ix = 1, size(core_x)
        write (unit, "(a,i0,',',*(es24.16,:,','))") "core", ix, core_x(ix), core_y(ix)
    end do
    do ix = 1, outline_count
        write (unit, "(a,i0,',',*(es24.16,:,','))") "coil", ix, coil_x(ix), coil_y(ix)
    end do
    close (unit)

    ! Primary physical view: magnitude, directions, and supplied-array layout.
    call figure(figsize=[8.6_dp, 6.4_dp])
    call contourf(x_grid, y_grid, bmag, cmap="viridis", show_colorbar=.true.)
    call colorbar(label="|B| (manufactured magnetostatic field)")
    call quiver(arrow_x, arrow_y, arrow_u, arrow_v, scale=1.0_dp, &
        scale_units="xy", angles="xy", color=[0.0_dp, 0.0_dp, 0.0_dp], &
        width=0.004_dp, headwidth=4.0_dp)
    ! Keep a line-segment fallback for minimal FortPlot backends that do not
    ! render quiver heads; the direction field must remain visible in the
    ! primary physical plot.
    do arrow = 1, arrow_count
        call plot([arrow_x(arrow), arrow_x(arrow) + arrow_u(arrow)], &
            [arrow_y(arrow), arrow_y(arrow) + arrow_v(arrow)], &
            color=[0.0_dp, 0.0_dp, 0.0_dp], linewidth=1.0_dp)
    end do
    call plot(core_x, core_y, color=[1.0_dp, 0.85_dp, 0.1_dp], linewidth=2.5_dp)
    call plot(coil_x, coil_y, color=[1.0_dp, 0.25_dp, 0.1_dp], linewidth=2.5_dp)
    call plot(gap_x, gap_y, color=[0.95_dp, 0.95_dp, 0.95_dp], linewidth=2.0_dp)
    call xlabel("x")
    call ylabel("y")
    call title("TEAM-3-shaped neutral C-core solution")
    call savefig(output_directory//"/solution.png")
    call record_gallery_stage("physical_solution")

    ! Lift the scalar field into z to expose a genuine 3-D solution surface.
    do iy = 1, ny
        do ix = 1, nx
            surface_x(iy, ix) = x_grid(ix)
            surface_y(iy, ix) = y_grid(iy)
            surface_z(iy, ix) = 0.22_dp*bmag(iy, ix)/scale
        end do
    end do
    call figure(figsize=[8.2_dp, 6.4_dp])
    call add_parametric_surface(surface_x, surface_y, surface_z, &
        color="royalblue", alpha=0.72_dp, linewidth=0.0_dp, filled=.true., &
        row_stride=surface_stride, column_stride=surface_stride, &
        label="|B| height surface")
    field_scale = maxval(bmag)
    do arrow = 1, arrow_count
        arrow_length = 0.14_dp*sqrt(arrow_u(arrow)**2 + arrow_v(arrow)**2)
        call add_3d_plot([arrow_x(arrow), arrow_x(arrow) + arrow_u(arrow)], &
            [arrow_y(arrow), arrow_y(arrow) + arrow_v(arrow)], &
            [0.22_dp*bmag_at(arrow_x(arrow), arrow_y(arrow))/scale, &
            0.22_dp*bmag_at(arrow_x(arrow) + arrow_u(arrow), &
            arrow_y(arrow) + arrow_v(arrow))/scale], color="black", linewidth=1.0_dp)
    end do
    call title("TEAM-3-shaped neutral solution: 3-D field surface")
    call savefig(output_directory//"/solution_3d.png")

    call figure(figsize=[8.0_dp, 4.8_dp])
    call plot(probe_x, probe_b, label="|B|", color=[0.12_dp, 0.35_dp, 0.75_dp], &
        linewidth=2.0_dp)
    call plot(probe_x, probe_j, label="|J_z|", color=[0.8_dp, 0.25_dp, 0.1_dp], &
        linewidth=2.0_dp)
    call xlabel("x at y=0 (air-gap cut)")
    call ylabel("manufactured amplitude")
    call title("TEAM-3 neutral probe response")
    call legend()
    call savefig(output_directory//"/probe.png")

    open (newunit=unit, file=output_directory//"/probe.csv", &
        status="replace", action="write")
    write (unit, "(a)") "x,y,B_magnitude,Jz_magnitude"
    do ix = 1, nx
        write (unit, "(*(es24.16,:,','))") probe_x(ix), 0.0_dp, probe_b(ix), probe_j(ix)
    end do
    close (unit)
    open (newunit=unit, file=output_directory//"/diagnostics.csv", &
        status="replace", action="write")
    write (unit, "(a,es24.16)") "mean_magnetic_energy,", energy
    write (unit, "(a,es24.16)") "mean_current_energy,", current_energy
    write (unit, "(a,es24.16)") "mean_source_energy,", source_energy
    write (unit, "(a,es24.16)") "divergence_proxy,", max_divergence
    write (unit, "(a,es24.16)") "max_relative_permeability,", maxval(mu_r_field)
    close (unit)
    call record_gallery_stage("diagnostics")

contains

    subroutine initialize_gallery_sequence()
        integer :: local_unit, local_status

        open (newunit=local_unit, file=output_directory//"/gallery_sequence.txt", &
            status="replace", action="write", iostat=local_status)
        if (local_status /= 0) error stop "cannot initialize gallery sequence"
        close (local_unit)
    end subroutine initialize_gallery_sequence

    subroutine record_gallery_stage(stage)
        character(*), intent(in) :: stage
        integer :: local_unit, local_status

        open (newunit=local_unit, file=output_directory//"/gallery_sequence.txt", &
            status="old", position="append", action="write", iostat=local_status)
        if (local_status /= 0) error stop "cannot record gallery sequence"
        write (local_unit, "(a)", iostat=local_status) stage
        close (local_unit)
        if (local_status /= 0) error stop "cannot write gallery sequence"
    end subroutine record_gallery_stage

    pure real(dp) function material_response(x, y)
        real(dp), intent(in) :: x, y
        real(dp) :: core_envelope

        core_envelope = exp(-3.5_dp*(max(abs(x + 0.25_dp), 0.0_dp))**2 - &
            1.8_dp*y**2)
        material_response = 1.0_dp + 9.0_dp*core_envelope
    end function material_response

    pure real(dp) function source_response(x, y)
        real(dp), intent(in) :: x, y

        source_response = 1.7_dp*exp(-35.0_dp*(x + 0.62_dp)**2 - 7.0_dp*y**2) &
            * (1.0_dp + 0.25_dp*cos(3.0_dp*pi*y))
    end function source_response

    pure subroutine manufactured_field(x, y, az, ax, ay, axx, ayy, bx, by, jz)
        real(dp), intent(in) :: x, y
        real(dp), intent(out) :: az, ax, ay, axx, ayy, bx, by, jz
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
        bx = ay
        by = -ax
        jz = -(axx + ayy)
    end subroutine manufactured_field

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

    pure real(dp) function bmag_at(x, y)
        real(dp), intent(in) :: x, y
        real(dp) :: az_local, ax_local, ay_local, axx_local, ayy_local
        real(dp) :: bx_local, by_local, jz_local

        call manufactured_field(x, y, az_local, ax_local, ay_local, &
            axx_local, ayy_local, bx_local, by_local, jz_local)
        bmag_at = sqrt(bx_local**2 + by_local**2)
    end function bmag_at

end program team3_neutral_benchmark
