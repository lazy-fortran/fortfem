---
title: free_boundary_port_gallery Example
---

# free_boundary_port_gallery Example

# Free-boundary port trace gallery

This small gallery fixture samples a manufactured vector trace on a circular
toroidal boundary and sends it through FortFEM's
`assemble_free_boundary_port_residual` contract.  The first image is the
physical boundary: the torus is shown in Cartesian coordinates, point colour
shows the supplied trace magnitude, and short black segments show its tangent
and normal components.

The residual is the positive weighted mismatch

\[
 r_q = w_q\,(t_q-g_q-k_q),
\]

where `t` is the physical trace, `g` is an external target, and `k` is a
manufactured sheet/jump target.  The program also checks the generated JVP by
a centred difference and the VJP by the real dot-product identity.  CSV files
contain the sampled geometry, traces, residual, and derivative diagnostics.

This is deliberately a **neutral numerical contract example**.  It does not
implement an equilibrium solver, coil model, vacuum/BEM/DTN operator, or
free-boundary physics.  Those choices belong to a caller-owned adapter; this
fixture only demonstrates the trace pairing and its differentiable residual.

The contract is documented in the source module
`src/operators/free_boundary_port_residual.f90` and is exposed through the
canonical `fortfem_boundary` facade.

## Usage

```bash
fpm run --example free_boundary_port_gallery
```

## Source Code

```fortran
program free_boundary_port_gallery
    !! Physical-first neutral free-boundary-port trace gallery.
    !!
    !! A manufactured vector trace is sampled on a torus and paired with a
    !! caller-owned external target and sheet-current target.  No equilibrium,
    !! coil, BEM, DtN, or free-boundary physics is inferred here.
    use fortfem_boundary, only: &
        assemble_free_boundary_port_residual, &
        assemble_free_boundary_port_residual_jvp, &
        assemble_free_boundary_port_residual_vjp
    use fortfem_kinds, only: dp
    use fortplot, only: add_3d_plot, add_parametric_surface, add_scatter, &
        colorbar, figure, pcolormesh, plot, savefig, title, view_init, &
        xlabel, ylabel
    use fortsparse, only: fortsparse_status_t
    implicit none

    integer, parameter :: theta_count = 25, phi_count = 33
    integer, parameter :: theta_samples = theta_count - 1
    integer, parameter :: phi_samples = phi_count - 1
    integer, parameter :: sample_count = theta_samples*phi_samples
    integer, parameter :: arrow_theta_stride = 3, arrow_phi_stride = 4
    real(dp), parameter :: pi = acos(-1.0_dp)
    real(dp), parameter :: major_radius = 2.4_dp, minor_radius = 0.7_dp
    real(dp), parameter :: orange(3) = [230.0_dp, 159.0_dp, 0.0_dp]/255.0_dp
    real(dp), parameter :: fd_step = 2.0e-7_dp
    character(*), parameter :: output_directory = &
        "output/example/free_boundary_port_gallery"

    real(dp) :: theta(theta_count), phi(phi_count)
    real(dp) :: surface_x(theta_count, phi_count)
    real(dp) :: surface_y(theta_count, phi_count)
    real(dp) :: surface_z(theta_count, phi_count)
    real(dp) :: trace_x(sample_count), trace_y(sample_count)
    real(dp) :: trace_z(sample_count)
    real(dp) :: physical_trace(sample_count, 3)
    real(dp) :: external_target(sample_count, 3)
    real(dp) :: sheet_current(sample_count, 3)
    real(dp) :: physical_trace_dot(sample_count, 3)
    real(dp) :: external_target_dot(sample_count, 3)
    real(dp) :: sheet_current_dot(sample_count, 3)
    real(dp) :: weights(sample_count), weights_dot(sample_count)
    real(dp) :: residual(sample_count, 3), residual_dot(sample_count, 3)
    real(dp) :: residual_plus(sample_count, 3), residual_minus(sample_count, 3)
    real(dp) :: residual_bar(sample_count, 3)
    real(dp) :: physical_trace_bar(sample_count, 3)
    real(dp) :: external_target_bar(sample_count, 3)
    real(dp) :: sheet_current_bar(sample_count, 3)
    real(dp) :: weights_bar(sample_count)
    real(dp) :: trace_norm(sample_count), residual_norm(sample_count)
    real(dp) :: trace_norm_grid(theta_count, phi_count)
    real(dp) :: residual_error, vjp_error, geometry_error
    real(dp) :: fd_error, objective, objective_dot
    real(dp) :: lhs, rhs, area_weight, radius, amplitude
    real(dp) :: normal(3), poloidal_tangent(3), toroidal_tangent(3)
    real(dp) :: line_x(2), line_y(2), line_z(2), vector_scale
    real(dp) :: start_time, elapsed_seconds
    integer :: i, j, point, arrow, unit, command_status
    type(fortsparse_status_t) :: status

    command_status = 0
    call execute_command_line("mkdir -p "//output_directory, &
        exitstat=command_status)
    if (command_status /= 0) error stop "cannot create gallery output directory"
    call initialize_gallery_sequence()

    do j = 1, phi_count
        phi(j) = 2.0_dp*pi*real(j - 1, dp)/real(phi_count - 1, dp)
    end do
    do i = 1, theta_count
        theta(i) = 2.0_dp*pi*real(i - 1, dp)/real(theta_count - 1, dp)
    end do

    do j = 1, phi_count
        do i = 1, theta_count
            radius = major_radius + minor_radius*cos(theta(i))
            surface_x(i, j) = radius*cos(phi(j))
            surface_y(i, j) = radius*sin(phi(j))
            surface_z(i, j) = minor_radius*sin(theta(i))
        end do
    end do

    point = 0
    geometry_error = 0.0_dp
    trace_norm_grid = 0.0_dp
    do j = 1, phi_samples
        do i = 1, theta_samples
            point = point + 1
            radius = major_radius + minor_radius*cos(theta(i))
            trace_x(point) = radius*cos(phi(j))
            trace_y(point) = radius*sin(phi(j))
            trace_z(point) = minor_radius*sin(theta(i))
            normal = [cos(theta(i))*cos(phi(j)), &
                cos(theta(i))*sin(phi(j)), sin(theta(i))]
            poloidal_tangent = [-sin(theta(i))*cos(phi(j)), &
                -sin(theta(i))*sin(phi(j)), cos(theta(i))]
            toroidal_tangent = [-sin(phi(j)), cos(phi(j)), 0.0_dp]
            amplitude = 0.78_dp + 0.18_dp*cos(theta(i) - 2.0_dp*phi(j)) + &
                0.08_dp*sin(2.0_dp*theta(i) + phi(j))
            external_target(point, :) = 0.42_dp*normal + &
                0.04_dp*cos(phi(j))*poloidal_tangent
            sheet_current(point, :) = 0.11_dp*sin(theta(i) + phi(j))* &
                toroidal_tangent
            physical_trace(point, :) = external_target(point, :) + &
                sheet_current(point, :) + amplitude*normal + &
                0.06_dp*poloidal_tangent
            physical_trace_dot(point, :) = 0.04_dp*normal + &
                0.02_dp*poloidal_tangent
            external_target_dot(point, :) = 0.015_dp*normal
            sheet_current_dot(point, :) = 0.01_dp*toroidal_tangent
            area_weight = minor_radius*radius* &
                (2.0_dp*pi/real(theta_samples, dp))* &
                (2.0_dp*pi/real(phi_samples, dp))
            weights(point) = area_weight*(1.0_dp + 0.08_dp*cos(theta(i))**2)
            weights_dot(point) = 0.01_dp*weights(point)*sin(phi(j))
            trace_norm(point) = sqrt(sum(physical_trace(point, :)**2))
            trace_norm_grid(i, j) = trace_norm(point)
            geometry_error = max(geometry_error, abs( &
                sqrt(trace_x(point)**2 + trace_y(point)**2) - radius))
        end do
        trace_norm_grid(theta_count, j) = trace_norm_grid(1, j)
    end do
    trace_norm_grid(:, phi_count) = trace_norm_grid(:, 1)

    call cpu_time(start_time)
    call assemble_free_boundary_port_residual( &
        physical_trace, external_target, weights, residual, status, &
        sheet_current)
    if (status%code /= 0) error stop "free-boundary port residual failed"
    call assemble_free_boundary_port_residual_jvp( &
        physical_trace, external_target, weights, physical_trace_dot, &
        external_target_dot, weights_dot, residual_dot, status, sheet_current, &
        sheet_current_dot)
    if (status%code /= 0) error stop "free-boundary port JVP failed"
    call assemble_free_boundary_port_residual( &
        physical_trace + fd_step*physical_trace_dot, &
        external_target + fd_step*external_target_dot, &
        weights + fd_step*weights_dot, residual_plus, status, &
        sheet_current + fd_step*sheet_current_dot)
    call assemble_free_boundary_port_residual( &
        physical_trace - fd_step*physical_trace_dot, &
        external_target - fd_step*external_target_dot, &
        weights - fd_step*weights_dot, residual_minus, status, &
        sheet_current - fd_step*sheet_current_dot)
    if (status%code /= 0) error stop "free-boundary port finite difference failed"
    fd_error = maxval(abs(residual_dot - &
        (residual_plus - residual_minus)/(2.0_dp*fd_step)))

    residual_bar = 0.25_dp
    call assemble_free_boundary_port_residual_vjp( &
        physical_trace, external_target, weights, residual_bar, &
        physical_trace_bar, external_target_bar, weights_bar, status, &
        sheet_current, sheet_current_bar)
    if (status%code /= 0) error stop "free-boundary port VJP failed"
    lhs = sum(residual_bar*residual_dot)
    rhs = sum(physical_trace_bar*physical_trace_dot) + &
        sum(external_target_bar*external_target_dot) + &
        dot_product(weights_bar, weights_dot) + &
        sum(sheet_current_bar*sheet_current_dot)
    vjp_error = abs(lhs - rhs)
    residual_norm = sqrt(sum(residual**2, dim=2))
    residual_error = 0.0_dp
    do point = 1, sample_count
        i = 1 + mod(point - 1, theta_samples)
        j = 1 + (point - 1)/theta_samples
        normal = [cos(theta(i))*cos(phi(j)), &
            cos(theta(i))*sin(phi(j)), sin(theta(i))]
        poloidal_tangent = [-sin(theta(i))*cos(phi(j)), &
            -sin(theta(i))*sin(phi(j)), cos(theta(i))]
        amplitude = 0.78_dp + 0.18_dp*cos(theta(i) - 2.0_dp*phi(j)) + &
            0.08_dp*sin(2.0_dp*theta(i) + phi(j))
        residual_error = max(residual_error, maxval(abs(residual(point, :) - &
            weights(point)*(amplitude*normal + 0.06_dp*poloidal_tangent))))
    end do
    objective = 0.5_dp*sum(weights*residual_norm**2)
    objective_dot = sum(0.5_dp*weights_dot*residual_norm**2 + &
        weights*sum(residual*residual_dot, dim=2))
    elapsed_seconds = 0.0_dp
    call cpu_time(elapsed_seconds)
    elapsed_seconds = elapsed_seconds - start_time
    if (fd_error > 3.0e-8_dp .or. vjp_error > 1.0e-12_dp .or. &
        residual_error > 1.0e-13_dp .or. geometry_error > 2.0e-12_dp) &
        error stop "free-boundary port analytical oracle failed"

    ! The physical boundary solution is intentionally the first gallery plot.
    call figure(figsize=[9.0_dp, 7.0_dp])
    call add_parametric_surface(surface_x, surface_y, surface_z, &
        color="lightsteelblue", alpha=0.24_dp, linewidth=0.0_dp, &
        filled=.true., label="manufactured toroidal boundary")
    call add_scatter(trace_x, trace_y, trace_z, c=trace_norm, cmap="viridis", &
        marker=".", markersize=5.0_dp, label="physical boundary trace")
    vector_scale = 0.18_dp/max(1.0e-12_dp, maxval(trace_norm))
    arrow = 0
    do j = 1, phi_samples, arrow_phi_stride
        do i = 1, theta_samples, arrow_theta_stride
            arrow = arrow + 1
            point = 1 + (j - 1)*theta_samples + i - 1
            line_x = [trace_x(point), trace_x(point) + &
                vector_scale*physical_trace(point, 1)]
            line_y = [trace_y(point), trace_y(point) + &
                vector_scale*physical_trace(point, 2)]
            line_z = [trace_z(point), trace_z(point) + &
                vector_scale*physical_trace(point, 3)]
            call add_3d_plot(line_x, line_y, line_z, color="black", &
                linewidth=1.2_dp)
        end do
    end do
    call colorbar(label="|physical trace| [1]")
    call xlabel("x [m]")
    call ylabel("y [m]")
    call title("Free-boundary port: manufactured toroidal boundary trace")
    call view_init(elev=28.0_dp, azim=-55.0_dp)
    call savefig(output_directory//"/free_boundary_port_solution_3d.png")
    call record_gallery_stage("physical_solution")

    call figure(figsize=[8.5_dp, 5.8_dp])
    call pcolormesh(theta, phi, transpose(trace_norm_grid), cmap="viridis")
    call colorbar(label="|physical trace| [1]")
    call xlabel("poloidal angle theta [rad]")
    call ylabel("toroidal angle phi [rad]")
    call title("Boundary trace magnitude in toroidal coordinates")
    call savefig(output_directory//"/free_boundary_port_trace_2d.png")

    call figure(figsize=[8.0_dp, 5.2_dp])
    call plot(theta(1:theta_samples), residual_norm(1:theta_samples), &
        color=orange, linewidth=2.0_dp)
    call xlabel("poloidal angle theta [rad] at phi = 0")
    call ylabel("weighted port residual norm [1]")
    call title("Free-boundary port residual diagnostic")
    call savefig(output_directory//"/free_boundary_port_diagnostics_1d.png")
    call record_gallery_stage("diagnostics")

    call write_outputs()
    open (newunit=unit, file=output_directory//"/benchmark.json", &
        status="replace", action="write")
    write (unit, "(a)") "{"// &
        '"samples":'//trim(adjustl(itoa(sample_count)))//","// &
        '"fd_error":'//trim(adjustl(rtoa(fd_error)))//","// &
        '"vjp_error":'//trim(adjustl(rtoa(vjp_error)))//","// &
        '"residual_error":'//trim(adjustl(rtoa(residual_error)))//","// &
        '"geometry_error":'//trim(adjustl(rtoa(geometry_error)))//","// &
        '"objective":'//trim(adjustl(rtoa(objective)))//","// &
        '"objective_dot":'//trim(adjustl(rtoa(objective_dot)))//","// &
        '"elapsed_seconds":'//trim(adjustl(rtoa(elapsed_seconds)))//","// &
        '"provenance":"analytic-torus-free-boundary-port-v1",'// &
        '"primary_plot":"free_boundary_port_solution_3d.png",'// &
        '"closure":"neutral-caller-owned-trace"}'
    close (unit)

contains

    subroutine initialize_gallery_sequence()
        integer :: sequence_unit
        open (newunit=sequence_unit, &
            file=output_directory//"/gallery_sequence.txt", status="replace", &
            action="write")
        close (sequence_unit)
    end subroutine initialize_gallery_sequence

    subroutine record_gallery_stage(stage)
        character(len=*), intent(in) :: stage
        integer :: sequence_unit
        open (newunit=sequence_unit, &
            file=output_directory//"/gallery_sequence.txt", status="old", &
            position="append", action="write")
        write (sequence_unit, "(a)") trim(stage)
        close (sequence_unit)
    end subroutine record_gallery_stage

    subroutine write_outputs()
        integer :: csv_unit, local_i, local_j, local_point
        open (newunit=csv_unit, file=output_directory//"/boundary_trace.csv", &
            status="replace", action="write")
        write (csv_unit, "(a)") &
            "theta,phi,x,y,z,trace_x,trace_y,trace_z,target_x,target_y,"// &
            "target_z,sheet_x,sheet_y,sheet_z,residual_x,residual_y,"// &
            "residual_z,weight"
        local_point = 0
        do local_j = 1, phi_samples
            do local_i = 1, theta_samples
                local_point = local_point + 1
                write (csv_unit, "(*(es24.16,:,','))") theta(local_i), &
                    phi(local_j), trace_x(local_point), trace_y(local_point), &
                    trace_z(local_point), physical_trace(local_point, :), &
                    external_target(local_point, :), sheet_current(local_point, :), &
                    residual(local_point, :), weights(local_point)
            end do
        end do
        close (csv_unit)
        open (newunit=csv_unit, &
            file=output_directory//"/free_boundary_port_diagnostics.csv", &
            status="replace", action="write")
        write (csv_unit, "(a)") &
            "samples,objective,objective_dot,fd_error,vjp_error,residual_error,geometry_error,elapsed_seconds"
        write (csv_unit, "(*(es24.16,:,','))") real(sample_count, dp), objective, &
            objective_dot, fd_error, vjp_error, residual_error, geometry_error, &
            elapsed_seconds
        close (csv_unit)
    end subroutine write_outputs

    function itoa(value) result(text)
        integer, intent(in) :: value
        character(32) :: text
        write (text, "(i0)") value
    end function itoa

    function rtoa(value) result(text)
        real(dp), intent(in) :: value
        character(64) :: text
        write (text, "(es24.16e3)") value
    end function rtoa

end program free_boundary_port_gallery
```

## Generated Plots

### primary.png

![primary.png](../../media/examples/free_boundary_port_gallery/primary.png)

### free_boundary_port_solution_3d.png

![free_boundary_port_solution_3d.png](../../media/examples/free_boundary_port_gallery/free_boundary_port_solution_3d.png)

---

[← Back to all examples](../index.html)
