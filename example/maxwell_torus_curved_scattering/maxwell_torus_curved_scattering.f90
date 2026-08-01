program maxwell_torus_curved_scattering
    !! Exact-curved torus PEC scattering with a regularized Maxwell CFIE.
    use fortfem_api, only: &
        evaluate_maxwell_torus_curved_far_field_rwg_3d, &
        generate_torus_surface_mesh, &
        solve_maxwell_pec_torus_curved_regularized_cfie_rwg_multiple_3d
    use fortfem_kinds, only: dp
    use fortplot, only: &
        add_3d_plot, add_scatter, colorbar, figure, legend, pcolormesh, plot, &
        savefig, title, xlabel, ylabel
    implicit none

    integer, parameter :: azimuth_cells = 36, polar_cells = 18
    integer, parameter :: trace_points = 73
    real(dp), parameter :: major_radius = 2.0_dp, minor_radius = 0.6_dp
    real(dp), parameter :: wave_number = 0.45_dp, impedance = 1.7_dp
    real(dp), parameter :: pi = acos(-1.0_dp)
    character(*), parameter :: output_directory = &
        "output/example/maxwell_torus_curved_scattering"
    complex(dp), allocatable :: currents(:, :)
    complex(dp) :: far_field(3), first_reciprocal(3), second_reciprocal(3)
    complex(dp), parameter :: z_polarization(3) = [ &
        cmplx(0.0_dp, 0.0_dp, dp), cmplx(0.0_dp, 0.0_dp, dp), &
        cmplx(1.0_dp, 0.0_dp, dp)]
    integer, allocatable :: triangles(:, :)
    real(dp), allocatable :: parameters(:, :), vertices(:, :)
    real(dp) :: azimuth, direction(3), end_time, far_field_seconds
    real(dp) :: first_amplitude, polar, radius, reciprocity_error
    real(dp) :: solve_seconds, start_time
    real(dp) :: azimuth_edges(azimuth_cells + 1)
    real(dp) :: polar_edges(polar_cells + 1)
    real(dp) :: rcs_map(azimuth_cells, polar_cells)
    real(dp) :: trace_azimuth(trace_points), trace_rcs(trace_points)
    real(dp) :: x(azimuth_cells*polar_cells)
    real(dp) :: y(azimuth_cells*polar_cells)
    real(dp) :: z(azimuth_cells*polar_cells)
    integer :: azimuth_index, point, polar_index, status, unit

    call execute_command_line("mkdir -p "//output_directory)
    call generate_torus_surface_mesh( &
        major_radius, minor_radius, 3, 3, vertices, triangles, parameters)
    call cpu_time(start_time)
    call solve_maxwell_pec_torus_curved_regularized_cfie_rwg_multiple_3d( &
        vertices, triangles, parameters, major_radius, minor_radius, reshape([ &
        1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, -1.0_dp, 0.0_dp], [3, 2]), &
        reshape([z_polarization, z_polarization], [3, 2]), wave_number, &
        impedance, 3, 3.0e-4_dp, 1, 0.12_dp, currents, status)
    call cpu_time(end_time)
    if (status /= 0) error stop "exact-torus Maxwell gallery solve failed"
    solve_seconds = end_time - start_time

    call cpu_time(start_time)
    call build_far_field_map()
    call build_trace()
    call evaluate_reciprocity()
    call cpu_time(end_time)
    far_field_seconds = end_time - start_time
    call render_plots()
    call write_outputs()

contains

    subroutine build_far_field_map()
        point = 0
        do azimuth_index = 1, azimuth_cells + 1
            azimuth_edges(azimuth_index) = 2.0_dp*pi* &
                real(azimuth_index - 1, dp)/real(azimuth_cells, dp)
        end do
        do polar_index = 1, polar_cells + 1
            polar_edges(polar_index) = pi* &
                real(polar_index - 1, dp)/real(polar_cells, dp)
        end do
        do polar_index = 1, polar_cells
            polar = 0.5_dp*( &
                polar_edges(polar_index) + polar_edges(polar_index + 1))
            do azimuth_index = 1, azimuth_cells
                azimuth = 0.5_dp*( &
                    azimuth_edges(azimuth_index) + &
                    azimuth_edges(azimuth_index + 1))
                direction = [ &
                    sin(polar)*cos(azimuth), sin(polar)*sin(azimuth), &
                    cos(polar)]
                call evaluate_maxwell_torus_curved_far_field_rwg_3d( &
                    vertices, triangles, parameters, major_radius, &
                    minor_radius, currents(:, 1), direction, wave_number, &
                    impedance, 8, far_field, status)
                if (status /= 0) error stop "toroidal far-field map failed"
                rcs_map(azimuth_index, polar_index) = &
                    4.0_dp*pi*sum(abs(far_field)**2)
                point = point + 1
                radius = rcs_map(azimuth_index, polar_index)
                x(point) = radius*direction(1)
                y(point) = radius*direction(2)
                z(point) = radius*direction(3)
            end do
        end do
        radius = max(maxval(rcs_map), tiny(1.0_dp))
        x = x/radius
        y = y/radius
        z = z/radius
    end subroutine build_far_field_map

    subroutine build_trace()
        integer :: sample

        do sample = 1, trace_points
            trace_azimuth(sample) = 2.0_dp*pi*real(sample - 1, dp)/ &
                real(trace_points - 1, dp)
            direction = [ &
                cos(trace_azimuth(sample)), sin(trace_azimuth(sample)), &
                0.0_dp]
            call evaluate_maxwell_torus_curved_far_field_rwg_3d( &
                vertices, triangles, parameters, major_radius, minor_radius, &
                currents(:, 1), direction, wave_number, impedance, 8, &
                far_field, status)
            if (status /= 0) error stop "toroidal azimuth trace failed"
            trace_rcs(sample) = 4.0_dp*pi*sum(abs(far_field)**2)
        end do
    end subroutine build_trace

    subroutine evaluate_reciprocity()
        call evaluate_maxwell_torus_curved_far_field_rwg_3d( &
            vertices, triangles, parameters, major_radius, minor_radius, &
            currents(:, 1), [0.0_dp, 1.0_dp, 0.0_dp], wave_number, impedance, &
            8, first_reciprocal, status)
        if (status /= 0) error stop "first toroidal reciprocal field failed"
        call evaluate_maxwell_torus_curved_far_field_rwg_3d( &
            vertices, triangles, parameters, major_radius, minor_radius, &
            currents(:, 2), [-1.0_dp, 0.0_dp, 0.0_dp], wave_number, &
            impedance, 8, second_reciprocal, status)
        if (status /= 0) error stop "second toroidal reciprocal field failed"
        first_amplitude = max( &
            abs(first_reciprocal(3)), abs(second_reciprocal(3)))
        reciprocity_error = &
            abs(first_reciprocal(3) - second_reciprocal(3))/first_amplitude
    end subroutine evaluate_reciprocity

    subroutine render_plots()
        integer, parameter :: visual_azimuth = 36, visual_polar = 18
        real(dp) :: curve_x(visual_azimuth + 1)
        real(dp) :: curve_y(visual_azimuth + 1)
        real(dp) :: curve_z(visual_azimuth + 1)
        real(dp) :: edge_x(4), edge_y(4), edge_z(4)
        real(dp) :: curve_azimuth, curve_polar, curve_radius
        integer :: curve_index, local_vertex, point_index, triangle, vertex

        call figure(figsize=[7.5_dp, 6.5_dp])
        call add_scatter( &
            vertices(1, :), vertices(2, :), vertices(3, :), &
            label="curved torus mesh vertices", marker=".", markersize=5.0_dp)
        do curve_index = 0, visual_polar
            curve_polar = 2.0_dp*pi*real(curve_index, dp)/ &
                real(visual_polar, dp)
            do point_index = 0, visual_azimuth
                curve_azimuth = 2.0_dp*pi*real(point_index, dp)/ &
                    real(visual_azimuth, dp)
                curve_radius = major_radius + minor_radius*cos(curve_polar)
                curve_x(point_index + 1) = curve_radius*cos(curve_azimuth)
                curve_y(point_index + 1) = curve_radius*sin(curve_azimuth)
                curve_z(point_index + 1) = minor_radius*sin(curve_polar)
            end do
            call add_3d_plot(curve_x, curve_y, curve_z, &
                label="exact torus geometry", color="teal", linewidth=0.7_dp)
        end do
        do curve_index = 0, visual_azimuth
            curve_azimuth = 2.0_dp*pi*real(curve_index, dp)/ &
                real(visual_azimuth, dp)
            do point_index = 0, visual_polar
                curve_polar = 2.0_dp*pi*real(point_index, dp)/ &
                    real(visual_polar, dp)
                curve_radius = major_radius + minor_radius*cos(curve_polar)
                curve_x(point_index + 1) = curve_radius*cos(curve_azimuth)
                curve_y(point_index + 1) = curve_radius*sin(curve_azimuth)
                curve_z(point_index + 1) = minor_radius*sin(curve_polar)
            end do
            call add_3d_plot(curve_x(:visual_polar + 1), &
                curve_y(:visual_polar + 1), curve_z(:visual_polar + 1), &
                color="teal", linewidth=0.7_dp)
        end do
        do triangle = 1, size(triangles, 2)
            do local_vertex = 1, 3
                vertex = triangles(local_vertex, triangle)
                edge_x(local_vertex) = vertices(1, vertex)
                edge_y(local_vertex) = vertices(2, vertex)
                edge_z(local_vertex) = vertices(3, vertex)
            end do
            edge_x(4) = edge_x(1)
            edge_y(4) = edge_y(1)
            edge_z(4) = edge_z(1)
            if (triangle == 1) then
                call add_3d_plot(edge_x, edge_y, edge_z, &
                    label="curved torus triangle edges", color="navy", &
                    linewidth=1.0_dp)
            else
                call add_3d_plot(edge_x, edge_y, edge_z, color="navy", &
                    linewidth=1.0_dp)
            end if
        end do
        call title("Exact-curved torus surface mesh used by Maxwell CFIE")
        call savefig(output_directory//"/torus_geometry_3d.png")

        call figure(figsize=[9.0_dp, 5.5_dp])
        call plot(trace_azimuth, trace_rcs, &
            label="exact-curved CFIE bistatic RCS", linestyle="-")
        call xlabel("observation azimuth phi [rad]")
        call ylabel("4 pi |E_infinity|^2")
        call title("PEC torus Maxwell scattering: equatorial cut")
        call legend()
        call savefig(output_directory//"/maxwell_torus_rcs_1d.png")

        call figure(figsize=[9.0_dp, 5.5_dp])
        call pcolormesh( &
            polar_edges, azimuth_edges, rcs_map, cmap="viridis")
        call colorbar(label="4 pi |E_infinity|^2")
        call xlabel("observation polar angle theta [rad]")
        call ylabel("observation azimuth phi [rad]")
        call title("PEC torus Maxwell bistatic RCS")
        call savefig(output_directory//"/maxwell_torus_rcs_2d.png")

        call figure(figsize=[7.5_dp, 6.5_dp])
        call add_scatter(x, y, z, label="normalized RCS surface", marker=".")
        call title("PEC torus Maxwell three-dimensional radiation pattern")
        call savefig(output_directory//"/maxwell_torus_rcs_3d.png")
    end subroutine render_plots

    subroutine write_outputs()
        integer :: sample

        open (newunit=unit, file=output_directory//"/rcs_trace.csv", &
            status="replace", action="write")
        write (unit, "(a)") "azimuth,radar_cross_section"
        do sample = 1, trace_points
            write (unit, "(*(es24.16,:,','))") &
                trace_azimuth(sample), trace_rcs(sample)
        end do
        close (unit)
        open (newunit=unit, file=output_directory//"/benchmark.txt", &
            status="replace", action="write")
        write (unit, "(a,i0)") "surface triangles: ", size(triangles, 2)
        write (unit, "(a,i0)") "RWG unknowns: ", size(currents, 1)
        write (unit, "(a,i0)") "incident fields: ", size(currents, 2)
        write (unit, "(a,i0)") "CFIE quadrature degree: ", 3
        write (unit, "(a,es14.6)") "batched CFIE solve seconds: ", solve_seconds
        write (unit, "(a,es14.6)") &
            "far-field sampling seconds: ", far_field_seconds
        write (unit, "(a,es14.6)") &
            "Lorentz reciprocity relative error: ", reciprocity_error
        write (unit, "(a,es14.6)") "maximum bistatic RCS: ", maxval(rcs_map)
        close (unit)
        if (reciprocity_error >= 2.0e-1_dp) &
            error stop "toroidal Maxwell gallery reciprocity regression"
    end subroutine write_outputs

end program maxwell_torus_curved_scattering
