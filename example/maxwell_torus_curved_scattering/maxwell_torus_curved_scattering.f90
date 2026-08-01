program maxwell_torus_curved_scattering
    !! Exact-curved torus PEC scattering with a regularized Maxwell CFIE.
    use fortfem_api, only: &
        evaluate_maxwell_torus_curved_far_field_rwg_3d, &
        evaluate_maxwell_torus_curved_magnetic_field_rwg_3d, &
        generate_torus_surface_mesh, &
        assemble_maxwell_torus_curved_dtn_rwg_3d, &
        apply_maxwell_trace_to_flux_map, &
        solve_maxwell_pec_torus_curved_regularized_cfie_rwg_multiple_3d
    use fortfem_kinds, only: dp
    use fortplot, only: &
        add_3d_plot, add_scatter, colorbar, figure, legend, pcolormesh, plot, &
        savefig, title, xlabel, ylabel
    implicit none

    integer, parameter :: azimuth_cells = 36, polar_cells = 18
    integer, parameter :: field_cells = 20
    integer, parameter :: mesh_polar_nodes = 3, mesh_azimuth_nodes = 3
    integer, parameter :: trace_points = 73
    real(dp), parameter :: major_radius = 2.0_dp, minor_radius = 0.6_dp
    real(dp), parameter :: wave_number = 0.45_dp, impedance = 1.7_dp
    real(dp), parameter :: pi = acos(-1.0_dp)
    character(*), parameter :: output_directory = &
        "output/example/maxwell_torus_curved_scattering"
    complex(dp), allocatable :: currents(:, :)
    complex(dp), allocatable :: dtn_map(:, :), dtn_trace(:), dtn_flux(:)
    complex(dp) :: far_field(3), first_reciprocal(3), second_reciprocal(3)
    complex(dp), parameter :: z_polarization(3) = [ &
        cmplx(0.0_dp, 0.0_dp, dp), cmplx(0.0_dp, 0.0_dp, dp), &
        cmplx(1.0_dp, 0.0_dp, dp)]
    integer, allocatable :: triangles(:, :)
    real(dp), allocatable :: parameters(:, :), vertices(:, :)
    real(dp) :: azimuth, direction(3), end_time, far_field_seconds
    real(dp) :: magnetic_seconds
    real(dp) :: first_amplitude, polar, radius, reciprocity_error
    real(dp) :: dtn_seconds, solve_seconds, start_time
    real(dp) :: azimuth_edges(azimuth_cells + 1)
    real(dp) :: polar_edges(polar_cells + 1)
    real(dp) :: field_x_edges(field_cells + 1)
    real(dp) :: field_z_edges(field_cells + 1)
    real(dp) :: rcs_map(azimuth_cells, polar_cells)
    real(dp) :: magnetic_map(field_cells, field_cells)
    real(dp) :: trace_azimuth(trace_points), trace_rcs(trace_points)
    real(dp) :: x(azimuth_cells*polar_cells)
    real(dp) :: y(azimuth_cells*polar_cells)
    real(dp) :: z(azimuth_cells*polar_cells)
    integer :: azimuth_index, point, polar_index, status, unit

    call execute_command_line("mkdir -p "//output_directory)
    call generate_torus_surface_mesh( &
        major_radius, minor_radius, mesh_polar_nodes, mesh_azimuth_nodes, &
        vertices, triangles, parameters)
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
    call build_dtn_trace()
    call cpu_time(end_time)
    dtn_seconds = end_time - start_time
    call cpu_time(start_time)
    call build_far_field_map()
    call build_magnetic_field_map()
    call build_trace()
    call evaluate_reciprocity()
    call cpu_time(end_time)
    far_field_seconds = end_time - start_time
    call render_plots()
    call write_outputs()

contains

    subroutine build_dtn_trace()
        integer :: coefficient

        call assemble_maxwell_torus_curved_dtn_rwg_3d( &
            vertices, triangles, parameters, major_radius, minor_radius, &
            wave_number, impedance, 3, 3.0e-4_dp, 1, 0.12_dp, dtn_map, status)
        if (status /= 0) error stop "curved torus Maxwell DtN assembly failed"
        allocate(dtn_trace(size(dtn_map, 2)), dtn_flux(size(dtn_map, 1)))
        do coefficient = 1, size(dtn_trace)
            dtn_trace(coefficient) = cmplx( &
                cos(0.17_dp*real(coefficient, dp)), &
                sin(0.11_dp*real(coefficient, dp)), dp)
        end do
        call apply_maxwell_trace_to_flux_map( &
            dtn_map, dtn_trace, dtn_flux, status)
        if (status /= 0) error stop "curved torus Maxwell DtN action failed"
    end subroutine build_dtn_trace

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

    subroutine build_magnetic_field_map()
        integer :: x_index, z_index
        real(dp) :: x_value, z_value, target(3)
        complex(dp) :: magnetic_field(3)

        call cpu_time(start_time)
        do x_index = 1, field_cells + 1
            field_x_edges(x_index) = -1.0_dp + 2.0_dp*real(x_index - 1, dp)/ &
                real(field_cells, dp)
            field_z_edges(x_index) = field_x_edges(x_index)
        end do
        do z_index = 1, field_cells
            z_value = 0.5_dp*(field_z_edges(z_index) + &
                field_z_edges(z_index + 1))
            do x_index = 1, field_cells
                x_value = 0.5_dp*(field_x_edges(x_index) + &
                    field_x_edges(x_index + 1))
                target = [x_value, 0.0_dp, z_value]
                call evaluate_maxwell_torus_curved_magnetic_field_rwg_3d( &
                    vertices, triangles, parameters, major_radius, minor_radius, &
                    currents(:, 1), target, wave_number, 4, magnetic_field, &
                    status)
                if (status /= 0) &
                    error stop "curved torus magnetic field evaluation failed"
                magnetic_map(x_index, z_index) = &
                    sqrt(sum(abs(magnetic_field)**2))
            end do
        end do
        call cpu_time(end_time)
        magnetic_seconds = end_time - start_time
    end subroutine build_magnetic_field_map

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
        real(dp) :: mesh_curve_x(mesh_polar_nodes + 1)
        real(dp) :: mesh_curve_y(mesh_polar_nodes + 1)
        real(dp) :: mesh_curve_z(mesh_polar_nodes + 1)
        real(dp) :: mesh_ring_x(mesh_azimuth_nodes + 1)
        real(dp) :: mesh_ring_y(mesh_azimuth_nodes + 1)
        real(dp) :: mesh_ring_z(mesh_azimuth_nodes + 1)
        real(dp) :: curve_azimuth, curve_polar, curve_radius
        integer :: coefficient, curve_index, mesh_vertex, point_index

        call figure(figsize=[8.5_dp, 6.5_dp])
        call pcolormesh( &
            field_x_edges, field_z_edges, magnetic_map, cmap="inferno")
        call colorbar(label="|H_scattered|")
        call xlabel("x at y = 0")
        call ylabel("z")
        call title("Computed torus Maxwell scattered magnetic field slice")
        call savefig(output_directory//"/maxwell_torus_solution_2d.png")

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
        ! Draw the coarse solver mesh along its parameter lines.  Drawing
        ! every triangle edge projects hidden diagonals through the torus and
        ! makes a valid periodic mesh look self-intersecting.
        do curve_index = 0, mesh_azimuth_nodes - 1
            do point_index = 0, mesh_polar_nodes
                mesh_vertex = curve_index*mesh_polar_nodes + &
                    modulo(point_index, mesh_polar_nodes) + 1
                mesh_curve_x(point_index + 1) = vertices(1, mesh_vertex)
                mesh_curve_y(point_index + 1) = vertices(2, mesh_vertex)
                mesh_curve_z(point_index + 1) = vertices(3, mesh_vertex)
            end do
            if (curve_index == 0) then
                call add_3d_plot(mesh_curve_x, mesh_curve_y, mesh_curve_z, &
                    label="coarse CFIE mesh parameter lines", color="gray", &
                    linewidth=0.7_dp)
            else
                call add_3d_plot(mesh_curve_x, mesh_curve_y, mesh_curve_z, &
                    color="gray", linewidth=0.7_dp)
            end if
        end do
        do curve_index = 0, mesh_polar_nodes - 1
            do point_index = 0, mesh_azimuth_nodes
                mesh_vertex = modulo(point_index, mesh_azimuth_nodes) * &
                    mesh_polar_nodes + curve_index + 1
                mesh_ring_x(point_index + 1) = vertices(1, mesh_vertex)
                mesh_ring_y(point_index + 1) = vertices(2, mesh_vertex)
                mesh_ring_z(point_index + 1) = vertices(3, mesh_vertex)
            end do
            call add_3d_plot(mesh_ring_x, mesh_ring_y, mesh_ring_z, &
                color="gray", linewidth=0.7_dp)
        end do
        call title("Exact-curved torus with periodic CFIE mesh lines")
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

        call figure(figsize=[9.0_dp, 5.5_dp])
        call plot(real([(coefficient, coefficient=1, size(dtn_flux))], dp), &
            abs(dtn_flux), label="RWG/RBC weak DtN response")
        call xlabel("RBC flux degree of freedom")
        call ylabel("absolute weak tangential flux")
        call title("Exact-curved torus Maxwell weak DtN response")
        call legend()
        call savefig(output_directory//"/maxwell_torus_dtn_1d.png")
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
        write (unit, "(a,es14.6)") "weak DtN assembly seconds: ", dtn_seconds
        write (unit, "(a,es14.6)") "weak DtN response norm: ", &
            sqrt(sum(abs(dtn_flux)**2))
        write (unit, "(a,es14.6)") &
            "far-field sampling seconds: ", far_field_seconds
        write (unit, "(a,es14.6)") &
            "curved magnetic-field slice seconds: ", magnetic_seconds
        write (unit, "(a,es14.6)") &
            "maximum scattered magnetic-field magnitude: ", maxval(magnetic_map)
        write (unit, "(a,es14.6)") &
            "Lorentz reciprocity relative error: ", reciprocity_error
        write (unit, "(a,es14.6)") "maximum bistatic RCS: ", maxval(rcs_map)
        close (unit)
        if (reciprocity_error >= 2.0e-1_dp) &
            error stop "toroidal Maxwell gallery reciprocity regression"
    end subroutine write_outputs

end program maxwell_torus_curved_scattering
