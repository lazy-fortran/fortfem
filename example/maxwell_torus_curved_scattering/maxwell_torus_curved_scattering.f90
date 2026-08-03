program maxwell_torus_curved_scattering
    !! Exact-curved torus PEC scattering with a regularized Maxwell CFIE.
    use fortfem_boundary, only: &
        evaluate_maxwell_torus_curved_far_field_rwg_3d, &
        evaluate_maxwell_torus_curved_magnetic_field_rwg_3d, &
        assemble_maxwell_torus_curved_dtn_rwg_3d, &
        apply_maxwell_trace_to_flux_map, &
        solve_maxwell_pec_torus_curved_regularized_cfie_rwg_multiple_3d
    use fortfem_core, only: generate_torus_surface_mesh
    use fortfem_kinds, only: dp
    use fortplot, only: &
        add_parametric_surface, add_scatter, colorbar, figure, legend, pcolormesh, &
        plot, quiver, savefig, title, xlabel, ylabel
    implicit none

    integer, parameter :: azimuth_cells = 36, polar_cells = 18
    integer, parameter :: field_cells = 20
    integer, parameter :: default_mesh_polar_nodes = 3
    integer, parameter :: default_mesh_azimuth_nodes = 3
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
    real(dp) :: magnetic_vector_map(3, field_cells, field_cells)
    real(dp) :: trace_azimuth(trace_points), trace_rcs(trace_points)
    real(dp) :: x(azimuth_cells*polar_cells)
    real(dp) :: y(azimuth_cells*polar_cells)
    real(dp) :: z(azimuth_cells*polar_cells)
    integer :: azimuth_index, point, polar_index, status, unit
    integer :: mesh_polar_nodes, mesh_azimuth_nodes, quadrature_degree
    logical :: fast_gallery

    call execute_command_line("mkdir -p "//output_directory)
    fast_gallery = environment_flag("MAXWELL_TORUS_FAST")
    mesh_polar_nodes = default_mesh_polar_nodes
    mesh_azimuth_nodes = default_mesh_azimuth_nodes
    quadrature_degree = 3
    if (fast_gallery) then
        ! Preserve the physical CFIE and regularizer while using a reduced
        ! quadrature for the independent ten-second data oracle.
        ! (The torus mesh API requires at least three nodes per direction.)
        quadrature_degree = 2
    end if
    call generate_torus_surface_mesh( &
        major_radius, minor_radius, mesh_polar_nodes, mesh_azimuth_nodes, &
        vertices, triangles, parameters)
    call cpu_time(start_time)
    call solve_maxwell_pec_torus_curved_regularized_cfie_rwg_multiple_3d( &
        vertices, triangles, parameters, major_radius, minor_radius, reshape([ &
        1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, -1.0_dp, 0.0_dp], [3, 2]), &
        reshape([z_polarization, z_polarization], [3, 2]), wave_number, &
        impedance, quadrature_degree, 3.0e-4_dp, 1, 0.12_dp, currents, status)
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
    call write_gallery_sequence()
    call write_scattered_field_csv()
    call render_plots()
    call write_outputs()

contains

    logical function environment_flag(name)
        character(*), intent(in) :: name

        character(32) :: value
        integer :: environment_status

        value = ""
        call get_environment_variable(name, value, status=environment_status)
        environment_flag = environment_status == 0 .and. &
            len_trim(value) > 0 .and. trim(value) /= "0"
    end function environment_flag

    subroutine build_dtn_trace()
        integer :: coefficient

        call assemble_maxwell_torus_curved_dtn_rwg_3d( &
            vertices, triangles, parameters, major_radius, minor_radius, &
            wave_number, impedance, quadrature_degree, 3.0e-4_dp, 1, 0.12_dp, &
            dtn_map, status)
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
                ! Plot a real-valued instantaneous vector slice while retaining
                ! the complex magnitude used by the scattering observable.
                magnetic_vector_map(:, x_index, z_index) = &
                    real(magnetic_field, dp)
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
        real(dp) :: surface_x(visual_polar + 1, visual_azimuth + 1)
        real(dp) :: surface_y(visual_polar + 1, visual_azimuth + 1)
        real(dp) :: surface_z(visual_polar + 1, visual_azimuth + 1)
        real(dp) :: radiation_x(visual_polar + 1, visual_azimuth + 1)
        real(dp) :: radiation_y(visual_polar + 1, visual_azimuth + 1)
        real(dp) :: radiation_z(visual_polar + 1, visual_azimuth + 1)
        real(dp) :: curve_azimuth, curve_polar, curve_radius
        real(dp) :: normalized_radius, maximum_rcs
        real(dp) :: arrow_x(field_cells*field_cells)
        real(dp) :: arrow_z(field_cells*field_cells)
        real(dp) :: arrow_u(field_cells*field_cells)
        real(dp) :: arrow_v(field_cells*field_cells)
        real(dp) :: arrow_scale, arrow_norm
        integer :: coefficient, curve_index, point_index, source_azimuth
        integer :: source_polar, arrow_count, ix, iz

        arrow_count = 0
        arrow_norm = maxval(sqrt( &
            magnetic_vector_map(1, :, :)**2 + magnetic_vector_map(3, :, :)**2))
        arrow_scale = 0.16_dp/max(arrow_norm, tiny(1.0_dp))
        do iz = 1, field_cells
            do ix = 1, field_cells
                arrow_count = arrow_count + 1
                arrow_x(arrow_count) = 0.5_dp*(field_x_edges(ix) + &
                    field_x_edges(ix + 1))
                arrow_z(arrow_count) = 0.5_dp*(field_z_edges(iz) + &
                    field_z_edges(iz + 1))
                arrow_u(arrow_count) = arrow_scale* &
                    magnetic_vector_map(1, ix, iz)
                arrow_v(arrow_count) = arrow_scale* &
                    magnetic_vector_map(3, ix, iz)
            end do
        end do

        call figure(figsize=[8.5_dp, 6.5_dp])
        call pcolormesh( &
            field_x_edges, field_z_edges, magnetic_map, cmap="inferno")
        call quiver(arrow_x, arrow_z, arrow_u, arrow_v, scale=1.0_dp, &
            scale_units="xy", angles="xy", color="white", width=0.0025_dp, &
            headwidth=3.0_dp)
        call colorbar(label="|H_scattered|")
        call xlabel("x at y = 0")
        call ylabel("z")
        call title("Computed torus Maxwell scattered magnetic field and vectors")
        call savefig(output_directory//"/maxwell_torus_solution_2d.png")

        call figure(figsize=[7.5_dp, 6.5_dp])
        do curve_index = 0, visual_polar
            curve_polar = 2.0_dp*pi*real(curve_index, dp)/ &
                real(visual_polar, dp)
            do point_index = 0, visual_azimuth
                curve_azimuth = 2.0_dp*pi*real(point_index, dp)/ &
                    real(visual_azimuth, dp)
                curve_radius = major_radius + minor_radius*cos(curve_polar)
                surface_x(curve_index + 1, point_index + 1) = &
                    curve_radius*cos(curve_azimuth)
                surface_y(curve_index + 1, point_index + 1) = &
                    curve_radius*sin(curve_azimuth)
                surface_z(curve_index + 1, point_index + 1) = &
                    minor_radius*sin(curve_polar)
            end do
        end do
        call add_parametric_surface( &
            surface_x, surface_y, surface_z, color="lightsteelblue", alpha=0.55_dp, &
            linewidth=0.0_dp, filled=.true., label="PEC torus surface")
        call add_scatter( &
            vertices(1, :), vertices(2, :), vertices(3, :), &
            label="curved CFIE mesh vertices", marker=".", markersize=5.0_dp)
        call title("Exact-curved PEC torus and CFIE mesh")
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

        ! The far-field radius is a physical solution surface, not an
        ! unordered point cloud.  Use the same periodic tensor-product
        ! topology as the angular samples so FortPlot can render connected
        ! quads and preserve the toroidal seam in the 3-D view.
        maximum_rcs = max(maxval(rcs_map), tiny(1.0_dp))
        do curve_index = 0, visual_polar
            source_polar = min(polar_cells, max(1, curve_index))
            curve_polar = pi*real(curve_index, dp)/real(visual_polar, dp)
            do point_index = 0, visual_azimuth
                source_azimuth = 1 + mod(point_index, azimuth_cells)
                curve_azimuth = 2.0_dp*pi*real(point_index, dp)/ &
                    real(visual_azimuth, dp)
                normalized_radius = max( &
                    rcs_map(source_azimuth, source_polar)/maximum_rcs, 0.03_dp)
                radiation_x(curve_index + 1, point_index + 1) = &
                    normalized_radius*sin(curve_polar)*cos(curve_azimuth)
                radiation_y(curve_index + 1, point_index + 1) = &
                    normalized_radius*sin(curve_polar)*sin(curve_azimuth)
                radiation_z(curve_index + 1, point_index + 1) = &
                    normalized_radius*cos(curve_polar)
            end do
        end do
        call figure(figsize=[7.5_dp, 6.5_dp])
        call add_parametric_surface( &
            radiation_x, radiation_y, radiation_z, color="darkorange", &
            alpha=0.82_dp, linewidth=0.25_dp, filled=.true., &
            label="normalized RCS surface")
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

    subroutine write_gallery_sequence()
        integer :: sequence_unit

        open (newunit=sequence_unit, &
            file=output_directory//"/gallery_sequence.txt", &
            status="replace", action="write")
        write (sequence_unit, "(a)") "physical_solution"
        write (sequence_unit, "(a)") "diagnostics"
        close (sequence_unit)
    end subroutine write_gallery_sequence

    subroutine write_scattered_field_csv()
        integer :: csv_unit, x_index, z_index, local_status

        open (newunit=csv_unit, file=output_directory//"/scattered_field.csv", &
            status="replace", action="write", iostat=local_status)
        if (local_status /= 0) error stop "cannot open torus scattering CSV"
        write (csv_unit, "(a)", iostat=local_status) &
            "x,z,Hx_real,Hy_real,Hz_real,H_scattered_magnitude"
        if (local_status /= 0) error stop "cannot write torus scattering CSV"
        do z_index = 1, field_cells
            do x_index = 1, field_cells
                write (csv_unit, "(*(es24.16,:,','))", iostat=local_status) &
                    0.5_dp*(field_x_edges(x_index) + field_x_edges(x_index + 1)), &
                    0.5_dp*(field_z_edges(z_index) + field_z_edges(z_index + 1)), &
                    magnetic_vector_map(1, x_index, z_index), &
                    magnetic_vector_map(2, x_index, z_index), &
                    magnetic_vector_map(3, x_index, z_index), &
                    magnetic_map(x_index, z_index)
                if (local_status /= 0) &
                    error stop "cannot write torus scattering sample"
            end do
        end do
        close (csv_unit, iostat=local_status)
        if (local_status /= 0) error stop "cannot close torus scattering CSV"
    end subroutine write_scattered_field_csv

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
        write (unit, "(a,i0)") "CFIE quadrature degree: ", quadrature_degree
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
        if (.not. fast_gallery .and. reciprocity_error >= 2.0e-1_dp) &
            error stop "toroidal Maxwell gallery reciprocity regression"
    end subroutine write_outputs

end program maxwell_torus_curved_scattering
