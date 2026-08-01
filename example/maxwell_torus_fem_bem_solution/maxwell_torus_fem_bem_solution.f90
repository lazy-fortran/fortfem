program maxwell_torus_fem_bem_solution
    !! Manufactured solved-state check for toroidal Maxwell FEM--BEM coupling.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortfem_api, only: &
        assemble_maxwell_fem_bem_torus_curved_system_3d, &
        build_tetra_nedelec_dof_map, evaluate_tetra_nedelec_first_kind, &
        generate_solid_torus_tetra_mesh, initialize_tetra_nedelec_first_kind, &
        invert_tetra_affine_map, map_tetra_nedelec_covariant, &
        solve_maxwell_fem_bem_torus_curved_system_3d, &
        tetra_nedelec_first_kind_t
    use fortfem_kinds, only: dp
    use fortplot, only: &
        add_3d_plot, add_scatter, colorbar, figure, legend, pcolormesh, &
        plot, quiver, savefig, title, xlabel, ylabel
    implicit none

    real(dp), parameter :: major_radius = 2.0_dp
    real(dp), parameter :: minor_radius = 0.65_dp
    real(dp), parameter :: curl_coefficient = 1.1_dp
    real(dp), parameter :: mass_coefficient = -0.7_dp
    real(dp), parameter :: wave_number = 0.8_dp
    real(dp), parameter :: impedance = 1.4_dp
    real(dp), parameter :: electric_field(3) = [0.7_dp, -0.4_dp, 0.2_dp]
    real(dp), parameter :: blue(3) = [0.0_dp, 114.0_dp, 178.0_dp]/255.0_dp
    real(dp), parameter :: orange(3) = [230.0_dp, 159.0_dp, 0.0_dp]/255.0_dp
    real(dp), parameter :: white(3) = [1.0_dp, 1.0_dp, 1.0_dp]
    real(dp), parameter :: slice_y = 0.017_dp
    integer, parameter :: slice_nx = 42, slice_nz = 30
    integer, parameter :: arrow_stride = 4
    integer, parameter :: geometry_polar = 16, geometry_toroidal = 24
    character(*), parameter :: output_directory = &
        "output/example/maxwell_torus_fem_bem_solution"
    integer, allocatable :: boundary_triangles(:, :), edges(:, :), faces(:, :)
    integer, allocatable :: edge_orientations(:, :), global_dofs(:, :)
    integer, allocatable :: face_permutations(:, :, :), tetrahedra(:, :)
    real(dp), allocatable :: parameters(:, :), vertices(:, :)
    complex(dp), allocatable :: system_matrix(:, :), right_hand_side(:)
    complex(dp), allocatable :: exact_state(:), field(:), current(:)
    real(dp) :: x_edges(slice_nx + 1), z_edges(slice_nz + 1)
    real(dp) :: field_map(slice_nx, slice_nz), exact_map(slice_nx, slice_nz)
    real(dp) :: error_map(slice_nx, slice_nz)
    real(dp) :: arrow_x(slice_nx*slice_nz), arrow_z(slice_nx*slice_nz)
    real(dp) :: arrow_u(slice_nx*slice_nz), arrow_v(slice_nx*slice_nz)
    real(dp) :: centre_x(slice_nx), computed_centre(slice_nx)
    real(dp) :: exact_centre(slice_nx), field_error, current_error
    real(dp) :: start_time, elapsed, x, z, value, exact_value, arrow_norm
    integer :: arrow_count, field_count, info, ix, iz, status, unit
    real(dp) :: sample_vector(3)

    call execute_command_line("mkdir -p "//output_directory)
    call generate_solid_torus_tetra_mesh( &
        major_radius, minor_radius, 1, 6, 4, vertices, tetrahedra, &
        boundary_triangles, parameters)
    call build_tetra_nedelec_dof_map( &
        1, tetrahedra, edges, faces, global_dofs, edge_orientations, &
        face_permutations, status)
    if (status /= 0) error stop "toroidal FEM--BEM Nedelec map failed"
    field_count = maxval(global_dofs)
    call assemble_maxwell_fem_bem_torus_curved_system_3d( &
        vertices, tetrahedra, parameters, boundary_triangles, major_radius, &
        minor_radius, curl_coefficient, mass_coefficient, wave_number, &
        impedance, 3, 1.0e-5_dp, 1, system_matrix, status)
    if (status /= 0) error stop "toroidal FEM--BEM system assembly failed"
    allocate(exact_state(size(system_matrix, 1)), right_hand_side(size(system_matrix, 1)))
    exact_state = cmplx(0.0_dp, 0.0_dp, dp)
    do info = 1, size(edges, 2)
        exact_state(info) = cmplx(dot_product( &
            electric_field, vertices(:, edges(2, info)) - &
            vertices(:, edges(1, info))), 0.0_dp, dp)
    end do
    right_hand_side = matmul(system_matrix, exact_state)

    call cpu_time(start_time)
    call solve_maxwell_fem_bem_torus_curved_system_3d( &
        vertices, tetrahedra, parameters, boundary_triangles, major_radius, &
        minor_radius, curl_coefficient, mass_coefficient, wave_number, &
        impedance, 3, 1.0e-5_dp, 1, right_hand_side, field, current, status)
    call cpu_time(elapsed)
    elapsed = elapsed - start_time
    if (status /= 0) error stop "toroidal FEM--BEM manufactured solve failed"
    field_error = maxval(abs(field - exact_state(:field_count)))
    current_error = maxval(abs(current))
    if (field_error > 3.0e-8_dp .or. current_error > 3.0e-8_dp) &
        error stop "toroidal FEM--BEM manufactured solution oracle failed"

    do ix = 1, slice_nx + 1
        x_edges(ix) = 1.25_dp + 1.5_dp*real(ix - 1, dp)/real(slice_nx, dp)
    end do
    do iz = 1, slice_nz + 1
        z_edges(iz) = -0.75_dp + 1.5_dp*real(iz - 1, dp)/real(slice_nz, dp)
    end do
    field_map = 0.0_dp
    exact_map = 0.0_dp
    error_map = 0.0_dp
    arrow_count = 0
    do iz = 1, slice_nz
        z = 0.5_dp*(z_edges(iz) + z_edges(iz + 1))
        do ix = 1, slice_nx
            x = 0.5_dp*(x_edges(ix) + x_edges(ix + 1))
            call evaluate_field([x, slice_y, z], value, status)
            if (status == 0 .and. ieee_is_finite(value)) then
                exact_value = sqrt(sum(electric_field**2))
                field_map(ix, iz) = value
                exact_map(ix, iz) = exact_value
                if (mod(ix - 1, arrow_stride) == 0 .and. &
                    mod(iz - 1, arrow_stride) == 0) then
                    arrow_count = arrow_count + 1
                    arrow_x(arrow_count) = x
                    arrow_z(arrow_count) = z
                    call evaluate_vector([x, slice_y, z], sample_vector, status)
                    if (status /= 0) error stop "toroidal vector sample failed"
                    arrow_u(arrow_count) = sample_vector(1)
                    arrow_v(arrow_count) = sample_vector(3)
                end if
            end if
        end do
    end do
    ! The manufactured solve is at round-off, so a linear error scale would
    ! collapse to a visually uninformative zero-valued colour bar.  Plot the
    ! resolved error on a guarded log scale while retaining the raw maximum
    ! in benchmark.txt and the independent test oracle.
    error_map = log10(max(abs(field_map - exact_map), 1.0e-16_dp))
    if (arrow_count > 0) then
        arrow_norm = maxval(sqrt(arrow_u(:arrow_count)**2 + &
            arrow_v(:arrow_count)**2))
        arrow_norm = max(arrow_norm, epsilon(1.0_dp))
        arrow_u(:arrow_count) = 0.16_dp*arrow_u(:arrow_count)/arrow_norm
        arrow_v(:arrow_count) = 0.16_dp*arrow_v(:arrow_count)/arrow_norm
    end if
    do ix = 1, slice_nx
        centre_x(ix) = 1.4_dp + 1.2_dp*real(ix - 1, dp)/ &
            real(slice_nx - 1, dp)
        call evaluate_field([centre_x(ix), slice_y, 0.0_dp], &
            computed_centre(ix), status)
        if (status /= 0) error stop "toroidal centerline sample failed"
        exact_centre(ix) = sqrt(sum(electric_field**2))
    end do

    call render_plots()
    open (newunit=unit, file=output_directory//"/benchmark.txt", &
        status="replace", action="write")
    write (unit, "(a)") "quantity,value"
    write (unit, "(a,i0)") "volume_tetrahedra,", size(tetrahedra, 2)
    write (unit, "(a,i0)") "surface_triangles,", size(boundary_triangles, 2)
    write (unit, "(a,i0)") "nedelec_unknowns,", field_count
    write (unit, "(a,i0)") "rwg_unknowns,", size(current)
    write (unit, "(a,es24.16)") "maximum_field_error,", field_error
    write (unit, "(a,es24.16)") "maximum_current_error,", current_error
    write (unit, "(a,es24.16)") "solve_seconds,", elapsed
    close (unit)

contains

    subroutine evaluate_field(point, magnitude, local_status)
        real(dp), intent(in) :: point(3)
        real(dp), intent(out) :: magnitude
        integer, intent(out) :: local_status

        real(dp) :: vector(3)

        call evaluate_vector(point, vector, local_status)
        if (local_status == 0) then
            magnitude = sqrt(sum(vector**2))
        else
            magnitude = 0.0_dp
        end if
    end subroutine evaluate_field

    subroutine evaluate_vector(point, vector, local_status)
        real(dp), intent(in) :: point(3)
        real(dp), intent(out) :: vector(3)
        integer, intent(out) :: local_status

        type(tetra_nedelec_first_kind_t) :: basis
        real(dp) :: local_vertices(3, 4), reference(3), jacobian(3, 3)
        real(dp) :: reference_values(3, 6), reference_curls(3, 6)
        real(dp) :: physical_values(3, 6), physical_curls(3, 6)
        real(dp) :: local_dofs(6)
        integer :: basis_status, edge, map_status, tetrahedron

        vector = 0.0_dp
        local_status = 1
        call initialize_tetra_nedelec_first_kind(1, basis, basis_status)
        if (basis_status /= 0) return
        do tetrahedron = 1, size(tetrahedra, 2)
            local_vertices = vertices(:, tetrahedra(:, tetrahedron))
            call invert_tetra_affine_map( &
                local_vertices, point, reference, map_status)
            if (map_status /= 0 .or. any(reference < -1.0e-9_dp) .or. &
                sum(reference) > 1.0_dp + 1.0e-9_dp) cycle
            call evaluate_tetra_nedelec_first_kind( &
                basis, reference, reference_values, reference_curls, basis_status)
            if (basis_status /= 0) return
            do edge = 1, 6
                local_dofs(edge) = real(edge_orientations(edge, tetrahedron), dp)* &
                    real(field(global_dofs(edge, tetrahedron)), dp)
            end do
            jacobian(:, 1) = local_vertices(:, 2) - local_vertices(:, 1)
            jacobian(:, 2) = local_vertices(:, 3) - local_vertices(:, 1)
            jacobian(:, 3) = local_vertices(:, 4) - local_vertices(:, 1)
            call map_tetra_nedelec_covariant( &
                jacobian, reference_values, reference_curls, physical_values, &
                physical_curls, basis_status)
            if (basis_status /= 0) return
            vector = matmul(physical_values, local_dofs)
            local_status = 0
            return
        end do
    end subroutine evaluate_vector

    subroutine render_plots()
        real(dp) :: curve_x(geometry_toroidal + 1)
        real(dp) :: curve_y(geometry_toroidal + 1)
        real(dp) :: curve_z(geometry_toroidal + 1)
        real(dp) :: section_x(129), section_z(129)
        real(dp) :: sample_x(slice_nx*slice_nz)
        real(dp) :: sample_z(slice_nx*slice_nz)
        real(dp) :: sample_value(slice_nx*slice_nz)
        real(dp) :: theta, phi, radius
        integer :: curve, point_index, sample_count, ix_local, iz_local

        do point_index = 0, 128
            theta = 2.0_dp*acos(-1.0_dp)*real(point_index, dp)/128.0_dp
            section_x(point_index + 1) = major_radius + minor_radius*cos(theta)
            section_z(point_index + 1) = minor_radius*sin(theta)
        end do

        call figure(figsize=[8.5_dp, 6.5_dp])
        call pcolormesh(x_edges, z_edges, transpose(field_map), cmap="viridis")
        call colorbar(label="|E_h|")
        call plot(section_x, section_z, label="exact torus section", &
            color=white, linewidth=1.2_dp)
        call xlabel("x at y = 0.017")
        call ylabel("z")
        call title("Solved toroidal Maxwell FEM--BEM field magnitude")
        call legend()
        call savefig(output_directory//"/maxwell_torus_fem_bem_solution_2d.png")

        sample_count = 0
        do iz_local = 1, slice_nz
            do ix_local = 1, slice_nx
                sample_count = sample_count + 1
                sample_x(sample_count) = &
                    0.5_dp*(x_edges(ix_local) + x_edges(ix_local + 1))
                sample_z(sample_count) = &
                    0.5_dp*(z_edges(iz_local) + z_edges(iz_local + 1))
                sample_value(sample_count) = field_map(ix_local, iz_local)
            end do
        end do
        call figure(figsize=[8.5_dp, 6.5_dp])
        call add_scatter(sample_x(:sample_count), sample_z(:sample_count), &
            c=sample_value(:sample_count), cmap="viridis", marker=".", &
            markersize=10.0_dp)
        call plot(section_x, section_z, color=white, linewidth=1.2_dp)
        if (arrow_count > 0) call quiver( &
            arrow_x(:arrow_count), arrow_z(:arrow_count), &
            arrow_u(:arrow_count), arrow_v(:arrow_count), color="black")
        call colorbar(label="|E_h|")
        call xlabel("x at y = 0.017")
        call ylabel("z")
        call title("Toroidal FEM--BEM solved H(curl) vector field")
        call savefig(output_directory//"/maxwell_torus_fem_bem_vector_2d.png")

        call figure(figsize=[8.5_dp, 6.5_dp])
        call pcolormesh(x_edges, z_edges, transpose(error_map), cmap="magma")
        call colorbar(label="log10 absolute field error")
        call plot(section_x, section_z, color=white, linewidth=1.2_dp)
        call xlabel("x at y = 0.017")
        call ylabel("z")
        call title("Toroidal Maxwell FEM--BEM round-off error (log10)")
        call savefig(output_directory//"/maxwell_torus_fem_bem_error_2d.png")

        call figure(figsize=[8.5_dp, 5.5_dp])
        call add_scatter(centre_x, computed_centre, label="computed |E_h|", &
            color=blue, marker=".", markersize=5.0_dp)
        call add_scatter(centre_x, exact_centre, label="constant analytical field", &
            color=orange, marker="x", markersize=5.0_dp)
        call xlabel("x at y = 0.017, z = 0")
        call ylabel("field magnitude")
        call title("Toroidal FEM--BEM cross-section oracle")
        call legend()
        call savefig(output_directory//"/maxwell_torus_fem_bem_centerline_1d.png")

        call figure(figsize=[7.5_dp, 6.5_dp])
        call add_scatter(vertices(1, :), vertices(2, :), vertices(3, :), &
            label="solid-torus volume vertices", marker=".", markersize=4.0_dp)
        do curve = 0, geometry_polar
            theta = 2.0_dp*acos(-1.0_dp)*real(curve, dp)/real(geometry_polar, dp)
            do point_index = 0, geometry_toroidal
                phi = 2.0_dp*acos(-1.0_dp)*real(point_index, dp)/ &
                    real(geometry_toroidal, dp)
                radius = major_radius + minor_radius*cos(theta)
                curve_x(point_index + 1) = radius*cos(phi)
                curve_y(point_index + 1) = radius*sin(phi)
                curve_z(point_index + 1) = minor_radius*sin(theta)
            end do
            call add_3d_plot(curve_x, curve_y, curve_z, &
                label="exact curved torus", color="teal", linewidth=0.7_dp)
        end do
        call title("Curved torus volume and FEM--BEM boundary geometry")
        call savefig(output_directory//"/maxwell_torus_fem_bem_geometry_3d.png")
    end subroutine render_plots

end program maxwell_torus_fem_bem_solution
