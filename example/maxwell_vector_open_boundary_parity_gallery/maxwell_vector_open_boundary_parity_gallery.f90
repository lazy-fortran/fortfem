program maxwell_vector_open_boundary_parity_gallery
    !! Physical-first gallery for vector curl-curl open-boundary parity.
    use fortfem_boundary, only: &
        apply_maxwell_trace_to_flux_map, &
        assemble_maxwell_fem_bem_torus_curved_system_3d, &
        assemble_maxwell_torus_curved_dtn_rwg_3d, &
        solve_maxwell_fem_bem_torus_curved_system_3d
    use fortfem_core, only: &
        generate_solid_torus_tetra_mesh, generate_structured_tetra_box_mesh, &
        generate_torus_surface_mesh, invert_tetra_affine_map
    use fortfem_feec, only: &
        build_tetra_edge_dof_map, build_tetra_nedelec_dof_map, &
        evaluate_tetra_nedelec_interpolant_at_point, &
        initialize_tetra_nedelec_first_kind, solve_tetra_nedelec_pml, &
        tetra_nedelec_first_kind_t
    use fortfem_kinds, only: dp
    use fortplot, only: &
        add_parametric_surface, add_scatter, colorbar, figure, legend, pcolormesh, &
        plot, quiver, savefig, title, xlabel, ylabel
    implicit none

    integer, parameter :: slice_nx = 32, slice_nz = 24
    integer, parameter :: pml_nx = 24, pml_ny = 24
    integer, parameter :: visual_polar = 20, visual_toroidal = 28
    real(dp), parameter :: major_radius = 2.0_dp, minor_radius = 0.6_dp
    real(dp), parameter :: wave_number = 1.1_dp, pml_imaginary = 0.28_dp
    real(dp), parameter :: vector_real(3) = [0.48_dp, -0.31_dp, 0.22_dp]
    real(dp), parameter :: vector_imag(3) = [-0.17_dp, 0.26_dp, 0.11_dp]
    real(dp), parameter :: rotation(3) = [0.0_dp, 0.0_dp, 0.27_dp]
    real(dp), parameter :: black(3) = [0.0_dp, 0.0_dp, 0.0_dp]
    real(dp), parameter :: white(3) = [1.0_dp, 1.0_dp, 1.0_dp]
    character(*), parameter :: output_directory = &
        "output/example/maxwell_vector_open_boundary_parity_gallery"
    real(dp), allocatable :: torus_vertices(:, :), torus_parameters(:, :)
    integer, allocatable :: torus_tetrahedra(:, :), torus_boundary(:, :)
    integer, allocatable :: torus_edges(:, :), torus_faces(:, :)
    integer, allocatable :: torus_global_dofs(:, :), torus_orientations(:, :)
    integer, allocatable :: torus_permutations(:, :, :)
    complex(dp), allocatable :: torus_field(:), torus_current(:)
    complex(dp), allocatable :: dtn_map(:, :), dtn_trace(:), dtn_flux(:)
    real(dp), allocatable :: small_vertices(:, :), large_vertices(:, :)
    integer, allocatable :: small_tetrahedra(:, :), large_tetrahedra(:, :)
    integer, allocatable :: small_edges(:, :), large_edges(:, :)
    integer, allocatable :: small_dofs(:, :), large_dofs(:, :)
    integer, allocatable :: small_orientations(:, :), large_orientations(:, :)
    integer, allocatable :: small_boundary(:), large_boundary(:)
    complex(dp), allocatable :: small_solution(:), large_solution(:)
    real(dp) :: torus_x_edges(slice_nx + 1), torus_z_edges(slice_nz + 1)
    real(dp) :: torus_magnitude(slice_nx, slice_nz)
    real(dp) :: torus_arrow_x(slice_nx*slice_nz), torus_arrow_z(slice_nx*slice_nz)
    real(dp) :: torus_arrow_u(slice_nx*slice_nz), torus_arrow_v(slice_nx*slice_nz)
    real(dp) :: pml_x_edges(pml_nx + 1), pml_y_edges(pml_ny + 1)
    real(dp) :: pml_magnitude(pml_nx, pml_ny)
    real(dp) :: pml_arrow_x(pml_nx*pml_ny), pml_arrow_y(pml_nx*pml_ny)
    real(dp) :: pml_arrow_u(pml_nx*pml_ny), pml_arrow_v(pml_nx*pml_ny)
    real(dp) :: domain_x(17)
    complex(dp) :: domain_small(17), domain_large(17)
    real(dp) :: dtn_seconds, fem_bem_seconds, pml_seconds, start_time, end_time
    real(dp) :: field_error, current_norm, dtn_norm, domain_difference
    integer :: arrow_count, pml_arrow_count, status, unit
    integer :: torus_poloidal_count, torus_toroidal_count
    logical :: fast_gallery

    call execute_command_line("mkdir -p "//output_directory)
    fast_gallery = environment_flag("MAXWELL_VECTOR_OPEN_BOUNDARY_FAST")
    torus_poloidal_count = 6
    torus_toroidal_count = 4
    if (fast_gallery) then
        torus_poloidal_count = 3
        torus_toroidal_count = 3
    end if
    call solve_torus_fem_bem(fem_bem_seconds, field_error, current_norm)
    call assemble_dtn(dtn_seconds, dtn_norm)
    call solve_pml_boxes(pml_seconds, domain_difference)
    call build_torus_slice()
    call build_pml_slice()

    open (newunit=unit, file=output_directory//"/gallery_sequence.txt", &
        status="replace", action="write")
    close (unit)

    ! Keep the physical fields first in the generated sequence and on disk.
    call render_torus_solution()
    call render_pml_solution()
    call render_torus_geometry()
    call write_torus_field_csv()
    call write_pml_field_csv()
    call write_sequence("physical_solution")

    call render_diagnostics()
    call write_sequence("diagnostics")
    call write_outputs(fem_bem_seconds, dtn_seconds, pml_seconds, field_error, &
        current_norm, dtn_norm, domain_difference)

contains

    subroutine solve_torus_fem_bem(seconds, error, surface_norm)
        real(dp), intent(out) :: seconds, error, surface_norm
        complex(dp), allocatable :: matrix(:, :), rhs(:), exact(:), solved(:)
        real(dp) :: midpoint(3), delta(3)
        integer :: field_count, edge, local_status

        call generate_solid_torus_tetra_mesh( &
            major_radius, minor_radius, 1, torus_poloidal_count, torus_toroidal_count, &
            torus_vertices, torus_tetrahedra, &
            torus_boundary, torus_parameters)
        call build_tetra_nedelec_dof_map( &
            1, torus_tetrahedra, torus_edges, torus_faces, torus_global_dofs, &
            torus_orientations, torus_permutations, local_status)
        if (local_status /= 0) error stop "vector parity torus edge map failed"
        field_count = maxval(torus_global_dofs)
        call cpu_time(start_time)
        call assemble_maxwell_fem_bem_torus_curved_system_3d( &
            torus_vertices, torus_tetrahedra, torus_parameters, torus_boundary, &
            major_radius, minor_radius, 1.2_dp, -0.6_dp, 0.65_dp, 1.4_dp, &
            3, 1.0e-5_dp, 1, matrix, local_status)
        if (local_status /= 0) error stop "vector parity FEM--BEM assembly failed"
        allocate(exact(size(matrix, 1)), rhs(size(matrix, 1)))
        exact = cmplx(0.0_dp, 0.0_dp, dp)
        do edge = 1, size(torus_edges, 2)
            delta = torus_vertices(:, torus_edges(2, edge)) - &
                torus_vertices(:, torus_edges(1, edge))
            midpoint = 0.5_dp*(torus_vertices(:, torus_edges(2, edge)) + &
                torus_vertices(:, torus_edges(1, edge)))
            exact(edge) = dot_product(manufactured_vector(midpoint), delta)
        end do
        rhs = matmul(matrix, exact)
        call solve_maxwell_fem_bem_torus_curved_system_3d( &
            torus_vertices, torus_tetrahedra, torus_parameters, torus_boundary, &
            major_radius, minor_radius, 1.2_dp, -0.6_dp, 0.65_dp, 1.4_dp, &
            3, 1.0e-5_dp, 1, rhs, solved, torus_current, local_status)
        call cpu_time(end_time)
        if (local_status /= 0) error stop "vector parity FEM--BEM solve failed"
        torus_field = solved
        seconds = end_time - start_time
        error = maxval(abs(solved - exact(:field_count)))
        surface_norm = sqrt(sum(abs(torus_current)**2))
    end subroutine solve_torus_fem_bem

    subroutine assemble_dtn(seconds, response_norm)
        real(dp), intent(out) :: seconds, response_norm
        real(dp), allocatable :: surface_vertices(:, :), surface_parameters(:, :)
        integer, allocatable :: surface_triangles(:, :)
        integer :: coefficient, local_status

        call generate_torus_surface_mesh( &
            major_radius, minor_radius, 3, 3, surface_vertices, surface_triangles, &
            surface_parameters)
        call cpu_time(start_time)
        call assemble_maxwell_torus_curved_dtn_rwg_3d( &
            surface_vertices, surface_triangles, surface_parameters, major_radius, &
            minor_radius, 0.45_dp, 1.7_dp, 3, 3.0e-4_dp, 1, 0.12_dp, dtn_map, &
            local_status)
        if (local_status /= 0) error stop "vector parity DtN assembly failed"
        allocate(dtn_trace(size(dtn_map, 2)), dtn_flux(size(dtn_map, 1)))
        do coefficient = 1, size(dtn_trace)
            dtn_trace(coefficient) = cmplx( &
                cos(0.19_dp*real(coefficient, dp)), &
                sin(0.13_dp*real(coefficient, dp)), dp)
        end do
        call apply_maxwell_trace_to_flux_map( &
            dtn_map, dtn_trace, dtn_flux, local_status)
        call cpu_time(end_time)
        if (local_status /= 0) error stop "vector parity DtN action failed"
        seconds = end_time - start_time
        response_norm = sqrt(sum(abs(dtn_flux)**2))
    end subroutine assemble_dtn

    subroutine solve_pml_boxes(seconds, difference)
        real(dp), intent(out) :: seconds, difference
        real(dp) :: small_bounds(3, 2), large_bounds(3, 2)
        integer :: local_status, sample
        complex(dp) :: small_field(3), large_field(3), expected(3)
        real(dp) :: error

        small_bounds(:, 1) = 0.0_dp
        small_bounds(:, 2) = 1.0_dp
        large_bounds(:, 1) = 0.0_dp
        large_bounds(:, 2) = [1.0_dp, 1.0_dp, 1.5_dp]
        call cpu_time(start_time)
        call solve_pml_box(small_bounds, [2, 2, 2], small_vertices, &
            small_tetrahedra, small_edges, small_dofs, small_orientations, &
            small_boundary, small_solution)
        call solve_pml_box(large_bounds, [2, 2, 3], large_vertices, &
            large_tetrahedra, large_edges, large_dofs, large_orientations, &
            large_boundary, large_solution)
        call cpu_time(end_time)
        seconds = end_time - start_time
        error = 0.0_dp
        difference = 0.0_dp
        do sample = 1, 3
            call evaluate_complex_field( &
                small_vertices, small_tetrahedra, small_dofs, small_orientations, &
                small_solution, target_point(sample), small_field, local_status)
            call evaluate_complex_field( &
                large_vertices, large_tetrahedra, large_dofs, large_orientations, &
                large_solution, target_point(sample), large_field, local_status)
            expected = plane_wave_field(target_point(sample))
            error = max(error, sqrt(sum(abs(small_field - expected)**2)))
            error = max(error, sqrt(sum(abs(large_field - expected)**2)))
            difference = max(difference, sqrt(sum(abs(small_field - large_field)**2)))
        end do
        if (error > 3.0e-1_dp .or. difference > 1.5e-1_dp) &
            error stop "vector parity PML larger-domain regression"
    end subroutine solve_pml_boxes

    subroutine solve_pml_box(bounds, counts, box_vertices, box_tetrahedra, &
            box_edges, box_dofs, box_orientations, box_boundary, solution)
        real(dp), intent(in) :: bounds(3, 2)
        integer, intent(in) :: counts(3)
        real(dp), allocatable, intent(out) :: box_vertices(:, :)
        integer, allocatable, intent(out) :: box_tetrahedra(:, :), box_edges(:, :)
        integer, allocatable, intent(out) :: box_dofs(:, :), box_orientations(:, :)
        integer, allocatable, intent(out) :: box_boundary(:)
        complex(dp), allocatable, intent(out) :: solution(:)
        complex(dp), allocatable :: stretch(:, :), load(:), exact(:)
        integer :: edge, local_status

        call generate_structured_tetra_box_mesh( &
            bounds, counts, box_vertices, box_tetrahedra, local_status)
        if (local_status /= 0) error stop "vector parity PML mesh failed"
        call build_tetra_edge_dof_map( &
            box_tetrahedra, box_edges, box_dofs, box_orientations, local_status)
        if (local_status /= 0) error stop "vector parity PML edge map failed"
        allocate(exact(size(box_edges, 2)), load(size(box_edges, 2)))
        do edge = 1, size(box_edges, 2)
            exact(edge) = plane_wave_edge( &
                box_vertices(:, box_edges(1, edge)), &
                box_vertices(:, box_edges(2, edge)))
        end do
        call boundary_edges(box_vertices, box_edges, bounds, box_boundary)
        allocate(stretch(3, size(box_tetrahedra, 2)))
        stretch = cmplx(1.0_dp, 0.0_dp, dp)
        stretch(1, :) = cmplx(1.0_dp, pml_imaginary, dp)
        load = cmplx(0.0_dp, 0.0_dp, dp)
        call solve_tetra_nedelec_pml( &
            box_vertices, box_tetrahedra, 1, stretch, wave_number, load, &
            box_boundary, exact(box_boundary), solution, local_status)
        if (local_status /= 0) error stop "vector parity PML solve failed"
    end subroutine solve_pml_box

    subroutine boundary_edges(box_vertices, box_edges, bounds, boundary)
        real(dp), intent(in) :: box_vertices(:, :), bounds(3, 2)
        integer, intent(in) :: box_edges(:, :)
        integer, allocatable, intent(out) :: boundary(:)
        logical, allocatable :: marked(:)
        integer :: coordinate, edge

        allocate(marked(size(box_edges, 2)))
        marked = .false.
        do edge = 1, size(box_edges, 2)
            do coordinate = 1, 3
                if (all(abs(box_vertices(coordinate, box_edges(:, edge)) - &
                    bounds(coordinate, 1)) < 1.0e-13_dp) .or. &
                    all(abs(box_vertices(coordinate, box_edges(:, edge)) - &
                    bounds(coordinate, 2)) < 1.0e-13_dp)) marked(edge) = .true.
            end do
        end do
        boundary = pack([(edge, edge=1, size(box_edges, 2))], marked)
    end subroutine boundary_edges

    subroutine build_torus_slice()
        integer :: ix, iz, local_status
        real(dp) :: x, z, vector(3)

        do ix = 1, slice_nx + 1
            torus_x_edges(ix) = 1.35_dp + 1.3_dp*real(ix - 1, dp)/real(slice_nx, dp)
        end do
        do iz = 1, slice_nz + 1
            torus_z_edges(iz) = -0.65_dp + 1.3_dp*real(iz - 1, dp)/real(slice_nz, dp)
        end do
        arrow_count = 0
        do iz = 1, slice_nz
            z = 0.5_dp*(torus_z_edges(iz) + torus_z_edges(iz + 1))
            do ix = 1, slice_nx
                x = 0.5_dp*(torus_x_edges(ix) + torus_x_edges(ix + 1))
                call evaluate_torus_field([x, 0.017_dp, z], vector, local_status)
                torus_magnitude(ix, iz) = 0.0_dp
                if (local_status == 0) torus_magnitude(ix, iz) = sqrt(sum(vector**2))
                if (local_status == 0 .and. mod(ix - 1, 4) == 0 .and. &
                    mod(iz - 1, 4) == 0) then
                    arrow_count = arrow_count + 1
                    torus_arrow_x(arrow_count) = x
                    torus_arrow_z(arrow_count) = z
                    torus_arrow_u(arrow_count) = vector(1)
                    torus_arrow_v(arrow_count) = vector(3)
                end if
            end do
        end do
        if (arrow_count > 0) then
            vector = 0.0_dp
            do ix = 1, arrow_count
                vector(1) = max(vector(1), abs(torus_arrow_u(ix)))
                vector(3) = max(vector(3), abs(torus_arrow_v(ix)))
            end do
            vector(1) = max(vector(1), vector(3))
            if (vector(1) > tiny(1.0_dp)) then
                torus_arrow_u(:arrow_count) = 0.16_dp * &
                    torus_arrow_u(:arrow_count)/vector(1)
                torus_arrow_v(:arrow_count) = 0.16_dp * &
                    torus_arrow_v(:arrow_count)/vector(1)
            end if
        end if
    end subroutine build_torus_slice

    subroutine build_pml_slice()
        integer :: ix, iy, local_status
        real(dp) :: x, y, arrow_norm
        complex(dp) :: vector(3)

        do ix = 1, pml_nx + 1
            pml_x_edges(ix) = real(ix - 1, dp)/real(pml_nx, dp)
        end do
        do iy = 1, pml_ny + 1
            pml_y_edges(iy) = real(iy - 1, dp)/real(pml_ny, dp)
        end do
        pml_arrow_count = 0
        do iy = 1, pml_ny
            y = 0.5_dp*(pml_y_edges(iy) + pml_y_edges(iy + 1))
            do ix = 1, pml_nx
                x = 0.5_dp*(pml_x_edges(ix) + pml_x_edges(ix + 1))
                call evaluate_complex_field( &
                    small_vertices, small_tetrahedra, small_dofs, small_orientations, &
                    small_solution, [x, y, 0.5_dp], vector, local_status)
                pml_magnitude(ix, iy) = 0.0_dp
                if (local_status == 0) pml_magnitude(ix, iy) = sqrt(sum(abs(vector)**2))
                if (local_status == 0 .and. mod(ix - 1, 4) == 0 .and. &
                    mod(iy - 1, 4) == 0) then
                    pml_arrow_count = pml_arrow_count + 1
                    pml_arrow_x(pml_arrow_count) = x
                    pml_arrow_y(pml_arrow_count) = y
                    pml_arrow_u(pml_arrow_count) = real(vector(1), dp)
                    pml_arrow_v(pml_arrow_count) = real(vector(2), dp)
                end if
            end do
        end do
        if (pml_arrow_count > 0) then
            arrow_norm = 0.0_dp
            do ix = 1, pml_arrow_count
                arrow_norm = max(arrow_norm, abs(pml_arrow_u(ix)))
                arrow_norm = max(arrow_norm, abs(pml_arrow_v(ix)))
            end do
            if (arrow_norm > tiny(1.0_dp)) then
                pml_arrow_u(:pml_arrow_count) = 0.16_dp * &
                    pml_arrow_u(:pml_arrow_count)/arrow_norm
                pml_arrow_v(:pml_arrow_count) = 0.16_dp * &
                    pml_arrow_v(:pml_arrow_count)/arrow_norm
            end if
        end if
    end subroutine build_pml_slice

    subroutine evaluate_torus_field(point, vector, local_status)
        real(dp), intent(in) :: point(3)
        real(dp), intent(out) :: vector(3)
        integer, intent(out) :: local_status
        type(tetra_nedelec_first_kind_t) :: basis
        real(dp) :: local_vertices(3, 4), reference(3), local_dofs(6)
        real(dp) :: values(3), curls(3)
        integer :: basis_status, map_status, edge, tetrahedron

        vector = 0.0_dp
        local_status = 1
        call initialize_tetra_nedelec_first_kind(1, basis, basis_status)
        if (basis_status /= 0) return
        do tetrahedron = 1, size(torus_tetrahedra, 2)
            local_vertices = torus_vertices(:, torus_tetrahedra(:, tetrahedron))
            call invert_tetra_affine_map(local_vertices, point, reference, map_status)
            if (map_status /= 0 .or. any(reference < -1.0e-9_dp) .or. &
                sum(reference) > 1.0_dp + 1.0e-9_dp) cycle
            do edge = 1, 6
                local_dofs(edge) = real(torus_orientations(edge, tetrahedron), dp)* &
                    real(torus_field(torus_global_dofs(edge, tetrahedron)), dp)
            end do
            call evaluate_tetra_nedelec_interpolant_at_point( &
                local_vertices, basis, local_dofs, point, values, curls, basis_status)
            if (basis_status /= 0) return
            vector = values
            local_status = 0
            return
        end do
    end subroutine evaluate_torus_field

    subroutine evaluate_complex_field(box_vertices, box_tetrahedra, box_dofs, &
            box_orientations, solution, point, field, local_status)
        real(dp), intent(in) :: box_vertices(:, :), point(3)
        integer, intent(in) :: box_tetrahedra(:, :), box_dofs(:, :), box_orientations(:, :)
        complex(dp), intent(in) :: solution(:)
        complex(dp), intent(out) :: field(3)
        integer, intent(out) :: local_status
        type(tetra_nedelec_first_kind_t) :: basis
        real(dp) :: local_vertices(3, 4), reference(3), real_dofs(6), imag_dofs(6)
        real(dp) :: real_value(3), imag_value(3), real_curl(3), imag_curl(3)
        integer :: basis_status, map_status, edge, tetrahedron

        field = cmplx(0.0_dp, 0.0_dp, dp)
        local_status = 1
        call initialize_tetra_nedelec_first_kind(1, basis, basis_status)
        if (basis_status /= 0) return
        do tetrahedron = 1, size(box_tetrahedra, 2)
            local_vertices = box_vertices(:, box_tetrahedra(:, tetrahedron))
            call invert_tetra_affine_map(local_vertices, point, reference, map_status)
            if (map_status /= 0 .or. any(reference < -1.0e-9_dp) .or. &
                sum(reference) > 1.0_dp + 1.0e-9_dp) cycle
            do edge = 1, 6
                real_dofs(edge) = real(box_orientations(edge, tetrahedron), dp)* &
                    real(solution(box_dofs(edge, tetrahedron)), dp)
                imag_dofs(edge) = real(box_orientations(edge, tetrahedron), dp)* &
                    aimag(solution(box_dofs(edge, tetrahedron)))
            end do
            call evaluate_tetra_nedelec_interpolant_at_point( &
                local_vertices, basis, real_dofs, point, real_value, real_curl, &
                basis_status)
            call evaluate_tetra_nedelec_interpolant_at_point( &
                local_vertices, basis, imag_dofs, point, imag_value, imag_curl, &
                basis_status)
            if (basis_status /= 0) return
            field = cmplx(real_value, imag_value, dp)
            local_status = 0
            return
        end do
    end subroutine evaluate_complex_field

    subroutine render_torus_solution()
        real(dp) :: section_x(129), section_z(129), theta
        integer :: point

        do point = 0, 128
            theta = 2.0_dp*acos(-1.0_dp)*real(point, dp)/128.0_dp
            section_x(point + 1) = major_radius + minor_radius*cos(theta)
            section_z(point + 1) = minor_radius*sin(theta)
        end do
        call figure(figsize=[8.5_dp, 6.5_dp])
        call pcolormesh(torus_x_edges, torus_z_edges, transpose(torus_magnitude), &
            cmap="viridis")
        call colorbar(label="|E_h| from solved FEM--BEM state")
        call plot(section_x, section_z, label="torus section", color=white, linewidth=1.2_dp)
        if (arrow_count > 0) call quiver( &
            torus_arrow_x(:arrow_count), torus_arrow_z(:arrow_count), &
            torus_arrow_u(:arrow_count), torus_arrow_v(:arrow_count), color=black)
        call xlabel("x at y = 0.017")
        call ylabel("z")
        call title("Solved nonzero toroidal curl--curl vector field")
        call legend()
        call savefig(output_directory//"/maxwell_vector_solution_2d.png")
    end subroutine render_torus_solution

    subroutine render_pml_solution()
        call figure(figsize=[8.5_dp, 6.5_dp])
        call pcolormesh(pml_x_edges, pml_y_edges, transpose(pml_magnitude), cmap="inferno")
        call colorbar(label="|E_h| from Nedelec PML state")
        if (pml_arrow_count > 0) call quiver( &
            pml_arrow_x(:pml_arrow_count), pml_arrow_y(:pml_arrow_count), &
            pml_arrow_u(:pml_arrow_count), pml_arrow_v(:pml_arrow_count), color=white)
        call xlabel("x")
        call ylabel("y at z = 0.5")
        call title("Solved vector PML field and arrows")
        call savefig(output_directory//"/maxwell_vector_pml_solution_2d.png")
    end subroutine render_pml_solution

    subroutine render_torus_geometry()
        real(dp) :: surface_x(visual_polar + 1, visual_toroidal + 1)
        real(dp) :: surface_y(visual_polar + 1, visual_toroidal + 1)
        real(dp) :: surface_z(visual_polar + 1, visual_toroidal + 1)
        real(dp) :: theta, phi, radius
        integer :: polar, toroidal

        do polar = 0, visual_polar
            theta = 2.0_dp*acos(-1.0_dp)*real(polar, dp)/real(visual_polar, dp)
            do toroidal = 0, visual_toroidal
                phi = 2.0_dp*acos(-1.0_dp)*real(toroidal, dp)/real(visual_toroidal, dp)
                radius = major_radius + minor_radius*cos(theta)
                surface_x(polar + 1, toroidal + 1) = radius*cos(phi)
                surface_y(polar + 1, toroidal + 1) = radius*sin(phi)
                surface_z(polar + 1, toroidal + 1) = minor_radius*sin(theta)
            end do
        end do
        call figure(figsize=[8.0_dp, 6.5_dp])
        call add_parametric_surface(surface_x, surface_y, surface_z, &
            color="lightsteelblue", alpha=0.62_dp, linewidth=0.2_dp, filled=.true., &
            label="curved torus FEM--BEM surface")
        call add_scatter( &
            torus_vertices(1, :), torus_vertices(2, :), torus_vertices(3, :), &
            marker=".", markersize=4.0_dp, label="volume mesh vertices")
        call title("Curved surface used by vector FEM--BEM/DtN paths")
        call savefig(output_directory//"/maxwell_vector_torus_geometry_3d.png")
    end subroutine render_torus_geometry

    subroutine render_diagnostics()
        integer :: sample
        complex(dp) :: small_value(3), large_value(3)

        do sample = 1, size(domain_x)
            domain_x(sample) = real(sample - 1, dp)/real(size(domain_x) - 1, dp)
            call evaluate_complex_field( &
                small_vertices, small_tetrahedra, small_dofs, small_orientations, &
                small_solution, [domain_x(sample), 0.5_dp, 0.5_dp], &
                small_value, status)
            call evaluate_complex_field( &
                large_vertices, large_tetrahedra, large_dofs, large_orientations, &
                large_solution, [domain_x(sample), 0.5_dp, 0.5_dp], &
                large_value, status)
            domain_small(sample) = cmplx(sqrt(sum(abs(small_value)**2)), 0.0_dp, dp)
            domain_large(sample) = cmplx(sqrt(sum(abs(large_value)**2)), 0.0_dp, dp)
        end do
        call figure(figsize=[9.0_dp, 5.5_dp])
        call plot(real([(sample, sample=1, size(dtn_flux))], dp), abs(dtn_flux), &
            label="curved RWG/RBC DtN response")
        call xlabel("surface trace degree of freedom")
        call ylabel("absolute weak flux")
        call title("Vector Maxwell DtN diagnostic")
        call legend()
        call savefig(output_directory//"/maxwell_vector_dtn_diagnostic_1d.png")
        call figure(figsize=[9.0_dp, 5.5_dp])
        call plot(domain_x, abs(domain_small), label="base PML domain")
        call plot(domain_x, abs(domain_large), label="larger PML domain", linestyle="--")
        call xlabel("common interior x target")
        call ylabel("|E_h|")
        call title("Larger-domain vector PML comparison")
        call legend()
        call savefig(output_directory//"/maxwell_vector_domain_comparison_1d.png")
    end subroutine render_diagnostics

    subroutine write_torus_field_csv()
        integer :: output_unit, ix, iz, local_status
        real(dp) :: x, z, vector(3)

        open (newunit=output_unit, file=output_directory//"/torus_vector_field.csv", &
            status="replace", action="write")
        write (output_unit, "(a)") "x,z,Ex_real,Ey_real,Ez_real,E_magnitude"
        do iz = 1, slice_nz
            z = 0.5_dp*(torus_z_edges(iz) + torus_z_edges(iz + 1))
            do ix = 1, slice_nx
                x = 0.5_dp*(torus_x_edges(ix) + torus_x_edges(ix + 1))
                call evaluate_torus_field([x, 0.017_dp, z], vector, local_status)
                if (local_status /= 0) vector = 0.0_dp
                write (output_unit, "(6(es24.16,','))") x, z, vector(1), vector(2), vector(3), &
                    sqrt(sum(vector**2))
            end do
        end do
        close (output_unit)
    end subroutine write_torus_field_csv

    subroutine write_pml_field_csv()
        integer :: output_unit, ix, iy, local_status
        real(dp) :: x, y
        complex(dp) :: vector(3)

        open (newunit=output_unit, file=output_directory//"/pml_vector_field.csv", &
            status="replace", action="write")
        write (output_unit, "(a)") &
            "x,y,Ex_real,Ey_real,Ez_real,Ex_imag,Ey_imag,Ez_imag,E_magnitude"
        do iy = 1, pml_ny
            y = 0.5_dp*(pml_y_edges(iy) + pml_y_edges(iy + 1))
            do ix = 1, pml_nx
                x = 0.5_dp*(pml_x_edges(ix) + pml_x_edges(ix + 1))
                call evaluate_complex_field( &
                    small_vertices, small_tetrahedra, small_dofs, small_orientations, &
                    small_solution, [x, y, 0.5_dp], vector, local_status)
                if (local_status /= 0) vector = cmplx(0.0_dp, 0.0_dp, dp)
                write (output_unit, "(9(es24.16,','))") x, y, real(vector(1), dp), &
                    real(vector(2), dp), real(vector(3), dp), aimag(vector(1)), &
                    aimag(vector(2)), aimag(vector(3)), sqrt(sum(abs(vector)**2))
            end do
        end do
        close (output_unit)
    end subroutine write_pml_field_csv

    subroutine write_sequence(stage)
        character(*), intent(in) :: stage
        integer :: sequence_unit
        logical :: exists

        inquire(file=output_directory//"/gallery_sequence.txt", exist=exists)
        open (newunit=sequence_unit, file=output_directory//"/gallery_sequence.txt", &
            status="old", action="write", position="append")
        write (sequence_unit, "(a)") stage
        close (sequence_unit)
    end subroutine write_sequence

    subroutine write_outputs(fem_seconds, dtn_seconds_local, pml_seconds_local, &
            error, surface_norm, response_norm, domain_error)
        real(dp), intent(in) :: fem_seconds, dtn_seconds_local, pml_seconds_local
        real(dp), intent(in) :: error, surface_norm, response_norm, domain_error
        integer :: output_unit

        open (newunit=output_unit, file=output_directory//"/benchmark.txt", &
            status="replace", action="write")
        write (output_unit, "(a)") "quantity,value"
        write (output_unit, "(a,es24.16)") "fem_bem_seconds,", fem_seconds
        write (output_unit, "(a,es24.16)") "dtn_seconds,", dtn_seconds_local
        write (output_unit, "(a,es24.16)") "pml_seconds,", pml_seconds_local
        write (output_unit, "(a,es24.16)") "fem_bem_field_error,", error
        write (output_unit, "(a,es24.16)") "fem_bem_surface_current_norm,", surface_norm
        write (output_unit, "(a,es24.16)") "dtn_response_norm,", response_norm
        write (output_unit, "(a,es24.16)") "larger_domain_difference,", domain_error
        write (output_unit, "(a,l1)") "fast_gallery,", fast_gallery
        close (output_unit)
    end subroutine write_outputs

    logical function environment_flag(name)
        character(*), intent(in) :: name

        character(64) :: value
        integer :: environment_status

        value = ""
        call get_environment_variable(name, value, status=environment_status)
        environment_flag = environment_status == 0 .and. &
            len_trim(value) > 0 .and. trim(value) /= "0"
    end function environment_flag

    subroutine evaluate_complex_field_value(box_vertices, box_tetrahedra, box_dofs, &
            box_orientations, solution, point, field)
        real(dp), intent(in) :: box_vertices(:, :), point(3)
        integer, intent(in) :: box_tetrahedra(:, :), box_dofs(:, :), box_orientations(:, :)
        complex(dp), intent(in) :: solution(:)
        complex(dp), intent(out) :: field(3)
        integer :: local_status

        call evaluate_complex_field(box_vertices, box_tetrahedra, box_dofs, &
            box_orientations, solution, point, field, local_status)
        if (local_status /= 0) error stop "complex field target evaluation failed"
    end subroutine evaluate_complex_field_value

    pure function target_point(index) result(point)
        integer, intent(in) :: index
        real(dp) :: point(3)

        select case (index)
        case (1)
            point = [0.25_dp, 0.25_dp, 0.50_dp]
        case (2)
            point = [0.50_dp, 0.50_dp, 0.50_dp]
        case default
            point = [0.75_dp, 0.75_dp, 0.50_dp]
        end select
    end function target_point

    pure function manufactured_vector(point) result(value)
        real(dp), intent(in) :: point(3)
        complex(dp) :: value(3)
        real(dp) :: rotational(3)

        rotational = [ &
            rotation(2)*point(3) - rotation(3)*point(2), &
            rotation(3)*point(1) - rotation(1)*point(3), &
            rotation(1)*point(2) - rotation(2)*point(1)]
        value = cmplx(vector_real + rotational, vector_imag, dp)
    end function manufactured_vector

    pure function plane_wave_edge(first, second) result(value)
        real(dp), intent(in) :: first(3), second(3)
        complex(dp) :: value, phase, increment
        real(dp) :: delta_x, delta_y

        delta_x = second(1) - first(1)
        delta_y = second(2) - first(2)
        phase = exp(cmplx(0.0_dp, wave_number, dp)* &
            cmplx(1.0_dp, pml_imaginary, dp)*first(1))
        if (abs(delta_x) < 1.0e-14_dp) then
            value = delta_y*phase
        else
            increment = cmplx(0.0_dp, wave_number, dp)* &
                cmplx(1.0_dp, pml_imaginary, dp)*delta_x
            value = delta_y*phase*(exp(increment) - 1.0_dp)/increment
        end if
    end function plane_wave_edge

    pure function plane_wave_field(point) result(value)
        real(dp), intent(in) :: point(3)
        complex(dp) :: value(3)

        value = [cmplx(0.0_dp, 0.0_dp, dp), exp( &
            cmplx(0.0_dp, wave_number, dp)* &
            cmplx(1.0_dp, pml_imaginary, dp)*point(1)), &
            cmplx(0.0_dp, 0.0_dp, dp)]
    end function plane_wave_field

end program maxwell_vector_open_boundary_parity_gallery
