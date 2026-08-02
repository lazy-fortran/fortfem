program test_maxwell_vector_open_boundary_parity
    !! Physical-first parity oracle for vector curl-curl open boundaries.
    !!
    !! The three paths intentionally share only a manufactured vector target:
    !! a curved toroidal FEM--BEM solve, a curved RWG/RBC DtN action, and
    !! two Nedelec PML boxes.  No plasma-specific data or external readers are
    !! involved.  The small complex algebra check below is independent of the
    !! geometric assemblies and guards the value/JVP/VJP contract.
    use check, only: check_condition, check_summary
    use fortfem_api, only: &
        apply_maxwell_trace_to_flux, &
        apply_maxwell_trace_to_flux_jvp, apply_maxwell_trace_to_flux_vjp, &
        assemble_maxwell_fem_bem_torus_curved_system_3d, &
        assemble_maxwell_torus_curved_dtn_rwg_3d, &
        build_tetra_edge_dof_map, build_tetra_nedelec_dof_map, &
        evaluate_tetra_nedelec_interpolant_at_point, &
        generate_solid_torus_tetra_mesh, generate_structured_tetra_box_mesh, &
        generate_torus_surface_mesh, initialize_tetra_nedelec_first_kind, &
        invert_tetra_affine_map, solve_maxwell_fem_bem_torus_curved_system_3d, &
        solve_tetra_nedelec_pml, tetra_nedelec_first_kind_t
    use fortfem_kinds, only: dp
    implicit none

    real(dp), parameter :: major_radius = 2.0_dp, minor_radius = 0.6_dp
    real(dp), parameter :: wave_number = 1.1_dp, pml_imaginary = 0.28_dp
    real(dp), parameter :: target_points(3, 3) = reshape([ &
        0.25_dp, 0.25_dp, 0.50_dp, &
        0.50_dp, 0.50_dp, 0.50_dp, &
        0.75_dp, 0.75_dp, 0.50_dp], [3, 3])
    real(dp), parameter :: vector_real(3) = [0.48_dp, -0.31_dp, 0.22_dp]
    real(dp), parameter :: vector_imag(3) = [-0.17_dp, 0.26_dp, 0.11_dp]
    real(dp), parameter :: rotation(3) = [0.0_dp, 0.0_dp, 0.27_dp]
    integer, allocatable :: boundary_triangles(:, :), edges(:, :), faces(:, :)
    integer, allocatable :: edge_orientations(:, :), global_dofs(:, :)
    integer, allocatable :: face_permutations(:, :, :), tetrahedra(:, :)
    real(dp), allocatable :: parameters(:, :), vertices(:, :)
    complex(dp), allocatable :: matrix(:, :), rhs(:), exact(:), field(:), current(:)
    complex(dp), allocatable :: dtn_map(:, :), dtn_trace(:), dtn_flux(:)
    real(dp) :: field_error, current_norm
    real(dp) :: reciprocity_error, dtn_norm
    integer :: status, field_count, edge
    logical :: all_passed

    all_passed = .true.

    call run_toroidal_fem_bem(field_error, current_norm)
    call record_condition(status_is_finite(field_error) .and. &
        abs(field_error) < 3.0e-8_dp, &
        "curved FEM--BEM recovers the nonzero manufactured curl-curl field")
    call record_condition(abs(current_norm) < 3.0e-8_dp, &
        "curved FEM--BEM keeps the manufactured surface current at zero")

    call generate_torus_surface_mesh( &
        major_radius, minor_radius, 3, 3, vertices, boundary_triangles, parameters)
    call assemble_maxwell_torus_curved_dtn_rwg_3d( &
        vertices, boundary_triangles, parameters, major_radius, minor_radius, &
        0.45_dp, 1.7_dp, 3, 3.0e-4_dp, 1, 0.12_dp, dtn_map, status)
    call record_condition(status == 0, &
        "curved RWG/RBC Maxwell DtN map assembles on the torus")
    if (status == 0) then
        allocate(dtn_trace(size(dtn_map, 2)), dtn_flux(size(dtn_map, 1)))
        do edge = 1, size(dtn_trace)
            dtn_trace(edge) = cmplx( &
                cos(0.19_dp*real(edge, dp)), sin(0.13_dp*real(edge, dp)), dp)
        end do
        dtn_flux = matmul(dtn_map, dtn_trace)
        dtn_norm = sqrt(sum(abs(dtn_flux)**2))
        call record_condition(dtn_norm > 1.0e-10_dp, &
            "curved Maxwell DtN action is a nonzero physical trace response")
    end if

    call run_complex_trace_oracle(reciprocity_error)
    call record_condition(reciprocity_error < 1.0e-8_dp, &
        "complex Maxwell trace JVP and VJP satisfy the adjoint identity")

    call run_pml_larger_domain()
    call check_summary("Vector Maxwell FEM--BEM/DtN/PML parity")
    if (.not. all_passed) error stop 1

contains

    subroutine run_toroidal_fem_bem(error, surface_norm)
        real(dp), intent(out) :: error, surface_norm

        real(dp), allocatable :: torus_parameters(:, :), torus_vertices(:, :)
        integer, allocatable :: torus_boundary(:, :), torus_tetrahedra(:, :)
        integer, allocatable :: torus_edges(:, :), torus_faces(:, :)
        integer, allocatable :: torus_dofs(:, :), torus_orientations(:, :)
        integer, allocatable :: torus_permutations(:, :, :)
        complex(dp), allocatable :: state(:), load(:), solved(:), sheet(:)
        complex(dp) :: edge_value
        real(dp) :: midpoint(3), delta(3)
        integer :: local_status, local_field_count, local_edge

        call generate_solid_torus_tetra_mesh( &
            major_radius, minor_radius, 1, 3, 3, torus_vertices, &
            torus_tetrahedra, torus_boundary, torus_parameters)
        call build_tetra_nedelec_dof_map( &
            1, torus_tetrahedra, torus_edges, torus_faces, torus_dofs, &
            torus_orientations, torus_permutations, local_status)
        if (local_status /= 0) then
            error = huge(1.0_dp)
            surface_norm = huge(1.0_dp)
            return
        end if
        local_field_count = maxval(torus_dofs)
        call assemble_maxwell_fem_bem_torus_curved_system_3d( &
            torus_vertices, torus_tetrahedra, torus_parameters, torus_boundary, &
            major_radius, minor_radius, 1.2_dp, -0.6_dp, 0.65_dp, 1.4_dp, &
            3, 1.0e-5_dp, 1, matrix, local_status)
        if (local_status /= 0) then
            error = huge(1.0_dp)
            surface_norm = huge(1.0_dp)
            return
        end if
        allocate(state(size(matrix, 1)), load(size(matrix, 1)))
        state = cmplx(0.0_dp, 0.0_dp, dp)
        do local_edge = 1, size(torus_edges, 2)
            delta = torus_vertices(:, torus_edges(2, local_edge)) - &
                torus_vertices(:, torus_edges(1, local_edge))
            midpoint = 0.5_dp*(torus_vertices(:, torus_edges(2, local_edge)) + &
                torus_vertices(:, torus_edges(1, local_edge)))
            edge_value = dot_product(manufactured_vector(midpoint), delta)
            state(local_edge) = edge_value
        end do
        load = matmul(matrix, state)
        call solve_maxwell_fem_bem_torus_curved_system_3d( &
            torus_vertices, torus_tetrahedra, torus_parameters, torus_boundary, &
            major_radius, minor_radius, 1.2_dp, -0.6_dp, 0.65_dp, 1.4_dp, &
            3, 1.0e-5_dp, 1, load, solved, sheet, local_status)
        if (local_status /= 0 .or. size(solved) /= local_field_count) then
            error = huge(1.0_dp)
            surface_norm = huge(1.0_dp)
            return
        end if
        error = maxval(abs(solved - state(:local_field_count)))
        surface_norm = sqrt(sum(abs(sheet)**2))
    end subroutine run_toroidal_fem_bem

    subroutine run_complex_trace_oracle(error)
        real(dp), intent(out) :: error

        complex(dp) :: electric(3, 3), magnetic(3, 3), mass(3, 3)
        complex(dp) :: electric_dot(3, 3), magnetic_dot(3, 3), mass_dot(3, 3)
        complex(dp) :: trace(3), trace_dot(3), flux(3), flux_dot(3)
        complex(dp) :: flux_bar(3), trace_bar(3), flux_base(3)
        complex(dp) :: electric_bar(3, 3), magnetic_bar(3, 3), mass_bar(3, 3)
        complex(dp) :: flux_plus(3), flux_minus(3)
        real(dp) :: epsilon, lhs, rhs, central_error, adjoint_error
        integer :: local_status, row, column

        do column = 1, 3
            do row = 1, 3
                electric(row, column) = cmplx( &
                    1.8_dp + 0.21_dp*real(row + column, dp), &
                    0.07_dp*real(row - column, dp), dp)
                magnetic(row, column) = cmplx( &
                    0.4_dp + 0.13_dp*real(row + 2*column, dp), &
                    -0.04_dp*real(row - column, dp), dp)
                mass(row, column) = cmplx( &
                    1.2_dp + 0.09_dp*real(row + column, dp), &
                    0.03_dp*real(column - row, dp), dp)
                electric_dot(row, column) = cmplx(0.01_dp*real(row, dp), &
                    -0.02_dp*real(column, dp), dp)
                magnetic_dot(row, column) = cmplx(-0.03_dp*real(column, dp), &
                    0.01_dp*real(row, dp), dp)
                mass_dot(row, column) = cmplx(0.02_dp*real(row + column, dp), &
                    0.015_dp*real(row - column, dp), dp)
            end do
        end do
        trace = [cmplx(0.7_dp, -0.2_dp, dp), cmplx(-0.3_dp, 0.4_dp, dp), &
            cmplx(0.2_dp, 0.6_dp, dp)]
        trace_dot = [cmplx(-0.12_dp, 0.08_dp, dp), cmplx(0.2_dp, -0.1_dp, dp), &
            cmplx(0.05_dp, 0.11_dp, dp)]
        call apply_maxwell_trace_to_flux_jvp( &
            electric, magnetic, mass, trace, electric_dot, magnetic_dot, &
            mass_dot, trace_dot, flux_dot, local_status)
        if (local_status /= 0) then
            error = huge(1.0_dp)
            return
        end if
        epsilon = 2.0e-6_dp
        call apply_maxwell_trace_to_flux( &
            electric, magnetic, mass, trace, flux_base, local_status)
        call apply_maxwell_trace_to_flux( &
            electric + epsilon*electric_dot, magnetic + epsilon*magnetic_dot, &
            mass + epsilon*mass_dot, trace + epsilon*trace_dot, &
            flux_plus, local_status)
        call apply_maxwell_trace_to_flux( &
            electric - epsilon*electric_dot, magnetic - epsilon*magnetic_dot, &
            mass - epsilon*mass_dot, trace - epsilon*trace_dot, &
            flux_minus, local_status)
        central_error = maxval(abs(flux_dot - &
            (flux_plus - flux_minus)/(2.0_dp*epsilon)))
        error = central_error

        flux_bar = [cmplx(0.31_dp, -0.2_dp, dp), cmplx(-0.17_dp, 0.08_dp, dp), &
            cmplx(0.23_dp, 0.11_dp, dp)]
        call apply_maxwell_trace_to_flux_vjp( &
            electric, magnetic, mass, trace, flux_bar, trace_bar, electric_bar, &
            magnetic_bar, mass_bar, local_status)
        if (local_status /= 0) then
            error = huge(1.0_dp)
            return
        end if
        call apply_maxwell_trace_to_flux_jvp( &
            electric, magnetic, mass, trace, 0.6_dp*electric_dot, &
            0.6_dp*magnetic_dot, 0.6_dp*mass_dot, 0.6_dp*trace_dot, flux, &
            local_status)
        lhs = real(sum(conjg(flux_bar)*flux), dp)
        rhs = real(sum(conjg(electric_bar)*(0.6_dp*electric_dot)) + &
            sum(conjg(magnetic_bar)*(0.6_dp*magnetic_dot)) + &
            sum(conjg(mass_bar)*(0.6_dp*mass_dot)) + &
            sum(conjg(trace_bar)*(0.6_dp*trace_dot)), dp)
        adjoint_error = abs(lhs - rhs)
        error = max(error, adjoint_error)
    end subroutine run_complex_trace_oracle

    subroutine run_pml_larger_domain()
        real(dp) :: small_bounds(3, 2), large_bounds(3, 2)
        real(dp), allocatable :: small_vertices(:, :), large_vertices(:, :)
        integer, allocatable :: small_tetrahedra(:, :), large_tetrahedra(:, :)
        integer, allocatable :: small_edges(:, :), large_edges(:, :)
        integer, allocatable :: small_dofs(:, :), large_dofs(:, :)
        integer, allocatable :: small_orientations(:, :), large_orientations(:, :)
        integer, allocatable :: small_boundary(:), large_boundary(:)
        complex(dp), allocatable :: small_solution(:), large_solution(:)
        complex(dp) :: small_field(3), large_field(3), expected(3)
        integer :: local_status, sample
        real(dp) :: error, difference

        small_bounds(:, 1) = 0.0_dp
        small_bounds(:, 2) = 1.0_dp
        large_bounds(:, 1) = 0.0_dp
        large_bounds(:, 2) = [1.0_dp, 1.0_dp, 1.5_dp]
        call solve_pml_box(small_bounds, [2, 2, 2], small_vertices, &
            small_tetrahedra, small_edges, small_dofs, small_orientations, &
            small_boundary, small_solution)
        call solve_pml_box(large_bounds, [2, 2, 3], large_vertices, &
            large_tetrahedra, large_edges, large_dofs, large_orientations, &
            large_boundary, large_solution)
        error = 0.0_dp
        difference = 0.0_dp
        do sample = 1, size(target_points, 2)
            call evaluate_complex_field( &
                small_vertices, small_tetrahedra, small_dofs, small_orientations, &
                small_solution, target_points(:, sample), small_field, local_status)
            call evaluate_complex_field( &
                large_vertices, large_tetrahedra, large_dofs, large_orientations, &
                large_solution, target_points(:, sample), large_field, local_status)
            expected = plane_wave_field(target_points(:, sample))
            error = max(error, sqrt(sum(abs(small_field - expected)**2)))
            error = max(error, sqrt(sum(abs(large_field - expected)**2)))
            difference = max(difference, sqrt(sum(abs(small_field - large_field)**2)))
        end do
        call record_condition(error < 3.0e-1_dp, &
            "Nedelec PML reproduces the nonzero vector target in both domains")
        call record_condition(difference < 1.5e-1_dp, &
            "larger-domain PML leaves common interior vector targets stable")
    end subroutine run_pml_larger_domain

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
        integer :: box_status, edge

        call generate_structured_tetra_box_mesh( &
            bounds, counts, box_vertices, box_tetrahedra, box_status)
        if (box_status /= 0) error stop "vector parity PML mesh failed"
        call build_tetra_edge_dof_map( &
            box_tetrahedra, box_edges, box_dofs, box_orientations, box_status)
        if (box_status /= 0) error stop "vector parity PML edge map failed"
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
            box_boundary, exact(box_boundary), solution, box_status)
        if (box_status /= 0) error stop "vector parity PML solve failed"
    end subroutine solve_pml_box

    subroutine boundary_edges(box_vertices, box_edges, bounds, boundary)
        real(dp), intent(in) :: box_vertices(:, :), bounds(3, 2)
        integer, intent(in) :: box_edges(:, :)
        integer, allocatable, intent(out) :: boundary(:)
        logical, allocatable :: marked(:)
        integer :: edge, coordinate

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

    subroutine evaluate_complex_field(box_vertices, box_tetrahedra, box_dofs, &
            box_orientations, solution, point, field, local_status)
        real(dp), intent(in) :: box_vertices(:, :), point(3)
        integer, intent(in) :: box_tetrahedra(:, :), box_dofs(:, :)
        integer, intent(in) :: box_orientations(:, :)
        complex(dp), intent(in) :: solution(:)
        complex(dp), intent(out) :: field(3)
        integer, intent(out) :: local_status

        type(tetra_nedelec_first_kind_t) :: basis
        real(dp) :: local_vertices(3, 4), reference(3), local_real(6), local_imag(6)
        real(dp) :: real_value(3), imag_value(3), real_curl(3), imag_curl(3)
        integer :: basis_status, map_status, local_edge, tetrahedron

        field = cmplx(0.0_dp, 0.0_dp, dp)
        local_status = 1
        call initialize_tetra_nedelec_first_kind(1, basis, basis_status)
        if (basis_status /= 0) return
        do tetrahedron = 1, size(box_tetrahedra, 2)
            local_vertices = box_vertices(:, box_tetrahedra(:, tetrahedron))
            call invert_tetra_affine_map(local_vertices, point, reference, map_status)
            if (map_status /= 0 .or. any(reference < -1.0e-10_dp) .or. &
                sum(reference) > 1.0_dp + 1.0e-10_dp) cycle
            do local_edge = 1, 6
                local_real(local_edge) = real(box_orientations(local_edge, tetrahedron), dp)* &
                    real(solution(box_dofs(local_edge, tetrahedron)), dp)
                local_imag(local_edge) = real(box_orientations(local_edge, tetrahedron), dp)* &
                    aimag(solution(box_dofs(local_edge, tetrahedron)))
            end do
            call evaluate_tetra_nedelec_interpolant_at_point( &
                local_vertices, basis, local_real, point, real_value, real_curl, &
                basis_status)
            call evaluate_tetra_nedelec_interpolant_at_point( &
                local_vertices, basis, local_imag, point, imag_value, imag_curl, &
                basis_status)
            if (basis_status /= 0) return
            field = cmplx(real_value, imag_value, dp)
            local_status = 0
            return
        end do
    end subroutine evaluate_complex_field

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

    logical pure function status_is_finite(value)
        real(dp), intent(in) :: value
        status_is_finite = value < huge(1.0_dp)
    end function status_is_finite

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_maxwell_vector_open_boundary_parity
