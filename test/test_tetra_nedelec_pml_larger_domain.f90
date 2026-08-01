program test_tetra_nedelec_pml_larger_domain
    use check, only: check_condition, check_summary
    use fortfem_api, only: &
        build_tetra_edge_dof_map, evaluate_tetra_nedelec_interpolant_at_point, &
        generate_structured_tetra_box_mesh, invert_tetra_affine_map, &
        initialize_tetra_nedelec_first_kind, solve_tetra_nedelec_pml, &
        tetra_nedelec_first_kind_t
    use fortfem_kinds, only: dp
    implicit none

    real(dp), parameter :: wave_number = 1.4_dp
    real(dp), parameter :: stretch_imaginary = 0.35_dp
    real(dp), parameter :: sample_points(3, 3) = reshape([ &
        0.25_dp, 0.25_dp, 0.50_dp, &
        0.50_dp, 0.50_dp, 0.50_dp, &
        0.75_dp, 0.75_dp, 0.50_dp], [3, 3])
    real(dp), allocatable :: small_vertices(:, :), large_vertices(:, :)
    integer, allocatable :: small_tetrahedra(:, :), large_tetrahedra(:, :)
    integer, allocatable :: small_edges(:, :), large_edges(:, :)
    integer, allocatable :: small_global_dofs(:, :), large_global_dofs(:, :)
    integer, allocatable :: small_orientations(:, :), large_orientations(:, :)
    integer, allocatable :: small_boundary(:), large_boundary(:)
    complex(dp), allocatable :: small_solution(:), large_solution(:)
    real(dp) :: small_bounds(3, 2), large_bounds(3, 2)
    complex(dp) :: small_field(3), large_field(3), exact_field(3)
    complex(dp) :: small_curl(3), large_curl(3)
    real(dp) :: field_error, domain_difference
    integer :: status, sample
    logical :: all_passed

    all_passed = .true.
    small_bounds(:, 1) = 0.0_dp
    small_bounds(:, 2) = 1.0_dp
    large_bounds(:, 1) = 0.0_dp
    large_bounds(:, 2) = [1.0_dp, 1.0_dp, 1.5_dp]

    call solve_box(small_bounds, [3, 3, 3], small_vertices, &
        small_tetrahedra, small_edges, small_global_dofs, small_orientations, &
        small_boundary, small_solution)
    call solve_box(large_bounds, [3, 3, 5], large_vertices, &
        large_tetrahedra, large_edges, large_global_dofs, large_orientations, &
        large_boundary, large_solution)

    field_error = 0.0_dp
    domain_difference = 0.0_dp
    do sample = 1, size(sample_points, 2)
        call evaluate_field( &
            small_vertices, small_tetrahedra, small_global_dofs, &
            small_orientations, small_solution, sample_points(:, sample), &
            small_field, small_curl)
        call evaluate_field( &
            large_vertices, large_tetrahedra, large_global_dofs, &
            large_orientations, large_solution, sample_points(:, sample), &
            large_field, large_curl)
        exact_field = [cmplx(0.0_dp, 0.0_dp, dp), &
            exp(cmplx(0.0_dp, wave_number*sample_points(1, sample), dp)* &
                cmplx(1.0_dp, stretch_imaginary, dp)), &
            cmplx(0.0_dp, 0.0_dp, dp)]
        field_error = max(field_error, complex_norm(small_field - exact_field))
        field_error = max(field_error, complex_norm(large_field - exact_field))
        domain_difference = max(domain_difference, &
            complex_norm(small_field - large_field))
    end do

    write (*, '(a,es12.4)') 'PML exact-field error: ', field_error
    write (*, '(a,es12.4)') 'PML larger-domain difference: ', domain_difference
    call record_condition(field_error < 2.5e-1_dp, &
        'both PML boxes reproduce the plane-wave field')
    call record_condition(domain_difference < 1.0e-1_dp, &
        'moving the PML boundary preserves the interior field')
    call check_summary('Tetrahedral Nedelec PML larger-domain control')
    if (.not. all_passed) error stop 1

contains

    subroutine solve_box(bounds, counts, vertices, tetrahedra, edges, &
            global_dofs, orientations, boundary, solution)
        real(dp), intent(in) :: bounds(3, 2)
        integer, intent(in) :: counts(3)
        real(dp), allocatable, intent(out) :: vertices(:, :)
        integer, allocatable, intent(out) :: tetrahedra(:, :), edges(:, :)
        integer, allocatable, intent(out) :: global_dofs(:, :), orientations(:, :)
        integer, allocatable, intent(out) :: boundary(:)
        complex(dp), allocatable, intent(out) :: solution(:)

        complex(dp), allocatable :: stretch(:, :), load(:), exact(:)
        integer :: edge, status

        call generate_structured_tetra_box_mesh( &
            bounds, counts, vertices, tetrahedra, status)
        if (status /= 0) error stop 'PML test mesh construction failed'
        call build_tetra_edge_dof_map( &
            tetrahedra, edges, global_dofs, orientations, status)
        if (status /= 0) error stop 'PML test edge map construction failed'
        allocate(exact(size(edges, 2)), load(size(edges, 2)))
        do edge = 1, size(edges, 2)
            exact(edge) = plane_wave_edge( &
                vertices(:, edges(1, edge)), vertices(:, edges(2, edge)))
        end do
        call find_boundary_edges(vertices, edges, bounds, boundary)
        allocate(stretch(3, size(tetrahedra, 2)))
        stretch(1, :) = cmplx(1.0_dp, stretch_imaginary, dp)
        stretch(2:3, :) = cmplx(1.0_dp, 0.0_dp, dp)
        load = cmplx(0.0_dp, 0.0_dp, dp)
        call solve_tetra_nedelec_pml( &
            vertices, tetrahedra, 1, stretch, wave_number, load, boundary, &
            exact(boundary), solution, status)
        if (status /= 0) error stop 'PML test solve failed'
    end subroutine solve_box

    subroutine find_boundary_edges(vertices, edges, bounds, boundary)
        real(dp), intent(in) :: vertices(:, :), bounds(3, 2)
        integer, intent(in) :: edges(:, :)
        integer, allocatable, intent(out) :: boundary(:)

        logical, allocatable :: is_boundary(:)
        integer :: edge, coordinate

        allocate(is_boundary(size(edges, 2)))
        is_boundary = .false.
        do edge = 1, size(edges, 2)
            do coordinate = 1, 3
                if (all(abs(vertices(coordinate, edges(:, edge)) - &
                        bounds(coordinate, 1)) < 1.0e-13_dp) .or. &
                    all(abs(vertices(coordinate, edges(:, edge)) - &
                        bounds(coordinate, 2)) < 1.0e-13_dp)) &
                    is_boundary(edge) = .true.
            end do
        end do
        boundary = pack([(edge, edge=1, size(edges, 2))], is_boundary)
    end subroutine find_boundary_edges

    subroutine evaluate_field(vertices, tetrahedra, global_dofs, orientations, &
            solution, point, field, curl)
        real(dp), intent(in) :: vertices(:, :), point(3)
        integer, intent(in) :: tetrahedra(:, :), global_dofs(:, :)
        integer, intent(in) :: orientations(:, :)
        complex(dp), intent(in) :: solution(:)
        complex(dp), intent(out) :: field(3), curl(3)

        type(tetra_nedelec_first_kind_t) :: basis
        real(dp) :: local_vertices(3, 4), local_dofs(6), reference(3)
        real(dp) :: imaginary_dofs(6), real_field(3), real_curl(3)
        real(dp) :: imaginary_field(3), imaginary_curl(3)
        integer :: basis_status, map_status, local_edge, local_status, tetrahedron

        call initialize_tetra_nedelec_first_kind(1, basis, basis_status)
        if (basis_status /= 0) error stop 'PML test basis construction failed'
        field = cmplx(0.0_dp, 0.0_dp, dp)
        curl = cmplx(0.0_dp, 0.0_dp, dp)
        do tetrahedron = 1, size(tetrahedra, 2)
            local_vertices = vertices(:, tetrahedra(:, tetrahedron))
            call invert_tetra_affine_map( &
                local_vertices, point, reference, map_status)
            if (map_status /= 0 .or. any(reference < -1.0e-10_dp) .or. &
                sum(reference) > 1.0_dp + 1.0e-10_dp) cycle
            do local_edge = 1, 6
                local_dofs(local_edge) = &
                    real(orientations(local_edge, tetrahedron), dp)* &
                    real(solution(global_dofs(local_edge, tetrahedron)), dp)
                imaginary_dofs(local_edge) = &
                    real(orientations(local_edge, tetrahedron), dp)* &
                    aimag(solution(global_dofs(local_edge, tetrahedron)))
            end do
            call evaluate_tetra_nedelec_interpolant_at_point( &
                local_vertices, basis, local_dofs, point, real_field, real_curl, &
                local_status)
            if (local_status /= 0) error stop 'PML test real evaluation failed'
            call evaluate_tetra_nedelec_interpolant_at_point( &
                local_vertices, basis, imaginary_dofs, point, &
                imaginary_field, imaginary_curl, local_status)
            if (local_status /= 0) &
                error stop 'PML test imaginary evaluation failed'
            field = cmplx(real_field, imaginary_field, dp)
            curl = cmplx(real_curl, imaginary_curl, dp)
            exit
        end do
    end subroutine evaluate_field

    pure complex(dp) function plane_wave_edge(first, second) result(value)
        real(dp), intent(in) :: first(3), second(3)

        complex(dp), parameter :: x_stretch = cmplx(1.0_dp, 0.35_dp, dp)
        complex(dp) :: phase, phase_increment
        real(dp) :: delta_x, delta_y

        delta_x = second(1) - first(1)
        delta_y = second(2) - first(2)
        phase = exp(cmplx(0.0_dp, wave_number, dp)*x_stretch*first(1))
        if (abs(delta_x) < 1.0e-14_dp) then
            value = delta_y*phase
        else
            phase_increment = &
                cmplx(0.0_dp, wave_number, dp)*x_stretch*delta_x
            value = delta_y*phase*(exp(phase_increment) - 1.0_dp)/ &
                phase_increment
        end if
    end function plane_wave_edge

    pure real(dp) function complex_norm(value) result(norm_value)
        complex(dp), intent(in) :: value(:)

        norm_value = sqrt(sum(abs(value)**2))
    end function complex_norm

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_tetra_nedelec_pml_larger_domain
