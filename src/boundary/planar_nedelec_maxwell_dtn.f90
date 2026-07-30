module fortfem_planar_nedelec_maxwell_dtn
    !! Pull back the planar Maxwell capacity map to tetrahedral Nedelec traces.
    use fortfem_kinds, only: dp
    use fortfem_planar_maxwell_dtn, only: assemble_planar_maxwell_dtn_form
    use fortfem_tetra_nedelec_arbitrary_order, only: &
        evaluate_tetra_nedelec_first_kind, &
        initialize_tetra_nedelec_first_kind, tetra_nedelec_dof_count, &
        tetra_nedelec_first_kind_t
    use fortfem_tetra_nedelec_global_dof_map, only: &
        build_tetra_nedelec_basis_transform, build_tetra_nedelec_dof_map
    use fortfem_tetra_piola_maps, only: map_tetra_nedelec_covariant
    use fortnum_linalg, only: det3, inv3
    implicit none
    private

    public :: assemble_planar_nedelec_maxwell_dtn_form
    public :: build_planar_nedelec_trace_sampling

contains

    subroutine assemble_planar_nedelec_maxwell_dtn_form( &
            vertices, tetrahedra, order, origin, periods, nx, ny, &
            wave_number, boundary_dofs, form, status)
        real(dp), intent(in) :: vertices(:, :), origin(3), periods(3, 2)
        integer, intent(in) :: tetrahedra(:, :), order, nx, ny
        real(dp), intent(in) :: wave_number
        integer, allocatable, intent(out) :: boundary_dofs(:)
        complex(dp), allocatable, intent(out) :: form(:, :)
        integer, intent(out) :: status

        complex(dp), allocatable :: capacity(:, :), sampling(:, :)
        real(dp) :: length_x, length_y

        status = 1
        if (allocated(form)) deallocate(form)
        call build_planar_nedelec_trace_sampling( &
            vertices, tetrahedra, order, origin, periods, nx, ny, &
            boundary_dofs, sampling, status)
        if (status /= 0) return
        length_x = norm2(periods(:, 1))
        length_y = norm2(periods(:, 2))
        call assemble_planar_maxwell_dtn_form( &
            nx, ny, wave_number, length_x, length_y, capacity, status)
        if (status /= 0) return
        allocate(form(size(boundary_dofs), size(boundary_dofs)))
        form = matmul(transpose(sampling), matmul(capacity, sampling))
        status = 0
    end subroutine assemble_planar_nedelec_maxwell_dtn_form

    subroutine build_planar_nedelec_trace_sampling( &
            vertices, tetrahedra, order, origin, periods, nx, ny, &
            boundary_dofs, sampling, status)
        real(dp), intent(in) :: vertices(:, :), origin(3), periods(3, 2)
        integer, intent(in) :: tetrahedra(:, :), order, nx, ny
        integer, allocatable, intent(out) :: boundary_dofs(:)
        complex(dp), allocatable, intent(out) :: sampling(:, :)
        integer, intent(out) :: status

        type(tetra_nedelec_first_kind_t) :: basis
        integer, allocatable :: edge_orientations(:, :), edges(:, :)
        integer, allocatable :: face_permutations(:, :, :), faces(:, :)
        integer, allocatable :: global_dofs(:, :)
        logical, allocatable :: face_selected(:), selected(:)
        real(dp), allocatable :: basis_transform(:, :)
        real(dp), allocatable :: physical_curls(:, :), physical_values(:, :)
        real(dp), allocatable :: reference_curls(:, :)
        real(dp), allocatable :: reference_values(:, :)
        real(dp), allocatable :: traced_basis(:, :)
        real(dp) :: axis_x(3), axis_y(3), determinant, inverse(3, 3)
        real(dp) :: jacobian(3, 3), normal(3), point(3), reference_point(3)
        real(dp) :: scale, tolerance
        integer :: boundary_count, column, dof, dof_count, edge
        integer :: face, face_dof_count, i, incidence, j, sample, tetrahedron

        status = 1
        if (allocated(boundary_dofs)) deallocate(boundary_dofs)
        if (allocated(sampling)) deallocate(sampling)
        if (.not. valid_input(vertices, tetrahedra, order, periods, nx, ny)) &
            return
        axis_x = periods(:, 1)/norm2(periods(:, 1))
        axis_y = periods(:, 2)/norm2(periods(:, 2))
        normal = cross_product(axis_x, axis_y)
        if (norm2(normal) <= 64.0_dp*epsilon(1.0_dp)) return
        normal = normal/norm2(normal)
        if (abs(dot_product(axis_x, axis_y)) > &
            256.0_dp*epsilon(1.0_dp)) return
        scale = max(1.0_dp, maxval(abs(vertices)), maxval(abs(periods)))
        tolerance = 512.0_dp*epsilon(1.0_dp)*scale

        call build_tetra_nedelec_dof_map( &
            order, tetrahedra, edges, faces, global_dofs, edge_orientations, &
            face_permutations, status)
        if (status /= 0) return
        allocate(selected(maxval(global_dofs)), source=.false.)
        allocate(face_selected(size(faces, 2)), source=.false.)
        do face = 1, size(faces, 2)
            incidence = 0
            do tetrahedron = 1, size(tetrahedra, 2)
                if (all_vertices_in_cell( &
                    faces(:, face), tetrahedra(:, tetrahedron))) then
                    incidence = incidence + 1
                end if
            end do
            face_selected(face) = incidence == 1 .and. all(abs( &
                projected_distances( &
                vertices(:, faces(:, face)), origin, normal)) <= &
                tolerance) .and. face_inside_patch( &
                vertices(:, faces(:, face)), origin, axis_x, axis_y, &
                norm2(periods(:, 1)), norm2(periods(:, 2)), tolerance)
        end do
        do edge = 1, size(edges, 2)
            if (.not. any([( &
                face_selected(face) .and. &
                all_vertices_in_cell(edges(:, edge), faces(:, face)), &
                face=1, size(faces, 2))])) cycle
            do dof = 1, order
                selected((edge - 1)*order + dof) = .true.
            end do
        end do
        face_dof_count = order*(order - 1)
        do face = 1, size(faces, 2)
            if (.not. face_selected(face)) cycle
            do dof = 1, face_dof_count
                selected(order*size(edges, 2) + &
                    (face - 1)*face_dof_count + dof) = .true.
            end do
        end do
        boundary_count = count(selected)
        if (boundary_count == 0) return
        allocate(boundary_dofs(boundary_count))
        boundary_dofs = pack([(dof, dof=1, size(selected))], selected)

        call initialize_tetra_nedelec_first_kind(order, basis, status)
        if (status /= 0) return
        dof_count = tetra_nedelec_dof_count(basis)
        allocate(basis_transform(dof_count, dof_count))
        allocate(reference_values(3, dof_count), reference_curls(3, dof_count))
        allocate(physical_values(3, dof_count), physical_curls(3, dof_count))
        allocate(traced_basis(2, dof_count))
        allocate(sampling(2*nx*ny, boundary_count))
        sampling = cmplx(0.0_dp, 0.0_dp, dp)
        sample = 0
        do j = 1, ny
            do i = 1, nx
                point = origin + &
                    real(i - 1, dp)/real(nx, dp)*periods(:, 1) + &
                    real(j - 1, dp)/real(ny, dp)*periods(:, 2)
                if (.not. point_on_selected_surface( &
                    vertices, faces, face_selected, point, tolerance)) &
                    return
                call locate_tetrahedron( &
                    vertices, tetrahedra, point, tolerance, tetrahedron, &
                    reference_point, jacobian, inverse, determinant, status)
                if (status /= 0) return
                call evaluate_tetra_nedelec_first_kind( &
                    basis, reference_point, reference_values, &
                    reference_curls, status)
                if (status /= 0) return
                call map_tetra_nedelec_covariant( &
                    jacobian, reference_values, reference_curls, &
                    physical_values, physical_curls, status)
                if (status /= 0) return
                call build_tetra_nedelec_basis_transform( &
                    order, edge_orientations(:, tetrahedron), &
                    face_permutations(:, :, tetrahedron), basis_transform, &
                    status)
                if (status /= 0) return
                traced_basis(1, :) = matmul(axis_x, &
                    matmul(physical_values, basis_transform))
                traced_basis(2, :) = matmul(axis_y, &
                    matmul(physical_values, basis_transform))
                sample = sample + 1
                do dof = 1, dof_count
                    column = find_index( &
                        boundary_dofs, global_dofs(dof, tetrahedron))
                    if (column == 0) cycle
                    sampling(2*sample - 1:2*sample, column) = &
                        sampling(2*sample - 1:2*sample, column) + &
                        cmplx(traced_basis(:, dof), 0.0_dp, dp)
                end do
            end do
        end do
        status = 0
    end subroutine build_planar_nedelec_trace_sampling

    pure logical function valid_input( &
            vertices, tetrahedra, order, periods, nx, ny) result(valid)
        real(dp), intent(in) :: vertices(:, :), periods(3, 2)
        integer, intent(in) :: tetrahedra(:, :), order, nx, ny

        valid = size(vertices, 1) == 3 .and. size(tetrahedra, 1) == 4 .and. &
            size(vertices, 2) >= 4 .and. size(tetrahedra, 2) >= 1 .and. &
            all(tetrahedra >= 1) .and. all(tetrahedra <= size(vertices, 2)) &
            .and. order >= 1 .and. nx >= 2 .and. ny >= 2 &
            .and. modulo(nx, 2) == 1 .and. modulo(ny, 2) == 1 .and. &
            norm2(periods(:, 1)) > tiny(1.0_dp) .and. &
            norm2(periods(:, 2)) > tiny(1.0_dp)
    end function valid_input

    subroutine locate_tetrahedron( &
            vertices, tetrahedra, point, tolerance, tetrahedron, reference, &
            jacobian, inverse, determinant, status)
        real(dp), intent(in) :: vertices(:, :), point(3), tolerance
        integer, intent(in) :: tetrahedra(:, :)
        integer, intent(out) :: tetrahedron
        real(dp), intent(out) :: reference(3), jacobian(3, 3)
        real(dp), intent(out) :: inverse(3, 3), determinant
        integer, intent(out) :: status

        integer :: candidate, inverse_status

        status = 1
        tetrahedron = 0
        do candidate = 1, size(tetrahedra, 2)
            jacobian(:, 1) = &
                vertices(:, tetrahedra(2, candidate)) - &
                vertices(:, tetrahedra(1, candidate))
            jacobian(:, 2) = &
                vertices(:, tetrahedra(3, candidate)) - &
                vertices(:, tetrahedra(1, candidate))
            jacobian(:, 3) = &
                vertices(:, tetrahedra(4, candidate)) - &
                vertices(:, tetrahedra(1, candidate))
            determinant = det3(jacobian)
            if (determinant <= tiny(1.0_dp)) cycle
            call inv3(jacobian, inverse, inverse_status)
            if (inverse_status /= 0) cycle
            reference = matmul( &
                inverse, point - vertices(:, tetrahedra(1, candidate)))
            if (any(reference < -tolerance)) cycle
            if (sum(reference) > 1.0_dp + tolerance) cycle
            tetrahedron = candidate
            status = 0
            return
        end do
    end subroutine locate_tetrahedron

    pure function projected_distances(points, origin, normal) result(distances)
        real(dp), intent(in) :: points(:, :), origin(3), normal(3)
        real(dp) :: distances(size(points, 2))
        integer :: point

        do point = 1, size(points, 2)
            distances(point) = dot_product(points(:, point) - origin, normal)
        end do
    end function projected_distances

    pure logical function face_inside_patch( &
            points, origin, axis_x, axis_y, length_x, length_y, tolerance) &
            result(inside)
        real(dp), intent(in) :: points(:, :), origin(3), axis_x(3), axis_y(3)
        real(dp), intent(in) :: length_x, length_y, tolerance
        real(dp) :: coordinate_x, coordinate_y
        integer :: point

        inside = .true.
        do point = 1, size(points, 2)
            coordinate_x = dot_product(points(:, point) - origin, axis_x)
            coordinate_y = dot_product(points(:, point) - origin, axis_y)
            if (coordinate_x < -tolerance .or. &
                coordinate_x > length_x + tolerance .or. &
                coordinate_y < -tolerance .or. &
                coordinate_y > length_y + tolerance) then
                inside = .false.
                return
            end if
        end do
    end function face_inside_patch

    pure logical function point_on_selected_surface( &
            vertices, faces, selected, point, tolerance) result(on_surface)
        real(dp), intent(in) :: vertices(:, :), point(3), tolerance
        integer, intent(in) :: faces(:, :)
        logical, intent(in) :: selected(:)
        real(dp) :: displacement(3), dot00, dot01, dot02, dot11, dot12
        real(dp) :: denominator, first_edge(3), second_edge(3), u, v
        integer :: face

        on_surface = .false.
        do face = 1, size(faces, 2)
            if (.not. selected(face)) cycle
            first_edge = &
                vertices(:, faces(2, face)) - vertices(:, faces(1, face))
            second_edge = &
                vertices(:, faces(3, face)) - vertices(:, faces(1, face))
            displacement = point - vertices(:, faces(1, face))
            dot00 = dot_product(first_edge, first_edge)
            dot01 = dot_product(first_edge, second_edge)
            dot02 = dot_product(first_edge, displacement)
            dot11 = dot_product(second_edge, second_edge)
            dot12 = dot_product(second_edge, displacement)
            denominator = dot00*dot11 - dot01*dot01
            if (denominator <= tiny(1.0_dp)) cycle
            u = (dot11*dot02 - dot01*dot12)/denominator
            v = (dot00*dot12 - dot01*dot02)/denominator
            if (u >= -tolerance .and. v >= -tolerance .and. &
                u + v <= 1.0_dp + tolerance) then
                on_surface = .true.
                return
            end if
        end do
    end function point_on_selected_surface

    pure integer function find_index(values, sought) result(index)
        integer, intent(in) :: values(:), sought

        do index = 1, size(values)
            if (values(index) == sought) return
        end do
        index = 0
    end function find_index

    pure logical function all_vertices_in_cell(vertices, cell) result(found)
        integer, intent(in) :: vertices(:), cell(:)
        integer :: vertex

        found = .true.
        do vertex = 1, size(vertices)
            if (.not. any(cell == vertices(vertex))) then
                found = .false.
                return
            end if
        end do
    end function all_vertices_in_cell

    pure function cross_product(first, second) result(product)
        real(dp), intent(in) :: first(3), second(3)
        real(dp) :: product(3)

        product = [ &
            first(2)*second(3) - first(3)*second(2), &
            first(3)*second(1) - first(1)*second(3), &
            first(1)*second(2) - first(2)*second(1)]
    end function cross_product

end module fortfem_planar_nedelec_maxwell_dtn
