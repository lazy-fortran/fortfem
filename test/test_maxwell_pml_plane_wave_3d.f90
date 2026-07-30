program test_maxwell_pml_plane_wave_3d
    use check, only: check_condition, check_summary
    use fortfem_api, only: build_tetra_edge_dof_map, &
        solve_tetra_nedelec_pml
    use fortfem_kinds, only: dp
    implicit none

    real(dp), parameter :: wave_number = 1.4_dp
    complex(dp), parameter :: x_stretch = cmplx(1.0_dp, 0.35_dp, dp)
    integer, parameter :: counts(2) = [3, 5]
    complex(dp), allocatable :: boundary_values(:), exact(:), load(:)
    complex(dp), allocatable :: solution(:), stretch(:, :)
    integer, allocatable :: boundary_dofs(:), edges(:, :), global_dofs(:, :)
    integer, allocatable :: orientations(:, :), tetrahedra(:, :)
    real(dp), allocatable :: vertices(:, :)
    real(dp) :: errors(2)
    integer :: edge, level, status
    logical :: all_passed

    all_passed = .true.
    do level = 1, 2
        call build_box_mesh(counts(level), vertices, tetrahedra)
        call build_tetra_edge_dof_map( &
            tetrahedra, edges, global_dofs, orientations, status)
        if (status /= 0) error stop "Maxwell PML edge map failed"
        allocate(exact(size(edges, 2)))
        do edge = 1, size(edges, 2)
            exact(edge) = exact_edge_integral( &
                vertices(:, edges(1, edge)), vertices(:, edges(2, edge)))
        end do
        call find_boundary_edges(vertices, edges, boundary_dofs)
        allocate(boundary_values(size(boundary_dofs)))
        boundary_values = exact(boundary_dofs)
        allocate(load(size(edges, 2)))
        load = cmplx(0.0_dp, 0.0_dp, dp)
        allocate(stretch(3, size(tetrahedra, 2)))
        stretch(1, :) = x_stretch
        stretch(2:3, :) = cmplx(1.0_dp, 0.0_dp, dp)
        call solve_tetra_nedelec_pml( &
            vertices, tetrahedra, 1, stretch, wave_number, load, &
            boundary_dofs, boundary_values, solution, status)
        if (status /= 0) error stop "Maxwell PML plane-wave solve failed"
        errors(level) = sqrt(sum(abs(solution - exact)**2))/ &
            sqrt(sum(abs(exact)**2))
        deallocate( &
            vertices, tetrahedra, edges, global_dofs, orientations, exact, &
            boundary_dofs, boundary_values, load, stretch, solution)
    end do

    write (*, '(a,2(es12.4,1x))') "Maxwell PML edge errors: ", errors
    call record_condition(errors(2) < 0.55_dp*errors(1), &
        "Maxwell PML plane wave converges under box refinement")
    call record_condition(errors(2) < 2.0e-2_dp, &
        "Maxwell PML matches the analytical complex-stretched plane wave")
    call check_summary("Physical Maxwell PML plane wave")
    if (.not. all_passed) error stop 1

contains

    subroutine build_box_mesh(count, vertices, tetrahedra)
        integer, intent(in) :: count
        real(dp), allocatable, intent(out) :: vertices(:, :)
        integer, allocatable, intent(out) :: tetrahedra(:, :)

        integer :: cell_vertices(8), i, j, k, local_tetrahedra(4, 6)
        integer :: local_node, tetrahedron, vertex

        allocate(vertices(3, count**3))
        vertex = 0
        do k = 0, count - 1
            do j = 0, count - 1
                do i = 0, count - 1
                    vertex = vertex + 1
                    vertices(:, vertex) = real([i, j, k], dp)/ &
                        real(count - 1, dp)
                end do
            end do
        end do
        allocate(tetrahedra(4, 6*(count - 1)**3))
        tetrahedron = 0
        do k = 0, count - 2
            do j = 0, count - 2
                do i = 0, count - 2
                    cell_vertices = [ &
                        box_node(i, j, k, count), &
                        box_node(i + 1, j, k, count), &
                        box_node(i, j + 1, k, count), &
                        box_node(i + 1, j + 1, k, count), &
                        box_node(i, j, k + 1, count), &
                        box_node(i + 1, j, k + 1, count), &
                        box_node(i, j + 1, k + 1, count), &
                        box_node(i + 1, j + 1, k + 1, count)]
                    local_tetrahedra = reshape([ &
                        1, 2, 4, 8, 1, 4, 3, 8, 1, 3, 7, 8, &
                        1, 7, 5, 8, 1, 5, 6, 8, 1, 6, 2, 8], [4, 6])
                    do vertex = 1, 6
                        do local_node = 1, 4
                            local_tetrahedra(local_node, vertex) = &
                                cell_vertices( &
                                local_tetrahedra(local_node, vertex))
                        end do
                    end do
                    do vertex = 1, 6
                        tetrahedron = tetrahedron + 1
                        tetrahedra(:, tetrahedron) = &
                            local_tetrahedra(:, vertex)
                        call orient_positive( &
                            vertices, tetrahedra(:, tetrahedron))
                    end do
                end do
            end do
        end do

    end subroutine build_box_mesh

    pure integer function box_node( &
            x_index, y_index, z_index, count) result(index)
        integer, intent(in) :: x_index, y_index, z_index, count

        index = 1 + x_index + count*(y_index + count*z_index)
    end function box_node

    subroutine find_boundary_edges(vertices, edges, boundary_dofs)
        real(dp), intent(in) :: vertices(:, :)
        integer, intent(in) :: edges(:, :)
        integer, allocatable, intent(out) :: boundary_dofs(:)

        logical, allocatable :: boundary(:)
        integer :: coordinate, edge

        allocate(boundary(size(edges, 2)))
        boundary = .false.
        do edge = 1, size(edges, 2)
            do coordinate = 1, 3
                if (all(abs(vertices(coordinate, edges(:, edge))) < &
                    1.0e-14_dp)) boundary(edge) = .true.
                if (all(abs(vertices(coordinate, edges(:, edge)) - 1.0_dp) < &
                    1.0e-14_dp)) boundary(edge) = .true.
            end do
        end do
        allocate(boundary_dofs(count(boundary)))
        boundary_dofs = pack([(edge, edge=1, size(edges, 2))], boundary)
    end subroutine find_boundary_edges

    pure complex(dp) function exact_edge_integral(first, second) result(value)
        real(dp), intent(in) :: first(3), second(3)

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
    end function exact_edge_integral

    pure subroutine orient_positive(vertices, tetrahedron)
        real(dp), intent(in) :: vertices(:, :)
        integer, intent(inout) :: tetrahedron(4)

        real(dp) :: determinant
        integer :: temporary

        determinant = dot_product( &
            vertices(:, tetrahedron(2)) - vertices(:, tetrahedron(1)), &
            vector_cross_product( &
            vertices(:, tetrahedron(3)) - vertices(:, tetrahedron(1)), &
            vertices(:, tetrahedron(4)) - vertices(:, tetrahedron(1))))
        if (determinant < 0.0_dp) then
            temporary = tetrahedron(3)
            tetrahedron(3) = tetrahedron(4)
            tetrahedron(4) = temporary
        end if
    end subroutine orient_positive

    pure function vector_cross_product(first, second) result(product)
        real(dp), intent(in) :: first(3), second(3)
        real(dp) :: product(3)

        product = [ &
            first(2)*second(3) - first(3)*second(2), &
            first(3)*second(1) - first(1)*second(3), &
            first(1)*second(2) - first(2)*second(1)]
    end function vector_cross_product

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_maxwell_pml_plane_wave_3d
