program test_scalar_helmholtz_pml_3d
    use check, only: check_condition, check_summary
    use fortfem_api, only: solve_scalar_helmholtz_pml_p1_3d
    use fortfem_kinds, only: dp
    implicit none

    real(dp), parameter :: wave_number = 1.6_dp
    complex(dp), parameter :: x_stretch = cmplx(1.0_dp, 0.4_dp, dp)
    integer, parameter :: counts(3) = [3, 5, 7]
    complex(dp), allocatable :: boundary_values(:), exact(:), load(:)
    complex(dp), allocatable :: solution(:), stretch(:, :)
    integer, allocatable :: boundary_nodes(:), tetrahedra(:, :)
    real(dp), allocatable :: vertices(:, :)
    real(dp) :: errors(3)
    integer :: level, status
    logical :: all_passed

    all_passed = .true.
    do level = 1, 3
        call build_box_mesh(counts(level), vertices, tetrahedra)
        call find_boundary_nodes(vertices, boundary_nodes)
        allocate(exact(size(vertices, 2)))
        exact = exp(cmplx(0.0_dp, wave_number, dp)* &
            x_stretch*vertices(1, :))
        allocate(boundary_values(size(boundary_nodes)))
        boundary_values = exact(boundary_nodes)
        allocate(load(size(vertices, 2)))
        load = cmplx(0.0_dp, 0.0_dp, dp)
        allocate(stretch(3, size(tetrahedra, 2)))
        stretch(1, :) = x_stretch
        stretch(2:3, :) = cmplx(1.0_dp, 0.0_dp, dp)
        call solve_scalar_helmholtz_pml_p1_3d( &
            vertices, tetrahedra, stretch, wave_number, load, boundary_nodes, &
            boundary_values, solution, status)
        if (status /= 0) error stop "3D scalar Helmholtz PML solve failed"
        errors(level) = sqrt(sum(abs(solution - exact)**2))/ &
            sqrt(sum(abs(exact)**2))
        deallocate( &
            vertices, tetrahedra, boundary_nodes, exact, boundary_values, &
            load, stretch, solution)
    end do

    write (*, '(a,3(es12.4,1x))') "3D scalar PML errors: ", errors
    call record_condition( &
        errors(2) < errors(1) .and. errors(3) < errors(2), &
        "3D scalar PML converges under tetrahedral refinement")
    call record_condition(errors(3) < 5.0e-3_dp, &
        "3D scalar PML matches a complex-stretched plane wave")
    call check_summary("Scalar Helmholtz PML in 3D")
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

    subroutine find_boundary_nodes(vertices, boundary_nodes)
        real(dp), intent(in) :: vertices(:, :)
        integer, allocatable, intent(out) :: boundary_nodes(:)

        logical, allocatable :: boundary(:)
        integer :: coordinate, node

        allocate(boundary(size(vertices, 2)))
        boundary = .false.
        do node = 1, size(vertices, 2)
            do coordinate = 1, 3
                if (abs(vertices(coordinate, node)) < 1.0e-14_dp) &
                    boundary(node) = .true.
                if (abs(vertices(coordinate, node) - 1.0_dp) < 1.0e-14_dp) &
                    boundary(node) = .true.
            end do
        end do
        allocate(boundary_nodes(count(boundary)))
        boundary_nodes = pack([(node, node=1, size(vertices, 2))], boundary)
    end subroutine find_boundary_nodes

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

end program test_scalar_helmholtz_pml_3d
