module fortfem_scalar_helmholtz_pml_2d
    !! P1 scalar Helmholtz assembly with elementwise diagonal complex stretch.
    use fortfem_cartesian_helmholtz_pml, only: &
        cartesian_scalar_helmholtz_pml_coefficients
    use fortfem_kinds, only: dp
    use fortsparse, only: csc_from_triplet, csc_z_t, fortsparse_status_t, &
        sparse_solve_once
    implicit none
    private

    public :: solve_scalar_helmholtz_pml_p1_2d

contains

    subroutine solve_scalar_helmholtz_pml_p1_2d( &
            vertices, triangles, stretch, wave_number, volume_load, &
            dirichlet_nodes, dirichlet_values, solution, status)
        real(dp), intent(in) :: vertices(:, :)
        integer, intent(in) :: triangles(:, :)
        complex(dp), intent(in) :: stretch(:, :)
        real(dp), intent(in) :: wave_number
        complex(dp), intent(in) :: volume_load(:)
        integer, intent(in) :: dirichlet_nodes(:)
        complex(dp), intent(in) :: dirichlet_values(:)
        complex(dp), allocatable, intent(out) :: solution(:)
        integer, intent(out) :: status

        type(csc_z_t) :: sparse_matrix
        type(fortsparse_status_t) :: sparse_status
        complex(dp), allocatable :: matrix(:, :), right_hand_side(:)
        complex(dp), allocatable :: triplet_values(:)
        complex(dp) :: gradient_coefficient(3), mass_coefficient
        complex(dp) :: stretch_3d(3)
        real(dp) :: area, determinant, gradients(2, 3), local_mass
        integer, allocatable :: columns(:), rows(:)
        integer :: column, element, entry, local_status, node
        integer :: row, local_nodes(3), vertex_count

        status = 1
        if (allocated(solution)) deallocate(solution)
        if (size(vertices, 1) /= 2) return
        vertex_count = size(vertices, 2)
        if (vertex_count < 3) return
        if (size(triangles, 1) /= 3) return
        if (size(triangles, 2) < 1) return
        if (size(stretch, 1) /= 2) return
        if (size(stretch, 2) /= size(triangles, 2)) return
        if (size(volume_load) /= vertex_count) return
        if (size(dirichlet_nodes) /= size(dirichlet_values)) return
        if (wave_number <= 0.0_dp) return
        if (any(triangles < 1) .or. any(triangles > vertex_count)) return
        if (any(dirichlet_nodes < 1) .or. &
            any(dirichlet_nodes > vertex_count)) return

        allocate(matrix(vertex_count, vertex_count))
        allocate(right_hand_side(vertex_count), solution(vertex_count))
        matrix = cmplx(0.0_dp, 0.0_dp, dp)
        right_hand_side = volume_load
        do element = 1, size(triangles, 2)
            local_nodes = triangles(:, element)
            determinant = &
                (vertices(1, local_nodes(2)) - vertices(1, local_nodes(1)))* &
                (vertices(2, local_nodes(3)) - vertices(2, local_nodes(1))) - &
                (vertices(1, local_nodes(3)) - vertices(1, local_nodes(1)))* &
                (vertices(2, local_nodes(2)) - vertices(2, local_nodes(1)))
            area = 0.5_dp*abs(determinant)
            if (area <= tiny(1.0_dp)) return
            gradients(:, 1) = [ &
                vertices(2, local_nodes(2)) - vertices(2, local_nodes(3)), &
                vertices(1, local_nodes(3)) - vertices(1, local_nodes(2))]/ &
                determinant
            gradients(:, 2) = [ &
                vertices(2, local_nodes(3)) - vertices(2, local_nodes(1)), &
                vertices(1, local_nodes(1)) - vertices(1, local_nodes(3))]/ &
                determinant
            gradients(:, 3) = [ &
                vertices(2, local_nodes(1)) - vertices(2, local_nodes(2)), &
                vertices(1, local_nodes(2)) - vertices(1, local_nodes(1))]/ &
                determinant
            stretch_3d = [ &
                stretch(1, element), stretch(2, element), &
                cmplx(1.0_dp, 0.0_dp, dp)]
            call cartesian_scalar_helmholtz_pml_coefficients( &
                stretch_3d, gradient_coefficient, mass_coefficient, &
                local_status)
            if (local_status /= 0) return
            do column = 1, 3
                do row = 1, 3
                    local_mass = area/12.0_dp
                    if (row == column) local_mass = area/6.0_dp
                    matrix(local_nodes(row), local_nodes(column)) = &
                        matrix(local_nodes(row), local_nodes(column)) + area*( &
                        gradient_coefficient(1)*gradients(1, row)* &
                        gradients(1, column) + gradient_coefficient(2)* &
                        gradients(2, row)*gradients(2, column)) - &
                        wave_number**2*mass_coefficient*local_mass
                end do
            end do
        end do

        do element = 1, size(dirichlet_nodes)
            node = dirichlet_nodes(element)
            right_hand_side = right_hand_side - &
                matrix(:, node)*dirichlet_values(element)
            matrix(:, node) = cmplx(0.0_dp, 0.0_dp, dp)
            matrix(node, :) = cmplx(0.0_dp, 0.0_dp, dp)
            matrix(node, node) = cmplx(1.0_dp, 0.0_dp, dp)
            right_hand_side(node) = dirichlet_values(element)
        end do
        allocate(rows(vertex_count**2), columns(vertex_count**2))
        allocate(triplet_values(vertex_count**2))
        entry = 0
        do column = 1, vertex_count
            do row = 1, vertex_count
                if (abs(matrix(row, column)) <= tiny(1.0_dp)) cycle
                entry = entry + 1
                rows(entry) = row
                columns(entry) = column
                triplet_values(entry) = matrix(row, column)
            end do
        end do
        call csc_from_triplet( &
            vertex_count, vertex_count, rows(:entry), columns(:entry), &
            triplet_values(:entry), sparse_matrix, sparse_status)
        if (sparse_status%code /= 0) return
        call sparse_solve_once( &
            sparse_matrix, right_hand_side, solution, sparse_status)
        if (sparse_status%code /= 0) return
        status = 0
    end subroutine solve_scalar_helmholtz_pml_p1_2d

end module fortfem_scalar_helmholtz_pml_2d
