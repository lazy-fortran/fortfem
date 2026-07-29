module fortfem_scalar_helmholtz_planar_dtn_2d
    use fortfem_kinds, only: dp
    use fortfem_planar_helmholtz_dtn, only: &
        assemble_planar_helmholtz_dtn_form
    use fortsparse, only: csc_from_triplet, csc_z_t, fortsparse_status_t, &
        sparse_solve_once
    implicit none
    private

    public :: solve_scalar_helmholtz_planar_dtn_p1

contains

    subroutine solve_scalar_helmholtz_planar_dtn_p1( &
            vertices, triangles, dtn_nodes, wavenumber, period, volume_load, &
            dirichlet_nodes, dirichlet_values, solution, status)
        real(dp), intent(in) :: vertices(:, :)
        integer, intent(in) :: triangles(:, :), dtn_nodes(:)
        real(dp), intent(in) :: wavenumber, period
        complex(dp), intent(in) :: volume_load(:)
        integer, intent(in) :: dirichlet_nodes(:)
        complex(dp), intent(in) :: dirichlet_values(:)
        complex(dp), intent(out) :: solution(:)
        integer, intent(out) :: status

        type(csc_z_t) :: sparse_matrix
        type(fortsparse_status_t) :: sparse_status
        complex(dp), allocatable :: dtn_form(:, :), matrix(:, :), rhs(:)
        complex(dp), allocatable :: triplet_values(:)
        integer, allocatable :: columns(:), rows(:)
        integer :: column, entry, node, row, vertex_count

        solution = cmplx(0.0_dp, 0.0_dp, dp)
        status = 1
        vertex_count = size(vertices, 2)
        if (.not. valid_inputs()) return

        allocate(matrix(vertex_count, vertex_count), rhs(vertex_count))
        call assemble_volume_matrix(matrix, status)
        if (status /= 0) return
        allocate(dtn_form(size(dtn_nodes), size(dtn_nodes)))
        call assemble_planar_helmholtz_dtn_form( &
            size(dtn_nodes), wavenumber, period, dtn_form, status)
        if (status /= 0) return
        matrix(dtn_nodes, dtn_nodes) = &
            matrix(dtn_nodes, dtn_nodes) - dtn_form
        rhs = volume_load

        do entry = 1, size(dirichlet_nodes)
            node = dirichlet_nodes(entry)
            rhs = rhs - matrix(:, node)*dirichlet_values(entry)
            matrix(:, node) = cmplx(0.0_dp, 0.0_dp, dp)
            matrix(node, :) = cmplx(0.0_dp, 0.0_dp, dp)
            matrix(node, node) = cmplx(1.0_dp, 0.0_dp, dp)
            rhs(node) = dirichlet_values(entry)
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
        call sparse_solve_once(sparse_matrix, rhs, solution, sparse_status)
        if (sparse_status%code /= 0) return
        status = 0

    contains

        logical function valid_inputs()
            valid_inputs = .false.
            if (size(vertices, 1) /= 2 .or. size(triangles, 1) /= 3) return
            if (vertex_count < 1 .or. size(triangles, 2) < 1) return
            if (wavenumber <= 0.0_dp .or. period <= 0.0_dp) return
            if (size(dtn_nodes) < 3) return
            if (size(volume_load) /= vertex_count .or. &
                size(solution) /= vertex_count) return
            if (size(dirichlet_nodes) /= size(dirichlet_values)) return
            if (any(triangles < 1) .or. any(triangles > vertex_count)) return
            if (any(dtn_nodes < 1) .or. any(dtn_nodes > vertex_count)) return
            if (any(dirichlet_nodes < 1) .or. &
                any(dirichlet_nodes > vertex_count)) return
            valid_inputs = .true.
        end function valid_inputs

        subroutine assemble_volume_matrix(volume_matrix, operator_status)
            complex(dp), intent(out) :: volume_matrix(:, :)
            integer, intent(out) :: operator_status

            real(dp) :: area, determinant
            real(dp) :: gradients(2, 3), local_mass(3, 3)
            real(dp) :: x1, x2, x3, y1, y2, y3
            integer :: first, second, element, local_nodes(3)

            volume_matrix = cmplx(0.0_dp, 0.0_dp, dp)
            operator_status = 1
            do element = 1, size(triangles, 2)
                local_nodes = triangles(:, element)
                x1 = vertices(1, local_nodes(1))
                y1 = vertices(2, local_nodes(1))
                x2 = vertices(1, local_nodes(2))
                y2 = vertices(2, local_nodes(2))
                x3 = vertices(1, local_nodes(3))
                y3 = vertices(2, local_nodes(3))
                determinant = (x2 - x1)*(y3 - y1) - &
                    (x3 - x1)*(y2 - y1)
                area = 0.5_dp*abs(determinant)
                if (area <= tiny(1.0_dp)) return
                gradients(:, 1) = [y2 - y3, x3 - x2]/determinant
                gradients(:, 2) = [y3 - y1, x1 - x3]/determinant
                gradients(:, 3) = [y1 - y2, x2 - x1]/determinant
                local_mass = area/12.0_dp
                local_mass(1, 1) = area/6.0_dp
                local_mass(2, 2) = area/6.0_dp
                local_mass(3, 3) = area/6.0_dp
                do first = 1, 3
                    do second = 1, 3
                        volume_matrix(local_nodes(first), local_nodes(second)) = &
                            volume_matrix(local_nodes(first), &
                            local_nodes(second)) + area*dot_product( &
                            gradients(:, first), gradients(:, second)) - &
                            wavenumber**2*local_mass(first, second)
                    end do
                end do
            end do
            operator_status = 0
        end subroutine assemble_volume_matrix
    end subroutine solve_scalar_helmholtz_planar_dtn_p1
end module fortfem_scalar_helmholtz_planar_dtn_2d
