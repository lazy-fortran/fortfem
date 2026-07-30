module fortfem_elasticity_p1_2d_core
    use fortfem_kinds, only: dp
    use fortsparse, only: csc_from_triplet, csc_z_t, fortsparse_status_t, &
        sparse_solve_once
    implicit none

    private

    public :: assemble_elasticity_dynamic_p1_2d
    public :: solve_dense_complex_dirichlet_fortsparse

contains

    subroutine assemble_elasticity_dynamic_p1_2d( &
            vertices, triangles, angular_frequency, young_modulus, &
            poisson_ratio, solid_density, matrix, status)
        real(dp), intent(in) :: vertices(:, :)
        integer, intent(in) :: triangles(:, :)
        real(dp), intent(in) :: angular_frequency
        complex(dp), intent(in) :: young_modulus
        real(dp), intent(in) :: poisson_ratio, solid_density
        complex(dp), intent(out) :: matrix(:, :)
        integer, intent(out) :: status

        real(dp) :: area, determinant
        complex(dp) :: b_matrix(3, 6), d_matrix(3, 3)
        complex(dp) :: lambda, local_matrix(6, 6), mu
        real(dp) :: local_mass(3, 3)
        real(dp) :: dx(3), dy(3), x1, x2, x3, y1, y2, y3
        integer :: element, first_local, local_dofs(6), local_nodes(3)
        integer :: second_local, vertex_count

        matrix = (0.0_dp, 0.0_dp)
        status = 1
        vertex_count = size(vertices, 2)
        if (size(vertices, 1) /= 2 .or. vertex_count < 1) return
        if (size(triangles, 1) /= 3 .or. size(triangles, 2) < 1) return
        if (size(matrix, 1) /= 2*vertex_count .or. &
            size(matrix, 2) /= 2*vertex_count) return
        if (angular_frequency <= 0.0_dp .or. abs(young_modulus) <= 0.0_dp .or. &
            poisson_ratio <= -1.0_dp .or. poisson_ratio >= 0.5_dp .or. &
            solid_density <= 0.0_dp) return
        if (any(triangles < 1) .or. any(triangles > vertex_count)) return

        mu = young_modulus/(2.0_dp*(1.0_dp + poisson_ratio))
        lambda = young_modulus*poisson_ratio/( &
            (1.0_dp + poisson_ratio)*(1.0_dp - 2.0_dp*poisson_ratio))
        d_matrix = (0.0_dp, 0.0_dp)
        d_matrix(1, 1) = lambda + 2.0_dp*mu
        d_matrix(1, 2) = lambda
        d_matrix(2, 1) = lambda
        d_matrix(2, 2) = lambda + 2.0_dp*mu
        d_matrix(3, 3) = mu
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
            dx = [y2 - y3, y3 - y1, y1 - y2]/determinant
            dy = [x3 - x2, x1 - x3, x2 - x1]/determinant
            b_matrix = (0.0_dp, 0.0_dp)
            do first_local = 1, 3
                b_matrix(1, 2*first_local - 1) = dx(first_local)
                b_matrix(2, 2*first_local) = dy(first_local)
                b_matrix(3, 2*first_local - 1) = dy(first_local)
                b_matrix(3, 2*first_local) = dx(first_local)
                local_dofs(2*first_local - 1) = &
                    2*local_nodes(first_local) - 1
                local_dofs(2*first_local) = 2*local_nodes(first_local)
            end do
            local_matrix = area*matmul( &
                transpose(b_matrix), matmul(d_matrix, b_matrix))
            local_mass = area/12.0_dp
            local_mass(1, 1) = area/6.0_dp
            local_mass(2, 2) = area/6.0_dp
            local_mass(3, 3) = area/6.0_dp
            do first_local = 1, 3
                do second_local = 1, 3
                    local_matrix(2*first_local - 1, 2*second_local - 1) = &
                        local_matrix(2*first_local - 1, &
                        2*second_local - 1) - solid_density* &
                        angular_frequency**2* &
                        local_mass(first_local, second_local)
                    local_matrix(2*first_local, 2*second_local) = &
                        local_matrix(2*first_local, 2*second_local) - &
                        solid_density*angular_frequency**2* &
                        local_mass(first_local, second_local)
                end do
            end do
            do first_local = 1, 6
                do second_local = 1, 6
                    matrix(local_dofs(first_local), &
                        local_dofs(second_local)) = &
                        matrix(local_dofs(first_local), &
                        local_dofs(second_local)) + &
                        local_matrix(first_local, second_local)
                end do
            end do
        end do
        status = 0
    end subroutine assemble_elasticity_dynamic_p1_2d

    subroutine solve_dense_complex_dirichlet_fortsparse( &
            matrix, right_hand_side, dirichlet_dofs, dirichlet_values, &
            solution, status)
        complex(dp), intent(in) :: matrix(:, :), right_hand_side(:)
        integer, intent(in) :: dirichlet_dofs(:)
        complex(dp), intent(in) :: dirichlet_values(:)
        complex(dp), intent(out) :: solution(:)
        integer, intent(out) :: status

        type(csc_z_t) :: sparse_matrix
        type(fortsparse_status_t) :: sparse_status
        complex(dp), allocatable :: constrained_matrix(:, :)
        complex(dp), allocatable :: constrained_rhs(:), values(:)
        integer, allocatable :: columns(:), rows(:)
        integer :: column, dof_count, entry, row

        solution = (0.0_dp, 0.0_dp)
        status = 1
        dof_count = size(right_hand_side)
        if (dof_count < 1 .or. size(matrix, 1) /= dof_count .or. &
            size(matrix, 2) /= dof_count .or. size(solution) /= dof_count) &
            return
        if (size(dirichlet_dofs) /= size(dirichlet_values)) return
        if (any(dirichlet_dofs < 1) .or. any(dirichlet_dofs > dof_count)) return

        allocate(constrained_matrix(dof_count, dof_count))
        allocate(constrained_rhs(dof_count))
        constrained_matrix = matrix
        constrained_rhs = right_hand_side
        do entry = 1, size(dirichlet_dofs)
            column = dirichlet_dofs(entry)
            constrained_rhs = constrained_rhs - &
                constrained_matrix(:, column)*dirichlet_values(entry)
            constrained_matrix(:, column) = (0.0_dp, 0.0_dp)
            constrained_matrix(column, :) = (0.0_dp, 0.0_dp)
            constrained_matrix(column, column) = (1.0_dp, 0.0_dp)
            constrained_rhs(column) = dirichlet_values(entry)
        end do

        allocate(rows(dof_count**2), columns(dof_count**2))
        allocate(values(dof_count**2))
        entry = 0
        do column = 1, dof_count
            do row = 1, dof_count
                if (abs(constrained_matrix(row, column)) <= tiny(1.0_dp)) &
                    cycle
                entry = entry + 1
                rows(entry) = row
                columns(entry) = column
                values(entry) = constrained_matrix(row, column)
            end do
        end do
        call csc_from_triplet( &
            dof_count, dof_count, rows(:entry), columns(:entry), &
            values(:entry), sparse_matrix, sparse_status)
        if (sparse_status%code /= 0) return
        call sparse_solve_once( &
            sparse_matrix, constrained_rhs, solution, sparse_status)
        if (sparse_status%code /= 0) return
        status = 0
    end subroutine solve_dense_complex_dirichlet_fortsparse

end module fortfem_elasticity_p1_2d_core
