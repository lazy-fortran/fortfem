module fortfem_elasticity_planar_acoustic_dtn_2d
    use fortfem_kinds, only: dp
    use fortfem_planar_acoustic_displacement_dtn, only: &
        assemble_planar_acoustic_displacement_dtn_form
    use fortsparse, only: csc_from_triplet, csc_z_t, fortsparse_status_t, &
        sparse_solve_once
    implicit none
    private

    public :: solve_elasticity_planar_acoustic_dtn_p1

contains

    subroutine solve_elasticity_planar_acoustic_dtn_p1( &
            vertices, triangles, interface_nodes, interface_normals, &
            angular_frequency, sound_speed, fluid_density, period, &
            maximum_mode, young_modulus, poisson_ratio, solid_density, &
            volume_load, incident_pressure, dirichlet_dofs, &
            dirichlet_values, solution, status)
        real(dp), intent(in) :: vertices(:, :)
        integer, intent(in) :: triangles(:, :), interface_nodes(:)
        real(dp), intent(in) :: interface_normals(:, :)
        real(dp), intent(in) :: angular_frequency, sound_speed, fluid_density
        real(dp), intent(in) :: period, young_modulus, poisson_ratio
        real(dp), intent(in) :: solid_density
        integer, intent(in) :: maximum_mode
        complex(dp), intent(in) :: volume_load(:), incident_pressure(:)
        integer, intent(in) :: dirichlet_dofs(:)
        complex(dp), intent(in) :: dirichlet_values(:)
        complex(dp), intent(out) :: solution(:)
        integer, intent(out) :: status

        type(csc_z_t) :: sparse_matrix
        type(fortsparse_status_t) :: sparse_status
        complex(dp), allocatable :: fluid_form(:, :), matrix(:, :), rhs(:)
        complex(dp), allocatable :: incident_weak(:), triplet_values(:)
        integer, allocatable :: columns(:), rows(:)
        integer :: column, component, dof_count, entry, first, node
        integer :: row, second, vertex_count

        solution = cmplx(0.0_dp, 0.0_dp, dp)
        status = 1
        vertex_count = size(vertices, 2)
        dof_count = 2*vertex_count
        if (.not. valid_inputs()) return

        allocate(matrix(dof_count, dof_count), rhs(dof_count))
        call assemble_elastic_volume_matrix(matrix, status)
        if (status /= 0) return
        allocate(fluid_form(size(interface_nodes), size(interface_nodes)))
        call assemble_planar_acoustic_displacement_dtn_form( &
            size(interface_nodes), angular_frequency, sound_speed, &
            fluid_density, period, maximum_mode, fluid_form, status)
        if (status /= 0) return
        do first = 1, size(interface_nodes)
            do second = 1, size(interface_nodes)
                do component = 1, 2
                    do column = 1, 2
                        matrix(2*interface_nodes(first) - 2 + component, &
                            2*interface_nodes(second) - 2 + column) = &
                            matrix(2*interface_nodes(first) - 2 + component, &
                            2*interface_nodes(second) - 2 + column) + &
                            interface_normals(component, first)* &
                            fluid_form(first, second)* &
                            interface_normals(column, second)
                    end do
                end do
            end do
        end do
        allocate(incident_weak(size(interface_nodes)))
        call apply_periodic_p1_mass(incident_pressure, period, incident_weak)
        rhs = volume_load
        do first = 1, size(interface_nodes)
            do component = 1, 2
                node = 2*interface_nodes(first) - 2 + component
                rhs(node) = rhs(node) - &
                    interface_normals(component, first)*incident_weak(first)
            end do
        end do

        do entry = 1, size(dirichlet_dofs)
            node = dirichlet_dofs(entry)
            rhs = rhs - matrix(:, node)*dirichlet_values(entry)
            matrix(:, node) = cmplx(0.0_dp, 0.0_dp, dp)
            matrix(node, :) = cmplx(0.0_dp, 0.0_dp, dp)
            matrix(node, node) = cmplx(1.0_dp, 0.0_dp, dp)
            rhs(node) = dirichlet_values(entry)
        end do

        allocate(rows(dof_count**2), columns(dof_count**2))
        allocate(triplet_values(dof_count**2))
        entry = 0
        do column = 1, dof_count
            do row = 1, dof_count
                if (abs(matrix(row, column)) <= tiny(1.0_dp)) cycle
                entry = entry + 1
                rows(entry) = row
                columns(entry) = column
                triplet_values(entry) = matrix(row, column)
            end do
        end do
        call csc_from_triplet( &
            dof_count, dof_count, rows(:entry), columns(:entry), &
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
            if (size(interface_nodes) < 3) return
            if (size(interface_normals, 1) /= 2 .or. &
                size(interface_normals, 2) /= size(interface_nodes)) return
            if (size(incident_pressure) /= size(interface_nodes)) return
            if (size(volume_load) /= dof_count .or. &
                size(solution) /= dof_count) return
            if (size(dirichlet_dofs) /= size(dirichlet_values)) return
            if (angular_frequency <= 0.0_dp .or. sound_speed <= 0.0_dp .or. &
                fluid_density <= 0.0_dp .or. solid_density <= 0.0_dp) return
            if (period <= 0.0_dp .or. young_modulus <= 0.0_dp) return
            if (poisson_ratio <= -1.0_dp .or. &
                poisson_ratio >= 0.5_dp) return
            if (any(triangles < 1) .or. any(triangles > vertex_count)) return
            if (any(interface_nodes < 1) .or. &
                any(interface_nodes > vertex_count)) return
            if (any(dirichlet_dofs < 1) .or. &
                any(dirichlet_dofs > dof_count)) return
            valid_inputs = .true.
        end function valid_inputs

        subroutine assemble_elastic_volume_matrix( &
                volume_matrix, operator_status)
            complex(dp), intent(out) :: volume_matrix(:, :)
            integer, intent(out) :: operator_status

            real(dp) :: area, determinant, lambda, mu
            real(dp) :: b_matrix(3, 6), d_matrix(3, 3)
            real(dp) :: local_matrix(6, 6), local_mass(3, 3)
            real(dp) :: dx(3), dy(3), x1, x2, x3, y1, y2, y3
            integer :: element, first_local, local_dofs(6), local_nodes(3)
            integer :: second_local

            volume_matrix = cmplx(0.0_dp, 0.0_dp, dp)
            operator_status = 1
            mu = young_modulus/(2.0_dp*(1.0_dp + poisson_ratio))
            lambda = young_modulus*poisson_ratio/( &
                (1.0_dp + poisson_ratio)*(1.0_dp - 2.0_dp*poisson_ratio))
            d_matrix = 0.0_dp
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
                b_matrix = 0.0_dp
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
                        local_matrix(2*first_local - 1, &
                            2*second_local - 1) = &
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
                        volume_matrix(local_dofs(first_local), &
                            local_dofs(second_local)) = &
                            volume_matrix(local_dofs(first_local), &
                            local_dofs(second_local)) + &
                            local_matrix(first_local, second_local)
                    end do
                end do
            end do
            operator_status = 0
        end subroutine assemble_elastic_volume_matrix
    end subroutine solve_elasticity_planar_acoustic_dtn_p1

    pure subroutine apply_periodic_p1_mass(values, period, weighted_values)
        complex(dp), intent(in) :: values(:)
        real(dp), intent(in) :: period
        complex(dp), intent(out) :: weighted_values(:)

        real(dp) :: spacing
        integer :: next, previous, row

        spacing = period/real(size(values), dp)
        do row = 1, size(values)
            previous = modulo(row - 2, size(values)) + 1
            next = modulo(row, size(values)) + 1
            weighted_values(row) = spacing/6.0_dp*( &
                values(previous) + 4.0_dp*values(row) + values(next))
        end do
    end subroutine apply_periodic_p1_mass
end module fortfem_elasticity_planar_acoustic_dtn_2d
