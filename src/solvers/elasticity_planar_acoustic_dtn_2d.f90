module fortfem_elasticity_planar_acoustic_dtn_2d
    use fortfem_elasticity_p1_2d_core, only: &
        assemble_elasticity_dynamic_p1_2d, &
        solve_dense_complex_dirichlet_fortsparse
    use fortfem_kinds, only: dp
    use fortfem_planar_acoustic_displacement_dtn, only: &
        assemble_planar_acoustic_displacement_dtn_form
    implicit none
    private

    public :: solve_elasticity_planar_acoustic_dtn_p1

    interface solve_elasticity_planar_acoustic_dtn_p1
        module procedure solve_elasticity_planar_acoustic_dtn_p1_real
        module procedure solve_elasticity_planar_acoustic_dtn_p1_complex
    end interface

contains

    subroutine solve_elasticity_planar_acoustic_dtn_p1_real( &
            vertices, triangles, interface_nodes, interface_normals, &
            angular_frequency, sound_speed, fluid_density, period, &
            maximum_mode, young_modulus, poisson_ratio, solid_density, &
            volume_load, incident_pressure, dirichlet_dofs, &
            dirichlet_values, solution, status, absorbing_edges, &
            absorbing_normals)
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
        integer, intent(in), optional :: absorbing_edges(:, :)
        real(dp), intent(in), optional :: absorbing_normals(:, :)

        if (present(absorbing_edges) .neqv. present(absorbing_normals)) then
            solution = cmplx(0.0_dp, 0.0_dp, dp)
            status = 1
        else if (present(absorbing_edges)) then
            call solve_elasticity_planar_acoustic_dtn_p1_complex( &
                vertices, triangles, interface_nodes, interface_normals, &
                angular_frequency, cmplx(sound_speed, 0.0_dp, dp), &
                fluid_density, period, maximum_mode, &
                cmplx(young_modulus, 0.0_dp, dp), poisson_ratio, &
                solid_density, volume_load, incident_pressure, &
                dirichlet_dofs, dirichlet_values, solution, status, &
                absorbing_edges, absorbing_normals)
        else
            call solve_elasticity_planar_acoustic_dtn_p1_complex( &
                vertices, triangles, interface_nodes, interface_normals, &
                angular_frequency, cmplx(sound_speed, 0.0_dp, dp), &
                fluid_density, period, maximum_mode, &
                cmplx(young_modulus, 0.0_dp, dp), poisson_ratio, &
                solid_density, volume_load, incident_pressure, &
                dirichlet_dofs, dirichlet_values, solution, status)
        end if
    end subroutine solve_elasticity_planar_acoustic_dtn_p1_real

    subroutine solve_elasticity_planar_acoustic_dtn_p1_complex( &
            vertices, triangles, interface_nodes, interface_normals, &
            angular_frequency, sound_speed, fluid_density, period, &
            maximum_mode, young_modulus, poisson_ratio, solid_density, &
            volume_load, incident_pressure, dirichlet_dofs, &
            dirichlet_values, solution, status, absorbing_edges, &
            absorbing_normals)
        real(dp), intent(in) :: vertices(:, :)
        integer, intent(in) :: triangles(:, :), interface_nodes(:)
        real(dp), intent(in) :: interface_normals(:, :)
        real(dp), intent(in) :: angular_frequency, fluid_density
        complex(dp), intent(in) :: sound_speed
        real(dp), intent(in) :: period, poisson_ratio
        complex(dp), intent(in) :: young_modulus
        real(dp), intent(in) :: solid_density
        integer, intent(in) :: maximum_mode
        complex(dp), intent(in) :: volume_load(:), incident_pressure(:)
        integer, intent(in) :: dirichlet_dofs(:)
        complex(dp), intent(in) :: dirichlet_values(:)
        complex(dp), intent(out) :: solution(:)
        integer, intent(out) :: status
        integer, intent(in), optional :: absorbing_edges(:, :)
        real(dp), intent(in), optional :: absorbing_normals(:, :)

        complex(dp), allocatable :: fluid_form(:, :), matrix(:, :), rhs(:)
        complex(dp), allocatable :: incident_weak(:)
        integer :: column, component, dof_count, first, node
        integer :: second, vertex_count

        solution = cmplx(0.0_dp, 0.0_dp, dp)
        status = 1
        vertex_count = size(vertices, 2)
        dof_count = 2*vertex_count
        if (.not. valid_inputs()) return

        allocate(matrix(dof_count, dof_count), rhs(dof_count))
        call assemble_elasticity_dynamic_p1_2d( &
            vertices, triangles, angular_frequency, young_modulus, &
            poisson_ratio, solid_density, matrix, status)
        if (status /= 0) return
        if (present(absorbing_edges)) then
            call assemble_absorbing_boundary(matrix, status)
            if (status /= 0) return
        end if
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

        call solve_dense_complex_dirichlet_fortsparse( &
            matrix, rhs, dirichlet_dofs, dirichlet_values, solution, status)

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
            if (present(absorbing_edges) .neqv. &
                present(absorbing_normals)) return
            if (present(absorbing_edges)) then
                if (size(absorbing_edges, 1) /= 2 .or. &
                    size(absorbing_normals, 1) /= 2) return
                if (size(absorbing_edges, 2) /= &
                    size(absorbing_normals, 2)) return
                if (any(absorbing_edges < 1) .or. &
                    any(absorbing_edges > vertex_count)) return
            end if
            if (angular_frequency <= 0.0_dp .or. abs(sound_speed) <= 0.0_dp .or. &
                fluid_density <= 0.0_dp .or. solid_density <= 0.0_dp) return
            if (period <= 0.0_dp .or. abs(young_modulus) <= 0.0_dp) return
            if (poisson_ratio <= -1.0_dp .or. &
                poisson_ratio >= 0.5_dp) return
            if (any(triangles < 1) .or. any(triangles > vertex_count)) return
            if (any(interface_nodes < 1) .or. &
                any(interface_nodes > vertex_count)) return
            if (any(dirichlet_dofs < 1) .or. &
                any(dirichlet_dofs > dof_count)) return
            valid_inputs = .true.
        end function valid_inputs

        subroutine assemble_absorbing_boundary( &
                volume_matrix, operator_status)
            complex(dp), intent(inout) :: volume_matrix(:, :)
            integer, intent(out) :: operator_status

            complex(dp) :: lambda, mu, p_impedance, s_impedance
            complex(dp) :: coefficient(2, 2), local(4, 4)
            real(dp) :: edge_mass(2, 2), length, normal(2)
            integer :: component, edge, first, first_dof, first_node
            integer :: other_component, second, second_dof, second_node

            operator_status = 1
            mu = young_modulus/(2.0_dp*(1.0_dp + poisson_ratio))
            lambda = young_modulus*poisson_ratio/( &
                (1.0_dp + poisson_ratio)*(1.0_dp - 2.0_dp*poisson_ratio))
            p_impedance = sqrt((lambda + 2.0_dp*mu)*solid_density)
            s_impedance = sqrt(mu*solid_density)
            do edge = 1, size(absorbing_edges, 2)
                first_node = absorbing_edges(1, edge)
                second_node = absorbing_edges(2, edge)
                length = norm2(vertices(:, second_node) - &
                    vertices(:, first_node))
                normal = absorbing_normals(:, edge)
                if (length <= tiny(1.0_dp) .or. &
                    abs(norm2(normal) - 1.0_dp) > 1.0e-10_dp) return
                edge_mass = length/6.0_dp*reshape( &
                    [2.0_dp, 1.0_dp, 1.0_dp, 2.0_dp], [2, 2])
                coefficient = cmplx(0.0_dp, 0.0_dp, dp)
                do component = 1, 2
                    coefficient(component, component) = &
                        cmplx(real(s_impedance, dp), 0.0_dp, dp)
                end do
                do first = 1, 2
                    do second = 1, 2
                        coefficient(first, second) = &
                            coefficient(first, second) + &
                            (real(p_impedance, dp) - &
                            real(s_impedance, dp))*normal(first)*normal(second)
                    end do
                end do
                local = cmplx(0.0_dp, 0.0_dp, dp)
                do first = 1, 2
                    do second = 1, 2
                        local(2*first - 1:2*first, &
                            2*second - 1:2*second) = &
                            cmplx(0.0_dp, angular_frequency* &
                            edge_mass(first, second), dp)*coefficient
                    end do
                end do
                do first = 1, 2
                    first_node = absorbing_edges(first, edge)
                    do second = 1, 2
                        second_node = absorbing_edges(second, edge)
                        do component = 1, 2
                            first_dof = 2*first_node - 2 + component
                            do other_component = 1, 2
                                second_dof = 2*second_node - 2 + &
                                    other_component
                                volume_matrix(first_dof, second_dof) = &
                                    volume_matrix(first_dof, second_dof) + &
                                    local(2*first - 2 + component, &
                                    2*second - 2 + other_component)
                            end do
                        end do
                    end do
                end do
            end do
            operator_status = 0
        end subroutine assemble_absorbing_boundary

    end subroutine solve_elasticity_planar_acoustic_dtn_p1_complex

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
