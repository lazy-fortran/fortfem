module fortfem_elasticity_curved_acoustic_ntd_2d
    use fortfem_curved_acoustic_displacement_ntd_2d, only: &
        assemble_curved_acoustic_displacement_ntd_form_2d
    use fortfem_elasticity_p1_2d_core, only: &
        assemble_elasticity_dynamic_p1_2d, &
        solve_dense_complex_dirichlet_fortsparse
    use fortfem_kinds, only: dp
    implicit none

    private

    public :: solve_elasticity_curved_acoustic_ntd_p1

contains

    subroutine solve_elasticity_curved_acoustic_ntd_p1( &
            vertices, triangles, interface_nodes, panel_nodes, &
            angular_frequency, sound_speed, fluid_density, young_modulus, &
            poisson_ratio, solid_density, volume_load, incident_pressure, &
            dirichlet_dofs, dirichlet_values, solution, status)
        real(dp), intent(in) :: vertices(:, :)
        integer, intent(in) :: triangles(:, :), interface_nodes(:)
        integer, intent(in) :: panel_nodes(:, :)
        real(dp), intent(in) :: angular_frequency, sound_speed, fluid_density
        complex(dp), intent(in) :: young_modulus
        real(dp), intent(in) :: poisson_ratio, solid_density
        complex(dp), intent(in) :: volume_load(:), incident_pressure(:)
        integer, intent(in) :: dirichlet_dofs(:)
        complex(dp), intent(in) :: dirichlet_values(:)
        complex(dp), intent(out) :: solution(:)
        integer, intent(out) :: status

        complex(dp), allocatable :: boundary_form(:, :), matrix(:, :), rhs(:)
        real(dp), allocatable :: panel_end(:, :), panel_start(:, :)
        real(dp) :: length, normal(2)
        integer :: component, dof_count, first, first_global, first_local
        integer :: panel, panel_count, second, second_global, second_local
        integer :: vertex_count

        solution = (0.0_dp, 0.0_dp)
        status = 1
        vertex_count = size(vertices, 2)
        dof_count = 2*vertex_count
        panel_count = size(interface_nodes)
        if (.not. valid_inputs()) return

        allocate(panel_start(2, panel_count), panel_end(2, panel_count))
        do panel = 1, panel_count
            panel_start(:, panel) = vertices(:, &
                interface_nodes(panel_nodes(1, panel)))
            panel_end(:, panel) = vertices(:, &
                interface_nodes(panel_nodes(2, panel)))
        end do

        allocate(matrix(dof_count, dof_count), rhs(dof_count))
        call assemble_elasticity_dynamic_p1_2d( &
            vertices, triangles, angular_frequency, young_modulus, &
            poisson_ratio, solid_density, matrix, status)
        if (status /= 0) return
        allocate(boundary_form(2*panel_count, 2*panel_count))
        call assemble_curved_acoustic_displacement_ntd_form_2d( &
            panel_start, panel_end, panel_nodes, &
            angular_frequency/sound_speed, fluid_density, angular_frequency, &
            16, boundary_form, status)
        if (status /= 0) return

        do first = 1, panel_count
            do component = 1, 2
                first_local = 2*first - 2 + component
                first_global = 2*interface_nodes(first) - 2 + component
                do second = 1, panel_count
                    do second_local = 1, 2
                        second_global = 2*interface_nodes(second) - 2 + &
                            second_local
                        matrix(first_global, second_global) = &
                            matrix(first_global, second_global) + &
                            boundary_form(first_local, &
                            2*second - 2 + second_local)
                    end do
                end do
            end do
        end do

        rhs = volume_load
        do panel = 1, panel_count
            length = norm2(panel_end(:, panel) - panel_start(:, panel))
            normal = [ &
                panel_end(2, panel) - panel_start(2, panel), &
                panel_start(1, panel) - panel_end(1, panel)]/length
            call subtract_incident_panel_load( &
                panel, length, normal, panel_nodes, interface_nodes, &
                incident_pressure, rhs)
        end do
        call solve_dense_complex_dirichlet_fortsparse( &
            matrix, rhs, dirichlet_dofs, dirichlet_values, solution, status)

    contains

        logical function valid_inputs()
            valid_inputs = .false.
            if (size(vertices, 1) /= 2 .or. vertex_count < 1) return
            if (size(triangles, 1) /= 3 .or. size(triangles, 2) < 1) return
            if (panel_count < 3 .or. size(panel_nodes, 1) /= 2 .or. &
                size(panel_nodes, 2) /= panel_count) return
            if (size(incident_pressure) /= panel_count) return
            if (size(volume_load) /= dof_count .or. &
                size(solution) /= dof_count) return
            if (size(dirichlet_dofs) /= size(dirichlet_values)) return
            if (any(triangles < 1) .or. any(triangles > vertex_count)) return
            if (any(interface_nodes < 1) .or. &
                any(interface_nodes > vertex_count)) return
            if (any(panel_nodes < 1) .or. any(panel_nodes > panel_count)) return
            if (any(dirichlet_dofs < 1) .or. &
                any(dirichlet_dofs > dof_count)) return
            if (angular_frequency <= 0.0_dp .or. sound_speed <= 0.0_dp .or. &
                fluid_density <= 0.0_dp .or. solid_density <= 0.0_dp) return
            if (abs(young_modulus) <= 0.0_dp .or. &
                poisson_ratio <= -1.0_dp .or. poisson_ratio >= 0.5_dp) return
            valid_inputs = .true.
        end function valid_inputs

    end subroutine solve_elasticity_curved_acoustic_ntd_p1

    pure subroutine subtract_incident_panel_load( &
            panel, length, normal, panel_nodes, interface_nodes, pressure, rhs)
        integer, intent(in) :: panel
        real(dp), intent(in) :: length, normal(2)
        integer, intent(in) :: panel_nodes(:, :), interface_nodes(:)
        complex(dp), intent(in) :: pressure(:)
        complex(dp), intent(inout) :: rhs(:)

        integer :: component, first, global_dof, second

        do first = 1, 2
            do second = 1, 2
                do component = 1, 2
                    global_dof = 2*interface_nodes( &
                        panel_nodes(first, panel)) - 2 + component
                    rhs(global_dof) = rhs(global_dof) - normal(component)* &
                        length*merge(2.0_dp, 1.0_dp, first == second)/6.0_dp* &
                        pressure(panel_nodes(second, panel))
                end do
            end do
        end do
    end subroutine subtract_incident_panel_load

end module fortfem_elasticity_curved_acoustic_ntd_2d
