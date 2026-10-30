module fortfem_elasticity_curved_acoustic_ntd_2d
    !! Vector P1 time-harmonic elasticity coupled to the exterior acoustic field
    !! of a closed polygonal interface through the curved displacement-to-pressure
    !! form (Burton-Miller BEM), with JVP and VJP.
    !!
    !! Time dependence is exp(-i omega t): the BEM kernel (i/4) H_0^(1)(k r) is
    !! outgoing, and a lossy solid has Im(E) < 0 (Im(E) > 0 is gain).
    !!
    !! One CSC pattern from fortfem_block_csc holds the triangles (generated
    !! elasticity cell) and one dense block over the 2m interface dofs (the BEM
    !! form, already contracted with the panel normals). Dirichlet values go to
    !! sparse_direct_solve_constrained. Storage is O(triangles + nnz + m**2).
    !!
    !! Derivatives cover the coordinates of vertices off the interface, fluid
    !! density, Young's modulus, Poisson ratio, solid density and incident
    !! pressure. The 2D Helmholtz boundary operators have no derivative in the
    !! panel geometry or the wavenumber, so the interface is held fixed and the
    !! angular frequency and sound speed are not differentiated.
    !! Cotangents follow fortfem_elasticity_acoustic_common_2d.
    use fortfem_block_csc, only: block_csc_add, block_csc_add_group, &
        block_csc_finalize, block_csc_gather, block_csc_init, &
        block_csc_new_matrix, block_csc_pattern_t
    use fortfem_curved_acoustic_displacement_ntd_2d, only: &
        assemble_curved_acoustic_displacement_ntd_form_2d
    use fortfem_elasticity_acoustic_common_2d, only: dirichlet_mask, &
        elastic_material, elastic_material_jvp, elastic_material_t, &
        elastic_material_vjp, uniform_cells_assemble, &
        uniform_cells_assemble_jvp, uniform_cells_assemble_vjp
    use fortfem_elasticity_p1_csc_2d, only: elasticity_p1_cell_dofs
    use fortfem_generated_ntd_panel_coupling, only: generated_ntd_panel_coupling
    use fortfem_kinds, only: dp
    use fortfem_sparse_direct, only: sparse_direct_solve_constrained, &
        sparse_direct_solve_constrained_jvp, sparse_direct_solve_constrained_vjp
    use fortsparse, only: csc_z_t
    implicit none

    private

    public :: solve_elasticity_curved_acoustic_ntd_p1
    public :: solve_elasticity_curved_acoustic_ntd_p1_jvp
    public :: solve_elasticity_curved_acoustic_ntd_p1_vjp

    type :: curved_system_t
        type(block_csc_pattern_t) :: pattern
        integer :: cell_group = 0
        integer :: interface_group = 0
        type(csc_z_t) :: matrix
        complex(dp), allocatable :: rhs(:), fixed(:)
        logical, allocatable :: constrained(:)
        type(elastic_material_t) :: material
        complex(dp), allocatable :: form(:, :, :) ! (2m, 2m, 1)
        real(dp), allocatable :: traction(:, :, :, :) ! (2, 2, 2, panel)
    end type curved_system_t

contains

    subroutine solve_elasticity_curved_acoustic_ntd_p1( &
            vertices, triangles, interface_nodes, panel_nodes, &
            angular_frequency, sound_speed, fluid_density, young_modulus, &
            poisson_ratio, solid_density, volume_load, incident_pressure, &
            dirichlet_dofs, dirichlet_values, solution, status)
        !! Time dependence exp(-i omega t): Im(E) < 0 is lossy.
        !! status: 0 success, 1 invalid input or panel, 2 failed BEM or sparse solve.
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

        type(curved_system_t) :: system

        solution = (0.0_dp, 0.0_dp)
        status = 1
        if (size(solution) /= 2*size(vertices, 2)) return
        call build_system(vertices, triangles, interface_nodes, panel_nodes, &
            angular_frequency, sound_speed, fluid_density, young_modulus, &
            poisson_ratio, solid_density, volume_load, incident_pressure, &
            dirichlet_dofs, dirichlet_values, system, status)
        if (status /= 0) return
        call sparse_direct_solve_constrained(system%matrix, system%rhs, &
            system%constrained, system%fixed, solution, status)
        if (status /= 0) then
            solution = (0.0_dp, 0.0_dp)
            status = 2
        end if
    end subroutine solve_elasticity_curved_acoustic_ntd_p1

    subroutine solve_elasticity_curved_acoustic_ntd_p1_jvp( &
            vertices, triangles, interface_nodes, panel_nodes, &
            angular_frequency, sound_speed, fluid_density, young_modulus, &
            poisson_ratio, solid_density, volume_load, incident_pressure, &
            dirichlet_dofs, dirichlet_values, vertices_dot, fluid_density_dot, &
            young_modulus_dot, poisson_ratio_dot, solid_density_dot, &
            incident_pressure_dot, solution_dot, status)
        !! Directional derivative of the solution. status 3: vertices_dot moves
        !! an interface vertex (not supported); otherwise as for the solve.
        real(dp), intent(in) :: vertices(:, :)
        integer, intent(in) :: triangles(:, :), interface_nodes(:)
        integer, intent(in) :: panel_nodes(:, :)
        real(dp), intent(in) :: angular_frequency, sound_speed, fluid_density
        complex(dp), intent(in) :: young_modulus
        real(dp), intent(in) :: poisson_ratio, solid_density
        complex(dp), intent(in) :: volume_load(:), incident_pressure(:)
        integer, intent(in) :: dirichlet_dofs(:)
        complex(dp), intent(in) :: dirichlet_values(:)
        real(dp), intent(in) :: vertices_dot(:, :), fluid_density_dot
        complex(dp), intent(in) :: young_modulus_dot
        real(dp), intent(in) :: poisson_ratio_dot, solid_density_dot
        complex(dp), intent(in) :: incident_pressure_dot(:)
        complex(dp), intent(out) :: solution_dot(:)
        integer, intent(out) :: status

        type(curved_system_t) :: system
        type(elastic_material_t) :: material_dot
        complex(dp), allocatable :: values_dot(:), rhs_dot(:)

        solution_dot = (0.0_dp, 0.0_dp)
        status = 1
        if (size(solution_dot) /= 2*size(vertices, 2)) return
        if (any(shape(vertices_dot) /= shape(vertices))) return
        if (size(incident_pressure_dot) /= size(incident_pressure)) return
        call build_system(vertices, triangles, interface_nodes, panel_nodes, &
            angular_frequency, sound_speed, fluid_density, young_modulus, &
            poisson_ratio, solid_density, volume_load, incident_pressure, &
            dirichlet_dofs, dirichlet_values, system, status)
        if (status /= 0) return
        if (any(vertices_dot(:, interface_nodes) /= 0.0_dp)) then
            status = 3
            return
        end if

        call elastic_material_jvp(young_modulus, poisson_ratio, solid_density, &
            angular_frequency, young_modulus_dot, poisson_ratio_dot, &
            solid_density_dot, 0.0_dp, material_dot)
        allocate(values_dot(system%pattern%nnz))
        values_dot = (0.0_dp, 0.0_dp)
        call uniform_cells_assemble_jvp(system%pattern, system%cell_group, &
            vertices, triangles, system%material, vertices_dot, material_dot, &
            values_dot, status)
        if (status /= 0) return
        ! The form is linear in the fluid density.
        call block_csc_add(system%pattern, system%interface_group, &
            system%form*(fluid_density_dot/fluid_density), values_dot)
        allocate(rhs_dot(size(system%rhs)))
        rhs_dot = (0.0_dp, 0.0_dp)
        call add_incident_load(system%traction, panel_nodes, interface_nodes, &
            incident_pressure_dot, rhs_dot)
        call sparse_direct_solve_constrained_jvp(system%matrix, system%rhs, &
            system%constrained, system%fixed, values_dot, rhs_dot, &
            0.0_dp*system%fixed, solution_dot, status)
        if (status /= 0) then
            solution_dot = (0.0_dp, 0.0_dp)
            status = 2
        end if
    end subroutine solve_elasticity_curved_acoustic_ntd_p1_jvp

    subroutine solve_elasticity_curved_acoustic_ntd_p1_vjp( &
            vertices, triangles, interface_nodes, panel_nodes, &
            angular_frequency, sound_speed, fluid_density, young_modulus, &
            poisson_ratio, solid_density, volume_load, incident_pressure, &
            dirichlet_dofs, dirichlet_values, solution_bar, vertices_bar, &
            fluid_density_bar, young_modulus_bar, poisson_ratio_bar, &
            solid_density_bar, incident_pressure_bar, status)
        !! Cotangents of the inputs of the JVP for solution_bar. The interface is
        !! held fixed, so the cotangents of interface vertices are zero.
        real(dp), intent(in) :: vertices(:, :)
        integer, intent(in) :: triangles(:, :), interface_nodes(:)
        integer, intent(in) :: panel_nodes(:, :)
        real(dp), intent(in) :: angular_frequency, sound_speed, fluid_density
        complex(dp), intent(in) :: young_modulus
        real(dp), intent(in) :: poisson_ratio, solid_density
        complex(dp), intent(in) :: volume_load(:), incident_pressure(:)
        integer, intent(in) :: dirichlet_dofs(:)
        complex(dp), intent(in) :: dirichlet_values(:), solution_bar(:)
        real(dp), intent(out) :: vertices_bar(:, :), fluid_density_bar
        complex(dp), intent(out) :: young_modulus_bar
        real(dp), intent(out) :: poisson_ratio_bar, solid_density_bar
        complex(dp), intent(out) :: incident_pressure_bar(:)
        integer, intent(out) :: status

        type(curved_system_t) :: system
        type(elastic_material_t) :: material_bar
        complex(dp), allocatable :: solution(:), values_bar(:), rhs_bar(:)
        complex(dp), allocatable :: fixed_bar(:), form_bar(:, :, :)
        real(dp) :: omega_bar
        integer :: c, dof, first, node, panel, second

        vertices_bar = 0.0_dp
        fluid_density_bar = 0.0_dp
        young_modulus_bar = (0.0_dp, 0.0_dp)
        poisson_ratio_bar = 0.0_dp
        solid_density_bar = 0.0_dp
        incident_pressure_bar = (0.0_dp, 0.0_dp)
        status = 1
        if (size(solution_bar) /= 2*size(vertices, 2)) return
        if (any(shape(vertices_bar) /= shape(vertices))) return
        if (size(incident_pressure_bar) /= size(incident_pressure)) return
        call build_system(vertices, triangles, interface_nodes, panel_nodes, &
            angular_frequency, sound_speed, fluid_density, young_modulus, &
            poisson_ratio, solid_density, volume_load, incident_pressure, &
            dirichlet_dofs, dirichlet_values, system, status)
        if (status /= 0) return
        allocate(solution(size(solution_bar)), values_bar(system%pattern%nnz))
        allocate(rhs_bar(size(solution_bar)), fixed_bar(size(solution_bar)))
        call sparse_direct_solve_constrained(system%matrix, system%rhs, &
            system%constrained, system%fixed, solution, status)
        if (status == 0) call sparse_direct_solve_constrained_vjp( &
            system%matrix, system%rhs, system%constrained, system%fixed, &
            solution, solution_bar, values_bar, rhs_bar, fixed_bar, status)
        if (status /= 0) then
            status = 2
            return
        end if

        call uniform_cells_assemble_vjp(system%pattern, system%cell_group, &
            vertices, triangles, system%material, values_bar, vertices_bar, &
            material_bar, status)
        if (status /= 0) return
        vertices_bar(:, interface_nodes) = 0.0_dp
        call elastic_material_vjp(young_modulus, poisson_ratio, solid_density, &
            angular_frequency, material_bar, young_modulus_bar, &
            poisson_ratio_bar, solid_density_bar, omega_bar)

        allocate(form_bar(size(system%form, 1), size(system%form, 2), 1))
        call block_csc_gather(system%pattern, system%interface_group, values_bar, &
            form_bar)
        fluid_density_bar = real(sum(conjg(system%form)*form_bar), dp)/fluid_density

        ! rhs = load - T p with real T, so p_bar = -T^T rhs_bar.
        do panel = 1, size(panel_nodes, 2)
            do first = 1, 2
                do second = 1, 2
                    do c = 1, 2
                        dof = 2*interface_nodes(panel_nodes(first, panel)) - 2 + c
                        node = panel_nodes(second, panel)
                        incident_pressure_bar(node) = incident_pressure_bar(node) - &
                            system%traction(c, first, second, panel)*rhs_bar(dof)
                    end do
                end do
            end do
        end do
    end subroutine solve_elasticity_curved_acoustic_ntd_p1_vjp

    subroutine build_system(vertices, triangles, interface_nodes, panel_nodes, &
            angular_frequency, sound_speed, fluid_density, young_modulus, &
            poisson_ratio, solid_density, volume_load, incident_pressure, &
            dirichlet_dofs, dirichlet_values, system, status)
        real(dp), intent(in) :: vertices(:, :)
        integer, intent(in) :: triangles(:, :), interface_nodes(:)
        integer, intent(in) :: panel_nodes(:, :)
        real(dp), intent(in) :: angular_frequency, sound_speed, fluid_density
        complex(dp), intent(in) :: young_modulus
        real(dp), intent(in) :: poisson_ratio, solid_density
        complex(dp), intent(in) :: volume_load(:), incident_pressure(:)
        integer, intent(in) :: dirichlet_dofs(:)
        complex(dp), intent(in) :: dirichlet_values(:)
        type(curved_system_t), intent(out) :: system
        integer, intent(out) :: status

        real(dp), allocatable :: panel_end(:, :), panel_start(:, :)
        integer, allocatable :: cell_dofs(:, :), interface_dofs(:, :)
        real(dp) :: delta(2), endpoint_mass, length, trace(2)
        integer :: dof_count, panel, panel_count

        status = 1
        dof_count = 2*size(vertices, 2)
        panel_count = size(interface_nodes)
        if (.not. valid_inputs()) return

        allocate(panel_start(2, panel_count), panel_end(2, panel_count))
        allocate(system%traction(2, 2, 2, panel_count))
        do panel = 1, panel_count
            panel_start(:, panel) = vertices(:, &
                interface_nodes(panel_nodes(1, panel)))
            panel_end(:, panel) = vertices(:, &
                interface_nodes(panel_nodes(2, panel)))
            delta = panel_end(:, panel) - panel_start(:, panel)
            call generated_ntd_panel_coupling(delta(1), delta(2), length, &
                endpoint_mass, trace, system%traction(:, :, :, panel))
        end do

        call elastic_material(young_modulus, poisson_ratio, solid_density, &
            angular_frequency, system%material)
        call elasticity_p1_cell_dofs(triangles, cell_dofs)
        call block_csc_init(system%pattern, dof_count)
        call block_csc_add_group(system%pattern, cell_dofs, system%cell_group)
        allocate(interface_dofs(2*panel_count, 1))
        interface_dofs(1::2, 1) = 2*interface_nodes - 1
        interface_dofs(2::2, 1) = 2*interface_nodes
        call block_csc_add_group(system%pattern, interface_dofs, &
            system%interface_group)
        call block_csc_finalize(system%pattern)
        call block_csc_new_matrix(system%pattern, system%matrix)
        call uniform_cells_assemble(system%pattern, system%cell_group, vertices, &
            triangles, system%material, system%matrix%val, status)
        if (status /= 0) return

        allocate(system%form(2*panel_count, 2*panel_count, 1))
        call assemble_curved_acoustic_displacement_ntd_form_2d( &
            panel_start, panel_end, panel_nodes, &
            angular_frequency/sound_speed, fluid_density, angular_frequency, &
            16, system%form(:, :, 1), status)
        if (status /= 0) return
        call block_csc_add(system%pattern, system%interface_group, system%form, &
            system%matrix%val)

        system%rhs = volume_load
        call add_incident_load(system%traction, panel_nodes, interface_nodes, &
            incident_pressure, system%rhs)
        call dirichlet_mask(dof_count, dirichlet_dofs, dirichlet_values, &
            system%constrained, system%fixed)
        status = 0

    contains

        logical function valid_inputs()
            valid_inputs = .false.
            if (size(vertices, 1) /= 2 .or. dof_count < 2) return
            if (size(triangles, 1) /= 3 .or. size(triangles, 2) < 1) return
            if (panel_count < 3 .or. size(panel_nodes, 1) /= 2 .or. &
                size(panel_nodes, 2) /= panel_count) return
            if (size(incident_pressure) /= panel_count) return
            if (size(volume_load) /= dof_count) return
            if (size(dirichlet_dofs) /= size(dirichlet_values)) return
            if (any(triangles < 1) .or. any(triangles > dof_count/2)) return
            if (any(interface_nodes < 1) .or. &
                any(interface_nodes > dof_count/2)) return
            if (any(panel_nodes < 1) .or. any(panel_nodes > panel_count)) return
            if (any(dirichlet_dofs < 1) .or. &
                any(dirichlet_dofs > dof_count)) return
            if (angular_frequency <= 0.0_dp .or. sound_speed <= 0.0_dp .or. &
                fluid_density <= 0.0_dp .or. solid_density <= 0.0_dp) return
            if (abs(young_modulus) <= 0.0_dp .or. &
                poisson_ratio <= -1.0_dp .or. poisson_ratio >= 0.5_dp) return
            valid_inputs = .true.
        end function valid_inputs

    end subroutine build_system

    pure subroutine add_incident_load(traction, panel_nodes, interface_nodes, &
            pressure, rhs)
        !! rhs -= P1 traction mass times the incident pressure (it pushes against
        !! the outward normal).
        real(dp), intent(in) :: traction(:, :, :, :)
        integer, intent(in) :: panel_nodes(:, :), interface_nodes(:)
        complex(dp), intent(in) :: pressure(:)
        complex(dp), intent(inout) :: rhs(:)

        integer :: c, dof, first, panel, second

        do panel = 1, size(panel_nodes, 2)
            do first = 1, 2
                do second = 1, 2
                    do c = 1, 2
                        dof = 2*interface_nodes(panel_nodes(first, panel)) - 2 + c
                        rhs(dof) = rhs(dof) - traction(c, first, second, panel)* &
                            pressure(panel_nodes(second, panel))
                    end do
                end do
            end do
        end do
    end subroutine add_incident_load

end module fortfem_elasticity_curved_acoustic_ntd_2d
