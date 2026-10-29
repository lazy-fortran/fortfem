module fortfem_elasticity_planar_acoustic_dtn_2d
    !! Vector P1 time-harmonic elasticity coupled to a periodic acoustic half-space
    !! through the planar displacement-to-pressure (DtN) form, with JVP and VJP.
    !!
    !! One CSC pattern from fortfem_block_csc holds three groups: the triangles
    !! (generated elasticity cell), the optional absorbing edges (generated
    !! impedance edge with the real parts of the P and S impedances) and one dense
    !! block over the 2m interface dofs, the m x m acoustic form contracted with
    !! the interface normals. Dirichlet values go to sparse_direct_solve_constrained.
    !! Storage is O(triangles + nnz + m**2).
    !!
    !! Derivatives cover vertex coordinates, fluid density, Young's modulus,
    !! Poisson ratio, solid density and incident pressure. The angular frequency
    !! and sound speed also enter the DtN mode multiplier, whose generated JVP and
    !! VJP (gen_planar_acoustic_dtn_products) await a FortSym emitter fix for
    !! the constant i, so they are not differentiated yet.
    !! Cotangents follow fortfem_elasticity_acoustic_common_2d.
    use fortfem_block_csc, only: block_csc_add, block_csc_add_group, &
        block_csc_finalize, block_csc_gather, block_csc_init, &
        block_csc_new_matrix, block_csc_pattern_t
    use fortfem_elasticity_acoustic_common_2d, only: dirichlet_mask, &
        elastic_material, elastic_material_jvp, elastic_material_t, &
        elastic_material_vjp, uniform_cells_assemble, &
        uniform_cells_assemble_jvp, uniform_cells_assemble_vjp
    use fortfem_elasticity_p1_csc_2d, only: elasticity_p1_cell_dofs
    use fortfem_generated_interface_normal_block, only: &
        generated_interface_normal_block
    use fortfem_generated_interface_normal_block_jvp, only: &
        generated_interface_normal_block_jvp
    use fortfem_generated_interface_normal_block_vjp, only: &
        generated_interface_normal_block_vjp
    use fortfem_generated_p1_elastic_impedance_edge, only: &
        generated_p1_elastic_impedance_edge
    use fortfem_generated_p1_elastic_impedance_edge_jvp, only: &
        generated_p1_elastic_impedance_edge_jvp
    use fortfem_generated_p1_elastic_impedance_edge_vjp, only: &
        generated_p1_elastic_impedance_edge_vjp
    use fortfem_generated_planar_incident_load, only: &
        generated_planar_incident_load
    use fortfem_generated_planar_incident_load_jvp, only: &
        generated_planar_incident_load_jvp
    use fortfem_generated_planar_incident_load_vjp, only: &
        generated_planar_incident_load_vjp
    use fortfem_generated_segment_length, only: generated_segment_length
    use fortfem_generated_segment_length_jvp, only: generated_segment_length_jvp
    use fortfem_generated_segment_length_vjp, only: generated_segment_length_vjp
    use fortfem_kinds, only: dp
    use fortfem_planar_acoustic_displacement_dtn, only: &
        assemble_planar_acoustic_displacement_dtn_form
    use fortfem_sparse_direct, only: sparse_direct_solve_constrained, &
        sparse_direct_solve_constrained_jvp, sparse_direct_solve_constrained_vjp
    use fortsparse, only: csc_z_t
    implicit none
    private

    public :: solve_elasticity_planar_acoustic_dtn_p1
    public :: solve_elasticity_planar_acoustic_dtn_p1_jvp
    public :: solve_elasticity_planar_acoustic_dtn_p1_vjp

    interface solve_elasticity_planar_acoustic_dtn_p1
        module procedure solve_elasticity_planar_acoustic_dtn_p1_real
        module procedure solve_elasticity_planar_acoustic_dtn_p1_complex
    end interface

    type :: planar_system_t
        type(block_csc_pattern_t) :: pattern
        integer :: cell_group = 0
        integer :: edge_group = 0
        integer :: interface_group = 0
        type(csc_z_t) :: matrix
        complex(dp), allocatable :: rhs(:), fixed(:)
        logical, allocatable :: constrained(:)
        type(elastic_material_t) :: material
        real(dp), allocatable :: lengths(:)
        complex(dp), allocatable :: form(:, :)
    end type planar_system_t

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

        call solve_elasticity_planar_acoustic_dtn_p1_complex( &
            vertices, triangles, interface_nodes, interface_normals, &
            angular_frequency, cmplx(sound_speed, 0.0_dp, dp), &
            fluid_density, period, maximum_mode, &
            cmplx(young_modulus, 0.0_dp, dp), poisson_ratio, &
            solid_density, volume_load, incident_pressure, &
            dirichlet_dofs, dirichlet_values, solution, status, &
            absorbing_edges, absorbing_normals)
    end subroutine solve_elasticity_planar_acoustic_dtn_p1_real

    subroutine solve_elasticity_planar_acoustic_dtn_p1_complex( &
            vertices, triangles, interface_nodes, interface_normals, &
            angular_frequency, sound_speed, fluid_density, period, &
            maximum_mode, young_modulus, poisson_ratio, solid_density, &
            volume_load, incident_pressure, dirichlet_dofs, &
            dirichlet_values, solution, status, absorbing_edges, &
            absorbing_normals)
        !! status: 0 success, 1 invalid input or failed form, 2 failed solve.
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

        type(planar_system_t) :: system

        solution = (0.0_dp, 0.0_dp)
        status = 1
        if (size(solution) /= 2*size(vertices, 2)) return
        call build_system(vertices, triangles, interface_nodes, &
            interface_normals, angular_frequency, sound_speed, fluid_density, &
            period, maximum_mode, young_modulus, poisson_ratio, solid_density, &
            volume_load, incident_pressure, dirichlet_dofs, dirichlet_values, &
            system, status, absorbing_edges, absorbing_normals)
        if (status /= 0) return
        call sparse_direct_solve_constrained(system%matrix, system%rhs, &
            system%constrained, system%fixed, solution, status)
        if (status /= 0) then
            solution = (0.0_dp, 0.0_dp)
            status = 2
        end if
    end subroutine solve_elasticity_planar_acoustic_dtn_p1_complex

    subroutine solve_elasticity_planar_acoustic_dtn_p1_jvp( &
            vertices, triangles, interface_nodes, interface_normals, &
            angular_frequency, sound_speed, fluid_density, period, &
            maximum_mode, young_modulus, poisson_ratio, solid_density, &
            volume_load, incident_pressure, dirichlet_dofs, dirichlet_values, &
            vertices_dot, fluid_density_dot, young_modulus_dot, &
            poisson_ratio_dot, solid_density_dot, incident_pressure_dot, &
            solution_dot, status, absorbing_edges, absorbing_normals)
        !! Directional derivative of the solution of the complex interface.
        !! status as for the solve.
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
        real(dp), intent(in) :: vertices_dot(:, :), fluid_density_dot
        complex(dp), intent(in) :: young_modulus_dot
        real(dp), intent(in) :: poisson_ratio_dot, solid_density_dot
        complex(dp), intent(in) :: incident_pressure_dot(:)
        complex(dp), intent(out) :: solution_dot(:)
        integer, intent(out) :: status
        integer, intent(in), optional :: absorbing_edges(:, :)
        real(dp), intent(in), optional :: absorbing_normals(:, :)

        type(planar_system_t) :: system
        type(elastic_material_t) :: material_dot
        complex(dp), allocatable :: values_dot(:), rhs_dot(:), edge_dot(:, :, :)
        complex(dp), allocatable :: block_dot(:, :, :), form_dot(:, :)
        real(dp) :: delta(2), delta_dot(2), length_dot
        integer :: a, b, e, m, previous, next

        solution_dot = (0.0_dp, 0.0_dp)
        status = 1
        if (size(solution_dot) /= 2*size(vertices, 2)) return
        if (any(shape(vertices_dot) /= shape(vertices))) return
        if (size(incident_pressure_dot) /= size(incident_pressure)) return
        call build_system(vertices, triangles, interface_nodes, &
            interface_normals, angular_frequency, sound_speed, fluid_density, &
            period, maximum_mode, young_modulus, poisson_ratio, solid_density, &
            volume_load, incident_pressure, dirichlet_dofs, dirichlet_values, &
            system, status, absorbing_edges, absorbing_normals)
        if (status /= 0) return
        m = size(interface_nodes)

        call elastic_material_jvp(young_modulus, poisson_ratio, solid_density, &
            angular_frequency, young_modulus_dot, poisson_ratio_dot, &
            solid_density_dot, 0.0_dp, material_dot)
        allocate(values_dot(system%pattern%nnz))
        values_dot = (0.0_dp, 0.0_dp)
        call uniform_cells_assemble_jvp(system%pattern, system%cell_group, &
            vertices, triangles, system%material, vertices_dot, material_dot, &
            values_dot, status)
        if (status /= 0) return
        if (present(absorbing_edges)) then
            allocate(edge_dot(4, 4, size(absorbing_edges, 2)))
            !$omp parallel do default(none) &
            !$omp shared(absorbing_edges, absorbing_normals, vertices) &
            !$omp shared(vertices_dot, system, material_dot, angular_frequency) &
            !$omp shared(edge_dot) private(e, delta, delta_dot, length_dot) &
            !$omp schedule(static)
            do e = 1, size(absorbing_edges, 2)
                delta = vertices(:, absorbing_edges(2, e)) - &
                    vertices(:, absorbing_edges(1, e))
                delta_dot = vertices_dot(:, absorbing_edges(2, e)) - &
                    vertices_dot(:, absorbing_edges(1, e))
                call generated_segment_length_jvp(delta(1), delta(2), &
                    delta_dot(1), delta_dot(2), length_dot)
                ! d Re(Z) = Re(dZ): the boundary uses the real impedances.
                call generated_p1_elastic_impedance_edge_jvp(system%lengths(e), &
                    absorbing_normals(1, e), absorbing_normals(2, e), &
                    real_part(system%material%p_impedance), &
                    real_part(system%material%s_impedance), angular_frequency, &
                    length_dot, 0.0_dp, 0.0_dp, &
                    real_part(material_dot%p_impedance), &
                    real_part(material_dot%s_impedance), 0.0_dp, &
                    edge_dot(:, :, e))
            end do
            !$omp end parallel do
            call block_csc_add(system%pattern, system%edge_group, edge_dot, &
                values_dot)
        end if
        ! The form is linear in the fluid density.
        form_dot = system%form*(fluid_density_dot/fluid_density)
        allocate(block_dot(2*m, 2*m, 1))
        !$omp parallel do default(none) &
        !$omp shared(m, interface_normals, system, form_dot, block_dot) &
        !$omp private(a, b) schedule(static)
        do b = 1, m
            do a = 1, m
                call generated_interface_normal_block_jvp( &
                    interface_normals(1, a), interface_normals(2, a), &
                    interface_normals(1, b), interface_normals(2, b), &
                    system%form(a, b), form_dot(a, b), &
                    block_dot(2*a - 1:2*a, 2*b - 1:2*b, 1))
            end do
        end do
        !$omp end parallel do
        call block_csc_add(system%pattern, system%interface_group, block_dot, &
            values_dot)

        allocate(rhs_dot(size(system%rhs)))
        rhs_dot = (0.0_dp, 0.0_dp)
        do a = 1, m
            previous = modulo(a - 2, m) + 1
            next = modulo(a, m) + 1
            call generated_planar_incident_load_jvp(period, real(m, dp), &
                incident_pressure(previous), incident_pressure(a), &
                incident_pressure(next), interface_normals(1, a), &
                interface_normals(2, a), 0.0_dp, &
                incident_pressure_dot(previous), incident_pressure_dot(a), &
                incident_pressure_dot(next), &
                rhs_dot(2*interface_nodes(a) - 1:2*interface_nodes(a)))
        end do
        call sparse_direct_solve_constrained_jvp(system%matrix, system%rhs, &
            system%constrained, system%fixed, values_dot, rhs_dot, &
            0.0_dp*system%fixed, solution_dot, status)
        if (status /= 0) then
            solution_dot = (0.0_dp, 0.0_dp)
            status = 2
        end if
    end subroutine solve_elasticity_planar_acoustic_dtn_p1_jvp

    subroutine solve_elasticity_planar_acoustic_dtn_p1_vjp( &
            vertices, triangles, interface_nodes, interface_normals, &
            angular_frequency, sound_speed, fluid_density, period, &
            maximum_mode, young_modulus, poisson_ratio, solid_density, &
            volume_load, incident_pressure, dirichlet_dofs, dirichlet_values, &
            solution_bar, vertices_bar, fluid_density_bar, young_modulus_bar, &
            poisson_ratio_bar, solid_density_bar, incident_pressure_bar, status, &
            absorbing_edges, absorbing_normals)
        !! Cotangents of the inputs of the JVP for the solution cotangent
        !! solution_bar; Re<solution_bar, solution_dot> equals their pairing
        !! with the tangents. The outputs are overwritten.
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
        complex(dp), intent(in) :: dirichlet_values(:), solution_bar(:)
        real(dp), intent(out) :: vertices_bar(:, :), fluid_density_bar
        complex(dp), intent(out) :: young_modulus_bar
        real(dp), intent(out) :: poisson_ratio_bar, solid_density_bar
        complex(dp), intent(out) :: incident_pressure_bar(:)
        integer, intent(out) :: status
        integer, intent(in), optional :: absorbing_edges(:, :)
        real(dp), intent(in), optional :: absorbing_normals(:, :)

        type(planar_system_t) :: system
        type(elastic_material_t) :: material_bar
        complex(dp), allocatable :: solution(:), values_bar(:), rhs_bar(:)
        complex(dp), allocatable :: fixed_bar(:), edge_bar(:, :, :)
        complex(dp), allocatable :: block_bar(:, :, :)
        complex(dp) :: length_c, n1_bar, n2_bar, zp_bar, zs_bar, omega_bar
        complex(dp) :: value_bar, period_bar, pressure_bar(3)
        real(dp) :: delta(2), delta_bar(2), length_bar, density_sum
        real(dp) :: omega_material_bar
        integer :: a, b, e, m, previous, next

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
        call build_system(vertices, triangles, interface_nodes, &
            interface_normals, angular_frequency, sound_speed, fluid_density, &
            period, maximum_mode, young_modulus, poisson_ratio, solid_density, &
            volume_load, incident_pressure, dirichlet_dofs, dirichlet_values, &
            system, status, absorbing_edges, absorbing_normals)
        if (status /= 0) return
        m = size(interface_nodes)
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
        if (present(absorbing_edges)) then
            allocate(edge_bar(4, 4, size(absorbing_edges, 2)))
            call block_csc_gather(system%pattern, system%edge_group, values_bar, &
                edge_bar)
            ! Serial: edges share vertices, and the boundary is small.
            do e = 1, size(absorbing_edges, 2)
                call generated_p1_elastic_impedance_edge_vjp(system%lengths(e), &
                    absorbing_normals(1, e), absorbing_normals(2, e), &
                    real_part(system%material%p_impedance), &
                    real_part(system%material%s_impedance), angular_frequency, &
                    conjg(edge_bar(:, :, e)), length_c, n1_bar, n2_bar, zp_bar, &
                    zs_bar, omega_bar)
                ! Re(Z) is real, so its cotangent is the real part, and the
                ! adjoint of Re passes it on unchanged.
                material_bar%p_impedance = material_bar%p_impedance + &
                    real(zp_bar, dp)
                material_bar%s_impedance = material_bar%s_impedance + &
                    real(zs_bar, dp)
                delta = vertices(:, absorbing_edges(2, e)) - &
                    vertices(:, absorbing_edges(1, e))
                length_bar = real(length_c, dp)
                call generated_segment_length_vjp(delta(1), delta(2), length_bar, &
                    delta_bar(1), delta_bar(2))
                vertices_bar(:, absorbing_edges(2, e)) = &
                    vertices_bar(:, absorbing_edges(2, e)) + delta_bar
                vertices_bar(:, absorbing_edges(1, e)) = &
                    vertices_bar(:, absorbing_edges(1, e)) - delta_bar
            end do
        end if
        call elastic_material_vjp(young_modulus, poisson_ratio, solid_density, &
            angular_frequency, material_bar, young_modulus_bar, &
            poisson_ratio_bar, solid_density_bar, omega_material_bar)

        allocate(block_bar(2*m, 2*m, 1))
        call block_csc_gather(system%pattern, system%interface_group, values_bar, &
            block_bar)
        density_sum = 0.0_dp
        do b = 1, m
            do a = 1, m
                call generated_interface_normal_block_vjp( &
                    interface_normals(1, a), interface_normals(2, a), &
                    interface_normals(1, b), interface_normals(2, b), &
                    system%form(a, b), conjg(block_bar(2*a - 1:2*a, 2*b - 1:2*b, 1)), &
                    value_bar)
                ! value_bar is conj(form_bar(a, b)); d form/d rho = form/rho.
                density_sum = density_sum + real(value_bar*system%form(a, b), dp)
            end do
        end do
        fluid_density_bar = density_sum/fluid_density

        do a = 1, m
            previous = modulo(a - 2, m) + 1
            next = modulo(a, m) + 1
            call generated_planar_incident_load_vjp(period, real(m, dp), &
                incident_pressure(previous), incident_pressure(a), &
                incident_pressure(next), interface_normals(1, a), &
                interface_normals(2, a), &
                conjg(rhs_bar(2*interface_nodes(a) - 1:2*interface_nodes(a))), &
                period_bar, pressure_bar(1), pressure_bar(2), pressure_bar(3))
            incident_pressure_bar(previous) = incident_pressure_bar(previous) + &
                conjg(pressure_bar(1))
            incident_pressure_bar(a) = incident_pressure_bar(a) + &
                conjg(pressure_bar(2))
            incident_pressure_bar(next) = incident_pressure_bar(next) + &
                conjg(pressure_bar(3))
        end do
    end subroutine solve_elasticity_planar_acoustic_dtn_p1_vjp

    subroutine build_system(vertices, triangles, interface_nodes, &
            interface_normals, angular_frequency, sound_speed, fluid_density, &
            period, maximum_mode, young_modulus, poisson_ratio, solid_density, &
            volume_load, incident_pressure, dirichlet_dofs, dirichlet_values, &
            system, status, absorbing_edges, absorbing_normals)
        !! Pattern, matrix values and right-hand side of the coupled problem.
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
        type(planar_system_t), intent(out) :: system
        integer, intent(out) :: status
        integer, intent(in), optional :: absorbing_edges(:, :)
        real(dp), intent(in), optional :: absorbing_normals(:, :)

        complex(dp), allocatable :: edge_blocks(:, :, :), blocks(:, :, :)
        complex(dp) :: load(2)
        integer, allocatable :: cell_dofs(:, :), edge_dofs(:, :), interface_dofs(:, :)
        real(dp) :: delta(2)
        integer :: a, b, e, m, dof_count, next, previous
        logical :: bad_edge

        status = 1
        m = size(interface_nodes)
        dof_count = 2*size(vertices, 2)
        if (.not. valid_inputs()) return

        call elastic_material(young_modulus, poisson_ratio, solid_density, &
            angular_frequency, system%material)
        call elasticity_p1_cell_dofs(triangles, cell_dofs)
        call block_csc_init(system%pattern, dof_count)
        call block_csc_add_group(system%pattern, cell_dofs, system%cell_group)
        if (present(absorbing_edges)) then
            allocate(edge_dofs(4, size(absorbing_edges, 2)))
            edge_dofs(1::2, :) = 2*absorbing_edges - 1
            edge_dofs(2::2, :) = 2*absorbing_edges
            call block_csc_add_group(system%pattern, edge_dofs, system%edge_group)
        end if
        allocate(interface_dofs(2*m, 1))
        interface_dofs(1::2, 1) = 2*interface_nodes - 1
        interface_dofs(2::2, 1) = 2*interface_nodes
        call block_csc_add_group(system%pattern, interface_dofs, &
            system%interface_group)
        call block_csc_finalize(system%pattern)
        call block_csc_new_matrix(system%pattern, system%matrix)

        call uniform_cells_assemble(system%pattern, system%cell_group, vertices, &
            triangles, system%material, system%matrix%val, status)
        if (status /= 0) return
        status = 1
        if (present(absorbing_edges)) then
            allocate(system%lengths(size(absorbing_edges, 2)))
            allocate(edge_blocks(4, 4, size(absorbing_edges, 2)))
            bad_edge = .false.
            !$omp parallel do default(none) &
            !$omp shared(absorbing_edges, absorbing_normals, vertices, system) &
            !$omp shared(angular_frequency, edge_blocks) private(e, delta) &
            !$omp reduction(.or.: bad_edge) schedule(static)
            do e = 1, size(absorbing_edges, 2)
                delta = vertices(:, absorbing_edges(2, e)) - &
                    vertices(:, absorbing_edges(1, e))
                call generated_segment_length(delta(1), delta(2), &
                    system%lengths(e))
                if (system%lengths(e) <= tiny(1.0_dp)) bad_edge = .true.
                if (abs(norm2(absorbing_normals(:, e)) - 1.0_dp) > 1.0e-10_dp) &
                    bad_edge = .true.
                call generated_p1_elastic_impedance_edge(system%lengths(e), &
                    absorbing_normals(1, e), absorbing_normals(2, e), &
                    real_part(system%material%p_impedance), &
                    real_part(system%material%s_impedance), angular_frequency, &
                    edge_blocks(:, :, e))
            end do
            !$omp end parallel do
            if (bad_edge) return
            call block_csc_add(system%pattern, system%edge_group, edge_blocks, &
                system%matrix%val)
        end if

        allocate(system%form(m, m))
        call assemble_planar_acoustic_displacement_dtn_form(m, angular_frequency, &
            sound_speed, fluid_density, period, maximum_mode, system%form, status)
        if (status /= 0) return
        allocate(blocks(2*m, 2*m, 1))
        !$omp parallel do default(none) &
        !$omp shared(m, interface_normals, system, blocks) private(a, b) &
        !$omp schedule(static)
        do b = 1, m
            do a = 1, m
                call generated_interface_normal_block( &
                    interface_normals(1, a), interface_normals(2, a), &
                    interface_normals(1, b), interface_normals(2, b), &
                    system%form(a, b), blocks(2*a - 1:2*a, 2*b - 1:2*b, 1))
            end do
        end do
        !$omp end parallel do
        call block_csc_add(system%pattern, system%interface_group, blocks, &
            system%matrix%val)

        system%rhs = volume_load
        do a = 1, m
            previous = modulo(a - 2, m) + 1
            next = modulo(a, m) + 1
            associate (node => interface_nodes(a))
                call generated_planar_incident_load(period, real(m, dp), &
                    incident_pressure(previous), incident_pressure(a), &
                    incident_pressure(next), interface_normals(1, a), &
                    interface_normals(2, a), load)
                system%rhs(2*node - 1:2*node) = system%rhs(2*node - 1:2*node) + &
                    load
            end associate
        end do
        call dirichlet_mask(dof_count, dirichlet_dofs, dirichlet_values, &
            system%constrained, system%fixed)
        status = 0

    contains

        logical function valid_inputs()
            valid_inputs = .false.
            if (size(vertices, 1) /= 2 .or. size(triangles, 1) /= 3) return
            if (dof_count < 2 .or. size(triangles, 2) < 1) return
            if (m < 3) return
            if (size(interface_normals, 1) /= 2 .or. &
                size(interface_normals, 2) /= m) return
            if (size(incident_pressure) /= m) return
            if (size(volume_load) /= dof_count) return
            if (size(dirichlet_dofs) /= size(dirichlet_values)) return
            if (present(absorbing_edges) .neqv. &
                present(absorbing_normals)) return
            if (present(absorbing_edges)) then
                if (size(absorbing_edges, 1) /= 2 .or. &
                    size(absorbing_normals, 1) /= 2) return
                if (size(absorbing_edges, 2) /= &
                    size(absorbing_normals, 2)) return
                if (any(absorbing_edges < 1) .or. &
                    any(absorbing_edges > dof_count/2)) return
            end if
            if (angular_frequency <= 0.0_dp .or. abs(sound_speed) <= 0.0_dp .or. &
                fluid_density <= 0.0_dp .or. solid_density <= 0.0_dp) return
            if (period <= 0.0_dp .or. abs(young_modulus) <= 0.0_dp) return
            if (poisson_ratio <= -1.0_dp .or. &
                poisson_ratio >= 0.5_dp) return
            if (any(triangles < 1) .or. any(triangles > dof_count/2)) return
            if (any(interface_nodes < 1) .or. &
                any(interface_nodes > dof_count/2)) return
            if (any(dirichlet_dofs < 1) .or. &
                any(dirichlet_dofs > dof_count)) return
            valid_inputs = .true.
        end function valid_inputs

    end subroutine build_system

    elemental function real_part(z) result(r)
        !! The absorbing condition uses the real parts of the impedances.
        complex(dp), intent(in) :: z
        complex(dp) :: r

        r = cmplx(real(z, dp), 0.0_dp, dp)
    end function real_part

end module fortfem_elasticity_planar_acoustic_dtn_2d
