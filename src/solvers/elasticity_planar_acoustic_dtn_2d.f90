module fortfem_elasticity_planar_acoustic_dtn_2d
    !! Vector P1 time-harmonic elasticity coupled to a periodic acoustic half-space
    !! through the planar displacement-to-pressure (DtN) form, with JVP and VJP.
    !!
    !! Time dependence is exp(-i omega t). The DtN mode multiplier
    !! -i rho_f omega**2/k_n with the principal root k_n is outgoing, absorbing
    !! edges add -i omega Z M, and a lossy material has Im(c) < 0 and Im(E) < 0
    !! (Im(c) > 0 or Im(E) > 0 is gain). The interface samples are m equispaced
    !! points of one period of the fluid half-space; the normals point from the
    !! solid into the fluid.
    !!
    !! One CSC pattern from fortfem_block_csc holds three groups: the triangles
    !! (generated elasticity cell), the optional absorbing edges (generated
    !! impedance edge with a complex (Z_p, Z_s) per edge, by default the real
    !! parts of the solid's P and S impedances) and the interface coupling.
    !! The DtN couples only the normal displacement component: when every
    !! interface normal is axis-aligned the interface group is the m x m block
    !! of the acoustic form on the normal-component dofs, contracted with the
    !! normals at assembly; otherwise it is the dense block over the 2m
    !! interface dofs. Dirichlet values go to sparse_direct_solve_constrained.
    !! Storage is O(triangles + nnz + (k*m)**2) with k the number of active
    !! normal components and a scratch of m**2, not (2m)**2, for the block.
    !!
    !! Derivatives cover vertex coordinates, angular frequency, complex sound
    !! speed, fluid density, Young's modulus, Poisson ratio, solid density,
    !! incident pressure and the per-edge impedances; the DtN form derivatives
    !! come from the generated mode multiplier (gen_planar_acoustic_dtn_products).
    !!
    !! An optional dtn_period_factor lfac > 1 extends the DtN period to lfac
    !! times the interface with a hard baffle (zero normal displacement) beyond
    !! it, removing the periodic coupling of the interface ends.
    !! Cotangents follow fortfem_elasticity_acoustic_common_2d.
    use, intrinsic :: iso_fortran_env, only: int64
    use fortfem_block_csc, only: block_csc_add, block_csc_add_group, &
        block_csc_finalize, block_csc_gather, block_csc_init, &
        block_csc_new_matrix, block_csc_pattern_t
    use fortfem_dirichlet_constraints, only: dirichlet_mask
    use fortfem_elasticity_acoustic_common_2d, only: elastic_material, &
        elastic_material_jvp, elastic_material_t, elastic_material_vjp, &
        uniform_cells_assemble, uniform_cells_assemble_jvp, &
        uniform_cells_assemble_vjp
    use fortfem_elasticity_p1_csc_2d, only: elasticity_p1_cell_dofs
    use fortfem_generated_interface_normal_block, only: &
        generated_interface_normal_block
    use fortfem_generated_interface_normal_block_jvp, only: &
        generated_interface_normal_block_jvp
    use fortfem_generated_interface_normal_block_vjp, only: &
        generated_interface_normal_block_vjp
    use fortfem_generated_interface_normal_component_block, only: &
        generated_interface_normal_component_block
    use fortfem_generated_interface_normal_component_block_jvp, only: &
        generated_interface_normal_component_block_jvp
    use fortfem_generated_interface_normal_component_block_vjp, only: &
        generated_interface_normal_component_block_vjp
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
        assemble_planar_acoustic_displacement_dtn_form, &
        assemble_planar_acoustic_displacement_dtn_form_jvp, &
        assemble_planar_acoustic_displacement_dtn_form_vjp
    use fortfem_sparse_direct, only: sparse_direct_solve_constrained, &
        sparse_direct_solve_constrained_jvp, sparse_direct_solve_constrained_vjp
    use fortsparse, only: csc_z_t
    implicit none
    private

    public :: solve_elasticity_planar_acoustic_dtn_p1
    public :: solve_elasticity_planar_acoustic_dtn_p1_jvp
    public :: solve_elasticity_planar_acoustic_dtn_p1_vjp
    public :: planar_dtn_interface_entries

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
        complex(dp), allocatable :: p_impedance(:), s_impedance(:)
        complex(dp), allocatable :: form(:, :)
        integer :: dtn_samples = 0, dtn_modes = 0
        real(dp) :: dtn_period = 0.0_dp
        logical :: baffled = .false.
        !! Interface coupling layout: normal_only with one active component
        !! stores the m x m form on the normal-component dofs; the dense
        !! layout keeps the 2m x 2m block over both components.
        logical :: normal_only = .false.
        integer :: normal_component = 0
    end type planar_system_t

contains

    subroutine solve_elasticity_planar_acoustic_dtn_p1_real( &
            vertices, triangles, interface_nodes, interface_normals, &
            angular_frequency, sound_speed, fluid_density, period, &
            maximum_mode, young_modulus, poisson_ratio, solid_density, &
            volume_load, incident_pressure, dirichlet_dofs, &
            dirichlet_values, solution, status, absorbing_edges, &
            absorbing_normals, absorbing_impedances, dtn_period_factor)
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
        complex(dp), intent(in), optional :: absorbing_impedances(:, :)
        real(dp), intent(in), optional :: dtn_period_factor

        call solve_elasticity_planar_acoustic_dtn_p1_complex( &
            vertices, triangles, interface_nodes, interface_normals, &
            angular_frequency, cmplx(sound_speed, 0.0_dp, dp), &
            fluid_density, period, maximum_mode, &
            cmplx(young_modulus, 0.0_dp, dp), poisson_ratio, &
            solid_density, volume_load, incident_pressure, &
            dirichlet_dofs, dirichlet_values, solution, status, &
            absorbing_edges, absorbing_normals, absorbing_impedances, &
            dtn_period_factor)
    end subroutine solve_elasticity_planar_acoustic_dtn_p1_real

    subroutine solve_elasticity_planar_acoustic_dtn_p1_complex( &
            vertices, triangles, interface_nodes, interface_normals, &
            angular_frequency, sound_speed, fluid_density, period, &
            maximum_mode, young_modulus, poisson_ratio, solid_density, &
            volume_load, incident_pressure, dirichlet_dofs, &
            dirichlet_values, solution, status, absorbing_edges, &
            absorbing_normals, absorbing_impedances, dtn_period_factor)
        !! Time dependence exp(-i omega t): Im(c) < 0 and Im(E) < 0 are lossy.
        !! status: 0 success, 1 invalid input or failed form, 2 failed solve.
        !! absorbing_impedances(:, e) = (Z_p, Z_s) of absorbing edge e, for
        !! example (rho_f c_f, 0) for a fluid half-space behind the face; without
        !! it every edge gets the real parts of the solid's own impedances.
        !! dtn_period_factor = lfac >= 1 (default 1, the periodic interface)
        !! sets the DtN period to lfac*period with lfac*m samples, the interface
        !! in the first m and a hard baffle (zero normal displacement) in the
        !! rest; lfac*m must be an integer. The mode cutoff scales with it, so
        !! the largest retained wavenumber stays 2 pi maximum_mode/period. The
        !! incident load extends the pressure constantly by one sample beyond
        !! each end instead of wrapping it.
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
        complex(dp), intent(in), optional :: absorbing_impedances(:, :)
        real(dp), intent(in), optional :: dtn_period_factor

        type(planar_system_t) :: system

        solution = (0.0_dp, 0.0_dp)
        status = 1
        if (size(solution) /= 2*size(vertices, 2)) return
        call build_system(vertices, triangles, interface_nodes, &
            interface_normals, angular_frequency, sound_speed, fluid_density, &
            period, maximum_mode, young_modulus, poisson_ratio, solid_density, &
            volume_load, incident_pressure, dirichlet_dofs, dirichlet_values, &
            system, status, absorbing_edges, absorbing_normals, &
            absorbing_impedances, dtn_period_factor)
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
            vertices_dot, angular_frequency_dot, sound_speed_dot, &
            fluid_density_dot, young_modulus_dot, poisson_ratio_dot, &
            solid_density_dot, incident_pressure_dot, solution_dot, status, &
            absorbing_edges, absorbing_normals, absorbing_impedances, &
            absorbing_impedances_dot, dtn_period_factor)
        !! Directional derivative of the solution of the complex interface.
        !! status as for the solve. absorbing_impedances_dot is the tangent of
        !! absorbing_impedances (zero when absent).
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
        real(dp), intent(in) :: vertices_dot(:, :), angular_frequency_dot
        complex(dp), intent(in) :: sound_speed_dot
        real(dp), intent(in) :: fluid_density_dot
        complex(dp), intent(in) :: young_modulus_dot
        real(dp), intent(in) :: poisson_ratio_dot, solid_density_dot
        complex(dp), intent(in) :: incident_pressure_dot(:)
        complex(dp), intent(out) :: solution_dot(:)
        integer, intent(out) :: status
        integer, intent(in), optional :: absorbing_edges(:, :)
        real(dp), intent(in), optional :: absorbing_normals(:, :)
        complex(dp), intent(in), optional :: absorbing_impedances(:, :)
        real(dp), intent(in), optional :: dtn_period_factor
        complex(dp), intent(in), optional :: absorbing_impedances_dot(:, :)

        type(planar_system_t) :: system
        type(elastic_material_t) :: material_dot
        complex(dp), allocatable :: values_dot(:), rhs_dot(:), edge_dot(:, :, :)
        complex(dp), allocatable :: block_dot(:, :, :), form_dot(:, :)
        complex(dp), allocatable :: zp_dot(:), zs_dot(:)
        integer :: form_status
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
            system, status, absorbing_edges, absorbing_normals, &
            absorbing_impedances, dtn_period_factor)
        if (status /= 0) return
        m = size(interface_nodes)

        call elastic_material_jvp(young_modulus, poisson_ratio, solid_density, &
            angular_frequency, young_modulus_dot, poisson_ratio_dot, &
            solid_density_dot, angular_frequency_dot, material_dot)
        allocate(values_dot(system%pattern%nnz))
        values_dot = (0.0_dp, 0.0_dp)
        call uniform_cells_assemble_jvp(system%pattern, system%cell_group, &
            vertices, triangles, system%material, vertices_dot, material_dot, &
            values_dot, status)
        if (status /= 0) return
        if (present(absorbing_edges)) then
            allocate(edge_dot(4, 4, size(absorbing_edges, 2)))
            allocate(zp_dot(size(absorbing_edges, 2)))
            allocate(zs_dot(size(absorbing_edges, 2)))
            if (present(absorbing_impedances)) then
                zp_dot = (0.0_dp, 0.0_dp)
                zs_dot = (0.0_dp, 0.0_dp)
                if (present(absorbing_impedances_dot)) then
                    if (any(shape(absorbing_impedances_dot) /= &
                        shape(absorbing_impedances))) then
                        status = 1
                        return
                    end if
                    zp_dot = absorbing_impedances_dot(1, :)
                    zs_dot = absorbing_impedances_dot(2, :)
                end if
            else
                ! d Re(Z) = Re(dZ): the default boundary uses the real impedances.
                zp_dot = real_part(material_dot%p_impedance)
                zs_dot = real_part(material_dot%s_impedance)
            end if
            !$omp parallel do default(none) &
            !$omp shared(absorbing_edges, absorbing_normals, vertices) &
            !$omp shared(vertices_dot, system, angular_frequency) &
            !$omp shared(angular_frequency_dot, zp_dot, zs_dot) &
            !$omp shared(edge_dot) private(e, delta, delta_dot, length_dot) &
            !$omp schedule(static)
            do e = 1, size(absorbing_edges, 2)
                delta = vertices(:, absorbing_edges(2, e)) - &
                    vertices(:, absorbing_edges(1, e))
                delta_dot = vertices_dot(:, absorbing_edges(2, e)) - &
                    vertices_dot(:, absorbing_edges(1, e))
                call generated_segment_length_jvp(delta(1), delta(2), &
                    delta_dot(1), delta_dot(2), length_dot)
                call generated_p1_elastic_impedance_edge_jvp(system%lengths(e), &
                    absorbing_normals(1, e), absorbing_normals(2, e), &
                    system%p_impedance(e), system%s_impedance(e), &
                    angular_frequency, length_dot, 0.0_dp, 0.0_dp, zp_dot(e), &
                    zs_dot(e), angular_frequency_dot, edge_dot(:, :, e))
            end do
            !$omp end parallel do
            call block_csc_add(system%pattern, system%edge_group, edge_dot, &
                values_dot)
        end if
        allocate(form_dot(m, m))
        call assemble_planar_acoustic_displacement_dtn_form_jvp(m, &
            angular_frequency, sound_speed, fluid_density, system%dtn_period, &
            system%dtn_modes, angular_frequency_dot, sound_speed_dot, &
            fluid_density_dot, form_dot, form_status, system%dtn_samples)
        if (form_status /= 0) then
            status = 1
            return
        end if
        if (system%normal_only) then
            allocate(block_dot(m, m, 1))
            !$omp parallel do default(none) &
            !$omp shared(m, interface_normals, system, form_dot, block_dot) &
            !$omp private(a, b) schedule(static)
            do b = 1, m
                do a = 1, m
                    call generated_interface_normal_component_block_jvp( &
                        interface_normals(system%normal_component, a), &
                        interface_normals(system%normal_component, b), &
                        system%form(a, b), form_dot(a, b), block_dot(a, b, 1))
                end do
            end do
            !$omp end parallel do
        else
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
        end if
        call block_csc_add(system%pattern, system%interface_group, block_dot, &
            values_dot)

        allocate(rhs_dot(size(system%rhs)))
        rhs_dot = (0.0_dp, 0.0_dp)
        do a = 1, m
            call load_neighbours(a, m, system%baffled, previous, next)
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
            solution_bar, vertices_bar, angular_frequency_bar, sound_speed_bar, &
            fluid_density_bar, young_modulus_bar, poisson_ratio_bar, &
            solid_density_bar, incident_pressure_bar, status, absorbing_edges, &
            absorbing_normals, absorbing_impedances, absorbing_impedances_bar, &
            dtn_period_factor)
        !! Cotangents of the inputs of the JVP for the solution cotangent
        !! solution_bar; Re<solution_bar, solution_dot> equals their pairing
        !! with the tangents. The outputs are overwritten; absorbing_impedances_bar
        !! pairs with absorbing_impedances_dot as Re(conj(z_bar) z_dot).
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
        real(dp), intent(out) :: vertices_bar(:, :), angular_frequency_bar
        complex(dp), intent(out) :: sound_speed_bar
        real(dp), intent(out) :: fluid_density_bar
        complex(dp), intent(out) :: young_modulus_bar
        real(dp), intent(out) :: poisson_ratio_bar, solid_density_bar
        complex(dp), intent(out) :: incident_pressure_bar(:)
        integer, intent(out) :: status
        integer, intent(in), optional :: absorbing_edges(:, :)
        real(dp), intent(in), optional :: absorbing_normals(:, :)
        complex(dp), intent(in), optional :: absorbing_impedances(:, :)
        real(dp), intent(in), optional :: dtn_period_factor
        complex(dp), intent(out), optional :: absorbing_impedances_bar(:, :)

        type(planar_system_t) :: system
        type(elastic_material_t) :: material_bar
        complex(dp), allocatable :: solution(:), values_bar(:), rhs_bar(:)
        complex(dp), allocatable :: fixed_bar(:), edge_bar(:, :, :)
        complex(dp), allocatable :: block_bar(:, :, :), form_bar(:, :)
        complex(dp) :: length_c, n1_bar, n2_bar, zp_bar, zs_bar, omega_bar
        complex(dp) :: value_bar, period_bar, pressure_bar(3)
        real(dp) :: delta(2), delta_bar(2), length_bar
        real(dp) :: omega_material_bar, omega_form_bar
        integer :: a, b, e, m, previous, next

        vertices_bar = 0.0_dp
        angular_frequency_bar = 0.0_dp
        sound_speed_bar = (0.0_dp, 0.0_dp)
        fluid_density_bar = 0.0_dp
        young_modulus_bar = (0.0_dp, 0.0_dp)
        poisson_ratio_bar = 0.0_dp
        solid_density_bar = 0.0_dp
        incident_pressure_bar = (0.0_dp, 0.0_dp)
        if (present(absorbing_impedances_bar)) &
            absorbing_impedances_bar = (0.0_dp, 0.0_dp)
        status = 1
        if (size(solution_bar) /= 2*size(vertices, 2)) return
        if (present(absorbing_impedances_bar)) then
            if (.not. present(absorbing_impedances)) return
            if (any(shape(absorbing_impedances_bar) /= &
                shape(absorbing_impedances))) return
        end if
        if (any(shape(vertices_bar) /= shape(vertices))) return
        if (size(incident_pressure_bar) /= size(incident_pressure)) return
        call build_system(vertices, triangles, interface_nodes, &
            interface_normals, angular_frequency, sound_speed, fluid_density, &
            period, maximum_mode, young_modulus, poisson_ratio, solid_density, &
            volume_load, incident_pressure, dirichlet_dofs, dirichlet_values, &
            system, status, absorbing_edges, absorbing_normals, &
            absorbing_impedances, dtn_period_factor)
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
                    system%p_impedance(e), system%s_impedance(e), &
                    angular_frequency, conjg(edge_bar(:, :, e)), length_c, &
                    n1_bar, n2_bar, zp_bar, zs_bar, omega_bar)
                if (present(absorbing_impedances)) then
                    if (present(absorbing_impedances_bar)) then
                        absorbing_impedances_bar(1, e) = conjg(zp_bar)
                        absorbing_impedances_bar(2, e) = conjg(zs_bar)
                    end if
                else
                    ! Re(Z) is real, so its cotangent is the real part, and the
                    ! adjoint of Re passes it on unchanged.
                    material_bar%p_impedance = material_bar%p_impedance + &
                        real(zp_bar, dp)
                    material_bar%s_impedance = material_bar%s_impedance + &
                        real(zs_bar, dp)
                end if
                angular_frequency_bar = angular_frequency_bar + real(omega_bar, dp)
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
        angular_frequency_bar = angular_frequency_bar + omega_material_bar

        if (system%normal_only) then
            allocate(block_bar(m, m, 1))
            call block_csc_gather(system%pattern, system%interface_group, &
                values_bar, block_bar)
            allocate(form_bar(m, m))
            do b = 1, m
                do a = 1, m
                    call generated_interface_normal_component_block_vjp( &
                        interface_normals(system%normal_component, a), &
                        interface_normals(system%normal_component, b), &
                        system%form(a, b), conjg(block_bar(a, b, 1)), value_bar)
                    form_bar(a, b) = conjg(value_bar)
                end do
            end do
        else
            call block_csc_gather(system%pattern, system%interface_group, &
                values_bar, block_bar)
            allocate(form_bar(m, m))
            do b = 1, m
                do a = 1, m
                    call generated_interface_normal_block_vjp( &
                        interface_normals(1, a), interface_normals(2, a), &
                        interface_normals(1, b), interface_normals(2, b), &
                        system%form(a, b), &
                        conjg(block_bar(2*a - 1:2*a, 2*b - 1:2*b, 1)), &
                        value_bar)
                    form_bar(a, b) = conjg(value_bar)
                end do
            end do
        end if
        call assemble_planar_acoustic_displacement_dtn_form_vjp(m, &
            angular_frequency, sound_speed, fluid_density, system%dtn_period, &
            system%dtn_modes, form_bar, omega_form_bar, sound_speed_bar, &
            fluid_density_bar, status, system%dtn_samples)
        if (status /= 0) then
            status = 1
            return
        end if
        angular_frequency_bar = angular_frequency_bar + omega_form_bar

        do a = 1, m
            call load_neighbours(a, m, system%baffled, previous, next)
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
            system, status, absorbing_edges, absorbing_normals, &
            absorbing_impedances, dtn_period_factor)
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
        complex(dp), intent(in), optional :: absorbing_impedances(:, :)
        real(dp), intent(in), optional :: dtn_period_factor

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
        if (.not. dtn_layout()) return
        call interface_layout(interface_normals, m, system%normal_only, &
            system%normal_component)

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
        if (system%normal_only) then
            allocate(interface_dofs(m, 1))
            if (system%normal_component == 1) then
                interface_dofs(:, 1) = 2*interface_nodes - 1
            else
                interface_dofs(:, 1) = 2*interface_nodes
            end if
        else
            allocate(interface_dofs(2*m, 1))
            interface_dofs(1::2, 1) = 2*interface_nodes - 1
            interface_dofs(2::2, 1) = 2*interface_nodes
        end if
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
            if (present(absorbing_impedances)) then
                system%p_impedance = absorbing_impedances(1, :)
                system%s_impedance = absorbing_impedances(2, :)
            else
                allocate(system%p_impedance(size(absorbing_edges, 2)), &
                    system%s_impedance(size(absorbing_edges, 2)))
                system%p_impedance = real_part(system%material%p_impedance)
                system%s_impedance = real_part(system%material%s_impedance)
            end if
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
                    system%p_impedance(e), system%s_impedance(e), &
                    angular_frequency, edge_blocks(:, :, e))
            end do
            !$omp end parallel do
            if (bad_edge) return
            call block_csc_add(system%pattern, system%edge_group, edge_blocks, &
                system%matrix%val)
        end if

        allocate(system%form(m, m))
        call assemble_planar_acoustic_displacement_dtn_form(m, angular_frequency, &
            sound_speed, fluid_density, system%dtn_period, system%dtn_modes, &
            system%form, status, system%dtn_samples)
        if (status /= 0) return
        if (system%normal_only) then
            !! The coupling lives on one scalar component: an m x m block of
            !! the form contracted with the normal components of its nodes.
            allocate(blocks(m, m, 1))
            !$omp parallel do default(none) &
            !$omp shared(m, interface_normals, system, blocks) private(a, b) &
            !$omp schedule(static)
            do b = 1, m
                do a = 1, m
                    call generated_interface_normal_component_block( &
                        interface_normals(system%normal_component, a), &
                        interface_normals(system%normal_component, b), &
                        system%form(a, b), blocks(a, b, 1))
                end do
            end do
            !$omp end parallel do
        else
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
        end if
        call block_csc_add(system%pattern, system%interface_group, blocks, &
            system%matrix%val)

        system%rhs = volume_load
        do a = 1, m
            call load_neighbours(a, m, system%baffled, previous, next)
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
            if (present(absorbing_impedances)) then
                if (.not. present(absorbing_edges)) return
                if (size(absorbing_impedances, 1) /= 2 .or. &
                    size(absorbing_impedances, 2) /= &
                    size(absorbing_edges, 2)) return
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

        logical function dtn_layout()
            !! DtN samples, period and mode cutoff of the (extended) period.
            real(dp) :: samples

            dtn_layout = .false.
            system%dtn_samples = m
            system%dtn_period = period
            system%dtn_modes = maximum_mode
            if (.not. present(dtn_period_factor)) then
                dtn_layout = .true.
                return
            end if
            if (dtn_period_factor < 1.0_dp) return
            samples = dtn_period_factor*real(m, dp)
            if (samples > real(huge(m)/2, dp)) return
            if (abs(samples - anint(samples)) > 1.0e-9_dp*samples) return
            system%dtn_samples = nint(samples)
            system%dtn_period = period*real(system%dtn_samples, dp)/real(m, dp)
            ! maximum_mode*M/m without overflow: same largest wavenumber.
            system%dtn_modes = int(real(maximum_mode, dp)*real(system%dtn_samples, &
                dp)/real(m, dp) + 1.0e-9_dp)
            system%baffled = system%dtn_samples > m
            dtn_layout = .true.
        end function dtn_layout

    end subroutine build_system

    pure subroutine interface_layout(normals, m, normal_only, component)
        !! The DtN couples only the normal component. When every interface
        !! normal is axis-aligned the coupling lives on that single scalar
        !! component and the interface block is m x m; otherwise both
        !! components are active and it is 2m x 2m.
        real(dp), intent(in) :: normals(:, :)
        integer, intent(in) :: m
        logical, intent(out) :: normal_only
        integer, intent(out) :: component

        logical :: active_x, active_y

        active_x = any(normals(1, 1:m) /= 0.0_dp)
        active_y = any(normals(2, 1:m) /= 0.0_dp)
        normal_only = active_x .neqv. active_y
        component = merge(1, 2, active_x)
    end subroutine interface_layout

    integer(int64) pure function planar_dtn_interface_entries(m, &
            interface_normals) result(entries)
        !! Number of interface-block entries the solver stores for m interface
        !! samples: m*m for axis-aligned normals (one active component),
        !! 4*m*m otherwise.
        integer, intent(in) :: m
        real(dp), intent(in) :: interface_normals(:, :)

        logical :: normal_only
        integer :: component
        integer(int64) :: width

        call interface_layout(interface_normals, m, normal_only, component)
        width = int(m, int64)
        if (.not. normal_only) width = 2*width
        entries = width*width
    end function planar_dtn_interface_entries

    pure subroutine load_neighbours(a, m, baffled, previous, next)
        !! Neighbour samples of the incident load: periodic, or, at the ends of
        !! a baffled interface, the end sample itself (constant extension).
        integer, intent(in) :: a, m
        logical, intent(in) :: baffled
        integer, intent(out) :: previous, next

        previous = modulo(a - 2, m) + 1
        next = modulo(a, m) + 1
        if (baffled) then
            previous = max(a - 1, 1)
            next = min(a + 1, m)
        end if
    end subroutine load_neighbours

    elemental function real_part(z) result(r)
        !! The absorbing condition uses the real parts of the impedances.
        complex(dp), intent(in) :: z
        complex(dp) :: r

        r = cmplx(real(z, dp), 0.0_dp, dp)
    end function real_part

end module fortfem_elasticity_planar_acoustic_dtn_2d
