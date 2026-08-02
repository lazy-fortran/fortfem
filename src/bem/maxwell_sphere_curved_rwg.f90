module fortfem_maxwell_sphere_curved_rwg
    !! Surface-Piola image of the affine RWG basis on a radial sphere panel.
    use fortfem_barycentric_surface_refinement, only: &
        barycentric_refine_sphere_surface_mesh_jvp, &
        barycentric_refine_sphere_surface_mesh_vjp
    use fortfem_kinds, only: dp
    use fortfem_maxwell_bc_surface, only: &
        build_maxwell_bc_transformation, &
        differentiate_maxwell_bc_transformation_jvp, &
        differentiate_maxwell_bc_transformation_vjp
    use fortfem_maxwell_efie_bc_3d, only: build_maxwell_bc_to_refined_rwg
    use fortfem_maxwell_rwg_surface, only: build_maxwell_rwg_surface_space
    use fortfem_sphere_curved_panel, only: &
        evaluate_sphere_curved_panel, evaluate_sphere_curved_panel_jvp, &
        evaluate_sphere_curved_panel_vjp, invert_sphere_curved_panel
    use fortfem_triangle_duffy_quadrature, only: triangle_duffy_quadrature
    use fortnum_linalg, only: dense_solve, inv3
    use fortnum_quadrature, only: gauss_legendre_ab
    implicit none
    private

    public :: evaluate_maxwell_sphere_curved_rwg_basis
    public :: assemble_maxwell_sphere_curved_rwg_mass_matrix
    public :: assemble_maxwell_sphere_curved_plane_wave_rhs_rwg_3d
    public :: assemble_maxwell_sphere_curved_plane_wave_rhs_bc_3d
    public :: assemble_maxwell_sphere_curved_plane_wave_rhs_bc_3d_jvp
    public :: assemble_maxwell_sphere_curved_plane_wave_rhs_bc_3d_vjp
    public :: evaluate_maxwell_sphere_curved_far_field_rwg_3d
    public :: integrate_maxwell_sphere_curved_coincident_rwg_pair_3d
    public :: integrate_maxwell_sphere_curved_adjacent_rwg_pair_3d
    public :: assemble_maxwell_sphere_curved_vector_potential_rwg_3d
    public :: assemble_maxwell_sphere_curved_potential_operators_rwg_3d
    public :: assemble_maxwell_sphere_curved_efie_rwg_3d
    public :: assemble_maxwell_sphere_curved_efie_imaginary_rwg_3d
    public :: assemble_maxwell_sphere_curved_efie_bc_imaginary_3d
    public :: assemble_maxwell_sphere_curved_regularized_cfie_rwg_3d
    public :: assemble_maxwell_sphere_curved_regularized_cfie_rhs_rwg_3d
    public :: solve_maxwell_pec_sphere_curved_regularized_cfie_rwg_3d
    public :: solve_maxwell_pec_sphere_curved_efie_rwg_3d
    public :: evaluate_maxwell_sphere_curved_magnetic_field_rwg_3d
    public :: evaluate_maxwell_sphere_curved_localized_rwg_basis
    public :: evaluate_maxwell_sphere_curved_localized_rwg_basis_jvp
    public :: evaluate_maxwell_sphere_curved_localized_rwg_basis_vjp
    public :: assemble_maxwell_sphere_curved_rwg_rbc_pairing
    public :: assemble_maxwell_sphere_curved_rwg_rbc_pairing_jvp
    public :: assemble_maxwell_sphere_curved_rwg_rbc_pairing_vjp
    public :: evaluate_maxwell_sphere_curved_rwg_basis_jvp
    public :: evaluate_maxwell_sphere_curved_rwg_basis_vjp
    public :: assemble_maxwell_sphere_curved_mfie_exterior_trace_rwg_rbc_3d
    public :: assemble_maxwell_sphere_curved_mfie_rwg_rbc_3d

contains

    subroutine solve_maxwell_pec_sphere_curved_regularized_cfie_rwg_3d( &
            vertices, triangles, radius, direction, polarization, wave_number, &
            impedance, quadrature_degree, tolerance, max_depth, mfie_offset, &
            density, status)
        real(dp), intent(in) :: vertices(:, :), radius, direction(3)
        complex(dp), intent(in) :: polarization(3)
        real(dp), intent(in) :: wave_number, impedance, tolerance, mfie_offset
        integer, intent(in) :: triangles(:, :), quadrature_degree, max_depth
        complex(dp), allocatable, intent(out) :: density(:)
        integer, intent(out) :: status

        complex(dp), allocatable :: cfie(:, :), efie(:, :), mfie(:, :)
        complex(dp), allocatable :: product(:, :), regularizer(:, :)
        complex(dp), allocatable :: right_hand_side(:)
        integer :: info

        call assemble_maxwell_sphere_curved_regularized_cfie_rwg_3d( &
            vertices, triangles, radius, wave_number, impedance, &
            quadrature_degree, tolerance, max_depth, mfie_offset, cfie, efie, &
            mfie, regularizer, product, status)
        if (status /= 0) return
        call assemble_maxwell_sphere_curved_regularized_cfie_rhs_rwg_3d( &
            vertices, triangles, radius, direction, polarization, wave_number, &
            impedance, quadrature_degree, regularizer, right_hand_side, status)
        if (status /= 0) return
        allocate(density(size(right_hand_side)))
        call dense_solve(cfie, right_hand_side, density, info)
        if (info /= 0) then
            status = 2
            return
        end if
        status = 0
    end subroutine &
        solve_maxwell_pec_sphere_curved_regularized_cfie_rwg_3d

    subroutine assemble_maxwell_sphere_curved_regularized_cfie_rhs_rwg_3d( &
            vertices, triangles, radius, direction, polarization, wave_number, &
            impedance, quadrature_degree, regularizer, right_hand_side, status)
        real(dp), intent(in) :: vertices(:, :), radius, direction(3)
        complex(dp), intent(in) :: polarization(3)
        real(dp), intent(in) :: wave_number, impedance
        integer, intent(in) :: triangles(:, :), quadrature_degree
        complex(dp), intent(in) :: regularizer(:, :)
        complex(dp), allocatable, intent(out) :: right_hand_side(:)
        integer, intent(out) :: status

        complex(dp), allocatable :: bc_rhs(:), efie_rhs(:), mass(:, :)
        complex(dp), allocatable :: mapped_rhs(:)
        real(dp), allocatable :: real_mass(:, :)
        integer :: info, system_size

        status = 1
        if (impedance <= 0.0_dp) return
        call assemble_maxwell_sphere_curved_plane_wave_rhs_rwg_3d( &
            vertices, triangles, radius, direction, polarization, wave_number, &
            quadrature_degree, efie_rhs, status)
        if (status /= 0) return
        call assemble_maxwell_sphere_curved_plane_wave_rhs_bc_3d( &
            vertices, triangles, radius, direction, polarization, wave_number, &
            quadrature_degree, bc_rhs, status)
        if (status /= 0) return
        call assemble_maxwell_sphere_curved_rwg_rbc_pairing( &
            vertices, triangles, radius, quadrature_degree, real_mass, status)
        if (status /= 0) return
        system_size = size(real_mass, 1)
        if (any(shape(regularizer) /= [system_size, system_size])) return
        allocate( &
            mass(system_size, system_size), mapped_rhs(system_size), &
            right_hand_side(system_size))
        mass = transpose(cmplx(real_mass, 0.0_dp, dp))
        call dense_solve(mass, efie_rhs, mapped_rhs, info)
        if (info /= 0) then
            status = 2
            return
        end if
        right_hand_side = bc_rhs - matmul(regularizer, mapped_rhs)
        status = 0
    end subroutine &
        assemble_maxwell_sphere_curved_regularized_cfie_rhs_rwg_3d

    subroutine assemble_maxwell_sphere_curved_plane_wave_rhs_bc_3d( &
            vertices, triangles, radius, direction, polarization, wave_number, &
            quadrature_degree, right_hand_side, status)
        real(dp), intent(in) :: vertices(:, :), radius, direction(3)
        complex(dp), intent(in) :: polarization(3)
        real(dp), intent(in) :: wave_number
        integer, intent(in) :: triangles(:, :), quadrature_degree
        complex(dp), allocatable, intent(out) :: right_hand_side(:)
        integer, intent(out) :: status

        integer, allocatable :: refined_triangles(:, :)
        real(dp), allocatable :: eta(:), refined_vertices(:, :)
        real(dp), allocatable :: transformation(:, :), weights(:), xi(:)
        real(dp) :: divergence, jacobian, local_value(3), point(3)
        complex(dp) :: incident_field(3)
        integer :: basis, local_edge, node, panel, row

        status = 1
        if (radius <= 0.0_dp .or. wave_number < 0.0_dp) return
        if (abs(norm2(direction) - 1.0_dp) > 128.0_dp*epsilon(1.0_dp)) return
        if (abs(sum(polarization*direction)) > &
            128.0_dp*epsilon(1.0_dp)*max(1.0_dp, &
            sqrt(sum(abs(polarization)**2)))) return
        call build_maxwell_bc_transformation( &
            vertices, triangles, refined_vertices, refined_triangles, &
            transformation, status, sphere_radius=radius)
        if (status /= 0) return
        call triangle_duffy_quadrature( &
            quadrature_degree, xi, eta, weights, status)
        if (status /= 0) return
        allocate(right_hand_side(size(transformation, 2)))
        right_hand_side = cmplx(0.0_dp, 0.0_dp, dp)
        do panel = 1, size(refined_triangles, 2)
            do node = 1, size(weights)
                do local_edge = 1, 3
                    call evaluate_maxwell_sphere_curved_localized_rwg_basis( &
                        refined_vertices, refined_triangles, panel, local_edge, &
                        radius, xi(node), eta(node), point, local_value, &
                        divergence, jacobian, status)
                    if (status /= 0) return
                    incident_field = polarization*exp(cmplx( &
                        0.0_dp, wave_number*dot_product(direction, point), dp))
                    row = 3*(panel - 1) + local_edge
                    do basis = 1, size(transformation, 2)
                        right_hand_side(basis) = right_hand_side(basis) - &
                            jacobian*weights(node)*transformation(row, basis)* &
                            sum(cmplx(local_value, 0.0_dp, dp)*incident_field)
                    end do
                end do
            end do
        end do
        status = 0
    end subroutine assemble_maxwell_sphere_curved_plane_wave_rhs_bc_3d

    subroutine assemble_maxwell_sphere_curved_plane_wave_rhs_bc_3d_jvp( &
            vertices, triangles, radius, direction, polarization, wave_number, &
            quadrature_degree, vertices_dot, radius_dot, direction_dot, &
            polarization_dot, wave_number_dot, right_hand_side, &
            right_hand_side_dot, status)
        real(dp), intent(in) :: vertices(:, :), radius, direction(3)
        complex(dp), intent(in) :: polarization(3), polarization_dot(3)
        real(dp), intent(in) :: wave_number
        integer, intent(in) :: triangles(:, :), quadrature_degree
        real(dp), intent(in) :: vertices_dot(:, :), radius_dot
        real(dp), intent(in) :: direction_dot(3), wave_number_dot
        complex(dp), allocatable, intent(out) :: right_hand_side(:)
        complex(dp), allocatable, intent(out) :: right_hand_side_dot(:)
        integer, intent(out) :: status

        integer, allocatable :: refined_triangles(:, :), refined_triangles_dot(:, :)
        real(dp), allocatable :: eta(:), refined_vertices(:, :)
        real(dp), allocatable :: refined_vertices_dot(:, :)
        real(dp), allocatable :: transformation(:, :), transformation_dot(:, :)
        real(dp), allocatable :: weights(:), xi(:)
        real(dp) :: divergence, divergence_dot, jacobian, jacobian_dot
        real(dp) :: local_value(3), local_value_dot(3), point(3), point_dot(3)
        real(dp) :: phase_argument_dot
        complex(dp) :: contraction, contraction_dot, incident_field(3)
        complex(dp) :: incident_field_dot(3), phase, phase_dot
        integer :: basis, local_edge, node, panel, row

        if (allocated(right_hand_side_dot)) deallocate(right_hand_side_dot)
        call assemble_maxwell_sphere_curved_plane_wave_rhs_bc_3d( &
            vertices, triangles, radius, direction, polarization, wave_number, &
            quadrature_degree, right_hand_side, status)
        if (status /= 0) return
        if (any(shape(vertices_dot) /= shape(vertices))) then
            status = 1
            return
        end if
        call build_maxwell_bc_transformation( &
            vertices, triangles, refined_vertices, refined_triangles, &
            transformation, status, sphere_radius=radius)
        if (status /= 0) return
        call barycentric_refine_sphere_surface_mesh_jvp( &
            vertices, triangles, radius, vertices_dot, radius_dot, &
            refined_vertices_dot, refined_triangles_dot, status)
        if (status /= 0) return
        if (any(refined_triangles_dot /= refined_triangles)) then
            status = 1
            return
        end if
        call differentiate_maxwell_bc_transformation_jvp( &
            vertices, triangles, refined_vertices, refined_triangles, &
            refined_vertices_dot, transformation_dot, status)
        if (status /= 0) return
        call triangle_duffy_quadrature( &
            quadrature_degree, xi, eta, weights, status)
        if (status /= 0) return
        allocate(right_hand_side_dot(size(transformation, 2)))
        right_hand_side_dot = cmplx(0.0_dp, 0.0_dp, dp)
        do panel = 1, size(refined_triangles, 2)
            do node = 1, size(weights)
                do local_edge = 1, 3
                    call evaluate_maxwell_sphere_curved_localized_rwg_basis_jvp( &
                        refined_vertices, refined_triangles, panel, local_edge, &
                        radius, xi(node), eta(node), refined_vertices_dot, radius_dot, &
                        point, local_value, divergence, jacobian, point_dot, &
                        local_value_dot, divergence_dot, jacobian_dot, status)
                    if (status /= 0) return
                    phase_argument_dot = wave_number_dot*dot_product( &
                        direction, point) + wave_number*( &
                        dot_product(direction_dot, point) + &
                        dot_product(direction, point_dot))
                    phase = exp(cmplx( &
                        0.0_dp, wave_number*dot_product(direction, point), dp))
                    phase_dot = cmplx(0.0_dp, phase_argument_dot, dp)*phase
                    incident_field = polarization*phase
                    incident_field_dot = polarization_dot*phase + &
                        polarization*phase_dot
                    contraction = sum(cmplx(local_value, 0.0_dp, dp)* &
                        incident_field)
                    contraction_dot = sum(cmplx(local_value_dot, 0.0_dp, dp)* &
                        incident_field + cmplx(local_value, 0.0_dp, dp)* &
                        incident_field_dot)
                    row = 3*(panel - 1) + local_edge
                    do basis = 1, size(transformation, 2)
                        right_hand_side_dot(basis) = &
                            right_hand_side_dot(basis) - weights(node)*( &
                            jacobian_dot*transformation(row, basis)*contraction + &
                            jacobian*(transformation_dot(row, basis)*contraction + &
                            transformation(row, basis)*contraction_dot))
                    end do
                end do
            end do
        end do
        status = 0
    end subroutine assemble_maxwell_sphere_curved_plane_wave_rhs_bc_3d_jvp

    subroutine assemble_maxwell_sphere_curved_plane_wave_rhs_bc_3d_vjp( &
            vertices, triangles, radius, direction, polarization, wave_number, &
            quadrature_degree, right_hand_side_bar, vertices_bar, radius_bar, &
            direction_bar, polarization_bar, wave_number_bar, status)
        real(dp), intent(in) :: vertices(:, :), radius, direction(3)
        complex(dp), intent(in) :: polarization(3), right_hand_side_bar(:)
        real(dp), intent(in) :: wave_number
        integer, intent(in) :: triangles(:, :), quadrature_degree
        real(dp), intent(out) :: vertices_bar(:, :), radius_bar
        real(dp), intent(out) :: direction_bar(3), wave_number_bar
        complex(dp), allocatable, intent(out) :: polarization_bar(:)
        integer, intent(out) :: status

        integer, allocatable :: refined_triangles(:, :)
        real(dp), allocatable :: eta(:), refined_vertices(:, :)
        real(dp), allocatable :: transformation(:, :), transformation_bar(:, :)
        real(dp), allocatable :: refined_vertices_bar(:, :)
        real(dp), allocatable :: local_vertices_bar(:, :)
        real(dp), allocatable :: weights(:), xi(:)
        real(dp) :: divergence, jacobian, jacobian_bar
        real(dp) :: local_value(3), local_value_bar(3), point(3), point_bar(3)
        real(dp) :: phase_argument_bar, local_radius_bar
        complex(dp) :: contraction, phase, phase_bar, seed
        integer :: basis, local_edge, node, panel, row

        vertices_bar = 0.0_dp
        radius_bar = 0.0_dp
        direction_bar = 0.0_dp
        wave_number_bar = 0.0_dp
        if (allocated(polarization_bar)) deallocate(polarization_bar)
        status = 1
        if (size(vertices, 1) /= 3 .or. radius <= 0.0_dp) return
        if (size(vertices_bar, 1) /= size(vertices, 1) .or. &
            size(vertices_bar, 2) /= size(vertices, 2)) return
        if (wave_number < 0.0_dp .or. quadrature_degree < 0) return
        if (abs(norm2(direction) - 1.0_dp) > 128.0_dp*epsilon(1.0_dp)) return
        if (abs(sum(polarization*direction)) > &
            128.0_dp*epsilon(1.0_dp)*max(1.0_dp, &
            sqrt(sum(abs(polarization)**2)))) return
        call build_maxwell_bc_transformation( &
            vertices, triangles, refined_vertices, refined_triangles, &
            transformation, status, sphere_radius=radius)
        if (status /= 0) return
        if (size(right_hand_side_bar) /= size(transformation, 2)) return
        call triangle_duffy_quadrature( &
            quadrature_degree, xi, eta, weights, status)
        if (status /= 0) return
        allocate( &
            polarization_bar(3), &
            transformation_bar(size(transformation, 1), size(transformation, 2)), &
            refined_vertices_bar(3, size(refined_vertices, 2)), &
            local_vertices_bar(3, size(refined_vertices, 2)))
        polarization_bar = cmplx(0.0_dp, 0.0_dp, dp)
        transformation_bar = 0.0_dp
        refined_vertices_bar = 0.0_dp
        do panel = 1, size(refined_triangles, 2)
            do node = 1, size(weights)
                do local_edge = 1, 3
                    call evaluate_maxwell_sphere_curved_localized_rwg_basis( &
                        refined_vertices, refined_triangles, panel, local_edge, &
                        radius, xi(node), eta(node), point, local_value, &
                        divergence, jacobian, status)
                    if (status /= 0) return
                    phase = exp(cmplx( &
                        0.0_dp, wave_number*dot_product(direction, point), dp))
                    contraction = sum(cmplx(local_value, 0.0_dp, dp)*polarization)
                    row = 3*(panel - 1) + local_edge
                    do basis = 1, size(transformation, 2)
                        seed = -weights(node)*jacobian* &
                            transformation(row, basis)*conjg(phase)* &
                            right_hand_side_bar(basis)
                        phase_bar = -weights(node)*jacobian* &
                            transformation(row, basis)*conjg(contraction)* &
                            right_hand_side_bar(basis)
                        jacobian_bar = real(conjg(right_hand_side_bar(basis))* &
                            (-weights(node)*transformation(row, basis)* &
                            phase*contraction), dp)
                        local_value_bar = real(conjg(seed)*polarization, dp)
                        polarization_bar = polarization_bar + &
                            seed*cmplx(local_value, 0.0_dp, dp)
                        phase_argument_bar = real(conjg(phase_bar)* &
                            cmplx(0.0_dp, 1.0_dp, dp)*phase, dp)
                        wave_number_bar = wave_number_bar + phase_argument_bar* &
                            dot_product(direction, point)
                        direction_bar = direction_bar + phase_argument_bar* &
                            wave_number*point
                        point_bar = phase_argument_bar*wave_number*direction
                        transformation_bar(row, basis) = &
                            transformation_bar(row, basis) + real( &
                            conjg(right_hand_side_bar(basis))* &
                            (-weights(node)*jacobian*phase*contraction), dp)
                        call evaluate_maxwell_sphere_curved_localized_rwg_basis_vjp( &
                            refined_vertices, refined_triangles, panel, local_edge, &
                            radius, xi(node), eta(node), point_bar, local_value_bar, &
                            0.0_dp, jacobian_bar, local_vertices_bar, &
                            local_radius_bar, status)
                        if (status /= 0) return
                        refined_vertices_bar = refined_vertices_bar + &
                            local_vertices_bar
                        radius_bar = radius_bar + local_radius_bar
                    end do
                end do
            end do
        end do
        call differentiate_maxwell_bc_transformation_vjp( &
            vertices, triangles, refined_vertices, refined_triangles, &
            transformation_bar, local_vertices_bar, status)
        if (status /= 0) return
        refined_vertices_bar = refined_vertices_bar + local_vertices_bar
        call barycentric_refine_sphere_surface_mesh_vjp( &
            vertices, triangles, radius, refined_vertices_bar, vertices_bar, &
            local_radius_bar, status)
        if (status /= 0) return
        radius_bar = radius_bar + local_radius_bar
        status = 0
    end subroutine assemble_maxwell_sphere_curved_plane_wave_rhs_bc_3d_vjp

    subroutine assemble_maxwell_sphere_curved_regularized_cfie_rwg_3d( &
            vertices, triangles, radius, wave_number, impedance, &
            quadrature_degree, tolerance, max_depth, mfie_offset, matrix, &
            efie, mfie, regularizer, regularized_efie, status)
        real(dp), intent(in) :: vertices(:, :), radius, wave_number, impedance
        real(dp), intent(in) :: tolerance, mfie_offset
        integer, intent(in) :: triangles(:, :), quadrature_degree, max_depth
        complex(dp), allocatable, intent(out) :: matrix(:, :), efie(:, :)
        complex(dp), allocatable, intent(out) :: mfie(:, :), regularizer(:, :)
        complex(dp), allocatable, intent(out) :: regularized_efie(:, :)
        integer, intent(out) :: status

        complex(dp), allocatable :: mass(:, :), mapped_efie(:, :)
        real(dp), allocatable :: real_mass(:, :)
        integer :: info, system_size

        status = 1
        if (radius <= 0.0_dp .or. wave_number <= 0.0_dp .or. &
            impedance <= 0.0_dp) return
        call assemble_maxwell_sphere_curved_efie_rwg_3d( &
            vertices, triangles, radius, wave_number, impedance, &
            quadrature_degree, tolerance, max_depth, efie, status)
        if (status /= 0) return
        call assemble_maxwell_sphere_curved_efie_bc_imaginary_3d( &
            vertices, triangles, radius, wave_number, impedance, &
            quadrature_degree, tolerance, max_depth, regularizer, status)
        if (status /= 0) return
        call assemble_maxwell_sphere_curved_mfie_rwg_rbc_3d( &
            vertices, triangles, radius, wave_number, quadrature_degree, &
            mfie_offset, mfie, status)
        if (status /= 0) return
        call assemble_maxwell_sphere_curved_rwg_rbc_pairing( &
            vertices, triangles, radius, quadrature_degree, real_mass, status)
        if (status /= 0) return
        system_size = size(real_mass, 1)
        if (size(real_mass, 2) /= system_size) return
        if (any(shape(efie) /= [system_size, system_size])) return
        if (any(shape(mfie) /= [system_size, system_size])) return
        if (any(shape(regularizer) /= [system_size, system_size])) return
        allocate( &
            mass(system_size, system_size), &
            mapped_efie(system_size, system_size))
        mass = transpose(cmplx(real_mass, 0.0_dp, dp))
        call dense_solve(mass, efie, mapped_efie, info)
        if (info /= 0) then
            status = 2
            return
        end if
        regularized_efie = matmul(regularizer, mapped_efie)
        matrix = mfie - regularized_efie
        status = 0
    end subroutine assemble_maxwell_sphere_curved_regularized_cfie_rwg_3d

    subroutine assemble_maxwell_sphere_curved_efie_bc_imaginary_3d( &
            vertices, triangles, radius, decay_rate, impedance, &
            quadrature_degree, tolerance, max_depth, matrix, status)
        real(dp), intent(in) :: vertices(:, :), radius, decay_rate, impedance
        integer, intent(in) :: triangles(:, :), quadrature_degree, max_depth
        real(dp), intent(in) :: tolerance
        complex(dp), allocatable, intent(out) :: matrix(:, :)
        integer, intent(out) :: status

        integer, allocatable :: refined_triangles(:, :)
        real(dp), allocatable :: refined_vertices(:, :), transformation(:, :)
        complex(dp), allocatable :: complex_transformation(:, :)
        complex(dp), allocatable :: refined_matrix(:, :)

        status = 1
        call build_maxwell_bc_to_refined_rwg( &
            vertices, triangles, refined_vertices, refined_triangles, &
            transformation, status, sphere_radius=radius)
        if (status /= 0) return
        call assemble_maxwell_sphere_curved_efie_imaginary_rwg_3d( &
            refined_vertices, refined_triangles, radius, decay_rate, impedance, &
            quadrature_degree, tolerance, max_depth, refined_matrix, status)
        if (status /= 0) return
        complex_transformation = cmplx(transformation, 0.0_dp, dp)
        matrix = matmul( &
            transpose(complex_transformation), &
            matmul(refined_matrix, complex_transformation))
        status = 0
    end subroutine assemble_maxwell_sphere_curved_efie_bc_imaginary_3d

    subroutine assemble_maxwell_sphere_curved_efie_imaginary_rwg_3d( &
            vertices, triangles, radius, decay_rate, impedance, &
            quadrature_degree, tolerance, max_depth, matrix, status)
        real(dp), intent(in) :: vertices(:, :), radius, decay_rate, impedance
        integer, intent(in) :: triangles(:, :), quadrature_degree, max_depth
        real(dp), intent(in) :: tolerance
        complex(dp), allocatable, intent(out) :: matrix(:, :)
        integer, intent(out) :: status

        complex(dp), allocatable :: scalar_potential(:, :)
        complex(dp), allocatable :: vector_potential(:, :)

        status = 1
        if (radius <= 0.0_dp .or. decay_rate <= 0.0_dp .or. &
            impedance <= 0.0_dp) return
        call assemble_maxwell_sphere_curved_potential_operators_rwg_3d( &
            vertices, triangles, radius, decay_rate, quadrature_degree, &
            tolerance, max_depth, vector_potential, scalar_potential, status, &
            decaying_kernel=.true.)
        if (status /= 0) return
        matrix = -impedance*( &
            decay_rate*vector_potential + scalar_potential/decay_rate)
        status = 0
    end subroutine assemble_maxwell_sphere_curved_efie_imaginary_rwg_3d

    subroutine assemble_maxwell_sphere_curved_mfie_rwg_rbc_3d( &
            vertices, triangles, radius, wave_number, quadrature_degree, &
            relative_offset, matrix, status)
        real(dp), intent(in) :: vertices(:, :), radius, wave_number
        real(dp), intent(in) :: relative_offset
        integer, intent(in) :: triangles(:, :), quadrature_degree
        complex(dp), allocatable, intent(out) :: matrix(:, :)
        integer, intent(out) :: status

        complex(dp), allocatable :: full_offset(:, :), half_offset(:, :)
        complex(dp), allocatable :: quarter_offset(:, :)

        status = 1
        call assemble_maxwell_sphere_curved_mfie_exterior_trace_rwg_rbc_3d( &
            vertices, triangles, radius, wave_number, quadrature_degree, &
            relative_offset, full_offset, status)
        if (status /= 0) return
        call assemble_maxwell_sphere_curved_mfie_exterior_trace_rwg_rbc_3d( &
            vertices, triangles, radius, wave_number, quadrature_degree, &
            relative_offset/2.0_dp, half_offset, status)
        if (status /= 0) return
        call assemble_maxwell_sphere_curved_mfie_exterior_trace_rwg_rbc_3d( &
            vertices, triangles, radius, wave_number, quadrature_degree, &
            relative_offset/4.0_dp, quarter_offset, status)
        if (status /= 0) return
        matrix = full_offset/3.0_dp - 2.0_dp*half_offset + &
            8.0_dp*quarter_offset/3.0_dp
        status = 0
    end subroutine assemble_maxwell_sphere_curved_mfie_rwg_rbc_3d

    subroutine assemble_maxwell_sphere_curved_mfie_exterior_trace_rwg_rbc_3d( &
            vertices, triangles, radius, wave_number, quadrature_degree, &
            relative_offset, matrix, status)
        real(dp), intent(in) :: vertices(:, :), radius, wave_number
        real(dp), intent(in) :: relative_offset
        integer, intent(in) :: triangles(:, :), quadrature_degree
        complex(dp), allocatable, intent(out) :: matrix(:, :)
        integer, intent(out) :: status

        integer, allocatable :: edge_triangles(:, :), edge_vertices(:, :)
        integer, allocatable :: refined_triangles(:, :)
        real(dp), allocatable :: eta(:), refined_vertices(:, :)
        real(dp), allocatable :: transformation(:, :), weights(:), xi(:)
        real(dp), allocatable :: bc_values(:, :)
        complex(dp), allocatable :: magnetic_fields(:, :)
        real(dp) :: divergence, jacobian, local_value(3), point(3), target(3)
        integer :: local_edge, node, refined_panel, row, test_basis, trial_basis

        status = 1
        if (radius <= 0.0_dp .or. wave_number < 0.0_dp .or. &
            relative_offset <= 0.0_dp) return
        call build_maxwell_bc_transformation( &
            vertices, triangles, refined_vertices, refined_triangles, &
            transformation, status, sphere_radius=radius)
        if (status /= 0) return
        call build_maxwell_rwg_surface_space( &
            vertices, triangles, edge_vertices, edge_triangles, status)
        if (status /= 0) return
        call triangle_duffy_quadrature( &
            quadrature_degree, xi, eta, weights, status)
        if (status /= 0) return
        allocate( &
            matrix(size(edge_vertices, 2), size(edge_vertices, 2)), &
            bc_values(3, size(edge_vertices, 2)), &
            magnetic_fields(3, size(edge_vertices, 2)))
        matrix = cmplx(0.0_dp, 0.0_dp, dp)
        do refined_panel = 1, size(refined_triangles, 2)
            do node = 1, size(weights)
                bc_values = 0.0_dp
                do local_edge = 1, 3
                    call evaluate_maxwell_sphere_curved_localized_rwg_basis( &
                        refined_vertices, refined_triangles, refined_panel, &
                        local_edge, radius, xi(node), eta(node), point, &
                        local_value, divergence, jacobian, status)
                    if (status /= 0) return
                    row = 3*(refined_panel - 1) + local_edge
                    do test_basis = 1, size(edge_vertices, 2)
                        bc_values(:, test_basis) = &
                            bc_values(:, test_basis) + &
                            transformation(row, test_basis)*local_value
                    end do
                end do
                target = (1.0_dp + relative_offset)*point
                call evaluate_all_curved_rwg_magnetic_fields( &
                    vertices, triangles, edge_vertices, edge_triangles, radius, &
                    target, wave_number, xi, eta, weights, magnetic_fields, &
                    status)
                if (status /= 0) return
                do test_basis = 1, size(edge_vertices, 2)
                    do trial_basis = 1, size(edge_vertices, 2)
                        matrix(test_basis, trial_basis) = &
                            matrix(test_basis, trial_basis) - &
                            jacobian*weights(node)*sum( &
                            bc_values(:, test_basis)* &
                            magnetic_fields(:, trial_basis))
                    end do
                end do
            end do
        end do
        status = 0
    end subroutine &
        assemble_maxwell_sphere_curved_mfie_exterior_trace_rwg_rbc_3d

    subroutine evaluate_all_curved_rwg_magnetic_fields( &
            vertices, triangles, edge_vertices, edge_triangles, radius, target, &
            wave_number, xi, eta, weights, magnetic_fields, status)
        real(dp), intent(in) :: vertices(:, :), radius, target(3), wave_number
        integer, intent(in) :: triangles(:, :), edge_vertices(:, :)
        integer, intent(in) :: edge_triangles(:, :)
        real(dp), intent(in) :: xi(:), eta(:), weights(:)
        complex(dp), intent(out) :: magnetic_fields(:, :)
        integer, intent(out) :: status

        real(dp) :: basis_value(3), displacement(3), divergence, jacobian
        real(dp) :: point(3), distance
        complex(dp) :: curl_integrand(3), gradient_green(3), green
        integer :: basis, node, panel

        magnetic_fields = cmplx(0.0_dp, 0.0_dp, dp)
        status = 1
        do panel = 1, size(triangles, 2)
            do node = 1, size(weights)
                do basis = 1, size(edge_vertices, 2)
                    if (.not. any(edge_triangles(:, basis) == panel)) cycle
                    call evaluate_maxwell_sphere_curved_rwg_basis( &
                        vertices, triangles, edge_vertices, edge_triangles, &
                        basis, panel, radius, xi(node), eta(node), point, &
                        basis_value, divergence, jacobian, status)
                    if (status /= 0) return
                    displacement = target - point
                    distance = norm2(displacement)
                    if (distance <= 128.0_dp*epsilon(1.0_dp)*radius) return
                    green = exp(cmplx(0.0_dp, wave_number*distance, dp))/ &
                        (4.0_dp*acos(-1.0_dp)*distance)
                    gradient_green = green* &
                        (cmplx(0.0_dp, wave_number, dp) - 1.0_dp/distance)* &
                        displacement/distance
                    curl_integrand = complex_cross_product( &
                        gradient_green, cmplx(basis_value, 0.0_dp, dp))
                    magnetic_fields(:, basis) = magnetic_fields(:, basis) + &
                        weights(node)*jacobian*curl_integrand
                end do
            end do
        end do
        status = 0
    end subroutine evaluate_all_curved_rwg_magnetic_fields

    subroutine assemble_maxwell_sphere_curved_rwg_rbc_pairing( &
            vertices, triangles, radius, quadrature_degree, matrix, status)
        real(dp), intent(in) :: vertices(:, :), radius
        integer, intent(in) :: triangles(:, :), quadrature_degree
        real(dp), allocatable, intent(out) :: matrix(:, :)
        integer, intent(out) :: status

        integer, allocatable :: edge_triangles(:, :), edge_vertices(:, :)
        integer, allocatable :: refined_triangles(:, :)
        real(dp), allocatable :: eta(:), refined_vertices(:, :)
        real(dp), allocatable :: transformation(:, :), weights(:), xi(:)
        real(dp), allocatable :: bc_values(:, :), rwg_values(:, :)
        real(dp) :: coarse_jacobian, coarse_vertices(3, 3), divergence
        real(dp) :: local_value(3), normal(3), point(3), refined_jacobian
        real(dp) :: rotated_bc(3)
        real(dp) :: coarse_eta, coarse_xi
        integer :: basis, local, local_edge, node, parent, refined_panel, row
        integer :: test_basis

        status = 1
        call build_maxwell_bc_transformation( &
            vertices, triangles, refined_vertices, refined_triangles, &
            transformation, status, sphere_radius=radius)
        if (status /= 0) return
        call build_maxwell_rwg_surface_space( &
            vertices, triangles, edge_vertices, edge_triangles, status)
        if (status /= 0) return
        call triangle_duffy_quadrature( &
            quadrature_degree, xi, eta, weights, status)
        if (status /= 0) return
        allocate( &
            matrix(size(edge_vertices, 2), size(edge_vertices, 2)), &
            bc_values(3, size(edge_vertices, 2)), &
            rwg_values(3, size(edge_vertices, 2)))
        matrix = 0.0_dp
        do refined_panel = 1, size(refined_triangles, 2)
            parent = (refined_panel - 1)/6 + 1
            do local = 1, 3
                coarse_vertices(:, local) = &
                    vertices(:, triangles(local, parent))
            end do
            do node = 1, size(weights)
                bc_values = 0.0_dp
                do local_edge = 1, 3
                    call evaluate_maxwell_sphere_curved_localized_rwg_basis( &
                        refined_vertices, refined_triangles, refined_panel, &
                        local_edge, radius, xi(node), eta(node), point, &
                        local_value, divergence, refined_jacobian, status)
                    if (status /= 0) return
                    row = 3*(refined_panel - 1) + local_edge
                    do test_basis = 1, size(edge_vertices, 2)
                        bc_values(:, test_basis) = &
                            bc_values(:, test_basis) + &
                            transformation(row, test_basis)*local_value
                    end do
                end do
                normal = point/radius
                call invert_sphere_curved_panel( &
                    coarse_vertices, radius, point, coarse_xi, coarse_eta, &
                    status)
                if (status /= 0) return
                rwg_values = 0.0_dp
                do basis = 1, size(edge_vertices, 2)
                    if (.not. any(edge_triangles(:, basis) == parent)) cycle
                    call evaluate_maxwell_sphere_curved_rwg_basis( &
                        vertices, triangles, edge_vertices, edge_triangles, &
                        basis, parent, radius, coarse_xi, coarse_eta, point, &
                        rwg_values(:, basis), divergence, coarse_jacobian, status)
                    if (status /= 0) return
                end do
                do test_basis = 1, size(edge_vertices, 2)
                    rotated_bc = real_cross_product( &
                        normal, bc_values(:, test_basis))
                    do basis = 1, size(edge_vertices, 2)
                        matrix(test_basis, basis) = matrix(test_basis, basis) + &
                            refined_jacobian*weights(node)*dot_product( &
                            rotated_bc, rwg_values(:, basis))
                    end do
                end do
            end do
        end do
        status = 0
    end subroutine assemble_maxwell_sphere_curved_rwg_rbc_pairing

    subroutine assemble_maxwell_sphere_curved_rwg_rbc_pairing_jvp( &
            vertices, triangles, radius, quadrature_degree, vertices_dot, &
            radius_dot, matrix, matrix_dot, status)
        real(dp), intent(in) :: vertices(:, :), radius
        integer, intent(in) :: triangles(:, :), quadrature_degree
        real(dp), intent(in) :: vertices_dot(:, :), radius_dot
        real(dp), allocatable, intent(out) :: matrix(:, :), matrix_dot(:, :)
        integer, intent(out) :: status

        integer, allocatable :: edge_triangles(:, :), edge_vertices(:, :)
        integer, allocatable :: refined_triangles(:, :), refined_triangles_dot(:, :)
        real(dp), allocatable :: eta(:), refined_vertices(:, :)
        real(dp), allocatable :: refined_vertices_dot(:, :)
        real(dp), allocatable :: transformation(:, :), transformation_dot(:, :)
        real(dp), allocatable :: weights(:), xi(:)
        real(dp), allocatable :: bc_values(:, :), bc_values_dot(:, :)
        real(dp), allocatable :: rwg_values(:, :), rwg_values_dot(:, :)
        real(dp) :: coarse_eta, coarse_eta_dot, coarse_xi, coarse_xi_dot
        real(dp) :: coarse_jacobian, coarse_jacobian_dot, divergence
        real(dp) :: divergence_dot, local_value(3), local_value_dot(3)
        real(dp) :: normal(3), normal_dot(3), point(3), point_dot(3)
        real(dp) :: refined_jacobian, refined_jacobian_dot
        real(dp) :: rotated_bc(3), rotated_bc_dot(3)
        integer :: basis, local_edge, node, parent, refined_panel, row
        integer :: test_basis

        status = 1
        if (allocated(matrix_dot)) deallocate(matrix_dot)
        call assemble_maxwell_sphere_curved_rwg_rbc_pairing( &
            vertices, triangles, radius, quadrature_degree, matrix, status)
        if (status /= 0) return
        if (radius <= 0.0_dp .or. any(shape(vertices_dot) /= shape(vertices))) return
        call build_maxwell_bc_transformation( &
            vertices, triangles, refined_vertices, refined_triangles, &
            transformation, status, sphere_radius=radius)
        if (status /= 0) return
        call barycentric_refine_sphere_surface_mesh_jvp( &
            vertices, triangles, radius, vertices_dot, radius_dot, &
            refined_vertices_dot, refined_triangles_dot, status)
        if (status /= 0) return
        if (any(refined_triangles_dot /= refined_triangles)) return
        call differentiate_maxwell_bc_transformation_jvp( &
            vertices, triangles, refined_vertices, refined_triangles, &
            refined_vertices_dot, transformation_dot, status)
        if (status /= 0) return
        call build_maxwell_rwg_surface_space( &
            vertices, triangles, edge_vertices, edge_triangles, status)
        if (status /= 0) return
        call triangle_duffy_quadrature( &
            quadrature_degree, xi, eta, weights, status)
        if (status /= 0) return
        allocate( &
            matrix_dot(size(edge_vertices, 2), size(edge_vertices, 2)), &
            bc_values(3, size(edge_vertices, 2)), &
            bc_values_dot(3, size(edge_vertices, 2)), &
            rwg_values(3, size(edge_vertices, 2)), &
            rwg_values_dot(3, size(edge_vertices, 2)))
        matrix_dot = 0.0_dp
        do refined_panel = 1, size(refined_triangles, 2)
            parent = (refined_panel - 1)/6 + 1
            do node = 1, size(weights)
                bc_values = 0.0_dp
                bc_values_dot = 0.0_dp
                do local_edge = 1, 3
                    call evaluate_maxwell_sphere_curved_localized_rwg_basis_jvp( &
                        refined_vertices, refined_triangles, refined_panel, &
                        local_edge, radius, xi(node), eta(node), &
                        refined_vertices_dot, radius_dot, point, local_value, &
                        divergence, refined_jacobian, point_dot, local_value_dot, &
                        divergence_dot, refined_jacobian_dot, status)
                    if (status /= 0) return
                    row = 3*(refined_panel - 1) + local_edge
                    do test_basis = 1, size(edge_vertices, 2)
                        bc_values(:, test_basis) = bc_values(:, test_basis) + &
                            transformation(row, test_basis)*local_value
                        bc_values_dot(:, test_basis) = &
                            bc_values_dot(:, test_basis) + &
                            transformation_dot(row, test_basis)*local_value + &
                            transformation(row, test_basis)*local_value_dot
                    end do
                end do
                call map_refined_sphere_point_to_parent_jvp( &
                    vertices(:, triangles(:, parent)), radius, point, &
                    vertices_dot(:, triangles(:, parent)), radius_dot, point_dot, &
                    coarse_xi, coarse_eta, coarse_xi_dot, coarse_eta_dot, status)
                if (status /= 0) return
                rwg_values = 0.0_dp
                rwg_values_dot = 0.0_dp
                do basis = 1, size(edge_vertices, 2)
                    if (.not. any(edge_triangles(:, basis) == parent)) cycle
                    call evaluate_maxwell_sphere_curved_rwg_basis_jvp( &
                        vertices, triangles, edge_vertices, edge_triangles, basis, &
                        parent, radius, coarse_xi, coarse_eta, vertices_dot, &
                        radius_dot, point, rwg_values(:, basis), divergence, &
                        coarse_jacobian, point_dot, local_value_dot, &
                        divergence_dot, coarse_jacobian_dot, status, coarse_xi_dot, &
                        coarse_eta_dot)
                    if (status /= 0) return
                    rwg_values_dot(:, basis) = local_value_dot
                end do
                normal = point/radius
                normal_dot = point_dot/radius - point*radius_dot/radius**2
                do test_basis = 1, size(edge_vertices, 2)
                    rotated_bc = real_cross_product( &
                        normal, bc_values(:, test_basis))
                    rotated_bc_dot = real_cross_product( &
                        normal_dot, bc_values(:, test_basis)) + &
                        real_cross_product(normal, bc_values_dot(:, test_basis))
                    do basis = 1, size(edge_vertices, 2)
                        matrix_dot(test_basis, basis) = &
                            matrix_dot(test_basis, basis) + &
                            refined_jacobian_dot*weights(node)*dot_product( &
                            rotated_bc, rwg_values(:, basis)) + &
                            refined_jacobian*weights(node)*( &
                            dot_product(rotated_bc_dot, rwg_values(:, basis)) + &
                            dot_product(rotated_bc, rwg_values_dot(:, basis)))
                    end do
                end do
            end do
        end do
        status = 0
    end subroutine assemble_maxwell_sphere_curved_rwg_rbc_pairing_jvp

    subroutine assemble_maxwell_sphere_curved_rwg_rbc_pairing_vjp( &
            vertices, triangles, radius, quadrature_degree, matrix_bar, &
            vertices_bar, radius_bar, status)
        real(dp), intent(in) :: vertices(:, :), radius
        integer, intent(in) :: triangles(:, :), quadrature_degree
        real(dp), intent(in) :: matrix_bar(:, :)
        real(dp), intent(out) :: vertices_bar(:, :), radius_bar
        integer, intent(out) :: status

        integer, allocatable :: edge_triangles(:, :), edge_vertices(:, :)
        integer, allocatable :: refined_triangles(:, :)
        real(dp), allocatable :: eta(:), refined_vertices(:, :)
        real(dp), allocatable :: transformation(:, :), transformation_bar(:, :)
        real(dp), allocatable :: weights(:), xi(:)
        real(dp), allocatable :: bc_values(:, :), bc_values_bar(:, :)
        real(dp), allocatable :: rwg_values(:, :), rwg_values_bar(:, :)
        real(dp), allocatable :: refined_vertices_bar(:, :)
        real(dp), allocatable :: local_refined_vertices_bar(:, :)
        real(dp), allocatable :: barycentric_vertices_bar(:, :)
        real(dp) :: coarse_eta, coarse_eta_bar, coarse_xi, coarse_xi_bar
        real(dp) :: coarse_jacobian, divergence, local_radius_bar
        real(dp) :: local_value(3), local_value_bar(3)
        real(dp) :: local_edge_vertices_bar(3, 2), local_panel_vertices_bar(3, 3)
        real(dp) :: local_point_bar(3), localized_values(3, 3)
        real(dp) :: map_point_bar(3), map_parent_vertices_bar(3, 3)
        real(dp) :: normal(3), normal_bar(3), point(3), point_bar(3)
        real(dp) :: refined_jacobian, refined_jacobian_bar
        real(dp) :: rotated_bc(3), rotated_bc_bar(3)
        real(dp) :: zero_divergence_bar, zero_jacobian_bar
        real(dp) :: local_xi_bar, local_eta_bar
        integer :: basis, local_edge, node, parent, refined_panel, row
        integer :: last_basis, test_basis, vertex

        vertices_bar = 0.0_dp
        radius_bar = 0.0_dp
        zero_divergence_bar = 0.0_dp
        zero_jacobian_bar = 0.0_dp
        status = 1
        if (size(vertices, 1) /= 3 .or. radius <= 0.0_dp) return
        call build_maxwell_bc_transformation( &
            vertices, triangles, refined_vertices, refined_triangles, &
            transformation, status, sphere_radius=radius)
        if (status /= 0) return
        call build_maxwell_rwg_surface_space( &
            vertices, triangles, edge_vertices, edge_triangles, status)
        if (status /= 0) return
        if (any(shape(matrix_bar) /= &
            [size(edge_vertices, 2), size(edge_vertices, 2)])) return
        call triangle_duffy_quadrature( &
            quadrature_degree, xi, eta, weights, status)
        if (status /= 0) return
        allocate( &
            transformation_bar(size(transformation, 1), size(transformation, 2)), &
            bc_values(3, size(edge_vertices, 2)), &
            bc_values_bar(3, size(edge_vertices, 2)), &
            rwg_values(3, size(edge_vertices, 2)), &
            rwg_values_bar(3, size(edge_vertices, 2)), &
            refined_vertices_bar(3, size(refined_vertices, 2)), &
            local_refined_vertices_bar(3, size(refined_vertices, 2)))
        transformation_bar = 0.0_dp
        refined_vertices_bar = 0.0_dp
        do refined_panel = 1, size(refined_triangles, 2)
            parent = (refined_panel - 1)/6 + 1
            do node = 1, size(weights)
                bc_values = 0.0_dp
                rwg_values = 0.0_dp
                do local_edge = 1, 3
                    call evaluate_maxwell_sphere_curved_localized_rwg_basis( &
                        refined_vertices, refined_triangles, refined_panel, &
                        local_edge, radius, xi(node), eta(node), point, &
                        local_value, divergence, refined_jacobian, status)
                    if (status /= 0) return
                    localized_values(:, local_edge) = local_value
                    row = 3*(refined_panel - 1) + local_edge
                    do test_basis = 1, size(edge_vertices, 2)
                        bc_values(:, test_basis) = &
                            bc_values(:, test_basis) + &
                            transformation(row, test_basis)*local_value
                    end do
                end do
                call map_refined_sphere_point_to_parent( &
                    vertices(:, triangles(:, parent)), radius, point, coarse_xi, &
                    coarse_eta, status)
                if (status /= 0) return
                rwg_values = 0.0_dp
                last_basis = 0
                do basis = 1, size(edge_vertices, 2)
                    if (.not. any(edge_triangles(:, basis) == parent)) cycle
                    last_basis = basis
                    call evaluate_maxwell_sphere_curved_rwg_basis( &
                        vertices, triangles, edge_vertices, edge_triangles, basis, &
                        parent, radius, coarse_xi, coarse_eta, point, &
                        rwg_values(:, basis), divergence, coarse_jacobian, status)
                    if (status /= 0) return
                end do
                normal = point/radius
                bc_values_bar = 0.0_dp
                rwg_values_bar = 0.0_dp
                normal_bar = 0.0_dp
                refined_jacobian_bar = 0.0_dp
                do test_basis = 1, size(edge_vertices, 2)
                    rotated_bc = real_cross_product( &
                        normal, bc_values(:, test_basis))
                    do basis = 1, size(edge_vertices, 2)
                        if (.not. any(edge_triangles(:, basis) == parent)) cycle
                        rotated_bc_bar = weights(node)*refined_jacobian* &
                            matrix_bar(test_basis, basis)*rwg_values(:, basis)
                        rwg_values_bar(:, basis) = rwg_values_bar(:, basis) + &
                            weights(node)*refined_jacobian* &
                            matrix_bar(test_basis, basis)*rotated_bc
                        normal_bar = normal_bar + real_cross_product( &
                            bc_values(:, test_basis), rotated_bc_bar)
                        bc_values_bar(:, test_basis) = &
                            bc_values_bar(:, test_basis) + &
                            real_cross_product(rotated_bc_bar, normal)
                        refined_jacobian_bar = refined_jacobian_bar + &
                            weights(node)*matrix_bar(test_basis, basis)* &
                            dot_product(rotated_bc, rwg_values(:, basis))
                    end do
                end do
                point_bar = normal_bar/radius
                local_radius_bar = -dot_product(normal_bar, point)/radius**2
                radius_bar = radius_bar + local_radius_bar
                coarse_xi_bar = 0.0_dp
                coarse_eta_bar = 0.0_dp
                do basis = 1, size(edge_vertices, 2)
                    if (.not. any(edge_triangles(:, basis) == parent)) cycle
                    local_xi_bar = 0.0_dp
                    local_eta_bar = 0.0_dp
                    local_point_bar = 0.0_dp
                    if (basis == last_basis) local_point_bar = point_bar
                    call evaluate_maxwell_sphere_curved_rwg_basis_vjp( &
                        vertices, triangles, edge_vertices, edge_triangles, basis, &
                        parent, radius, coarse_xi, coarse_eta, local_point_bar, &
                        rwg_values_bar(:, basis), zero_divergence_bar, &
                        zero_jacobian_bar, local_edge_vertices_bar, &
                        local_panel_vertices_bar, local_radius_bar, status, &
                        local_xi_bar, local_eta_bar)
                    if (status /= 0) return
                    coarse_xi_bar = coarse_xi_bar + local_xi_bar
                    coarse_eta_bar = coarse_eta_bar + local_eta_bar
                    do vertex = 1, 2
                        vertices_bar(:, edge_vertices(vertex, basis)) = &
                            vertices_bar(:, edge_vertices(vertex, basis)) + &
                            local_edge_vertices_bar(:, vertex)
                    end do
                    do vertex = 1, 3
                        vertices_bar(:, triangles(vertex, parent)) = &
                            vertices_bar(:, triangles(vertex, parent)) + &
                            local_panel_vertices_bar(:, vertex)
                    end do
                    radius_bar = radius_bar + local_radius_bar
                end do
                call map_refined_sphere_point_to_parent_vjp( &
                    vertices(:, triangles(:, parent)), radius, point, coarse_xi_bar, &
                    coarse_eta_bar, map_point_bar, map_parent_vertices_bar, &
                    local_radius_bar, status)
                if (status /= 0) return
                do vertex = 1, 3
                    vertices_bar(:, triangles(vertex, parent)) = &
                        vertices_bar(:, triangles(vertex, parent)) + &
                        map_parent_vertices_bar(:, vertex)
                end do
                radius_bar = radius_bar + local_radius_bar
                do local_edge = 1, 3
                    row = 3*(refined_panel - 1) + local_edge
                    local_value_bar = 0.0_dp
                    local_point_bar = 0.0_dp
                    if (local_edge == 3) local_point_bar = map_point_bar
                    do test_basis = 1, size(edge_vertices, 2)
                        local_value_bar = local_value_bar + &
                            transformation(row, test_basis)* &
                            bc_values_bar(:, test_basis)
                        transformation_bar(row, test_basis) = &
                            transformation_bar(row, test_basis) + &
                            dot_product(bc_values_bar(:, test_basis), &
                            localized_values(:, local_edge))
                    end do
                    call evaluate_maxwell_sphere_curved_localized_rwg_basis_vjp( &
                        refined_vertices, refined_triangles, refined_panel, &
                        local_edge, radius, xi(node), eta(node), local_point_bar, &
                        local_value_bar, zero_divergence_bar, &
                        merge(refined_jacobian_bar, zero_jacobian_bar, &
                        local_edge == 3), local_refined_vertices_bar, &
                        local_radius_bar, status)
                    if (status /= 0) return
                    refined_vertices_bar = refined_vertices_bar + &
                        local_refined_vertices_bar
                    radius_bar = radius_bar + local_radius_bar
                end do
            end do
        end do
        allocate( &
            barycentric_vertices_bar(size(vertices, 1), size(vertices, 2)))
        call differentiate_maxwell_bc_transformation_vjp( &
            vertices, triangles, refined_vertices, refined_triangles, &
            transformation_bar, local_refined_vertices_bar, status)
        if (status /= 0) return
        refined_vertices_bar = refined_vertices_bar + local_refined_vertices_bar
        call barycentric_refine_sphere_surface_mesh_vjp( &
            vertices, triangles, radius, refined_vertices_bar, &
            barycentric_vertices_bar, local_radius_bar, status)
        if (status /= 0) return
        vertices_bar = vertices_bar + barycentric_vertices_bar
        radius_bar = radius_bar + local_radius_bar
        status = 0
    end subroutine assemble_maxwell_sphere_curved_rwg_rbc_pairing_vjp

    pure subroutine evaluate_maxwell_sphere_curved_localized_rwg_basis( &
            vertices, triangles, panel, local_edge, radius, xi, eta, point, &
            value, surface_divergence, surface_jacobian, status)
        real(dp), intent(in) :: vertices(:, :), radius, xi, eta
        integer, intent(in) :: triangles(:, :), panel, local_edge
        real(dp), intent(out) :: point(3), value(3), surface_divergence
        real(dp), intent(out) :: surface_jacobian
        integer, intent(out) :: status

        real(dp) :: edge_length, opposite_coordinates(2), panel_vertices(3, 3)
        real(dp) :: tangent_eta(3), tangent_xi(3)
        integer :: edge_vertices(2), local, opposite

        point = 0.0_dp
        value = 0.0_dp
        surface_divergence = 0.0_dp
        surface_jacobian = 0.0_dp
        status = 1
        if (panel < 1 .or. panel > size(triangles, 2)) return
        if (local_edge < 1 .or. local_edge > 3) return
        select case (local_edge)
        case (1)
            edge_vertices = [1, 2]
            opposite = 3
        case (2)
            edge_vertices = [3, 1]
            opposite = 2
        case (3)
            edge_vertices = [2, 3]
            opposite = 1
        end select
        do local = 1, 3
            panel_vertices(:, local) = vertices(:, triangles(local, panel))
        end do
        call evaluate_sphere_curved_panel( &
            panel_vertices, radius, xi, eta, point, tangent_xi, tangent_eta, &
            surface_jacobian, status)
        if (status /= 0) return
        select case (opposite)
        case (1)
            opposite_coordinates = [0.0_dp, 0.0_dp]
        case (2)
            opposite_coordinates = [1.0_dp, 0.0_dp]
        case (3)
            opposite_coordinates = [0.0_dp, 1.0_dp]
        end select
        edge_length = norm2( &
            panel_vertices(:, edge_vertices(2)) - &
            panel_vertices(:, edge_vertices(1)))
        value = edge_length/surface_jacobian*( &
            (xi - opposite_coordinates(1))*tangent_xi + &
            (eta - opposite_coordinates(2))*tangent_eta)
        surface_divergence = 2.0_dp*edge_length/surface_jacobian
        status = 0
    end subroutine evaluate_maxwell_sphere_curved_localized_rwg_basis

    pure subroutine evaluate_maxwell_sphere_curved_localized_rwg_basis_jvp( &
            vertices, triangles, panel, local_edge, radius, xi, eta, &
            vertices_dot, radius_dot, point, value, surface_divergence, &
            surface_jacobian, point_dot, value_dot, surface_divergence_dot, &
            surface_jacobian_dot, status)
        real(dp), intent(in) :: vertices(:, :), radius, xi, eta
        integer, intent(in) :: triangles(:, :), panel, local_edge
        real(dp), intent(in) :: vertices_dot(:, :), radius_dot
        real(dp), intent(out) :: point(3), value(3), surface_divergence
        real(dp), intent(out) :: surface_jacobian, point_dot(3), value_dot(3)
        real(dp), intent(out) :: surface_divergence_dot, surface_jacobian_dot
        integer, intent(out) :: status

        real(dp) :: edge(3), edge_dot(3), edge_length, edge_length_dot
        real(dp) :: panel_vertices(3, 3), panel_vertices_dot(3, 3)
        real(dp) :: tangent_eta(3), tangent_eta_dot(3)
        real(dp) :: tangent_xi(3), tangent_xi_dot(3)
        real(dp) :: vector(3), vector_dot(3), coefficient, coefficient_dot
        real(dp) :: opposite_coordinates(2)
        integer :: edge_local_vertices(2), local, opposite

        point = 0.0_dp
        value = 0.0_dp
        surface_divergence = 0.0_dp
        surface_jacobian = 0.0_dp
        point_dot = 0.0_dp
        value_dot = 0.0_dp
        surface_divergence_dot = 0.0_dp
        surface_jacobian_dot = 0.0_dp
        status = 1
        if (size(vertices, 1) /= 3 .or. &
            any(shape(vertices_dot) /= shape(vertices))) return
        if (panel < 1 .or. panel > size(triangles, 2)) return
        if (local_edge < 1 .or. local_edge > 3) return
        select case (local_edge)
        case (1)
            edge_local_vertices = [1, 2]
            opposite = 3
        case (2)
            edge_local_vertices = [3, 1]
            opposite = 2
        case (3)
            edge_local_vertices = [2, 3]
            opposite = 1
        end select
        do local = 1, 3
            panel_vertices(:, local) = vertices(:, triangles(local, panel))
            panel_vertices_dot(:, local) = &
                vertices_dot(:, triangles(local, panel))
        end do
        call evaluate_sphere_curved_panel( &
            panel_vertices, radius, xi, eta, point, tangent_xi, tangent_eta, &
            surface_jacobian, status)
        if (status /= 0) return
        call evaluate_sphere_curved_panel_jvp( &
            panel_vertices, radius, xi, eta, panel_vertices_dot, radius_dot, &
            0.0_dp, 0.0_dp, point_dot, tangent_xi_dot, tangent_eta_dot, &
            surface_jacobian_dot, status)
        if (status /= 0) return
        select case (opposite)
        case (1)
            opposite_coordinates = [0.0_dp, 0.0_dp]
        case (2)
            opposite_coordinates = [1.0_dp, 0.0_dp]
        case (3)
            opposite_coordinates = [0.0_dp, 1.0_dp]
        end select
        edge = panel_vertices(:, edge_local_vertices(2)) - &
            panel_vertices(:, edge_local_vertices(1))
        edge_dot = panel_vertices_dot(:, edge_local_vertices(2)) - &
            panel_vertices_dot(:, edge_local_vertices(1))
        edge_length = norm2(edge)
        if (edge_length <= tiny(1.0_dp)) return
        edge_length_dot = dot_product(edge, edge_dot)/edge_length
        vector = (xi - opposite_coordinates(1))*tangent_xi + &
            (eta - opposite_coordinates(2))*tangent_eta
        vector_dot = (xi - opposite_coordinates(1))*tangent_xi_dot + &
            (eta - opposite_coordinates(2))*tangent_eta_dot
        coefficient = edge_length/surface_jacobian
        coefficient_dot = edge_length_dot/surface_jacobian - &
            edge_length*surface_jacobian_dot/surface_jacobian**2
        value = coefficient*vector
        value_dot = coefficient_dot*vector + coefficient*vector_dot
        surface_divergence = 2.0_dp*coefficient
        surface_divergence_dot = 2.0_dp*coefficient_dot
        status = 0
    end subroutine evaluate_maxwell_sphere_curved_localized_rwg_basis_jvp

    pure subroutine evaluate_maxwell_sphere_curved_localized_rwg_basis_vjp( &
            vertices, triangles, panel, local_edge, radius, xi, eta, point_bar, &
            value_bar, surface_divergence_bar, surface_jacobian_bar, &
            vertices_bar, radius_bar, status)
        real(dp), intent(in) :: vertices(:, :), radius, xi, eta
        integer, intent(in) :: triangles(:, :), panel, local_edge
        real(dp), intent(in) :: point_bar(3), value_bar(3)
        real(dp), intent(in) :: surface_divergence_bar, surface_jacobian_bar
        real(dp), intent(out) :: vertices_bar(:, :), radius_bar
        integer, intent(out) :: status

        real(dp) :: edge(3), edge_bar(3), edge_length, edge_length_bar
        real(dp) :: panel_vertices(3, 3), panel_vertices_bar(3, 3)
        real(dp) :: tangent_eta(3), tangent_eta_bar(3)
        real(dp) :: tangent_xi(3), tangent_xi_bar(3)
        real(dp) :: vector(3), vector_bar(3), coefficient
        real(dp) :: coefficient_bar, jacobian_bar
        real(dp) :: opposite_coordinates(2), surface_jacobian
        real(dp) :: point_local(3), xi_local_bar, eta_local_bar
        integer :: edge_local_vertices(2), local, opposite

        vertices_bar = 0.0_dp
        radius_bar = 0.0_dp
        status = 1
        if (size(vertices, 1) /= 3 .or. &
            size(vertices_bar, 1) /= size(vertices, 1) .or. &
            size(vertices_bar, 2) /= size(vertices, 2)) return
        if (panel < 1 .or. panel > size(triangles, 2)) return
        if (local_edge < 1 .or. local_edge > 3) return
        select case (local_edge)
        case (1)
            edge_local_vertices = [1, 2]
            opposite = 3
        case (2)
            edge_local_vertices = [3, 1]
            opposite = 2
        case (3)
            edge_local_vertices = [2, 3]
            opposite = 1
        end select
        do local = 1, 3
            panel_vertices(:, local) = vertices(:, triangles(local, panel))
        end do
        call evaluate_sphere_curved_panel( &
            panel_vertices, radius, xi, eta, point_local, tangent_xi, &
            tangent_eta, surface_jacobian, status)
        if (status /= 0) return
        select case (opposite)
        case (1)
            opposite_coordinates = [0.0_dp, 0.0_dp]
        case (2)
            opposite_coordinates = [1.0_dp, 0.0_dp]
        case (3)
            opposite_coordinates = [0.0_dp, 1.0_dp]
        end select
        edge = panel_vertices(:, edge_local_vertices(2)) - &
            panel_vertices(:, edge_local_vertices(1))
        edge_length = norm2(edge)
        if (edge_length <= tiny(1.0_dp)) return
        coefficient = edge_length/surface_jacobian
        vector = (xi - opposite_coordinates(1))*tangent_xi + &
            (eta - opposite_coordinates(2))*tangent_eta
        coefficient_bar = dot_product(value_bar, vector)
        vector_bar = coefficient*value_bar
        edge_length_bar = coefficient_bar/surface_jacobian + &
            2.0_dp*surface_divergence_bar/surface_jacobian
        jacobian_bar = surface_jacobian_bar - &
            edge_length*coefficient_bar/surface_jacobian**2 - &
            2.0_dp*edge_length*surface_divergence_bar/surface_jacobian**2
        edge_bar = edge_length_bar*edge/edge_length
        panel_vertices_bar = 0.0_dp
        tangent_xi_bar = (xi - opposite_coordinates(1))*vector_bar
        tangent_eta_bar = (eta - opposite_coordinates(2))*vector_bar
        xi_local_bar = dot_product(vector_bar, tangent_xi)
        eta_local_bar = dot_product(vector_bar, tangent_eta)
        call evaluate_sphere_curved_panel_vjp( &
            panel_vertices, radius, xi, eta, point_bar, tangent_xi_bar, &
            tangent_eta_bar, jacobian_bar, panel_vertices_bar, radius_bar, &
            xi_local_bar, eta_local_bar, status)
        if (status /= 0) return
        panel_vertices_bar(:, edge_local_vertices(1)) = &
            panel_vertices_bar(:, edge_local_vertices(1)) - edge_bar
        panel_vertices_bar(:, edge_local_vertices(2)) = &
            panel_vertices_bar(:, edge_local_vertices(2)) + edge_bar
        do local = 1, 3
            vertices_bar(:, triangles(local, panel)) = &
                vertices_bar(:, triangles(local, panel)) + panel_vertices_bar(:, local)
        end do
        status = 0
    end subroutine evaluate_maxwell_sphere_curved_localized_rwg_basis_vjp

    subroutine evaluate_maxwell_sphere_curved_magnetic_field_rwg_3d( &
            vertices, triangles, radius, coefficients, observation, &
            wave_number, quadrature_degree, magnetic_field, status)
        real(dp), intent(in) :: vertices(:, :), radius, observation(3)
        complex(dp), intent(in) :: coefficients(:)
        real(dp), intent(in) :: wave_number
        integer, intent(in) :: triangles(:, :), quadrature_degree
        complex(dp), intent(out) :: magnetic_field(3)
        integer, intent(out) :: status

        integer, allocatable :: edge_triangles(:, :), edge_vertices(:, :)
        real(dp), allocatable :: eta(:), weights(:), xi(:)
        real(dp) :: basis_value(3), displacement(3), divergence, jacobian
        real(dp) :: point(3), distance
        complex(dp) :: curl_integrand(3), gradient_green(3), green
        complex(dp) :: surface_current(3)
        integer :: basis, node, panel

        magnetic_field = cmplx(0.0_dp, 0.0_dp, dp)
        status = 1
        if (radius <= 0.0_dp .or. wave_number < 0.0_dp) return
        call build_maxwell_rwg_surface_space( &
            vertices, triangles, edge_vertices, edge_triangles, status)
        if (status /= 0 .or. size(coefficients) /= size(edge_vertices, 2)) return
        call triangle_duffy_quadrature( &
            quadrature_degree, xi, eta, weights, status)
        if (status /= 0) return
        do panel = 1, size(triangles, 2)
            do node = 1, size(weights)
                surface_current = cmplx(0.0_dp, 0.0_dp, dp)
                do basis = 1, size(edge_vertices, 2)
                    if (.not. any(edge_triangles(:, basis) == panel)) cycle
                    call evaluate_maxwell_sphere_curved_rwg_basis( &
                        vertices, triangles, edge_vertices, edge_triangles, &
                        basis, panel, radius, xi(node), eta(node), point, &
                        basis_value, divergence, jacobian, status)
                    if (status /= 0) return
                    surface_current = surface_current + &
                        coefficients(basis)*basis_value
                end do
                displacement = observation - point
                distance = norm2(displacement)
                if (distance <= 128.0_dp*epsilon(1.0_dp)*radius) return
                green = exp(cmplx(0.0_dp, wave_number*distance, dp))/ &
                    (4.0_dp*acos(-1.0_dp)*distance)
                gradient_green = green* &
                    (cmplx(0.0_dp, wave_number, dp) - 1.0_dp/distance)* &
                    displacement/distance
                curl_integrand = complex_cross_product( &
                    gradient_green, surface_current)
                magnetic_field = magnetic_field + &
                    weights(node)*jacobian*curl_integrand
            end do
        end do
        status = 0
    end subroutine evaluate_maxwell_sphere_curved_magnetic_field_rwg_3d

    subroutine solve_maxwell_pec_sphere_curved_efie_rwg_3d( &
            vertices, triangles, radius, direction, polarization, wave_number, &
            impedance, quadrature_degree, tolerance, max_depth, density, status)
        real(dp), intent(in) :: vertices(:, :), radius, direction(3)
        complex(dp), intent(in) :: polarization(3)
        real(dp), intent(in) :: wave_number, impedance, tolerance
        integer, intent(in) :: triangles(:, :), quadrature_degree, max_depth
        complex(dp), allocatable, intent(out) :: density(:)
        integer, intent(out) :: status

        complex(dp), allocatable :: matrix(:, :), right_hand_side(:)
        integer :: info

        status = 1
        call assemble_maxwell_sphere_curved_efie_rwg_3d( &
            vertices, triangles, radius, wave_number, impedance, &
            quadrature_degree, tolerance, max_depth, matrix, status)
        if (status /= 0) return
        call assemble_maxwell_sphere_curved_plane_wave_rhs_rwg_3d( &
            vertices, triangles, radius, direction, polarization, wave_number, &
            quadrature_degree, right_hand_side, status)
        if (status /= 0) return
        allocate(density(size(right_hand_side)))
        call dense_solve(matrix, right_hand_side, density, info)
        if (info /= 0) then
            status = 2
            return
        end if
        status = 0
    end subroutine solve_maxwell_pec_sphere_curved_efie_rwg_3d

    subroutine assemble_maxwell_sphere_curved_efie_rwg_3d( &
            vertices, triangles, radius, wave_number, impedance, &
            quadrature_degree, tolerance, max_depth, matrix, status)
        real(dp), intent(in) :: vertices(:, :), radius, wave_number, impedance
        integer, intent(in) :: triangles(:, :), quadrature_degree, max_depth
        real(dp), intent(in) :: tolerance
        complex(dp), allocatable, intent(out) :: matrix(:, :)
        integer, intent(out) :: status

        complex(dp), allocatable :: scalar_potential(:, :)
        complex(dp), allocatable :: vector_potential(:, :)

        status = 1
        if (radius <= 0.0_dp .or. wave_number <= 0.0_dp .or. &
            impedance <= 0.0_dp) return
        call assemble_maxwell_sphere_curved_potential_operators_rwg_3d( &
            vertices, triangles, radius, wave_number, quadrature_degree, &
            tolerance, max_depth, vector_potential, scalar_potential, status)
        if (status /= 0) return
        matrix = cmplx(0.0_dp, wave_number*impedance, dp)*vector_potential - &
            cmplx(0.0_dp, impedance/wave_number, dp)*scalar_potential
        status = 0
    end subroutine assemble_maxwell_sphere_curved_efie_rwg_3d

    subroutine assemble_maxwell_sphere_curved_potential_operators_rwg_3d( &
            vertices, triangles, radius, wave_number, quadrature_degree, &
            tolerance, max_depth, vector_potential, scalar_potential, status, &
            decaying_kernel)
        real(dp), intent(in) :: vertices(:, :), radius, wave_number, tolerance
        integer, intent(in) :: triangles(:, :), quadrature_degree, max_depth
        complex(dp), allocatable, intent(out) :: vector_potential(:, :)
        complex(dp), allocatable, intent(out) :: scalar_potential(:, :)
        integer, intent(out) :: status
        logical, optional, intent(in) :: decaying_kernel

        integer, allocatable :: edge_triangles(:, :), edge_vertices(:, :)
        real(dp), allocatable :: reference_divergence(:, :)
        complex(dp) :: contribution, scalar_green
        real(dp) :: divergence, jacobian, point(3), rwg_value(3)
        integer :: basis, first, first_panel, first_slot, second, second_panel
        integer :: second_slot
        logical :: use_decaying_kernel

        status = 1
        use_decaying_kernel = .false.
        if (present(decaying_kernel)) use_decaying_kernel = decaying_kernel
        call build_maxwell_rwg_surface_space( &
            vertices, triangles, edge_vertices, edge_triangles, status)
        if (status /= 0 .or. size(edge_vertices, 2) == 0) return
        allocate( &
            vector_potential(size(edge_vertices, 2), size(edge_vertices, 2)), &
            scalar_potential(size(edge_vertices, 2), size(edge_vertices, 2)), &
            reference_divergence(2, size(edge_vertices, 2)))
        vector_potential = cmplx(0.0_dp, 0.0_dp, dp)
        scalar_potential = cmplx(0.0_dp, 0.0_dp, dp)
        do basis = 1, size(edge_vertices, 2)
            do first_slot = 1, 2
                call evaluate_maxwell_sphere_curved_rwg_basis( &
                    vertices, triangles, edge_vertices, edge_triangles, basis, &
                    edge_triangles(first_slot, basis), radius, 1.0_dp/3.0_dp, &
                    1.0_dp/3.0_dp, point, rwg_value, divergence, jacobian, &
                    status)
                if (status /= 0) return
                reference_divergence(first_slot, basis) = divergence*jacobian
            end do
        end do
        do first = 1, size(edge_vertices, 2)
            do second = 1, first
                do first_slot = 1, 2
                    first_panel = edge_triangles(first_slot, first)
                    do second_slot = 1, 2
                        second_panel = edge_triangles(second_slot, second)
                        if (first_panel == second_panel) then
                            call &
                                integrate_maxwell_sphere_curved_coincident_rwg_pair_3d( &
                                vertices, triangles, edge_vertices, &
                                edge_triangles, first, first_panel, second, &
                                radius, wave_number, quadrature_degree, &
                                contribution, status, scalar_green, &
                                use_decaying_kernel)
                        else
                            call &
                                integrate_maxwell_sphere_curved_adjacent_rwg_pair_3d( &
                                vertices, triangles, edge_vertices, &
                                edge_triangles, first, first_panel, second, &
                                second_panel, radius, wave_number, &
                                quadrature_degree, tolerance, max_depth, &
                                contribution, status, scalar_green, &
                                use_decaying_kernel)
                        end if
                        if (status /= 0) return
                        vector_potential(first, second) = &
                            vector_potential(first, second) + contribution
                        scalar_potential(first, second) = &
                            scalar_potential(first, second) + &
                            reference_divergence(first_slot, first)* &
                            reference_divergence(second_slot, second)* &
                            scalar_green
                    end do
                end do
                vector_potential(second, first) = vector_potential(first, second)
                scalar_potential(second, first) = scalar_potential(first, second)
            end do
        end do
        status = 0
    end subroutine assemble_maxwell_sphere_curved_potential_operators_rwg_3d

    subroutine assemble_maxwell_sphere_curved_vector_potential_rwg_3d( &
            vertices, triangles, radius, wave_number, quadrature_degree, &
            tolerance, max_depth, matrix, status)
        real(dp), intent(in) :: vertices(:, :), radius, wave_number, tolerance
        integer, intent(in) :: triangles(:, :), quadrature_degree, max_depth
        complex(dp), allocatable, intent(out) :: matrix(:, :)
        integer, intent(out) :: status

        integer, allocatable :: edge_triangles(:, :), edge_vertices(:, :)
        complex(dp) :: contribution
        integer :: first, first_panel, first_slot, second, second_panel
        integer :: second_slot

        status = 1
        call build_maxwell_rwg_surface_space( &
            vertices, triangles, edge_vertices, edge_triangles, status)
        if (status /= 0 .or. size(edge_vertices, 2) == 0) return
        allocate(matrix(size(edge_vertices, 2), size(edge_vertices, 2)))
        matrix = cmplx(0.0_dp, 0.0_dp, dp)
        do first = 1, size(edge_vertices, 2)
            do second = 1, first
                do first_slot = 1, 2
                    first_panel = edge_triangles(first_slot, first)
                    do second_slot = 1, 2
                        second_panel = edge_triangles(second_slot, second)
                        if (first_panel == second_panel) then
                            call &
                                integrate_maxwell_sphere_curved_coincident_rwg_pair_3d( &
                                vertices, triangles, edge_vertices, &
                                edge_triangles, first, first_panel, second, &
                                radius, wave_number, quadrature_degree, &
                                contribution, status)
                        else
                            call &
                                integrate_maxwell_sphere_curved_adjacent_rwg_pair_3d( &
                                vertices, triangles, edge_vertices, &
                                edge_triangles, first, first_panel, second, &
                                second_panel, radius, wave_number, &
                                quadrature_degree, tolerance, max_depth, &
                                contribution, status)
                        end if
                        if (status /= 0) return
                        matrix(first, second) = &
                            matrix(first, second) + contribution
                    end do
                end do
                matrix(second, first) = matrix(first, second)
            end do
        end do
        status = 0
    end subroutine assemble_maxwell_sphere_curved_vector_potential_rwg_3d

    subroutine integrate_maxwell_sphere_curved_adjacent_rwg_pair_3d( &
            vertices, triangles, edge_vertices, edge_triangles, first_basis, &
            first_panel, second_basis, second_panel, radius, wave_number, &
            quadrature_degree, tolerance, max_depth, value, status, &
            scalar_value, decaying_kernel)
        real(dp), intent(in) :: vertices(:, :), radius, wave_number, tolerance
        integer, intent(in) :: triangles(:, :), edge_vertices(:, :)
        integer, intent(in) :: edge_triangles(:, :)
        integer, intent(in) :: first_basis, first_panel, second_basis
        integer, intent(in) :: second_panel, quadrature_degree, max_depth
        complex(dp), intent(out) :: value
        integer, intent(out) :: status
        complex(dp), optional, intent(out) :: scalar_value
        logical, optional, intent(in) :: decaying_kernel

        real(dp), allocatable :: eta(:), weights(:), xi(:)
        real(dp) :: reference_triangle(2, 3)
        complex(dp) :: scalar_integral
        logical :: use_decaying_kernel

        value = cmplx(0.0_dp, 0.0_dp, dp)
        if (present(scalar_value)) scalar_value = cmplx(0.0_dp, 0.0_dp, dp)
        use_decaying_kernel = .false.
        if (present(decaying_kernel)) use_decaying_kernel = decaying_kernel
        status = 1
        if (radius <= 0.0_dp .or. wave_number < 0.0_dp .or. &
            tolerance <= 0.0_dp .or. max_depth < 1) return
        if (first_panel == second_panel) return
        if (.not. any(edge_triangles(:, first_basis) == first_panel) .or. &
            .not. any(edge_triangles(:, second_basis) == second_panel)) return
        call triangle_duffy_quadrature( &
            quadrature_degree, xi, eta, weights, status)
        if (status /= 0) return
        reference_triangle = reshape( &
            [0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, 1.0_dp], [2, 3])
        call integrate_adaptive_curved_rwg_pair( &
            vertices, triangles, edge_vertices, edge_triangles, first_basis, &
            first_panel, second_basis, second_panel, radius, wave_number, &
            reference_triangle, reference_triangle, xi, eta, weights, &
            tolerance, present(scalar_value), use_decaying_kernel, 0, &
            max_depth, value, &
            scalar_integral, status)
        if (present(scalar_value)) scalar_value = scalar_integral
    end subroutine integrate_maxwell_sphere_curved_adjacent_rwg_pair_3d

    recursive subroutine integrate_adaptive_curved_rwg_pair( &
            vertices, triangles, edge_vertices, edge_triangles, first_basis, &
            first_panel, second_basis, second_panel, radius, wave_number, &
            first_reference, second_reference, xi, eta, weights, tolerance, &
            need_scalar, decaying_kernel, depth, max_depth, value, scalar_value, &
            status)
        real(dp), intent(in) :: vertices(:, :), radius, wave_number
        integer, intent(in) :: triangles(:, :), edge_vertices(:, :)
        integer, intent(in) :: edge_triangles(:, :)
        integer, intent(in) :: first_basis, first_panel, second_basis
        integer, intent(in) :: second_panel, depth, max_depth
        real(dp), intent(in) :: first_reference(2, 3)
        real(dp), intent(in) :: second_reference(2, 3)
        real(dp), intent(in) :: xi(:), eta(:), weights(:), tolerance
        logical, intent(in) :: need_scalar
        logical, intent(in) :: decaying_kernel
        complex(dp), intent(out) :: value, scalar_value
        integer, intent(out) :: status

        real(dp) :: first_children(2, 3, 4), second_children(2, 3, 4)
        complex(dp) :: coarse, coarse_scalar, contribution
        complex(dp) :: contribution_scalar, refined, refined_scalar
        integer :: first_child, second_child

        call integrate_regular_curved_rwg_pair( &
            vertices, triangles, edge_vertices, edge_triangles, first_basis, &
            first_panel, second_basis, second_panel, radius, wave_number, &
            first_reference, second_reference, xi, eta, weights, coarse, &
            coarse_scalar, decaying_kernel, status)
        if (status /= 0) return
        call subdivide_reference_triangle(first_reference, first_children)
        call subdivide_reference_triangle(second_reference, second_children)
        refined = cmplx(0.0_dp, 0.0_dp, dp)
        refined_scalar = cmplx(0.0_dp, 0.0_dp, dp)
        do first_child = 1, 4
            do second_child = 1, 4
                call integrate_regular_curved_rwg_pair( &
                    vertices, triangles, edge_vertices, edge_triangles, &
                    first_basis, first_panel, second_basis, second_panel, &
                    radius, wave_number, first_children(:, :, first_child), &
                    second_children(:, :, second_child), xi, eta, weights, &
                    contribution, contribution_scalar, decaying_kernel, status)
                if (status /= 0) return
                refined = refined + contribution
                refined_scalar = refined_scalar + contribution_scalar
            end do
        end do
        if (depth + 1 >= max_depth .or. (abs(refined - coarse) <= &
            tolerance*max(tiny(1.0_dp), abs(refined)) .and. &
            (.not. need_scalar .or. abs(refined_scalar - coarse_scalar) <= &
            tolerance*max(tiny(1.0_dp), abs(refined_scalar))))) then
            value = refined
            scalar_value = refined_scalar
            status = 0
            return
        end if
        value = cmplx(0.0_dp, 0.0_dp, dp)
        scalar_value = cmplx(0.0_dp, 0.0_dp, dp)
        do first_child = 1, 4
            do second_child = 1, 4
                if (curved_reference_children_touch( &
                    vertices, triangles, first_panel, second_panel, radius, &
                    first_children(:, :, first_child), &
                    second_children(:, :, second_child))) then
                    call integrate_adaptive_curved_rwg_pair( &
                        vertices, triangles, edge_vertices, edge_triangles, &
                        first_basis, first_panel, second_basis, second_panel, &
                        radius, wave_number, first_children(:, :, first_child), &
                        second_children(:, :, second_child), xi, eta, weights, &
                        tolerance, need_scalar, decaying_kernel, depth + 1, &
                        max_depth, &
                        contribution, &
                        contribution_scalar, status)
                else
                    call integrate_regular_curved_rwg_pair( &
                        vertices, triangles, edge_vertices, edge_triangles, &
                        first_basis, first_panel, second_basis, second_panel, &
                        radius, wave_number, first_children(:, :, first_child), &
                        second_children(:, :, second_child), xi, eta, weights, &
                        contribution, contribution_scalar, decaying_kernel, &
                        status)
                end if
                if (status /= 0) return
                value = value + contribution
                scalar_value = scalar_value + contribution_scalar
            end do
        end do
        status = 0
    end subroutine integrate_adaptive_curved_rwg_pair

    subroutine integrate_regular_curved_rwg_pair( &
            vertices, triangles, edge_vertices, edge_triangles, first_basis, &
            first_panel, second_basis, second_panel, radius, wave_number, &
            first_reference, second_reference, xi, eta, weights, value, &
            scalar_value, decaying_kernel, status)
        real(dp), intent(in) :: vertices(:, :), radius, wave_number
        integer, intent(in) :: triangles(:, :), edge_vertices(:, :)
        integer, intent(in) :: edge_triangles(:, :)
        integer, intent(in) :: first_basis, first_panel, second_basis
        integer, intent(in) :: second_panel
        real(dp), intent(in) :: first_reference(2, 3)
        real(dp), intent(in) :: second_reference(2, 3)
        real(dp), intent(in) :: xi(:), eta(:), weights(:)
        complex(dp), intent(out) :: value, scalar_value
        logical, intent(in) :: decaying_kernel
        integer, intent(out) :: status

        real(dp) :: first_divergence, first_jacobian, first_point(3)
        real(dp) :: first_reference_jacobian, first_value(3), first_xi_eta(2)
        real(dp) :: physical_distance, second_divergence, second_jacobian
        real(dp) :: second_point(3), second_reference_jacobian
        real(dp) :: second_value(3), second_xi_eta(2)
        complex(dp) :: green
        integer :: first_node, second_node

        value = cmplx(0.0_dp, 0.0_dp, dp)
        scalar_value = cmplx(0.0_dp, 0.0_dp, dp)
        first_reference_jacobian = reference_triangle_jacobian(first_reference)
        second_reference_jacobian = &
            reference_triangle_jacobian(second_reference)
        do first_node = 1, size(weights)
            first_xi_eta = reference_point( &
                first_reference, xi(first_node), eta(first_node))
            call evaluate_maxwell_sphere_curved_rwg_basis( &
                vertices, triangles, edge_vertices, edge_triangles, &
                first_basis, first_panel, radius, first_xi_eta(1), &
                first_xi_eta(2), first_point, first_value, first_divergence, &
                first_jacobian, status)
            if (status /= 0) return
            do second_node = 1, size(weights)
                second_xi_eta = reference_point( &
                    second_reference, xi(second_node), eta(second_node))
                call evaluate_maxwell_sphere_curved_rwg_basis( &
                    vertices, triangles, edge_vertices, edge_triangles, &
                    second_basis, second_panel, radius, second_xi_eta(1), &
                    second_xi_eta(2), second_point, second_value, &
                    second_divergence, second_jacobian, status)
                if (status /= 0) return
                physical_distance = norm2(first_point - second_point)
                green = boundary_green( &
                    wave_number, physical_distance, decaying_kernel)
                value = value + first_reference_jacobian* &
                    second_reference_jacobian*weights(first_node)* &
                    weights(second_node)*first_jacobian*second_jacobian*green* &
                    dot_product(first_value, second_value)
                scalar_value = scalar_value + first_reference_jacobian* &
                    second_reference_jacobian*weights(first_node)* &
                    weights(second_node)*green
            end do
        end do
        status = 0
    end subroutine integrate_regular_curved_rwg_pair

    subroutine integrate_maxwell_sphere_curved_coincident_rwg_pair_3d( &
            vertices, triangles, edge_vertices, edge_triangles, first_basis, &
            panel, second_basis, radius, wave_number, quadrature_degree, &
            value, status, scalar_value, decaying_kernel)
        real(dp), intent(in) :: vertices(:, :), radius, wave_number
        integer, intent(in) :: triangles(:, :), edge_vertices(:, :)
        integer, intent(in) :: edge_triangles(:, :)
        integer, intent(in) :: first_basis, panel, second_basis
        integer, intent(in) :: quadrature_degree
        complex(dp), intent(out) :: value
        integer, intent(out) :: status
        complex(dp), optional, intent(out) :: scalar_value
        logical, optional, intent(in) :: decaying_kernel

        real(dp), allocatable :: eta(:), line_nodes(:), line_weights(:)
        real(dp), allocatable :: weights(:), xi(:)
        real(dp) :: direction(2), first_divergence, first_jacobian
        real(dp) :: first_point(3), first_value(3), reference_vertices(2, 3)
        real(dp) :: physical_distance, rho, second_divergence, second_eta
        real(dp) :: second_jacobian, second_point(3), second_value(3)
        real(dp) :: second_xi, t, wedge_first(2), wedge_jacobian
        real(dp) :: wedge_second(2)
        complex(dp) :: green, scalar_integral
        logical :: use_decaying_kernel
        integer :: first_node, line_count, radial_node, tangential_node, wedge

        value = cmplx(0.0_dp, 0.0_dp, dp)
        scalar_integral = cmplx(0.0_dp, 0.0_dp, dp)
        if (present(scalar_value)) scalar_value = scalar_integral
        use_decaying_kernel = .false.
        if (present(decaying_kernel)) use_decaying_kernel = decaying_kernel
        status = 1
        if (radius <= 0.0_dp .or. wave_number < 0.0_dp) return
        call triangle_duffy_quadrature( &
            quadrature_degree, xi, eta, weights, status)
        if (status /= 0) return
        line_count = max(2, quadrature_degree)
        allocate(line_nodes(line_count), line_weights(line_count))
        call gauss_legendre_ab( &
            line_count, 0.0_dp, 1.0_dp, line_nodes, line_weights)
        reference_vertices = reshape( &
            [0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, 1.0_dp], [2, 3])
        do first_node = 1, size(weights)
            call evaluate_maxwell_sphere_curved_rwg_basis( &
                vertices, triangles, edge_vertices, edge_triangles, &
                first_basis, panel, radius, xi(first_node), eta(first_node), &
                first_point, first_value, first_divergence, first_jacobian, &
                status)
            if (status /= 0) return
            do wedge = 1, 3
                wedge_first(1) = reference_vertices(1, wedge) - xi(first_node)
                wedge_first(2) = reference_vertices(2, wedge) - eta(first_node)
                wedge_second(1) = &
                    reference_vertices(1, modulo(wedge, 3) + 1) - xi(first_node)
                wedge_second(2) = &
                    reference_vertices(2, modulo(wedge, 3) + 1) - eta(first_node)
                wedge_jacobian = abs( &
                    wedge_first(1)*wedge_second(2) - &
                    wedge_first(2)*wedge_second(1))
                do radial_node = 1, line_count
                    rho = line_nodes(radial_node)
                    do tangential_node = 1, line_count
                        t = line_nodes(tangential_node)
                        direction = (1.0_dp - t)*wedge_first + t*wedge_second
                        second_xi = xi(first_node) + rho*direction(1)
                        second_eta = eta(first_node) + rho*direction(2)
                        call evaluate_maxwell_sphere_curved_rwg_basis( &
                            vertices, triangles, edge_vertices, edge_triangles, &
                            second_basis, panel, radius, second_xi, second_eta, &
                            second_point, second_value, second_divergence, &
                            second_jacobian, status)
                        if (status /= 0) return
                        physical_distance = norm2(first_point - second_point)
                        green = boundary_green( &
                            wave_number, physical_distance, use_decaying_kernel)
                        value = value + weights(first_node)*first_jacobian* &
                            line_weights(radial_node)* &
                            line_weights(tangential_node)*rho*wedge_jacobian* &
                            second_jacobian*green* &
                            dot_product(first_value, second_value)
                        scalar_integral = scalar_integral + &
                            weights(first_node)*line_weights(radial_node)* &
                            line_weights(tangential_node)*rho*wedge_jacobian* &
                            green
                    end do
                end do
            end do
        end do
        if (present(scalar_value)) scalar_value = scalar_integral
        status = 0
    end subroutine integrate_maxwell_sphere_curved_coincident_rwg_pair_3d

    subroutine evaluate_maxwell_sphere_curved_far_field_rwg_3d( &
            vertices, triangles, radius, coefficients, direction, wave_number, &
            impedance, quadrature_degree, far_field, status)
        real(dp), intent(in) :: vertices(:, :), radius, direction(3)
        complex(dp), intent(in) :: coefficients(:)
        real(dp), intent(in) :: wave_number, impedance
        integer, intent(in) :: triangles(:, :), quadrature_degree
        complex(dp), intent(out) :: far_field(3)
        integer, intent(out) :: status

        integer, allocatable :: edge_triangles(:, :), edge_vertices(:, :)
        real(dp), allocatable :: eta(:), weights(:), xi(:)
        real(dp) :: basis_value(3), divergence, jacobian, point(3)
        complex(dp) :: phase, surface_current(3), transverse_current(3)
        integer :: basis, node, panel

        far_field = cmplx(0.0_dp, 0.0_dp, dp)
        status = 1
        if (radius <= 0.0_dp .or. wave_number <= 0.0_dp .or. &
            impedance <= 0.0_dp) return
        if (abs(norm2(direction) - 1.0_dp) > 128.0_dp*epsilon(1.0_dp)) return
        call build_maxwell_rwg_surface_space( &
            vertices, triangles, edge_vertices, edge_triangles, status)
        if (status /= 0 .or. size(coefficients) /= size(edge_vertices, 2)) return
        call triangle_duffy_quadrature( &
            quadrature_degree, xi, eta, weights, status)
        if (status /= 0) return
        do panel = 1, size(triangles, 2)
            do node = 1, size(weights)
                surface_current = cmplx(0.0_dp, 0.0_dp, dp)
                do basis = 1, size(edge_vertices, 2)
                    if (.not. any(edge_triangles(:, basis) == panel)) cycle
                    call evaluate_maxwell_sphere_curved_rwg_basis( &
                        vertices, triangles, edge_vertices, edge_triangles, &
                        basis, panel, radius, xi(node), eta(node), point, &
                        basis_value, divergence, jacobian, status)
                    if (status /= 0) return
                    surface_current = surface_current + &
                        coefficients(basis)*basis_value
                end do
                transverse_current = surface_current - &
                    direction*sum(direction*surface_current)
                phase = exp(cmplx( &
                    0.0_dp, -wave_number*dot_product(direction, point), dp))
                far_field = far_field + weights(node)*jacobian*phase* &
                    transverse_current
            end do
        end do
        far_field = cmplx(0.0_dp, wave_number*impedance/ &
            (4.0_dp*acos(-1.0_dp)), dp)*far_field
        status = 0
    end subroutine evaluate_maxwell_sphere_curved_far_field_rwg_3d

    subroutine assemble_maxwell_sphere_curved_plane_wave_rhs_rwg_3d( &
            vertices, triangles, radius, direction, polarization, wave_number, &
            quadrature_degree, right_hand_side, status)
        real(dp), intent(in) :: vertices(:, :), radius, direction(3)
        complex(dp), intent(in) :: polarization(3)
        real(dp), intent(in) :: wave_number
        integer, intent(in) :: triangles(:, :), quadrature_degree
        complex(dp), allocatable, intent(out) :: right_hand_side(:)
        integer, intent(out) :: status

        integer, allocatable :: edge_triangles(:, :), edge_vertices(:, :)
        real(dp), allocatable :: eta(:), weights(:), xi(:)
        real(dp) :: basis_value(3), divergence, jacobian, point(3)
        complex(dp) :: incident_field(3)
        integer :: basis, node, panel

        status = 1
        if (radius <= 0.0_dp .or. wave_number < 0.0_dp) return
        if (abs(norm2(direction) - 1.0_dp) > 128.0_dp*epsilon(1.0_dp)) return
        if (abs(sum(polarization*direction)) > &
            128.0_dp*epsilon(1.0_dp)*max(1.0_dp, &
            sqrt(sum(abs(polarization)**2)))) return
        if (sqrt(sum(abs(polarization)**2)) <= tiny(1.0_dp)) return
        call build_maxwell_rwg_surface_space( &
            vertices, triangles, edge_vertices, edge_triangles, status)
        if (status /= 0 .or. size(edge_vertices, 2) == 0) return
        call triangle_duffy_quadrature( &
            quadrature_degree, xi, eta, weights, status)
        if (status /= 0) return
        allocate(right_hand_side(size(edge_vertices, 2)))
        right_hand_side = cmplx(0.0_dp, 0.0_dp, dp)
        do basis = 1, size(edge_vertices, 2)
            do panel = 1, 2
                do node = 1, size(weights)
                    call evaluate_maxwell_sphere_curved_rwg_basis( &
                        vertices, triangles, edge_vertices, edge_triangles, &
                        basis, edge_triangles(panel, basis), radius, xi(node), &
                        eta(node), point, basis_value, divergence, jacobian, &
                        status)
                    if (status /= 0) return
                    incident_field = polarization*exp(cmplx( &
                        0.0_dp, wave_number*dot_product(direction, point), dp))
                    right_hand_side(basis) = right_hand_side(basis) - &
                        weights(node)*jacobian* &
                        sum(cmplx(basis_value, 0.0_dp, dp)*incident_field)
                end do
            end do
        end do
        status = 0
    end subroutine assemble_maxwell_sphere_curved_plane_wave_rhs_rwg_3d

    subroutine assemble_maxwell_sphere_curved_rwg_mass_matrix( &
            vertices, triangles, radius, quadrature_degree, matrix, status)
        real(dp), intent(in) :: vertices(:, :), radius
        integer, intent(in) :: triangles(:, :), quadrature_degree
        real(dp), allocatable, intent(out) :: matrix(:, :)
        integer, intent(out) :: status

        integer, allocatable :: edge_triangles(:, :), edge_vertices(:, :)
        real(dp), allocatable :: eta(:), weights(:), xi(:)
        real(dp), allocatable :: values(:, :)
        real(dp) :: divergence, jacobian, point(3)
        integer :: basis, node, panel, test_basis

        status = 1
        call build_maxwell_rwg_surface_space( &
            vertices, triangles, edge_vertices, edge_triangles, status)
        if (status /= 0 .or. size(edge_vertices, 2) == 0) return
        call triangle_duffy_quadrature( &
            quadrature_degree, xi, eta, weights, status)
        if (status /= 0) return
        allocate( &
            matrix(size(edge_vertices, 2), size(edge_vertices, 2)), &
            values(3, size(edge_vertices, 2)))
        matrix = 0.0_dp
        do panel = 1, size(triangles, 2)
            do node = 1, size(weights)
                values = 0.0_dp
                do basis = 1, size(edge_vertices, 2)
                    if (.not. any(edge_triangles(:, basis) == panel)) cycle
                    call evaluate_maxwell_sphere_curved_rwg_basis( &
                        vertices, triangles, edge_vertices, edge_triangles, &
                        basis, panel, radius, xi(node), eta(node), point, &
                        values(:, basis), divergence, jacobian, status)
                    if (status /= 0) return
                end do
                do basis = 1, size(edge_vertices, 2)
                    if (.not. any(edge_triangles(:, basis) == panel)) cycle
                    do test_basis = 1, size(edge_vertices, 2)
                        if (.not. any( &
                            edge_triangles(:, test_basis) == panel)) cycle
                        matrix(test_basis, basis) = matrix(test_basis, basis) + &
                            weights(node)*jacobian*dot_product( &
                            values(:, test_basis), values(:, basis))
                    end do
                end do
            end do
        end do
        status = 0
    end subroutine assemble_maxwell_sphere_curved_rwg_mass_matrix

    pure subroutine evaluate_maxwell_sphere_curved_rwg_basis( &
            vertices, triangles, edge_vertices, edge_triangles, basis, panel, &
            radius, xi, eta, point, value, surface_divergence, &
            surface_jacobian, status)
        real(dp), intent(in) :: vertices(:, :), radius, xi, eta
        integer, intent(in) :: triangles(:, :), edge_vertices(:, :)
        integer, intent(in) :: edge_triangles(:, :), basis, panel
        real(dp), intent(out) :: point(3), value(3), surface_divergence
        real(dp), intent(out) :: surface_jacobian
        integer, intent(out) :: status

        real(dp) :: edge_length, opposite_coordinates(2), panel_vertices(3, 3)
        real(dp) :: tangent_eta(3), tangent_xi(3)
        integer :: local, next, opposite, orientation

        point = 0.0_dp
        value = 0.0_dp
        surface_divergence = 0.0_dp
        surface_jacobian = 0.0_dp
        status = 1
        if (basis < 1 .or. basis > size(edge_vertices, 2)) return
        if (panel < 1 .or. panel > size(triangles, 2)) return
        if (.not. any(edge_triangles(:, basis) == panel)) return
        orientation = 0
        opposite = 0
        do local = 1, 3
            next = modulo(local, 3) + 1
            if (triangles(local, panel) == edge_vertices(1, basis) .and. &
                triangles(next, panel) == edge_vertices(2, basis)) then
                orientation = 1
                opposite = modulo(next, 3) + 1
                exit
            end if
            if (triangles(local, panel) == edge_vertices(2, basis) .and. &
                triangles(next, panel) == edge_vertices(1, basis)) then
                orientation = -1
                opposite = modulo(next, 3) + 1
                exit
            end if
        end do
        if (orientation == 0) return
        do local = 1, 3
            panel_vertices(:, local) = vertices(:, triangles(local, panel))
        end do
        call evaluate_sphere_curved_panel( &
            panel_vertices, radius, xi, eta, point, &
            tangent_xi, tangent_eta, surface_jacobian, status)
        if (status /= 0) return
        select case (opposite)
        case (1)
            opposite_coordinates = [0.0_dp, 0.0_dp]
        case (2)
            opposite_coordinates = [1.0_dp, 0.0_dp]
        case (3)
            opposite_coordinates = [0.0_dp, 1.0_dp]
        end select
        edge_length = norm2( &
            vertices(:, edge_vertices(2, basis)) - &
            vertices(:, edge_vertices(1, basis)))
        value = real(orientation, dp)*edge_length/surface_jacobian*( &
            (xi - opposite_coordinates(1))*tangent_xi + &
            (eta - opposite_coordinates(2))*tangent_eta)
        surface_divergence = &
            2.0_dp*real(orientation, dp)*edge_length/surface_jacobian
        status = 0
    end subroutine evaluate_maxwell_sphere_curved_rwg_basis

    pure subroutine evaluate_maxwell_sphere_curved_rwg_basis_jvp( &
            vertices, triangles, edge_vertices, edge_triangles, basis, panel, &
            radius, xi, eta, vertices_dot, radius_dot, point, value, &
            surface_divergence, surface_jacobian, point_dot, value_dot, &
            surface_divergence_dot, surface_jacobian_dot, status, xi_dot, &
            eta_dot)
        real(dp), intent(in) :: vertices(:, :), radius, xi, eta
        integer, intent(in) :: triangles(:, :), edge_vertices(:, :)
        integer, intent(in) :: edge_triangles(:, :), basis, panel
        real(dp), intent(in) :: vertices_dot(:, :), radius_dot
        real(dp), intent(out) :: point(3), value(3), surface_divergence
        real(dp), intent(out) :: surface_jacobian, point_dot(3), value_dot(3)
        real(dp), intent(out) :: surface_divergence_dot, surface_jacobian_dot
        integer, intent(out) :: status
        real(dp), optional, intent(in) :: xi_dot, eta_dot

        real(dp) :: edge(3), edge_dot(3), edge_length, edge_length_dot
        real(dp) :: panel_vertices(3, 3), panel_vertices_dot(3, 3)
        real(dp) :: tangent_eta(3), tangent_eta_dot(3)
        real(dp) :: tangent_xi(3), tangent_xi_dot(3)
        real(dp) :: vector(3), vector_dot(3), coefficient, coefficient_dot
        real(dp) :: opposite_coordinates(2), local_xi_dot, local_eta_dot
        integer :: local, next, opposite, orientation

        point = 0.0_dp
        value = 0.0_dp
        surface_divergence = 0.0_dp
        surface_jacobian = 0.0_dp
        point_dot = 0.0_dp
        value_dot = 0.0_dp
        surface_divergence_dot = 0.0_dp
        surface_jacobian_dot = 0.0_dp
        status = 1
        call evaluate_maxwell_sphere_curved_rwg_basis( &
            vertices, triangles, edge_vertices, edge_triangles, basis, panel, &
            radius, xi, eta, point, value, surface_divergence, &
            surface_jacobian, status)
        if (status /= 0) return
        if (any(shape(vertices_dot) /= shape(vertices))) return
        local_xi_dot = 0.0_dp
        local_eta_dot = 0.0_dp
        if (present(xi_dot)) local_xi_dot = xi_dot
        if (present(eta_dot)) local_eta_dot = eta_dot
        orientation = 0
        opposite = 0
        do local = 1, 3
            next = modulo(local, 3) + 1
            if (triangles(local, panel) == edge_vertices(1, basis) .and. &
                triangles(next, panel) == edge_vertices(2, basis)) then
                orientation = 1
                opposite = modulo(next, 3) + 1
                exit
            end if
            if (triangles(local, panel) == edge_vertices(2, basis) .and. &
                triangles(next, panel) == edge_vertices(1, basis)) then
                orientation = -1
                opposite = modulo(next, 3) + 1
                exit
            end if
        end do
        if (orientation == 0) return
        do local = 1, 3
            panel_vertices(:, local) = vertices(:, triangles(local, panel))
            panel_vertices_dot(:, local) = &
                vertices_dot(:, triangles(local, panel))
        end do
        call evaluate_sphere_curved_panel( &
            panel_vertices, radius, xi, eta, point, tangent_xi, tangent_eta, &
            surface_jacobian, status)
        if (status /= 0) return
        call evaluate_sphere_curved_panel_jvp( &
            panel_vertices, radius, xi, eta, panel_vertices_dot, radius_dot, &
            local_xi_dot, local_eta_dot, point_dot, tangent_xi_dot, &
            tangent_eta_dot, surface_jacobian_dot, status)
        if (status /= 0) return
        select case (opposite)
        case (1)
            opposite_coordinates = [0.0_dp, 0.0_dp]
        case (2)
            opposite_coordinates = [1.0_dp, 0.0_dp]
        case (3)
            opposite_coordinates = [0.0_dp, 1.0_dp]
        end select
        edge = vertices(:, edge_vertices(2, basis)) - &
            vertices(:, edge_vertices(1, basis))
        edge_dot = vertices_dot(:, edge_vertices(2, basis)) - &
            vertices_dot(:, edge_vertices(1, basis))
        edge_length = norm2(edge)
        if (edge_length <= tiny(1.0_dp)) return
        edge_length_dot = dot_product(edge, edge_dot)/edge_length
        vector = (xi - opposite_coordinates(1))*tangent_xi + &
            (eta - opposite_coordinates(2))*tangent_eta
        vector_dot = local_xi_dot*tangent_xi + local_eta_dot*tangent_eta + &
            (xi - opposite_coordinates(1))*tangent_xi_dot + &
            (eta - opposite_coordinates(2))*tangent_eta_dot
        coefficient = real(orientation, dp)*edge_length/surface_jacobian
        coefficient_dot = real(orientation, dp)*( &
            edge_length_dot/surface_jacobian - &
            edge_length*surface_jacobian_dot/surface_jacobian**2)
        value = coefficient*vector
        value_dot = coefficient_dot*vector + coefficient*vector_dot
        surface_divergence = 2.0_dp*coefficient
        surface_divergence_dot = 2.0_dp*coefficient_dot
        status = 0
    end subroutine evaluate_maxwell_sphere_curved_rwg_basis_jvp

    pure subroutine evaluate_maxwell_sphere_curved_rwg_basis_vjp( &
            vertices, triangles, edge_vertices, edge_triangles, basis, panel, &
            radius, xi, eta, point_bar, value_bar, surface_divergence_bar, &
            surface_jacobian_bar, edge_vertices_bar, panel_vertices_bar, &
            radius_bar, status, xi_bar, eta_bar)
        real(dp), intent(in) :: vertices(:, :), radius, xi, eta
        integer, intent(in) :: triangles(:, :), edge_vertices(:, :)
        integer, intent(in) :: edge_triangles(:, :), basis, panel
        real(dp), intent(in) :: point_bar(3), value_bar(3)
        real(dp), intent(in) :: surface_divergence_bar, surface_jacobian_bar
        real(dp), intent(out) :: edge_vertices_bar(3, 2), panel_vertices_bar(3, 3)
        real(dp), intent(out) :: radius_bar
        integer, intent(out) :: status
        real(dp), optional, intent(out) :: xi_bar, eta_bar

        real(dp) :: edge(3), edge_bar(3), edge_length, edge_length_bar
        real(dp) :: panel_vertices(3, 3), tangent_eta(3), tangent_eta_bar(3)
        real(dp) :: tangent_xi(3), tangent_xi_bar(3)
        real(dp) :: vector(3), vector_bar(3), coefficient
        real(dp) :: coefficient_bar, jacobian_bar
        real(dp) :: opposite_coordinates(2), surface_jacobian
        real(dp) :: dummy_point(3), dummy_value(3), dummy_divergence
        real(dp) :: xi_local_bar, eta_local_bar
        real(dp) :: panel_xi_bar, panel_eta_bar
        integer :: local, next, opposite, orientation

        edge_vertices_bar = 0.0_dp
        panel_vertices_bar = 0.0_dp
        radius_bar = 0.0_dp
        if (present(xi_bar)) xi_bar = 0.0_dp
        if (present(eta_bar)) eta_bar = 0.0_dp
        status = 1
        call evaluate_maxwell_sphere_curved_rwg_basis( &
            vertices, triangles, edge_vertices, edge_triangles, basis, panel, &
            radius, xi, eta, dummy_point, dummy_value, dummy_divergence, &
            surface_jacobian, status)
        if (status /= 0) return
        orientation = 0
        opposite = 0
        do local = 1, 3
            next = modulo(local, 3) + 1
            if (triangles(local, panel) == edge_vertices(1, basis) .and. &
                triangles(next, panel) == edge_vertices(2, basis)) then
                orientation = 1
                opposite = modulo(next, 3) + 1
                exit
            end if
            if (triangles(local, panel) == edge_vertices(2, basis) .and. &
                triangles(next, panel) == edge_vertices(1, basis)) then
                orientation = -1
                opposite = modulo(next, 3) + 1
                exit
            end if
        end do
        if (orientation == 0) return
        do local = 1, 3
            panel_vertices(:, local) = vertices(:, triangles(local, panel))
        end do
        call evaluate_sphere_curved_panel( &
            panel_vertices, radius, xi, eta, dummy_point, tangent_xi, &
            tangent_eta, surface_jacobian, status)
        if (status /= 0) return
        select case (opposite)
        case (1)
            opposite_coordinates = [0.0_dp, 0.0_dp]
        case (2)
            opposite_coordinates = [1.0_dp, 0.0_dp]
        case (3)
            opposite_coordinates = [0.0_dp, 1.0_dp]
        end select
        edge = vertices(:, edge_vertices(2, basis)) - &
            vertices(:, edge_vertices(1, basis))
        edge_length = norm2(edge)
        if (edge_length <= tiny(1.0_dp)) return
        coefficient = real(orientation, dp)*edge_length/surface_jacobian
        vector = (xi - opposite_coordinates(1))*tangent_xi + &
            (eta - opposite_coordinates(2))*tangent_eta
        coefficient_bar = dot_product(value_bar, vector)
        vector_bar = coefficient*value_bar
        edge_length_bar = real(orientation, dp)*( &
            coefficient_bar/surface_jacobian + &
            2.0_dp*surface_divergence_bar/surface_jacobian)
        jacobian_bar = surface_jacobian_bar - &
            real(orientation, dp)*edge_length*coefficient_bar/ &
            surface_jacobian**2 - &
            2.0_dp*real(orientation, dp)*edge_length* &
            surface_divergence_bar/surface_jacobian**2
        edge_bar = edge_length_bar*edge/edge_length
        edge_vertices_bar(:, 1) = -edge_bar
        edge_vertices_bar(:, 2) = edge_bar
        panel_vertices_bar = 0.0_dp
        tangent_xi_bar = (xi - opposite_coordinates(1))*vector_bar
        tangent_eta_bar = (eta - opposite_coordinates(2))*vector_bar
        xi_local_bar = dot_product(vector_bar, tangent_xi)
        eta_local_bar = dot_product(vector_bar, tangent_eta)
        call evaluate_sphere_curved_panel_vjp( &
            panel_vertices, radius, xi, eta, point_bar, tangent_xi_bar, &
            tangent_eta_bar, jacobian_bar, panel_vertices_bar, radius_bar, &
            panel_xi_bar, panel_eta_bar, status)
        xi_local_bar = xi_local_bar + panel_xi_bar
        eta_local_bar = eta_local_bar + panel_eta_bar
        if (present(xi_bar)) xi_bar = xi_local_bar
        if (present(eta_bar)) eta_bar = eta_local_bar
    end subroutine evaluate_maxwell_sphere_curved_rwg_basis_vjp

    pure subroutine map_refined_sphere_point_to_parent( &
            parent_vertices, radius, point, parent_xi, parent_eta, status)
        real(dp), intent(in) :: parent_vertices(3, 3), radius, point(3)
        real(dp), intent(out) :: parent_xi, parent_eta
        integer, intent(out) :: status

        real(dp) :: lambda

        parent_xi = 0.0_dp
        parent_eta = 0.0_dp
        call invert_sphere_curved_panel( &
            parent_vertices, radius, point, parent_xi, parent_eta, status)
        if (status /= 0) return
        call sphere_parent_lambda( &
            parent_vertices, radius, point, parent_xi, parent_eta, lambda, status)
    end subroutine map_refined_sphere_point_to_parent

    pure subroutine map_refined_sphere_point_to_parent_jvp( &
            parent_vertices, radius, point, parent_vertices_dot, radius_dot, &
            point_dot, parent_xi, parent_eta, parent_xi_dot, parent_eta_dot, &
            status)
        real(dp), intent(in) :: parent_vertices(3, 3), radius, point(3)
        real(dp), intent(in) :: parent_vertices_dot(3, 3), radius_dot
        real(dp), intent(in) :: point_dot(3)
        real(dp), intent(out) :: parent_xi, parent_eta, parent_xi_dot
        real(dp), intent(out) :: parent_eta_dot
        integer, intent(out) :: status

        real(dp) :: inverse_matrix(3, 3), matrix(3, 3), matrix_dot(3, 3)
        real(dp) :: right_hand_side(3), right_hand_side_dot(3), solution(3)
        real(dp) :: solution_dot(3), direction(3), direction_dot(3)
        integer :: inverse_status

        parent_xi = 0.0_dp
        parent_eta = 0.0_dp
        parent_xi_dot = 0.0_dp
        parent_eta_dot = 0.0_dp
        status = 1
        if (radius <= 0.0_dp) return
        call sphere_parent_system( &
            parent_vertices, radius, point, matrix, right_hand_side, status)
        if (status /= 0) return
        call inv3(matrix, inverse_matrix, inverse_status)
        if (inverse_status /= 0) return
        solution = matmul(inverse_matrix, right_hand_side)
        if (solution(3) <= 0.0_dp) return
        parent_xi = solution(1)
        parent_eta = solution(2)
        direction = point/radius
        direction_dot = point_dot/radius - point*radius_dot/radius**2
        matrix_dot(:, 1) = parent_vertices_dot(:, 2) - parent_vertices_dot(:, 1)
        matrix_dot(:, 2) = parent_vertices_dot(:, 3) - parent_vertices_dot(:, 1)
        matrix_dot(:, 3) = -direction_dot
        right_hand_side_dot = -parent_vertices_dot(:, 1)
        solution_dot = matmul(inverse_matrix, &
            right_hand_side_dot - matmul(matrix_dot, solution))
        parent_xi_dot = solution_dot(1)
        parent_eta_dot = solution_dot(2)
        if (parent_xi < -512.0_dp*epsilon(1.0_dp) .or. &
            parent_eta < -512.0_dp*epsilon(1.0_dp) .or. &
            parent_xi + parent_eta > 1.0_dp + &
            512.0_dp*epsilon(1.0_dp)) return
        parent_xi = max(0.0_dp, parent_xi)
        parent_eta = max(0.0_dp, parent_eta)
        status = 0
    end subroutine map_refined_sphere_point_to_parent_jvp

    pure subroutine map_refined_sphere_point_to_parent_vjp( &
            parent_vertices, radius, point, parent_xi_bar, parent_eta_bar, &
            point_bar, parent_vertices_bar, radius_bar, status)
        real(dp), intent(in) :: parent_vertices(3, 3), radius, point(3)
        real(dp), intent(in) :: parent_xi_bar, parent_eta_bar
        real(dp), intent(out) :: point_bar(3), parent_vertices_bar(3, 3)
        real(dp), intent(out) :: radius_bar
        integer, intent(out) :: status

        real(dp) :: inverse_matrix(3, 3), matrix(3, 3), right_hand_side(3)
        real(dp) :: solution(3), solution_bar(3), right_hand_side_bar(3)
        real(dp) :: matrix_bar(3, 3), edge1_bar(3), edge2_bar(3)
        real(dp) :: direction_bar(3)
        integer :: inverse_status

        point_bar = 0.0_dp
        parent_vertices_bar = 0.0_dp
        radius_bar = 0.0_dp
        status = 1
        if (radius <= 0.0_dp) return
        call sphere_parent_system( &
            parent_vertices, radius, point, matrix, right_hand_side, status)
        if (status /= 0) return
        call inv3(matrix, inverse_matrix, inverse_status)
        if (inverse_status /= 0) then
            status = 1
            return
        end if
        solution = matmul(inverse_matrix, right_hand_side)
        if (solution(3) <= 0.0_dp) return
        solution_bar = [parent_xi_bar, parent_eta_bar, 0.0_dp]
        right_hand_side_bar = matmul(transpose(inverse_matrix), solution_bar)
        matrix_bar = -spread(right_hand_side_bar, 2, 3)* &
            spread(solution, 1, 3)
        parent_vertices_bar(:, 1) = -right_hand_side_bar
        edge1_bar = matrix_bar(:, 1)
        edge2_bar = matrix_bar(:, 2)
        parent_vertices_bar(:, 1) = parent_vertices_bar(:, 1) - &
            edge1_bar - edge2_bar
        parent_vertices_bar(:, 2) = edge1_bar
        parent_vertices_bar(:, 3) = edge2_bar
        direction_bar = -matrix_bar(:, 3)
        point_bar = direction_bar/radius
        radius_bar = -dot_product(direction_bar, point)/radius**2
        if (parent_xi_bar == 0.0_dp .and. parent_eta_bar == 0.0_dp) then
            status = 0
            return
        end if
        status = 0
    end subroutine map_refined_sphere_point_to_parent_vjp

    pure subroutine sphere_parent_system( &
            parent_vertices, radius, point, matrix, right_hand_side, status)
        real(dp), intent(in) :: parent_vertices(3, 3), radius, point(3)
        real(dp), intent(out) :: matrix(3, 3), right_hand_side(3)
        integer, intent(out) :: status

        real(dp) :: direction(3)

        matrix = 0.0_dp
        right_hand_side = 0.0_dp
        status = 1
        if (radius <= 0.0_dp) return
        if (abs(norm2(point) - radius) > &
            1024.0_dp*epsilon(1.0_dp)*max(1.0_dp, radius)) return
        direction = point/radius
        matrix(:, 1) = parent_vertices(:, 2) - parent_vertices(:, 1)
        matrix(:, 2) = parent_vertices(:, 3) - parent_vertices(:, 1)
        matrix(:, 3) = -direction
        right_hand_side = -parent_vertices(:, 1)
        status = 0
    end subroutine sphere_parent_system

    pure subroutine sphere_parent_lambda( &
            parent_vertices, radius, point, parent_xi, parent_eta, lambda, status)
        real(dp), intent(in) :: parent_vertices(3, 3), radius, point(3)
        real(dp), intent(in) :: parent_xi, parent_eta
        real(dp), intent(out) :: lambda
        integer, intent(out) :: status

        real(dp) :: affine_point(3), direction(3)

        lambda = 0.0_dp
        status = 1
        if (radius <= 0.0_dp) return
        direction = point/radius
        affine_point = parent_vertices(:, 1) + &
            parent_xi*(parent_vertices(:, 2) - parent_vertices(:, 1)) + &
            parent_eta*(parent_vertices(:, 3) - parent_vertices(:, 1))
        lambda = dot_product(affine_point, direction)
        if (lambda <= 0.0_dp) then
            status = 1
            return
        end if
        status = 0
    end subroutine sphere_parent_lambda

    pure function helmholtz_green(wave_number, radius) result(value)
        real(dp), intent(in) :: wave_number, radius
        complex(dp) :: value

        value = exp(cmplx(0.0_dp, wave_number*radius, dp))/ &
            (4.0_dp*acos(-1.0_dp)*radius)
    end function helmholtz_green

    pure function boundary_green( &
            wave_number, radius, decaying_kernel) result(value)
        real(dp), intent(in) :: wave_number, radius
        logical, intent(in) :: decaying_kernel
        complex(dp) :: value

        if (decaying_kernel) then
            value = cmplx(exp(-wave_number*radius)/ &
                (4.0_dp*acos(-1.0_dp)*radius), 0.0_dp, dp)
        else
            value = helmholtz_green(wave_number, radius)
        end if
    end function boundary_green

    pure function complex_cross_product(first, second) result(product)
        complex(dp), intent(in) :: first(3), second(3)
        complex(dp) :: product(3)

        product(1) = first(2)*second(3) - first(3)*second(2)
        product(2) = first(3)*second(1) - first(1)*second(3)
        product(3) = first(1)*second(2) - first(2)*second(1)
    end function complex_cross_product

    pure function real_cross_product(first, second) result(product)
        real(dp), intent(in) :: first(3), second(3)
        real(dp) :: product(3)

        product(1) = first(2)*second(3) - first(3)*second(2)
        product(2) = first(3)*second(1) - first(1)*second(3)
        product(3) = first(1)*second(2) - first(2)*second(1)
    end function real_cross_product

    pure logical function curved_reference_children_touch( &
            vertices, triangles, first_panel, second_panel, radius, &
            first_reference, second_reference) result(touch)
        real(dp), intent(in) :: vertices(:, :), radius
        integer, intent(in) :: triangles(:, :), first_panel, second_panel
        real(dp), intent(in) :: first_reference(2, 3)
        real(dp), intent(in) :: second_reference(2, 3)

        real(dp) :: first_panel_vertices(3, 3), first_point(3), jacobian
        real(dp) :: scale, second_panel_vertices(3, 3), second_point(3)
        real(dp) :: tangent_eta(3), tangent_xi(3)
        integer :: first_vertex, second_vertex, status

        scale = max(1.0_dp, radius)
        touch = .false.
        do first_vertex = 1, 3
            first_panel_vertices(:, first_vertex) = &
                vertices(:, triangles(first_vertex, first_panel))
            second_panel_vertices(:, first_vertex) = &
                vertices(:, triangles(first_vertex, second_panel))
        end do
        do first_vertex = 1, 3
            call evaluate_sphere_curved_panel( &
                first_panel_vertices, radius, &
                first_reference(1, first_vertex), &
                first_reference(2, first_vertex), first_point, tangent_xi, &
                tangent_eta, jacobian, status)
            if (status /= 0) return
            do second_vertex = 1, 3
                call evaluate_sphere_curved_panel( &
                    second_panel_vertices, radius, &
                    second_reference(1, second_vertex), &
                    second_reference(2, second_vertex), second_point, &
                    tangent_xi, tangent_eta, jacobian, status)
                if (status /= 0) return
                if (norm2(first_point - second_point) <= &
                    512.0_dp*epsilon(1.0_dp)*scale) then
                    touch = .true.
                    return
                end if
            end do
        end do
    end function curved_reference_children_touch

    pure function reference_point(vertices, xi, eta) result(point)
        real(dp), intent(in) :: vertices(2, 3), xi, eta
        real(dp) :: point(2)

        point = vertices(:, 1) + xi*(vertices(:, 2) - vertices(:, 1)) + &
            eta*(vertices(:, 3) - vertices(:, 1))
    end function reference_point

    pure function reference_triangle_jacobian(vertices) result(jacobian)
        real(dp), intent(in) :: vertices(2, 3)
        real(dp) :: jacobian

        jacobian = abs( &
            (vertices(1, 2) - vertices(1, 1))* &
            (vertices(2, 3) - vertices(2, 1)) - &
            (vertices(2, 2) - vertices(2, 1))* &
            (vertices(1, 3) - vertices(1, 1)))
    end function reference_triangle_jacobian

    pure subroutine subdivide_reference_triangle(vertices, children)
        real(dp), intent(in) :: vertices(2, 3)
        real(dp), intent(out) :: children(2, 3, 4)

        real(dp) :: midpoint_12(2), midpoint_23(2), midpoint_31(2)

        midpoint_12 = 0.5_dp*(vertices(:, 1) + vertices(:, 2))
        midpoint_23 = 0.5_dp*(vertices(:, 2) + vertices(:, 3))
        midpoint_31 = 0.5_dp*(vertices(:, 3) + vertices(:, 1))
        children(:, 1, 1) = vertices(:, 1)
        children(:, 2, 1) = midpoint_12
        children(:, 3, 1) = midpoint_31
        children(:, 1, 2) = midpoint_12
        children(:, 2, 2) = vertices(:, 2)
        children(:, 3, 2) = midpoint_23
        children(:, 1, 3) = midpoint_31
        children(:, 2, 3) = midpoint_23
        children(:, 3, 3) = vertices(:, 3)
        children(:, 1, 4) = midpoint_12
        children(:, 2, 4) = midpoint_23
        children(:, 3, 4) = midpoint_31
    end subroutine subdivide_reference_triangle

end module fortfem_maxwell_sphere_curved_rwg
