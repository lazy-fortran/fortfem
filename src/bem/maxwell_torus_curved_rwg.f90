module fortfem_maxwell_torus_curved_rwg
    !! Surface-Piola RWG basis on exact parametric torus panels.
    use fortfem_kinds, only: dp
    use fortfem_barycentric_surface_refinement, only: &
        barycentric_refine_torus_surface_mesh_jvp, &
        barycentric_refine_torus_surface_mesh_vjp
    use fortfem_maxwell_bc_surface, only: &
        build_maxwell_bc_transformation, &
        differentiate_maxwell_bc_transformation_jvp, &
        differentiate_maxwell_bc_transformation_vjp
    use fortfem_maxwell_efie_bc_3d, only: build_maxwell_bc_to_refined_rwg
    use fortfem_maxwell_rwg_surface, only: build_maxwell_rwg_surface_space
    use fortfem_torus_curved_panel, only: &
        evaluate_torus_curved_panel, evaluate_torus_curved_panel_jvp, &
        evaluate_torus_curved_panel_vjp
    use fortfem_triangle_duffy_quadrature, only: triangle_duffy_quadrature
    use fortnum_linalg, only: dense_solve
    use fortnum_quadrature, only: gauss_legendre_ab
    implicit none
    private

    public :: assemble_maxwell_torus_curved_rwg_mass_matrix
    public :: assemble_maxwell_torus_curved_rwg_mass_matrix_jvp
    public :: assemble_maxwell_torus_curved_rwg_mass_matrix_vjp
    public :: assemble_maxwell_torus_curved_rwg_rbc_pairing
    public :: assemble_maxwell_torus_curved_rwg_rbc_pairing_jvp
    public :: assemble_maxwell_torus_curved_rwg_rbc_pairing_vjp
    public :: assemble_maxwell_torus_curved_plane_wave_rhs_rwg_3d
    public :: assemble_maxwell_torus_curved_plane_wave_rhs_rwg_3d_jvp
    public :: assemble_maxwell_torus_curved_plane_wave_rhs_rwg_3d_vjp
    public :: assemble_maxwell_torus_curved_efie_rwg_3d
    public :: assemble_maxwell_torus_curved_efie_imaginary_rwg_3d
    public :: assemble_maxwell_torus_curved_efie_bc_imaginary_3d
    public :: assemble_maxwell_torus_curved_mfie_offset_trace_rwg_rbc_3d
    public :: assemble_maxwell_torus_curved_mfie_rwg_rbc_3d
    public :: assemble_maxwell_torus_curved_potential_operators_rwg_3d
    public :: assemble_maxwell_torus_curved_regularized_cfie_rwg_3d
    public :: assemble_maxwell_torus_curved_plane_wave_rhs_bc_3d
    public :: evaluate_maxwell_torus_curved_far_field_rwg_3d
    public :: evaluate_maxwell_torus_curved_far_field_rwg_3d_jvp
    public :: evaluate_maxwell_torus_curved_far_field_rwg_3d_vjp
    public :: evaluate_maxwell_torus_curved_magnetic_field_rwg_3d
    public :: evaluate_maxwell_torus_curved_magnetic_field_rwg_3d_jvp
    public :: evaluate_maxwell_torus_magnetic_geometry_jvp
    public :: evaluate_maxwell_torus_magnetic_geometry_vjp
    public :: evaluate_maxwell_torus_curved_magnetic_field_rwg_3d_vjp
    public :: evaluate_maxwell_torus_curved_rwg_basis_jvp
    public :: evaluate_maxwell_torus_curved_rwg_basis_vjp
    public :: evaluate_maxwell_torus_curved_localized_rwg_basis
    public :: evaluate_maxwell_torus_curved_localized_rwg_basis_jvp
    public :: evaluate_maxwell_torus_curved_localized_rwg_basis_vjp
    public :: evaluate_maxwell_torus_curved_rwg_basis
    public :: integrate_maxwell_torus_curved_adjacent_rwg_pair_3d
    public :: integrate_maxwell_torus_curved_coincident_rwg_pair_3d
    public :: solve_maxwell_pec_torus_curved_regularized_cfie_rwg_3d
    public :: solve_maxwell_pec_torus_curved_regularized_cfie_rwg_multiple_3d

contains

    subroutine evaluate_maxwell_torus_curved_magnetic_field_rwg_3d( &
            vertices, triangles, parameters, major_radius, minor_radius, &
            coefficients, observation, wave_number, quadrature_degree, &
            magnetic_field, status)
        real(dp), intent(in) :: vertices(:, :), parameters(:, :)
        real(dp), intent(in) :: major_radius, minor_radius, observation(3)
        real(dp), intent(in) :: wave_number
        integer, intent(in) :: triangles(:, :), quadrature_degree
        complex(dp), intent(in) :: coefficients(:)
        complex(dp), intent(out) :: magnetic_field(3)
        integer, intent(out) :: status

        complex(dp), allocatable :: basis_fields(:, :)

        call evaluate_torus_curved_magnetic_basis( &
            vertices, triangles, parameters, major_radius, minor_radius, &
            observation, wave_number, quadrature_degree, basis_fields, status)
        if (status /= 0) return
        if (size(coefficients) /= size(basis_fields, 2)) then
            status = 1
            magnetic_field = cmplx(0.0_dp, 0.0_dp, dp)
            return
        end if
        magnetic_field = matmul(basis_fields, coefficients)
        status = 0
    end subroutine evaluate_maxwell_torus_curved_magnetic_field_rwg_3d

    subroutine evaluate_maxwell_torus_curved_magnetic_field_rwg_3d_jvp( &
            vertices, triangles, parameters, major_radius, minor_radius, &
            coefficients_dot, observation, wave_number, quadrature_degree, &
            magnetic_field_dot, status)
        real(dp), intent(in) :: vertices(:, :), parameters(:, :)
        real(dp), intent(in) :: major_radius, minor_radius, observation(3)
        real(dp), intent(in) :: wave_number
        integer, intent(in) :: triangles(:, :), quadrature_degree
        complex(dp), intent(in) :: coefficients_dot(:)
        complex(dp), intent(out) :: magnetic_field_dot(3)
        integer, intent(out) :: status

        call evaluate_maxwell_torus_curved_magnetic_field_rwg_3d( &
            vertices, triangles, parameters, major_radius, minor_radius, &
            coefficients_dot, observation, wave_number, quadrature_degree, &
            magnetic_field_dot, status)
    end subroutine &
        evaluate_maxwell_torus_curved_magnetic_field_rwg_3d_jvp

    subroutine evaluate_maxwell_torus_magnetic_geometry_jvp( &
            vertices, triangles, parameters, major_radius, minor_radius, &
            coefficients, observation, wave_number, quadrature_degree, &
            vertices_dot, parameters_dot, major_radius_dot, minor_radius_dot, &
            coefficients_dot, observation_dot, wave_number_dot, &
            magnetic_field_dot, status)
        real(dp), intent(in) :: vertices(:, :), parameters(:, :)
        real(dp), intent(in) :: major_radius, minor_radius, observation(3)
        real(dp), intent(in) :: wave_number
        integer, intent(in) :: triangles(:, :), quadrature_degree
        complex(dp), intent(in) :: coefficients(:), coefficients_dot(:)
        real(dp), intent(in) :: vertices_dot(:, :), parameters_dot(:, :)
        real(dp), intent(in) :: major_radius_dot, minor_radius_dot
        real(dp), intent(in) :: observation_dot(3), wave_number_dot
        complex(dp), intent(out) :: magnetic_field_dot(3)
        complex(dp), allocatable :: basis_fields(:, :), basis_fields_dot(:, :)
        integer, allocatable :: edge_triangles(:, :), edge_vertices(:, :)
        real(dp), allocatable :: eta(:), weights(:), xi(:)
        integer :: status

        magnetic_field_dot = cmplx(0.0_dp, 0.0_dp, dp)
        status = 1
        if (size(vertices, 1) /= 3 .or. size(parameters, 1) /= 2) return
        if (size(parameters, 2) /= size(vertices, 2)) return
        if (any(shape(vertices_dot) /= shape(vertices))) return
        if (any(shape(parameters_dot) /= shape(parameters))) return
        if (size(coefficients_dot) /= size(coefficients)) return
        if (major_radius <= minor_radius .or. minor_radius <= 0.0_dp) return
        if (wave_number < 0.0_dp .or. quadrature_degree < 0) return
        call build_maxwell_rwg_surface_space( &
            vertices, triangles, edge_vertices, edge_triangles, status)
        if (status /= 0) return
        if (size(edge_vertices, 2) == 0) return
        if (size(coefficients) /= size(edge_vertices, 2)) return
        call triangle_duffy_quadrature( &
            quadrature_degree, xi, eta, weights, status)
        if (status /= 0) return
        allocate( &
            basis_fields(3, size(edge_vertices, 2)), &
            basis_fields_dot(3, size(edge_vertices, 2)))
        call evaluate_all_torus_curved_rwg_magnetic_fields_jvp( &
            vertices, triangles, parameters, edge_vertices, edge_triangles, &
            major_radius, minor_radius, observation, wave_number, xi, eta, &
            weights, vertices_dot, parameters_dot, major_radius_dot, &
            minor_radius_dot, observation_dot, wave_number_dot, basis_fields, &
            basis_fields_dot, status)
        if (status /= 0) return
        magnetic_field_dot = matmul(basis_fields_dot, coefficients) + &
            matmul(basis_fields, coefficients_dot)
        status = 0
    end subroutine &
        evaluate_maxwell_torus_magnetic_geometry_jvp

    subroutine evaluate_maxwell_torus_magnetic_geometry_vjp( &
            vertices, triangles, parameters, major_radius, minor_radius, &
            coefficients, observation, wave_number, quadrature_degree, &
            magnetic_field_bar, vertices_bar, parameters_bar, &
            major_radius_bar, minor_radius_bar, coefficients_bar, &
            observation_bar, wave_number_bar, status)
        real(dp), intent(in) :: vertices(:, :), parameters(:, :)
        real(dp), intent(in) :: major_radius, minor_radius, observation(3)
        real(dp), intent(in) :: wave_number
        integer, intent(in) :: triangles(:, :), quadrature_degree
        complex(dp), intent(in) :: coefficients(:), magnetic_field_bar(3)
        real(dp), intent(out) :: vertices_bar(:, :), parameters_bar(:, :)
        real(dp), intent(out) :: major_radius_bar, minor_radius_bar
        complex(dp), allocatable, intent(out) :: coefficients_bar(:)
        real(dp), intent(out) :: observation_bar(3), wave_number_bar
        complex(dp), allocatable :: magnetic_fields(:, :)
        integer, allocatable :: edge_triangles(:, :), edge_vertices(:, :)
        real(dp), allocatable :: eta(:), weights(:), xi(:)
        integer :: status

        vertices_bar = 0.0_dp
        parameters_bar = 0.0_dp
        major_radius_bar = 0.0_dp
        minor_radius_bar = 0.0_dp
        observation_bar = 0.0_dp
        wave_number_bar = 0.0_dp
        if (allocated(coefficients_bar)) deallocate(coefficients_bar)
        status = 1
        if (size(vertices, 1) /= 3 .or. size(parameters, 1) /= 2) return
        if (size(parameters, 2) /= size(vertices, 2)) return
        if (size(vertices_bar, 1) /= size(vertices, 1)) return
        if (size(vertices_bar, 2) /= size(vertices, 2)) return
        if (size(parameters_bar, 1) /= size(parameters, 1)) return
        if (size(parameters_bar, 2) /= size(parameters, 2)) return
        if (major_radius <= minor_radius .or. minor_radius <= 0.0_dp) return
        if (wave_number < 0.0_dp .or. quadrature_degree < 0) return
        call build_maxwell_rwg_surface_space( &
            vertices, triangles, edge_vertices, edge_triangles, status)
        if (status /= 0) return
        if (size(edge_vertices, 2) == 0) return
        if (size(coefficients) /= size(edge_vertices, 2)) return
        call triangle_duffy_quadrature( &
            quadrature_degree, xi, eta, weights, status)
        if (status /= 0) return
        allocate(coefficients_bar(size(edge_vertices, 2)))
        allocate(magnetic_fields(3, size(edge_vertices, 2)))
        call evaluate_all_torus_curved_rwg_magnetic_fields_vjp( &
            vertices, triangles, parameters, edge_vertices, edge_triangles, &
            major_radius, minor_radius, observation, wave_number, xi, eta, &
            weights, coefficients, magnetic_field_bar, vertices_bar, parameters_bar, &
            major_radius_bar, minor_radius_bar, observation_bar, &
            wave_number_bar, magnetic_fields, status)
        if (status /= 0) return
        coefficients_bar = matmul( &
            conjg(transpose(magnetic_fields)), magnetic_field_bar)
        status = 0
    end subroutine evaluate_maxwell_torus_magnetic_geometry_vjp

    subroutine evaluate_maxwell_torus_curved_magnetic_field_rwg_3d_vjp( &
            vertices, triangles, parameters, major_radius, minor_radius, &
            observation, wave_number, quadrature_degree, magnetic_field_bar, &
            coefficients_bar, status)
        real(dp), intent(in) :: vertices(:, :), parameters(:, :)
        real(dp), intent(in) :: major_radius, minor_radius, observation(3)
        real(dp), intent(in) :: wave_number
        integer, intent(in) :: triangles(:, :), quadrature_degree
        complex(dp), intent(in) :: magnetic_field_bar(3)
        complex(dp), allocatable, intent(out) :: coefficients_bar(:)
        integer, intent(out) :: status

        complex(dp), allocatable :: basis_fields(:, :)

        if (allocated(coefficients_bar)) deallocate(coefficients_bar)
        call evaluate_torus_curved_magnetic_basis( &
            vertices, triangles, parameters, major_radius, minor_radius, &
            observation, wave_number, quadrature_degree, basis_fields, status)
        if (status /= 0) return
        allocate(coefficients_bar(size(basis_fields, 2)))
        coefficients_bar = matmul(conjg(transpose(basis_fields)), &
            magnetic_field_bar)
        status = 0
    end subroutine &
        evaluate_maxwell_torus_curved_magnetic_field_rwg_3d_vjp

    subroutine evaluate_torus_curved_magnetic_basis( &
            vertices, triangles, parameters, major_radius, minor_radius, &
            observation, wave_number, quadrature_degree, basis_fields, status)
        real(dp), intent(in) :: vertices(:, :), parameters(:, :)
        real(dp), intent(in) :: major_radius, minor_radius, observation(3)
        real(dp), intent(in) :: wave_number
        integer, intent(in) :: triangles(:, :), quadrature_degree
        complex(dp), allocatable, intent(out) :: basis_fields(:, :)
        integer, intent(out) :: status

        integer, allocatable :: edge_triangles(:, :), edge_vertices(:, :)
        real(dp), allocatable :: eta(:), weights(:), xi(:)

        if (allocated(basis_fields)) deallocate(basis_fields)
        status = 1
        if (size(vertices, 1) /= 3 .or. size(parameters, 1) /= 2) return
        if (size(parameters, 2) /= size(vertices, 2)) return
        if (major_radius <= minor_radius .or. minor_radius <= 0.0_dp) return
        if (wave_number < 0.0_dp) return
        call build_maxwell_rwg_surface_space( &
            vertices, triangles, edge_vertices, edge_triangles, status)
        if (status /= 0) return
        call triangle_duffy_quadrature( &
            quadrature_degree, xi, eta, weights, status)
        if (status /= 0) return
        allocate(basis_fields(3, size(edge_vertices, 2)))
        call evaluate_all_torus_curved_rwg_magnetic_fields( &
            vertices, triangles, parameters, edge_vertices, edge_triangles, &
            major_radius, minor_radius, observation, wave_number, xi, eta, &
            weights, basis_fields, status)
    end subroutine evaluate_torus_curved_magnetic_basis

    subroutine solve_maxwell_pec_torus_curved_regularized_cfie_rwg_3d( &
            vertices, triangles, parameters, major_radius, minor_radius, &
            direction, polarization, wave_number, impedance, quadrature_degree, &
            tolerance, max_depth, mfie_offset, density, status)
        real(dp), intent(in) :: vertices(:, :), parameters(:, :)
        real(dp), intent(in) :: major_radius, minor_radius, direction(3)
        complex(dp), intent(in) :: polarization(3)
        real(dp), intent(in) :: wave_number, impedance, tolerance, mfie_offset
        integer, intent(in) :: triangles(:, :), quadrature_degree, max_depth
        complex(dp), allocatable, intent(out) :: density(:)
        integer, intent(out) :: status

        complex(dp), allocatable :: densities(:, :)
        real(dp) :: directions(3, 1)
        complex(dp) :: polarizations(3, 1)

        directions(:, 1) = direction
        polarizations(:, 1) = polarization
        call &
            solve_maxwell_pec_torus_curved_regularized_cfie_rwg_multiple_3d( &
            vertices, triangles, parameters, major_radius, minor_radius, &
            directions, polarizations, wave_number, impedance, &
            quadrature_degree, tolerance, max_depth, mfie_offset, densities, &
            status)
        if (status /= 0) return
        allocate(density(size(densities, 1)))
        density = densities(:, 1)
    end subroutine solve_maxwell_pec_torus_curved_regularized_cfie_rwg_3d

    subroutine &
            solve_maxwell_pec_torus_curved_regularized_cfie_rwg_multiple_3d( &
            vertices, triangles, parameters, major_radius, minor_radius, &
            directions, polarizations, wave_number, impedance, &
            quadrature_degree, tolerance, max_depth, mfie_offset, densities, &
            status)
        real(dp), intent(in) :: vertices(:, :), parameters(:, :)
        real(dp), intent(in) :: major_radius, minor_radius, directions(:, :)
        complex(dp), intent(in) :: polarizations(:, :)
        real(dp), intent(in) :: wave_number, impedance, tolerance, mfie_offset
        integer, intent(in) :: triangles(:, :), quadrature_degree, max_depth
        complex(dp), allocatable, intent(out) :: densities(:, :)
        integer, intent(out) :: status

        complex(dp), allocatable :: bc_rhs(:), cfie(:, :), efie(:, :)
        complex(dp), allocatable :: efie_rhs(:), mapped_rhs(:, :)
        complex(dp), allocatable :: mapped_solution(:, :)
        complex(dp), allocatable :: mfie(:, :), product(:, :)
        complex(dp), allocatable :: regularizer(:, :), right_hand_side(:, :)
        complex(dp), allocatable :: mass(:, :)
        real(dp), allocatable :: real_mass(:, :)
        integer :: incidence, incidence_count, info, system_size

        status = 1
        if (allocated(densities)) deallocate(densities)
        if (size(directions, 1) /= 3) return
        if (size(polarizations, 1) /= 3) return
        incidence_count = size(directions, 2)
        if (incidence_count < 1) return
        if (size(polarizations, 2) /= incidence_count) return
        call assemble_maxwell_torus_curved_regularized_cfie_rwg_3d( &
            vertices, triangles, parameters, major_radius, minor_radius, &
            wave_number, impedance, quadrature_degree, tolerance, max_depth, &
            mfie_offset, cfie, efie, mfie, regularizer, product, status)
        if (status /= 0) return
        call assemble_maxwell_torus_curved_rwg_rbc_pairing( &
            vertices, triangles, parameters, major_radius, minor_radius, &
            quadrature_degree, real_mass, status)
        if (status /= 0) return
        system_size = size(real_mass, 1)
        allocate( &
            mass(system_size, system_size), &
            mapped_rhs(system_size, incidence_count), &
            mapped_solution(system_size, incidence_count), &
            right_hand_side(system_size, incidence_count), &
            densities(system_size, incidence_count))
        do incidence = 1, incidence_count
            call assemble_maxwell_torus_curved_plane_wave_rhs_rwg_3d( &
                vertices, triangles, parameters, major_radius, minor_radius, &
                directions(:, incidence), polarizations(:, incidence), &
                wave_number, quadrature_degree, efie_rhs, status)
            if (status /= 0) return
            call assemble_maxwell_torus_curved_plane_wave_rhs_bc_3d( &
                vertices, triangles, parameters, major_radius, minor_radius, &
                directions(:, incidence), polarizations(:, incidence), &
                wave_number, quadrature_degree, bc_rhs, status)
            if (status /= 0) return
            mapped_rhs(:, incidence) = efie_rhs
            right_hand_side(:, incidence) = bc_rhs
        end do
        mass = transpose(cmplx(real_mass, 0.0_dp, dp))
        call dense_solve(mass, mapped_rhs, mapped_solution, info)
        if (info /= 0) then
            status = 2
            return
        end if
        right_hand_side = &
            right_hand_side - matmul(regularizer, mapped_solution)
        call dense_solve(cfie, right_hand_side, densities, info)
        if (info /= 0) then
            status = 3
            return
        end if
        status = 0
    end subroutine &
        solve_maxwell_pec_torus_curved_regularized_cfie_rwg_multiple_3d

    subroutine assemble_maxwell_torus_curved_plane_wave_rhs_bc_3d( &
            vertices, triangles, parameters, major_radius, minor_radius, &
            direction, polarization, wave_number, quadrature_degree, &
            right_hand_side, status)
        real(dp), intent(in) :: vertices(:, :), parameters(:, :)
        real(dp), intent(in) :: major_radius, minor_radius, direction(3)
        complex(dp), intent(in) :: polarization(3)
        real(dp), intent(in) :: wave_number
        integer, intent(in) :: triangles(:, :), quadrature_degree
        complex(dp), allocatable, intent(out) :: right_hand_side(:)
        integer, intent(out) :: status

        integer, allocatable :: refined_triangles(:, :)
        real(dp), allocatable :: eta(:), refined_parameters(:, :)
        real(dp), allocatable :: refined_vertices(:, :)
        real(dp), allocatable :: transformation(:, :), weights(:), xi(:)
        real(dp) :: divergence, jacobian, local_value(3), point(3)
        complex(dp) :: incident_field(3)
        integer :: basis, local_edge, node, refined_panel, row

        status = 1
        if (allocated(right_hand_side)) deallocate(right_hand_side)
        if (wave_number < 0.0_dp) return
        if (abs(norm2(direction) - 1.0_dp) > &
            128.0_dp*epsilon(1.0_dp)) return
        if (abs(sum(polarization*direction)) > &
            128.0_dp*epsilon(1.0_dp)*max(1.0_dp, &
            sqrt(sum(abs(polarization)**2)))) return
        if (sqrt(sum(abs(polarization)**2)) <= tiny(1.0_dp)) return
        call build_maxwell_bc_transformation( &
            vertices, triangles, refined_vertices, refined_triangles, &
            transformation, status, torus_parameters=parameters, &
            torus_major_radius=major_radius, torus_minor_radius=minor_radius, &
            refined_torus_parameters=refined_parameters)
        if (status /= 0) return
        call triangle_duffy_quadrature( &
            quadrature_degree, xi, eta, weights, status)
        if (status /= 0) return
        allocate(right_hand_side(size(transformation, 2)))
        right_hand_side = cmplx(0.0_dp, 0.0_dp, dp)
        do refined_panel = 1, size(refined_triangles, 2)
            do node = 1, size(weights)
                incident_field = cmplx(0.0_dp, 0.0_dp, dp)
                do local_edge = 1, 3
                    call evaluate_maxwell_torus_curved_localized_rwg_basis( &
                        refined_vertices, refined_triangles, &
                        refined_parameters, refined_panel, local_edge, &
                        major_radius, minor_radius, xi(node), eta(node), point, &
                        local_value, divergence, jacobian, status)
                    if (status /= 0) return
                    incident_field = polarization*exp(cmplx( &
                        0.0_dp, wave_number*dot_product(direction, point), dp))
                    row = 3*(refined_panel - 1) + local_edge
                    do basis = 1, size(transformation, 2)
                        right_hand_side(basis) = right_hand_side(basis) - &
                            jacobian*weights(node)*transformation(row, basis)* &
                            sum(cmplx(local_value, 0.0_dp, dp)*incident_field)
                    end do
                end do
            end do
        end do
        status = 0
    end subroutine assemble_maxwell_torus_curved_plane_wave_rhs_bc_3d

    subroutine assemble_maxwell_torus_curved_regularized_cfie_rwg_3d( &
            vertices, triangles, parameters, major_radius, minor_radius, &
            wave_number, impedance, quadrature_degree, tolerance, max_depth, &
            mfie_offset, matrix, efie, mfie, regularizer, regularized_efie, &
            status)
        real(dp), intent(in) :: vertices(:, :), parameters(:, :)
        real(dp), intent(in) :: major_radius, minor_radius, wave_number
        real(dp), intent(in) :: impedance, tolerance, mfie_offset
        integer, intent(in) :: triangles(:, :), quadrature_degree, max_depth
        complex(dp), allocatable, intent(out) :: matrix(:, :), efie(:, :)
        complex(dp), allocatable, intent(out) :: mfie(:, :), regularizer(:, :)
        complex(dp), allocatable, intent(out) :: regularized_efie(:, :)
        integer, intent(out) :: status

        complex(dp), allocatable :: mass(:, :), mapped_efie(:, :)
        real(dp), allocatable :: real_mass(:, :)
        integer :: info, system_size

        status = 1
        if (wave_number <= 0.0_dp .or. impedance <= 0.0_dp) return
        call assemble_maxwell_torus_curved_efie_rwg_3d( &
            vertices, triangles, parameters, major_radius, minor_radius, &
            wave_number, impedance, quadrature_degree, tolerance, max_depth, &
            efie, status)
        if (status /= 0) return
        call assemble_maxwell_torus_curved_efie_bc_imaginary_3d( &
            vertices, triangles, parameters, major_radius, minor_radius, &
            wave_number, impedance, quadrature_degree, tolerance, max_depth, &
            regularizer, status)
        if (status /= 0) return
        call assemble_maxwell_torus_curved_mfie_rwg_rbc_3d( &
            vertices, triangles, parameters, major_radius, minor_radius, &
            wave_number, quadrature_degree, mfie_offset, mfie, status)
        if (status /= 0) return
        call assemble_maxwell_torus_curved_rwg_rbc_pairing( &
            vertices, triangles, parameters, major_radius, minor_radius, &
            quadrature_degree, real_mass, status)
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
    end subroutine assemble_maxwell_torus_curved_regularized_cfie_rwg_3d

    subroutine assemble_maxwell_torus_curved_efie_bc_imaginary_3d( &
            vertices, triangles, parameters, major_radius, minor_radius, &
            decay_rate, impedance, quadrature_degree, tolerance, max_depth, &
            matrix, status)
        real(dp), intent(in) :: vertices(:, :), parameters(:, :)
        real(dp), intent(in) :: major_radius, minor_radius, decay_rate
        real(dp), intent(in) :: impedance, tolerance
        integer, intent(in) :: triangles(:, :), quadrature_degree, max_depth
        complex(dp), allocatable, intent(out) :: matrix(:, :)
        integer, intent(out) :: status

        integer, allocatable :: refined_triangles(:, :)
        real(dp), allocatable :: refined_parameters(:, :)
        real(dp), allocatable :: refined_vertices(:, :), transformation(:, :)
        complex(dp), allocatable :: complex_transformation(:, :)
        complex(dp), allocatable :: refined_matrix(:, :)

        status = 1
        if (allocated(matrix)) deallocate(matrix)
        call build_maxwell_bc_to_refined_rwg( &
            vertices, triangles, refined_vertices, refined_triangles, &
            transformation, status, torus_parameters=parameters, &
            torus_major_radius=major_radius, torus_minor_radius=minor_radius, &
            refined_torus_parameters=refined_parameters)
        if (status /= 0) return
        call assemble_maxwell_torus_curved_efie_imaginary_rwg_3d( &
            refined_vertices, refined_triangles, refined_parameters, &
            major_radius, minor_radius, decay_rate, impedance, &
            quadrature_degree, tolerance, max_depth, refined_matrix, status)
        if (status /= 0) return
        complex_transformation = cmplx(transformation, 0.0_dp, dp)
        matrix = matmul( &
            transpose(complex_transformation), &
            matmul(refined_matrix, complex_transformation))
        status = 0
    end subroutine assemble_maxwell_torus_curved_efie_bc_imaginary_3d

    subroutine assemble_maxwell_torus_curved_efie_imaginary_rwg_3d( &
            vertices, triangles, parameters, major_radius, minor_radius, &
            decay_rate, impedance, quadrature_degree, tolerance, max_depth, &
            matrix, status)
        real(dp), intent(in) :: vertices(:, :), parameters(:, :)
        real(dp), intent(in) :: major_radius, minor_radius, decay_rate
        real(dp), intent(in) :: impedance, tolerance
        integer, intent(in) :: triangles(:, :), quadrature_degree, max_depth
        complex(dp), allocatable, intent(out) :: matrix(:, :)
        integer, intent(out) :: status

        complex(dp), allocatable :: scalar_potential(:, :)
        complex(dp), allocatable :: vector_potential(:, :)

        status = 1
        if (allocated(matrix)) deallocate(matrix)
        if (decay_rate <= 0.0_dp .or. impedance <= 0.0_dp) return
        call assemble_maxwell_torus_curved_potential_operators_rwg_3d( &
            vertices, triangles, parameters, major_radius, minor_radius, &
            decay_rate, quadrature_degree, tolerance, max_depth, &
            vector_potential, scalar_potential, status, decaying_kernel=.true.)
        if (status /= 0) return
        matrix = -impedance*( &
            decay_rate*vector_potential + scalar_potential/decay_rate)
        status = 0
    end subroutine assemble_maxwell_torus_curved_efie_imaginary_rwg_3d

    subroutine assemble_maxwell_torus_curved_mfie_rwg_rbc_3d( &
            vertices, triangles, parameters, major_radius, minor_radius, &
            wave_number, quadrature_degree, relative_offset, matrix, status)
        real(dp), intent(in) :: vertices(:, :), parameters(:, :)
        real(dp), intent(in) :: major_radius, minor_radius, wave_number
        real(dp), intent(in) :: relative_offset
        integer, intent(in) :: triangles(:, :), quadrature_degree
        complex(dp), allocatable, intent(out) :: matrix(:, :)
        integer, intent(out) :: status

        complex(dp), allocatable :: full_offset(:, :), half_offset(:, :)
        complex(dp), allocatable :: quarter_offset(:, :)

        status = 1
        if (allocated(matrix)) deallocate(matrix)
        if (relative_offset <= 0.0_dp) return
        call assemble_maxwell_torus_curved_mfie_offset_trace_rwg_rbc_3d( &
            vertices, triangles, parameters, major_radius, minor_radius, &
            wave_number, quadrature_degree, relative_offset, full_offset, &
            status)
        if (status /= 0) return
        call assemble_maxwell_torus_curved_mfie_offset_trace_rwg_rbc_3d( &
            vertices, triangles, parameters, major_radius, minor_radius, &
            wave_number, quadrature_degree, relative_offset/2.0_dp, &
            half_offset, status)
        if (status /= 0) return
        call assemble_maxwell_torus_curved_mfie_offset_trace_rwg_rbc_3d( &
            vertices, triangles, parameters, major_radius, minor_radius, &
            wave_number, quadrature_degree, relative_offset/4.0_dp, &
            quarter_offset, status)
        if (status /= 0) return
        matrix = full_offset/3.0_dp - 2.0_dp*half_offset + &
            8.0_dp*quarter_offset/3.0_dp
        status = 0
    end subroutine assemble_maxwell_torus_curved_mfie_rwg_rbc_3d

    subroutine assemble_maxwell_torus_curved_mfie_offset_trace_rwg_rbc_3d( &
            vertices, triangles, parameters, major_radius, minor_radius, &
            wave_number, quadrature_degree, relative_offset, matrix, status)
        real(dp), intent(in) :: vertices(:, :), parameters(:, :)
        real(dp), intent(in) :: major_radius, minor_radius, wave_number
        real(dp), intent(in) :: relative_offset
        integer, intent(in) :: triangles(:, :), quadrature_degree
        complex(dp), allocatable, intent(out) :: matrix(:, :)
        integer, intent(out) :: status

        integer, allocatable :: edge_triangles(:, :), edge_vertices(:, :)
        integer, allocatable :: refined_triangles(:, :)
        real(dp), allocatable :: eta(:), refined_parameters(:, :)
        real(dp), allocatable :: refined_vertices(:, :)
        real(dp), allocatable :: transformation(:, :), weights(:), xi(:)
        real(dp), allocatable :: bc_values(:, :)
        complex(dp), allocatable :: magnetic_fields(:, :)
        real(dp) :: divergence, jacobian, local_value(3), normal(3), point(3)
        real(dp) :: target(3)
        integer :: local_edge, node, refined_panel, row, test_basis, trial_basis

        status = 1
        if (allocated(matrix)) deallocate(matrix)
        if (wave_number < 0.0_dp) return
        if (abs(relative_offset) <= tiny(1.0_dp)) return
        call build_maxwell_bc_transformation( &
            vertices, triangles, refined_vertices, refined_triangles, &
            transformation, status, torus_parameters=parameters, &
            torus_major_radius=major_radius, torus_minor_radius=minor_radius, &
            refined_torus_parameters=refined_parameters)
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
                    call evaluate_maxwell_torus_curved_localized_rwg_basis( &
                        refined_vertices, refined_triangles, &
                        refined_parameters, refined_panel, local_edge, &
                        major_radius, minor_radius, xi(node), eta(node), point, &
                        local_value, divergence, jacobian, status)
                    if (status /= 0) return
                    row = 3*(refined_panel - 1) + local_edge
                    do test_basis = 1, size(edge_vertices, 2)
                        bc_values(:, test_basis) = &
                            bc_values(:, test_basis) + &
                            transformation(row, test_basis)*local_value
                    end do
                end do
                normal = torus_unit_normal(point, major_radius)
                target = point + relative_offset*minor_radius*normal
                call evaluate_all_torus_curved_rwg_magnetic_fields( &
                    vertices, triangles, parameters, edge_vertices, &
                    edge_triangles, major_radius, minor_radius, target, &
                    wave_number, xi, eta, weights, magnetic_fields, status)
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
        assemble_maxwell_torus_curved_mfie_offset_trace_rwg_rbc_3d

    subroutine evaluate_all_torus_curved_rwg_magnetic_fields( &
            vertices, triangles, parameters, edge_vertices, edge_triangles, &
            major_radius, minor_radius, target, wave_number, xi, eta, weights, &
            magnetic_fields, status)
        real(dp), intent(in) :: vertices(:, :), parameters(:, :)
        real(dp), intent(in) :: major_radius, minor_radius, target(3)
        real(dp), intent(in) :: wave_number, xi(:), eta(:), weights(:)
        integer, intent(in) :: triangles(:, :), edge_vertices(:, :)
        integer, intent(in) :: edge_triangles(:, :)
        complex(dp), intent(out) :: magnetic_fields(:, :)
        integer, intent(out) :: status

        real(dp) :: basis_value(3), displacement(3), distance, divergence
        real(dp) :: jacobian, point(3)
        complex(dp) :: curl_integrand(3), gradient_green(3), green
        integer :: basis, node, panel

        magnetic_fields = cmplx(0.0_dp, 0.0_dp, dp)
        status = 1
        do panel = 1, size(triangles, 2)
            do node = 1, size(weights)
                do basis = 1, size(edge_vertices, 2)
                    if (.not. any(edge_triangles(:, basis) == panel)) cycle
                    call evaluate_maxwell_torus_curved_rwg_basis( &
                        vertices, triangles, parameters, edge_vertices, &
                        edge_triangles, basis, panel, major_radius, &
                        minor_radius, xi(node), eta(node), point, basis_value, &
                        divergence, jacobian, status)
                    if (status /= 0) return
                    displacement = target - point
                    distance = norm2(displacement)
                    if (distance <= 128.0_dp*epsilon(1.0_dp)*minor_radius) return
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
    end subroutine evaluate_all_torus_curved_rwg_magnetic_fields

    subroutine evaluate_all_torus_curved_rwg_magnetic_fields_jvp( &
            vertices, triangles, parameters, edge_vertices, edge_triangles, &
            major_radius, minor_radius, target, wave_number, xi, eta, weights, &
            vertices_dot, parameters_dot, major_radius_dot, minor_radius_dot, &
            target_dot, wave_number_dot, magnetic_fields, magnetic_fields_dot, &
            status)
        real(dp), intent(in) :: vertices(:, :), parameters(:, :)
        real(dp), intent(in) :: major_radius, minor_radius, target(3)
        real(dp), intent(in) :: wave_number, xi(:), eta(:), weights(:)
        integer, intent(in) :: triangles(:, :), edge_vertices(:, :)
        integer, intent(in) :: edge_triangles(:, :)
        real(dp), intent(in) :: vertices_dot(:, :), parameters_dot(:, :)
        real(dp), intent(in) :: major_radius_dot, minor_radius_dot
        real(dp), intent(in) :: target_dot(3), wave_number_dot
        complex(dp), intent(out) :: magnetic_fields(:, :), magnetic_fields_dot(:, :)
        real(dp) :: basis_value(3), basis_value_dot(3)
        real(dp) :: displacement(3), displacement_dot(3), distance
        real(dp) :: distance_dot, divergence, divergence_dot, jacobian
        real(dp) :: jacobian_dot, point(3), point_dot(3), unit(3), unit_dot(3)
        complex(dp) :: curl_integrand(3), curl_integrand_dot(3)
        complex(dp) :: gradient_green(3), gradient_green_dot(3), green, green_dot
        complex(dp) :: q, q_dot
        integer :: basis, node, panel, status

        magnetic_fields = cmplx(0.0_dp, 0.0_dp, dp)
        magnetic_fields_dot = cmplx(0.0_dp, 0.0_dp, dp)
        status = 1
        do panel = 1, size(triangles, 2)
            do node = 1, size(weights)
                do basis = 1, size(edge_vertices, 2)
                    if (.not. any(edge_triangles(:, basis) == panel)) cycle
                    call evaluate_maxwell_torus_curved_rwg_basis_jvp( &
                        vertices, triangles, parameters, edge_vertices, &
                        edge_triangles, basis, panel, major_radius, minor_radius, &
                        xi(node), eta(node), vertices_dot, parameters_dot, &
                        major_radius_dot, minor_radius_dot, point, basis_value, &
                        divergence, jacobian, point_dot, basis_value_dot, &
                        divergence_dot, jacobian_dot, status)
                    if (status /= 0) return
                    displacement = target - point
                    displacement_dot = target_dot - point_dot
                    distance = norm2(displacement)
                    if (distance <= &
                        128.0_dp*epsilon(1.0_dp)*minor_radius) return
                    distance_dot = dot_product(displacement, displacement_dot)/ &
                        distance
                    green = exp(cmplx(0.0_dp, wave_number*distance, dp))/ &
                        (4.0_dp*acos(-1.0_dp)*distance)
                    green_dot = green*(cmplx( &
                        0.0_dp, wave_number_dot*distance + &
                        wave_number*distance_dot, dp) - distance_dot/distance)
                    unit = displacement/distance
                    unit_dot = (displacement_dot - unit*distance_dot)/distance
                    q = cmplx(0.0_dp, wave_number, dp) - 1.0_dp/distance
                    q_dot = cmplx(0.0_dp, wave_number_dot, dp) + &
                        distance_dot/distance**2
                    gradient_green = green*q*unit
                    gradient_green_dot = (green_dot*q + green*q_dot)*unit + &
                        green*q*unit_dot
                    curl_integrand = complex_cross_product( &
                        gradient_green, cmplx(basis_value, 0.0_dp, dp))
                    curl_integrand_dot = complex_cross_product( &
                        gradient_green_dot, cmplx(basis_value, 0.0_dp, dp)) + &
                        complex_cross_product( &
                        gradient_green, cmplx(basis_value_dot, 0.0_dp, dp))
                    magnetic_fields(:, basis) = magnetic_fields(:, basis) + &
                        weights(node)*jacobian*curl_integrand
                    magnetic_fields_dot(:, basis) = &
                        magnetic_fields_dot(:, basis) + weights(node)*( &
                        jacobian_dot*curl_integrand + &
                        jacobian*curl_integrand_dot)
                end do
            end do
        end do
        status = 0
    end subroutine evaluate_all_torus_curved_rwg_magnetic_fields_jvp

    subroutine evaluate_all_torus_curved_rwg_magnetic_fields_vjp( &
            vertices, triangles, parameters, edge_vertices, edge_triangles, &
            major_radius, minor_radius, target, wave_number, xi, eta, weights, &
            coefficients, magnetic_field_bar, vertices_bar, parameters_bar, &
            major_radius_bar, &
            minor_radius_bar, target_bar, wave_number_bar, magnetic_fields, &
            status)
        real(dp), intent(in) :: vertices(:, :), parameters(:, :)
        real(dp), intent(in) :: major_radius, minor_radius, target(3)
        real(dp), intent(in) :: wave_number, xi(:), eta(:), weights(:)
        integer, intent(in) :: triangles(:, :), edge_vertices(:, :)
        integer, intent(in) :: edge_triangles(:, :)
        complex(dp), intent(in) :: coefficients(:), magnetic_field_bar(3)
        real(dp), intent(out) :: vertices_bar(:, :), parameters_bar(:, :)
        real(dp), intent(out) :: major_radius_bar, minor_radius_bar
        real(dp), intent(out) :: target_bar(3), wave_number_bar
        complex(dp), intent(out) :: magnetic_fields(:, :)
        real(dp) :: basis_value(3), displacement(3), distance, divergence
        real(dp) :: jacobian, point(3), unit(3), unit_bar(3)
        real(dp) :: distance_bar, displacement_bar(3)
        real(dp) :: point_bar(3), value_bar(3), jacobian_bar
        real(dp) :: local_major_bar, local_minor_bar
        real(dp) :: local_vertices_bar(3, 2), local_parameters_bar(2, 3)
        complex(dp) :: curl_integrand(3), curl_bar(3), gradient_green(3)
        complex(dp) :: green, q, factor
        complex(dp) :: gradient_bar(3), wave_factor, unit_factor
        complex(dp) :: seed(3)
        integer :: basis, local, node, panel, status_local
        integer :: vertex, status

        vertices_bar = 0.0_dp
        parameters_bar = 0.0_dp
        major_radius_bar = 0.0_dp
        minor_radius_bar = 0.0_dp
        target_bar = 0.0_dp
        wave_number_bar = 0.0_dp
        magnetic_fields = cmplx(0.0_dp, 0.0_dp, dp)
        status = 1
        do panel = 1, size(triangles, 2)
            do node = 1, size(weights)
                do basis = 1, size(edge_vertices, 2)
                    if (.not. any(edge_triangles(:, basis) == panel)) cycle
                    call evaluate_maxwell_torus_curved_rwg_basis( &
                        vertices, triangles, parameters, edge_vertices, &
                        edge_triangles, basis, panel, major_radius, minor_radius, &
                        xi(node), eta(node), point, basis_value, divergence, &
                        jacobian, status_local)
                    if (status_local /= 0) return
                    displacement = target - point
                    distance = norm2(displacement)
                    if (distance <= &
                        128.0_dp*epsilon(1.0_dp)*minor_radius) return
                    unit = displacement/distance
                    green = exp(cmplx(0.0_dp, wave_number*distance, dp))/ &
                        (4.0_dp*acos(-1.0_dp)*distance)
                    q = cmplx(0.0_dp, wave_number, dp) - 1.0_dp/distance
                    factor = green*q
                    gradient_green = factor*unit
                    curl_integrand = complex_cross_product( &
                        gradient_green, cmplx(basis_value, 0.0_dp, dp))
                    magnetic_fields(:, basis) = magnetic_fields(:, basis) + &
                        weights(node)*jacobian*curl_integrand

                    seed = magnetic_field_bar*conjg(coefficients(basis))
                    curl_bar = weights(node)*jacobian*seed
                    jacobian_bar = weights(node)*real(sum(conjg( &
                        seed)*curl_integrand), dp)
                    gradient_bar = complex_cross_product( &
                        cmplx(basis_value, 0.0_dp, dp), curl_bar)
                    value_bar = real(complex_cross_product( &
                        curl_bar, conjg(gradient_green)), dp)
                    unit_bar = real(conjg(gradient_bar)*factor, dp)
                    wave_factor = cmplx(0.0_dp, 1.0_dp, dp)
                    unit_factor = green*(q**2 + 1.0_dp/distance**2)
                    distance_bar = real(sum(conjg(gradient_bar)* &
                        unit_factor*unit), dp)
                    wave_number_bar = wave_number_bar + real(sum( &
                        conjg(gradient_bar)*green*wave_factor*( &
                        distance*q + 1.0_dp)*unit), dp)
                    displacement_bar = distance_bar*unit + &
                        (unit_bar - unit*dot_product(unit, unit_bar))/distance
                    target_bar = target_bar + displacement_bar
                    point_bar = -displacement_bar
                    call evaluate_maxwell_torus_curved_rwg_basis_vjp( &
                        vertices, triangles, parameters, edge_vertices, &
                        edge_triangles, basis, panel, major_radius, minor_radius, &
                        xi(node), eta(node), point_bar, value_bar, 0.0_dp, &
                        jacobian_bar, local_vertices_bar, local_parameters_bar, &
                        local_major_bar, local_minor_bar, status_local)
                    if (status_local /= 0) return
                    do vertex = 1, 2
                        vertices_bar(:, edge_vertices(vertex, basis)) = &
                            vertices_bar(:, edge_vertices(vertex, basis)) + &
                            local_vertices_bar(:, vertex)
                    end do
                    do local = 1, 3
                        parameters_bar(:, triangles(local, panel)) = &
                            parameters_bar(:, triangles(local, panel)) + &
                            local_parameters_bar(:, local)
                    end do
                    major_radius_bar = major_radius_bar + local_major_bar
                    minor_radius_bar = minor_radius_bar + local_minor_bar
                end do
            end do
        end do
        status = 0
    end subroutine evaluate_all_torus_curved_rwg_magnetic_fields_vjp

    subroutine assemble_maxwell_torus_curved_rwg_rbc_pairing( &
            vertices, triangles, parameters, major_radius, minor_radius, &
            quadrature_degree, matrix, status)
        real(dp), intent(in) :: vertices(:, :), parameters(:, :)
        real(dp), intent(in) :: major_radius, minor_radius
        integer, intent(in) :: triangles(:, :), quadrature_degree
        real(dp), allocatable, intent(out) :: matrix(:, :)
        integer, intent(out) :: status

        integer, allocatable :: edge_triangles(:, :), edge_vertices(:, :)
        integer, allocatable :: refined_triangles(:, :)
        real(dp), allocatable :: eta(:), refined_parameters(:, :)
        real(dp), allocatable :: refined_vertices(:, :)
        real(dp), allocatable :: transformation(:, :), weights(:), xi(:)
        real(dp), allocatable :: bc_values(:, :), rwg_values(:, :)
        real(dp) :: coarse_eta, coarse_jacobian, coarse_xi, divergence
        real(dp) :: local_value(3)
        real(dp) :: normal(3), point(3), refined_jacobian, rotated_bc(3)
        integer :: basis, local_edge, node, parent, refined_panel, row
        integer :: test_basis

        status = 1
        if (allocated(matrix)) deallocate(matrix)
        call build_maxwell_bc_transformation( &
            vertices, triangles, refined_vertices, refined_triangles, &
            transformation, status, torus_parameters=parameters, &
            torus_major_radius=major_radius, torus_minor_radius=minor_radius, &
            refined_torus_parameters=refined_parameters)
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
            do node = 1, size(weights)
                bc_values = 0.0_dp
                do local_edge = 1, 3
                    call evaluate_maxwell_torus_curved_localized_rwg_basis( &
                        refined_vertices, refined_triangles, &
                        refined_parameters, refined_panel, local_edge, &
                        major_radius, minor_radius, xi(node), eta(node), point, &
                        local_value, divergence, refined_jacobian, status)
                    if (status /= 0) return
                    row = 3*(refined_panel - 1) + local_edge
                    do test_basis = 1, size(edge_vertices, 2)
                        bc_values(:, test_basis) = &
                            bc_values(:, test_basis) + &
                            transformation(row, test_basis)*local_value
                    end do
                end do
                call map_refined_torus_point_to_parent( &
                    refined_parameters(:, refined_triangles(:, refined_panel)), &
                    parameters(:, triangles(:, parent)), xi(node), eta(node), &
                    coarse_xi, coarse_eta, status)
                if (status /= 0) return
                rwg_values = 0.0_dp
                do basis = 1, size(edge_vertices, 2)
                    if (.not. any(edge_triangles(:, basis) == parent)) cycle
                    call evaluate_maxwell_torus_curved_rwg_basis( &
                        vertices, triangles, parameters, edge_vertices, &
                        edge_triangles, basis, parent, major_radius, &
                        minor_radius, coarse_xi, coarse_eta, point, &
                        rwg_values(:, basis), divergence, coarse_jacobian, status)
                    if (status /= 0) return
                end do
                normal = torus_unit_normal(point, major_radius)
                do test_basis = 1, size(edge_vertices, 2)
                    rotated_bc = real_cross_product( &
                        normal, bc_values(:, test_basis))
                    do basis = 1, size(edge_vertices, 2)
                        matrix(test_basis, basis) = &
                            matrix(test_basis, basis) + &
                            refined_jacobian*weights(node)*dot_product( &
                            rotated_bc, rwg_values(:, basis))
                    end do
                end do
            end do
        end do
        status = 0
    end subroutine assemble_maxwell_torus_curved_rwg_rbc_pairing

    subroutine assemble_maxwell_torus_curved_rwg_rbc_pairing_jvp( &
            vertices, triangles, parameters, major_radius, minor_radius, &
            quadrature_degree, vertices_dot, parameters_dot, major_radius_dot, &
            minor_radius_dot, matrix, matrix_dot, status)
        real(dp), intent(in) :: vertices(:, :), parameters(:, :)
        real(dp), intent(in) :: major_radius, minor_radius
        integer, intent(in) :: triangles(:, :), quadrature_degree
        real(dp), intent(in) :: vertices_dot(:, :), parameters_dot(:, :)
        real(dp), intent(in) :: major_radius_dot, minor_radius_dot
        real(dp), allocatable, intent(out) :: matrix(:, :), matrix_dot(:, :)
        integer, intent(out) :: status

        integer, allocatable :: edge_triangles(:, :), edge_vertices(:, :)
        integer, allocatable :: refined_triangles(:, :)
        integer, allocatable :: refined_triangles_dot(:, :)
        real(dp), allocatable :: eta(:), refined_parameters(:, :)
        real(dp), allocatable :: refined_parameters_dot(:, :)
        real(dp), allocatable :: refined_vertices(:, :), refined_vertices_dot(:, :)
        real(dp), allocatable :: transformation(:, :), transformation_dot(:, :)
        real(dp), allocatable :: weights(:), xi(:)
        real(dp), allocatable :: bc_values(:, :), bc_values_dot(:, :)
        real(dp), allocatable :: rwg_values(:, :), rwg_values_dot(:, :)
        real(dp) :: coarse_eta, coarse_eta_dot, coarse_jacobian
        real(dp) :: coarse_jacobian_dot, coarse_xi, coarse_xi_dot
        real(dp) :: divergence, divergence_dot
        real(dp) :: local_value(3), local_value_dot(3)
        real(dp) :: normal(3), normal_dot(3), point(3), point_dot(3)
        real(dp) :: refined_jacobian, refined_jacobian_dot
        real(dp) :: rotated_bc(3), rotated_bc_dot(3)
        integer :: basis, local_edge, node, parent, refined_panel, row
        integer :: test_basis

        status = 1
        if (allocated(matrix_dot)) deallocate(matrix_dot)
        call assemble_maxwell_torus_curved_rwg_rbc_pairing( &
            vertices, triangles, parameters, major_radius, minor_radius, &
            quadrature_degree, matrix, status)
        if (status /= 0) return
        if (any(shape(vertices_dot) /= shape(vertices))) return
        if (any(shape(parameters_dot) /= shape(parameters))) return
        call build_maxwell_bc_transformation( &
            vertices, triangles, refined_vertices, refined_triangles, &
            transformation, status, torus_parameters=parameters, &
            torus_major_radius=major_radius, torus_minor_radius=minor_radius, &
            refined_torus_parameters=refined_parameters)
        if (status /= 0) return
        call barycentric_refine_torus_surface_mesh_jvp( &
            vertices, triangles, parameters, major_radius, minor_radius, &
            vertices_dot, parameters_dot, major_radius_dot, minor_radius_dot, &
            refined_vertices_dot, refined_triangles_dot, refined_parameters_dot, &
            status)
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
                    call evaluate_maxwell_torus_curved_localized_rwg_basis_jvp( &
                        refined_vertices, refined_triangles, refined_parameters, &
                        refined_panel, local_edge, major_radius, minor_radius, &
                        xi(node), eta(node), refined_vertices_dot, &
                        refined_parameters_dot, major_radius_dot, minor_radius_dot, &
                        point, local_value, divergence, refined_jacobian, &
                        point_dot, local_value_dot, divergence_dot, &
                        refined_jacobian_dot, status)
                    if (status /= 0) return
                    row = 3*(refined_panel - 1) + local_edge
                    do test_basis = 1, size(edge_vertices, 2)
                        bc_values(:, test_basis) = &
                            bc_values(:, test_basis) + &
                            transformation(row, test_basis)*local_value
                        bc_values_dot(:, test_basis) = &
                            bc_values_dot(:, test_basis) + &
                            transformation_dot(row, test_basis)*local_value + &
                            transformation(row, test_basis)*local_value_dot
                    end do
                end do
                call map_refined_torus_point_to_parent_jvp( &
                    refined_parameters(:, refined_triangles(:, refined_panel)), &
                    parameters(:, triangles(:, parent)), xi(node), eta(node), &
                    refined_parameters_dot(:, refined_triangles(:, refined_panel)), &
                    parameters_dot(:, triangles(:, parent)), 0.0_dp, 0.0_dp, &
                    coarse_xi, coarse_eta, coarse_xi_dot, coarse_eta_dot, status)
                if (status /= 0) return
                rwg_values = 0.0_dp
                rwg_values_dot = 0.0_dp
                do basis = 1, size(edge_vertices, 2)
                    if (.not. any(edge_triangles(:, basis) == parent)) cycle
                    call evaluate_maxwell_torus_curved_rwg_basis_jvp( &
                        vertices, triangles, parameters, edge_vertices, &
                        edge_triangles, basis, parent, major_radius, minor_radius, &
                        coarse_xi, coarse_eta, vertices_dot, parameters_dot, &
                        major_radius_dot, minor_radius_dot, point, &
                        rwg_values(:, basis), divergence, coarse_jacobian, &
                        point_dot, local_value_dot, divergence_dot, &
                        coarse_jacobian_dot, status, coarse_xi_dot, coarse_eta_dot)
                    if (status /= 0) return
                    rwg_values_dot(:, basis) = local_value_dot
                end do
                normal = torus_unit_normal(point, major_radius)
                call torus_unit_normal_jvp( &
                    point, major_radius, point_dot, major_radius_dot, normal_dot)
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
    end subroutine assemble_maxwell_torus_curved_rwg_rbc_pairing_jvp

    subroutine assemble_maxwell_torus_curved_rwg_rbc_pairing_vjp( &
            vertices, triangles, parameters, major_radius, minor_radius, &
            quadrature_degree, matrix_bar, vertices_bar, parameters_bar, &
            major_radius_bar, minor_radius_bar, status)
        real(dp), intent(in) :: vertices(:, :), parameters(:, :)
        real(dp), intent(in) :: major_radius, minor_radius
        integer, intent(in) :: triangles(:, :), quadrature_degree
        real(dp), intent(in) :: matrix_bar(:, :)
        real(dp), intent(out) :: vertices_bar(:, :), parameters_bar(:, :)
        real(dp), intent(out) :: major_radius_bar, minor_radius_bar
        integer, intent(out) :: status

        integer, allocatable :: edge_triangles(:, :), edge_vertices(:, :)
        integer, allocatable :: refined_triangles(:, :)
        real(dp), allocatable :: eta(:), refined_parameters(:, :)
        real(dp), allocatable :: refined_vertices(:, :)
        real(dp), allocatable :: transformation(:, :), transformation_bar(:, :)
        real(dp), allocatable :: weights(:), xi(:)
        real(dp), allocatable :: bc_values(:, :), bc_values_bar(:, :)
        real(dp), allocatable :: rwg_values(:, :), rwg_values_bar(:, :)
        real(dp), allocatable :: refined_vertices_bar(:, :)
        real(dp), allocatable :: refined_parameters_bar(:, :)
        real(dp), allocatable :: refined_vertices_bar_bc(:, :)
        real(dp), allocatable :: local_refined_vertices_bar(:, :)
        real(dp), allocatable :: local_refined_parameters_bar(:, :)
        real(dp), allocatable :: barycentric_vertices_bar(:, :)
        real(dp), allocatable :: barycentric_parameters_bar(:, :)
        real(dp) :: coarse_eta, coarse_xi, coarse_eta_bar, coarse_xi_bar
        real(dp) :: local_eta_bar, local_xi_bar
        real(dp) :: coarse_jacobian, divergence
        real(dp) :: local_major_bar, local_minor_bar
        real(dp) :: local_value(3), local_value_bar(3)
        real(dp) :: local_vertices_bar(3, 2), local_parameters_bar(2, 3)
        real(dp) :: local_point_bar(3), localized_values(3, 3)
        real(dp) :: local_map_refined_parameters_bar(2, 3)
        real(dp) :: local_map_parent_parameters_bar(2, 3)
        real(dp) :: normal(3), normal_bar(3), point(3), point_bar(3)
        real(dp) :: refined_jacobian, refined_jacobian_bar
        real(dp) :: rotated_bc(3), rotated_bc_bar(3)
        real(dp) :: zero_divergence_bar
        real(dp) :: zero_jacobian_bar
        integer :: basis, local_edge, node, parent, refined_panel, row
        integer :: last_basis, test_basis, vertex

        vertices_bar = 0.0_dp
        parameters_bar = 0.0_dp
        major_radius_bar = 0.0_dp
        minor_radius_bar = 0.0_dp
        zero_divergence_bar = 0.0_dp
        zero_jacobian_bar = 0.0_dp
        status = 1
        if (size(vertices, 1) /= 3 .or. size(parameters, 1) /= 2) return
        if (size(parameters, 2) /= size(vertices, 2)) return
        if (major_radius <= minor_radius .or. minor_radius <= 0.0_dp) return
        call build_maxwell_bc_transformation( &
            vertices, triangles, refined_vertices, refined_triangles, &
            transformation, status, torus_parameters=parameters, &
            torus_major_radius=major_radius, torus_minor_radius=minor_radius, &
            refined_torus_parameters=refined_parameters)
        if (status /= 0) return
        call build_maxwell_rwg_surface_space( &
            vertices, triangles, edge_vertices, edge_triangles, status)
        if (status /= 0) return
        if (size(matrix_bar, 1) /= size(edge_vertices, 2)) return
        if (size(matrix_bar, 2) /= size(edge_vertices, 2)) return
        call triangle_duffy_quadrature( &
            quadrature_degree, xi, eta, weights, status)
        if (status /= 0) return
        allocate( &
            transformation_bar(size(transformation, 1), size(transformation, 2)), &
            bc_values(3, size(edge_vertices, 2)), &
            bc_values_bar(3, size(edge_vertices, 2)), &
            rwg_values(3, size(edge_vertices, 2)), &
            rwg_values_bar(3, size(edge_vertices, 2)), &
            refined_vertices_bar(size(refined_vertices, 1), &
            size(refined_vertices, 2)), &
            refined_parameters_bar(size(refined_parameters, 1), &
            size(refined_parameters, 2)), &
            local_refined_vertices_bar(size(refined_vertices, 1), &
            size(refined_vertices, 2)), &
            local_refined_parameters_bar(size(refined_parameters, 1), &
            size(refined_parameters, 2)))
        transformation_bar = 0.0_dp
        refined_vertices_bar = 0.0_dp
        refined_parameters_bar = 0.0_dp
        do refined_panel = 1, size(refined_triangles, 2)
            parent = (refined_panel - 1)/6 + 1
            do node = 1, size(weights)
                bc_values = 0.0_dp
                rwg_values = 0.0_dp
                do local_edge = 1, 3
                    call evaluate_maxwell_torus_curved_localized_rwg_basis( &
                        refined_vertices, refined_triangles, &
                        refined_parameters, refined_panel, local_edge, &
                        major_radius, minor_radius, xi(node), eta(node), point, &
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
                call map_refined_torus_point_to_parent( &
                    refined_parameters(:, refined_triangles(:, refined_panel)), &
                    parameters(:, triangles(:, parent)), xi(node), eta(node), &
                    coarse_xi, coarse_eta, status)
                if (status /= 0) return
                rwg_values = 0.0_dp
                last_basis = 0
                do basis = 1, size(edge_vertices, 2)
                    if (.not. any(edge_triangles(:, basis) == parent)) cycle
                    last_basis = basis
                    call evaluate_maxwell_torus_curved_rwg_basis( &
                        vertices, triangles, parameters, edge_vertices, &
                        edge_triangles, basis, parent, major_radius, &
                        minor_radius, coarse_xi, coarse_eta, point, &
                        rwg_values(:, basis), divergence, coarse_jacobian, status)
                    if (status /= 0) return
                end do
                normal = torus_unit_normal(point, major_radius)
                bc_values_bar = 0.0_dp
                rwg_values_bar = 0.0_dp
                normal_bar = 0.0_dp
                refined_jacobian_bar = 0.0_dp
                do test_basis = 1, size(edge_vertices, 2)
                    rotated_bc = real_cross_product( &
                        normal, bc_values(:, test_basis))
                    do basis = 1, size(edge_vertices, 2)
                        if (.not. any( &
                            edge_triangles(:, basis) == parent)) cycle
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
                call torus_unit_normal_vjp( &
                    point, major_radius, normal_bar, point_bar, local_major_bar)
                major_radius_bar = major_radius_bar + local_major_bar
                coarse_xi_bar = 0.0_dp
                coarse_eta_bar = 0.0_dp
                do basis = 1, size(edge_vertices, 2)
                    if (.not. any(edge_triangles(:, basis) == parent)) cycle
                    local_xi_bar = 0.0_dp
                    local_eta_bar = 0.0_dp
                    local_point_bar = 0.0_dp
                    if (basis == last_basis) local_point_bar = point_bar
                    call evaluate_maxwell_torus_curved_rwg_basis_vjp( &
                        vertices, triangles, parameters, edge_vertices, &
                        edge_triangles, basis, parent, major_radius, &
                        minor_radius, coarse_xi, coarse_eta, local_point_bar, &
                        rwg_values_bar(:, basis), zero_divergence_bar, &
                        zero_jacobian_bar, local_vertices_bar, &
                        local_parameters_bar, local_major_bar, local_minor_bar, &
                        status, local_xi_bar, local_eta_bar)
                    if (status /= 0) return
                    coarse_xi_bar = coarse_xi_bar + local_xi_bar
                    coarse_eta_bar = coarse_eta_bar + local_eta_bar
                    do vertex = 1, 2
                        vertices_bar(:, edge_vertices(vertex, basis)) = &
                            vertices_bar(:, edge_vertices(vertex, basis)) + &
                            local_vertices_bar(:, vertex)
                    end do
                    do vertex = 1, 3
                        parameters_bar(:, triangles(vertex, parent)) = &
                            parameters_bar(:, triangles(vertex, parent)) + &
                            local_parameters_bar(:, vertex)
                    end do
                    major_radius_bar = major_radius_bar + local_major_bar
                    minor_radius_bar = minor_radius_bar + local_minor_bar
                end do
                call map_refined_torus_point_to_parent_vjp( &
                    refined_parameters(:, refined_triangles(:, refined_panel)), &
                    parameters(:, triangles(:, parent)), xi(node), eta(node), &
                    coarse_xi_bar, coarse_eta_bar, &
                    local_map_refined_parameters_bar, &
                    local_map_parent_parameters_bar, status)
                if (status /= 0) return
                do vertex = 1, 3
                    refined_parameters_bar(:, refined_triangles(vertex, &
                        refined_panel)) = refined_parameters_bar(:, &
                        refined_triangles(vertex, refined_panel)) + &
                        local_map_refined_parameters_bar(:, vertex)
                    parameters_bar(:, triangles(vertex, parent)) = &
                        parameters_bar(:, triangles(vertex, parent)) + &
                        local_map_parent_parameters_bar(:, vertex)
                end do
                do local_edge = 1, 3
                    row = 3*(refined_panel - 1) + local_edge
                    local_value_bar = 0.0_dp
                    do test_basis = 1, size(edge_vertices, 2)
                        local_value_bar = local_value_bar + &
                            transformation(row, test_basis)* &
                            bc_values_bar(:, test_basis)
                        transformation_bar(row, test_basis) = &
                            transformation_bar(row, test_basis) + &
                            dot_product(bc_values_bar(:, test_basis), &
                            localized_values(:, local_edge))
                    end do
                    local_point_bar = 0.0_dp
                    if (status /= 0) return
                    call evaluate_maxwell_torus_curved_localized_rwg_basis_vjp( &
                        refined_vertices, refined_triangles, refined_parameters, &
                        refined_panel, local_edge, major_radius, minor_radius, &
                        xi(node), eta(node), local_point_bar, local_value_bar, &
                        zero_divergence_bar, &
                        merge(refined_jacobian_bar, zero_jacobian_bar, &
                        local_edge == 3), local_refined_vertices_bar, &
                        local_refined_parameters_bar, local_major_bar, &
                        local_minor_bar, &
                        status)
                    if (status /= 0) return
                    refined_vertices_bar = refined_vertices_bar + &
                        local_refined_vertices_bar
                    refined_parameters_bar = refined_parameters_bar + &
                        local_refined_parameters_bar
                    major_radius_bar = major_radius_bar + local_major_bar
                    minor_radius_bar = minor_radius_bar + local_minor_bar
                end do
            end do
        end do
        allocate( &
            refined_vertices_bar_bc(size(refined_vertices, 1), &
            size(refined_vertices, 2)), &
            barycentric_vertices_bar(size(vertices, 1), size(vertices, 2)), &
            barycentric_parameters_bar(size(parameters, 1), size(parameters, 2)))
        call differentiate_maxwell_bc_transformation_vjp( &
            vertices, triangles, refined_vertices, refined_triangles, &
            transformation_bar, refined_vertices_bar_bc, status)
        if (status /= 0) return
        refined_vertices_bar = refined_vertices_bar + refined_vertices_bar_bc
        call barycentric_refine_torus_surface_mesh_vjp( &
            vertices, triangles, parameters, major_radius, minor_radius, &
            refined_vertices_bar, refined_parameters_bar, &
            barycentric_vertices_bar, barycentric_parameters_bar, &
            local_major_bar, local_minor_bar, status)
        if (status /= 0) return
        vertices_bar = vertices_bar + barycentric_vertices_bar
        parameters_bar = parameters_bar + barycentric_parameters_bar
        major_radius_bar = major_radius_bar + local_major_bar
        minor_radius_bar = minor_radius_bar + local_minor_bar
        status = 0
    end subroutine assemble_maxwell_torus_curved_rwg_rbc_pairing_vjp

    pure subroutine evaluate_maxwell_torus_curved_localized_rwg_basis( &
            vertices, triangles, parameters, panel, local_edge, major_radius, &
            minor_radius, xi, eta, point, value, surface_divergence, &
            surface_jacobian, status)
        real(dp), intent(in) :: vertices(:, :), parameters(:, :)
        real(dp), intent(in) :: major_radius, minor_radius, xi, eta
        integer, intent(in) :: triangles(:, :), panel, local_edge
        real(dp), intent(out) :: point(3), value(3), surface_divergence
        real(dp), intent(out) :: surface_jacobian
        integer, intent(out) :: status

        real(dp) :: edge_length, opposite_coordinates(2)
        real(dp) :: panel_parameters(2, 3), tangent_eta(3), tangent_xi(3)
        integer :: edge_local_vertices(2), local, opposite

        point = 0.0_dp
        value = 0.0_dp
        surface_divergence = 0.0_dp
        surface_jacobian = 0.0_dp
        status = 1
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
            panel_parameters(:, local) = &
                parameters(:, triangles(local, panel))
        end do
        call evaluate_torus_curved_panel( &
            panel_parameters, major_radius, minor_radius, xi, eta, point, &
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
            vertices(:, triangles(edge_local_vertices(2), panel)) - &
            vertices(:, triangles(edge_local_vertices(1), panel)))
        value = edge_length/surface_jacobian*( &
            (xi - opposite_coordinates(1))*tangent_xi + &
            (eta - opposite_coordinates(2))*tangent_eta)
        surface_divergence = 2.0_dp*edge_length/surface_jacobian
        status = 0
    end subroutine evaluate_maxwell_torus_curved_localized_rwg_basis

    pure subroutine evaluate_maxwell_torus_curved_localized_rwg_basis_jvp( &
            vertices, triangles, parameters, panel, local_edge, major_radius, &
            minor_radius, xi, eta, vertices_dot, parameters_dot, &
            major_radius_dot, minor_radius_dot, point, value, surface_divergence, &
            surface_jacobian, point_dot, value_dot, surface_divergence_dot, &
            surface_jacobian_dot, status)
        real(dp), intent(in) :: vertices(:, :), parameters(:, :)
        integer, intent(in) :: triangles(:, :), panel, local_edge
        real(dp), intent(in) :: major_radius, minor_radius, xi, eta
        real(dp), intent(in) :: vertices_dot(:, :), parameters_dot(:, :)
        real(dp), intent(in) :: major_radius_dot, minor_radius_dot
        real(dp), intent(out) :: point(3), value(3), surface_divergence
        real(dp), intent(out) :: surface_jacobian, point_dot(3), value_dot(3)
        real(dp), intent(out) :: surface_divergence_dot, surface_jacobian_dot
        integer, intent(out) :: status

        real(dp) :: edge(3), edge_dot(3), edge_length, edge_length_dot
        real(dp) :: panel_parameters(2, 3), panel_parameters_dot(2, 3)
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
        if (size(vertices, 1) /= 3) return
        if (size(parameters, 1) /= 2) return
        if (size(parameters, 2) /= size(vertices, 2)) return
        if (any(shape(vertices_dot) /= shape(vertices))) return
        if (any(shape(parameters_dot) /= shape(parameters))) return
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
            panel_parameters(:, local) = parameters(:, triangles(local, panel))
            panel_parameters_dot(:, local) = &
                parameters_dot(:, triangles(local, panel))
        end do
        call evaluate_torus_curved_panel( &
            panel_parameters, major_radius, minor_radius, xi, eta, point, &
            tangent_xi, tangent_eta, surface_jacobian, status)
        if (status /= 0) return
        call evaluate_torus_curved_panel_jvp( &
            panel_parameters, major_radius, minor_radius, xi, eta, &
            panel_parameters_dot, major_radius_dot, minor_radius_dot, 0.0_dp, &
            0.0_dp, point_dot, tangent_xi_dot, tangent_eta_dot, &
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
        edge = vertices(:, triangles(edge_local_vertices(2), panel)) - &
            vertices(:, triangles(edge_local_vertices(1), panel))
        edge_dot = vertices_dot(:, triangles(edge_local_vertices(2), panel)) - &
            vertices_dot(:, triangles(edge_local_vertices(1), panel))
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
    end subroutine evaluate_maxwell_torus_curved_localized_rwg_basis_jvp

    pure subroutine evaluate_maxwell_torus_curved_localized_rwg_basis_vjp( &
            vertices, triangles, parameters, panel, local_edge, major_radius, &
            minor_radius, xi, eta, point_bar, value_bar, &
            surface_divergence_bar, surface_jacobian_bar, vertices_bar, &
            parameters_bar, major_radius_bar, minor_radius_bar, status)
        real(dp), intent(in) :: vertices(:, :), parameters(:, :)
        integer, intent(in) :: triangles(:, :), panel, local_edge
        real(dp), intent(in) :: major_radius, minor_radius, xi, eta
        real(dp), intent(in) :: point_bar(3), value_bar(3)
        real(dp), intent(in) :: surface_divergence_bar, surface_jacobian_bar
        real(dp), intent(out) :: vertices_bar(:, :), parameters_bar(:, :)
        real(dp), intent(out) :: major_radius_bar, minor_radius_bar
        integer, intent(out) :: status

        real(dp) :: edge(3), edge_bar(3), edge_length, edge_length_bar
        real(dp) :: panel_parameters(2, 3), panel_parameters_bar(2, 3)
        real(dp) :: tangent_eta(3), tangent_eta_bar(3)
        real(dp) :: tangent_xi(3), tangent_xi_bar(3)
        real(dp) :: vector(3), vector_bar(3), coefficient, coefficient_bar
        real(dp) :: local_jacobian_bar, local_major_bar, local_minor_bar
        real(dp) :: opposite_coordinates(2), surface_jacobian
        real(dp) :: dummy_point(3), dummy_tangent_xi(3), dummy_tangent_eta(3)
        integer :: edge_local_vertices(2), local, opposite

        vertices_bar = 0.0_dp
        parameters_bar = 0.0_dp
        major_radius_bar = 0.0_dp
        minor_radius_bar = 0.0_dp
        status = 1
        if (size(vertices, 1) /= 3) return
        if (size(parameters, 1) /= 2) return
        if (size(parameters, 2) /= size(vertices, 2)) return
        if (size(vertices_bar, 1) /= size(vertices, 1)) return
        if (size(vertices_bar, 2) /= size(vertices, 2)) return
        if (size(parameters_bar, 1) /= size(parameters, 1)) return
        if (size(parameters_bar, 2) /= size(parameters, 2)) return
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
            panel_parameters(:, local) = parameters(:, triangles(local, panel))
        end do
        call evaluate_torus_curved_panel( &
            panel_parameters, major_radius, minor_radius, xi, eta, dummy_point, &
            dummy_tangent_xi, dummy_tangent_eta, surface_jacobian, status)
        if (status /= 0) return
        select case (opposite)
        case (1)
            opposite_coordinates = [0.0_dp, 0.0_dp]
        case (2)
            opposite_coordinates = [1.0_dp, 0.0_dp]
        case (3)
            opposite_coordinates = [0.0_dp, 1.0_dp]
        end select
        edge = vertices(:, triangles(edge_local_vertices(2), panel)) - &
            vertices(:, triangles(edge_local_vertices(1), panel))
        edge_length = norm2(edge)
        if (edge_length <= tiny(1.0_dp)) return
        tangent_xi = dummy_tangent_xi
        tangent_eta = dummy_tangent_eta
        coefficient = edge_length/surface_jacobian
        vector = (xi - opposite_coordinates(1))*tangent_xi + &
            (eta - opposite_coordinates(2))*tangent_eta
        coefficient_bar = dot_product(value_bar, vector)
        vector_bar = coefficient*value_bar
        edge_length_bar = coefficient_bar/surface_jacobian + &
            2.0_dp*surface_divergence_bar/surface_jacobian
        local_jacobian_bar = surface_jacobian_bar - &
            edge_length*coefficient_bar/surface_jacobian**2 - &
            2.0_dp*edge_length*surface_divergence_bar/surface_jacobian**2
        edge_bar = edge_length_bar*edge/edge_length
        vertices_bar(:, triangles(edge_local_vertices(1), panel)) = &
            vertices_bar(:, triangles(edge_local_vertices(1), panel)) - edge_bar
        vertices_bar(:, triangles(edge_local_vertices(2), panel)) = &
            vertices_bar(:, triangles(edge_local_vertices(2), panel)) + edge_bar
        tangent_xi_bar = (xi - opposite_coordinates(1))*vector_bar
        tangent_eta_bar = (eta - opposite_coordinates(2))*vector_bar
        call evaluate_torus_curved_panel_vjp( &
            panel_parameters, major_radius, minor_radius, xi, eta, point_bar, &
            tangent_xi_bar, tangent_eta_bar, local_jacobian_bar, &
            panel_parameters_bar, local_major_bar, local_minor_bar, &
            dummy_point(1), dummy_point(2), status)
        if (status /= 0) return
        do local = 1, 3
            parameters_bar(:, triangles(local, panel)) = &
                parameters_bar(:, triangles(local, panel)) + &
                panel_parameters_bar(:, local)
        end do
        major_radius_bar = local_major_bar
        minor_radius_bar = local_minor_bar
        status = 0
    end subroutine evaluate_maxwell_torus_curved_localized_rwg_basis_vjp

    subroutine assemble_maxwell_torus_curved_efie_rwg_3d( &
            vertices, triangles, parameters, major_radius, minor_radius, &
            wave_number, impedance, quadrature_degree, tolerance, max_depth, &
            matrix, status)
        real(dp), intent(in) :: vertices(:, :), parameters(:, :)
        real(dp), intent(in) :: major_radius, minor_radius, wave_number
        real(dp), intent(in) :: impedance, tolerance
        integer, intent(in) :: triangles(:, :), quadrature_degree, max_depth
        complex(dp), allocatable, intent(out) :: matrix(:, :)
        integer, intent(out) :: status

        complex(dp), allocatable :: scalar_potential(:, :)
        complex(dp), allocatable :: vector_potential(:, :)

        status = 1
        if (allocated(matrix)) deallocate(matrix)
        if (wave_number <= 0.0_dp .or. impedance <= 0.0_dp) return
        call assemble_maxwell_torus_curved_potential_operators_rwg_3d( &
            vertices, triangles, parameters, major_radius, minor_radius, &
            wave_number, quadrature_degree, tolerance, max_depth, &
            vector_potential, scalar_potential, status)
        if (status /= 0) return
        matrix = cmplx(0.0_dp, wave_number*impedance, dp)*vector_potential - &
            cmplx(0.0_dp, impedance/wave_number, dp)*scalar_potential
        status = 0
    end subroutine assemble_maxwell_torus_curved_efie_rwg_3d

    subroutine assemble_maxwell_torus_curved_potential_operators_rwg_3d( &
            vertices, triangles, parameters, major_radius, minor_radius, &
            wave_number, quadrature_degree, tolerance, max_depth, &
            vector_potential, scalar_potential, status, decaying_kernel)
        real(dp), intent(in) :: vertices(:, :), parameters(:, :)
        real(dp), intent(in) :: major_radius, minor_radius, wave_number
        real(dp), intent(in) :: tolerance
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
        if (allocated(vector_potential)) deallocate(vector_potential)
        if (allocated(scalar_potential)) deallocate(scalar_potential)
        use_decaying_kernel = .false.
        if (present(decaying_kernel)) use_decaying_kernel = decaying_kernel
        call build_maxwell_rwg_surface_space( &
            vertices, triangles, edge_vertices, edge_triangles, status)
        if (status /= 0) return
        if (size(edge_vertices, 2) == 0) then
            status = 1
            return
        end if
        allocate( &
            vector_potential(size(edge_vertices, 2), size(edge_vertices, 2)), &
            scalar_potential(size(edge_vertices, 2), size(edge_vertices, 2)), &
            reference_divergence(2, size(edge_vertices, 2)))
        vector_potential = cmplx(0.0_dp, 0.0_dp, dp)
        scalar_potential = cmplx(0.0_dp, 0.0_dp, dp)
        do basis = 1, size(edge_vertices, 2)
            do first_slot = 1, 2
                call evaluate_maxwell_torus_curved_rwg_basis( &
                    vertices, triangles, parameters, edge_vertices, &
                    edge_triangles, basis, edge_triangles(first_slot, basis), &
                    major_radius, minor_radius, 1.0_dp/3.0_dp, 1.0_dp/3.0_dp, &
                    point, rwg_value, divergence, jacobian, status)
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
                                integrate_maxwell_torus_curved_coincident_rwg_pair_3d( &
                                vertices, triangles, parameters, edge_vertices, &
                                edge_triangles, first, first_panel, second, &
                                major_radius, minor_radius, wave_number, &
                                quadrature_degree, contribution, status, &
                                scalar_green, use_decaying_kernel)
                        else
                            call &
                                integrate_maxwell_torus_curved_adjacent_rwg_pair_3d( &
                                vertices, triangles, parameters, edge_vertices, &
                                edge_triangles, first, first_panel, second, &
                                second_panel, major_radius, minor_radius, &
                                wave_number, quadrature_degree, tolerance, &
                                max_depth, contribution, status, scalar_green, &
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
    end subroutine assemble_maxwell_torus_curved_potential_operators_rwg_3d

    subroutine integrate_maxwell_torus_curved_adjacent_rwg_pair_3d( &
            vertices, triangles, parameters, edge_vertices, edge_triangles, &
            first_basis, first_panel, second_basis, second_panel, major_radius, &
            minor_radius, wave_number, quadrature_degree, tolerance, max_depth, &
            value, status, scalar_value, decaying_kernel)
        real(dp), intent(in) :: vertices(:, :), parameters(:, :)
        real(dp), intent(in) :: major_radius, minor_radius, wave_number
        real(dp), intent(in) :: tolerance
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
        scalar_integral = cmplx(0.0_dp, 0.0_dp, dp)
        if (present(scalar_value)) scalar_value = scalar_integral
        use_decaying_kernel = .false.
        if (present(decaying_kernel)) use_decaying_kernel = decaying_kernel
        status = 1
        if (major_radius <= minor_radius .or. minor_radius <= 0.0_dp) return
        if (wave_number < 0.0_dp .or. tolerance <= 0.0_dp) return
        if (max_depth < 1 .or. first_panel == second_panel) return
        if (.not. any(edge_triangles(:, first_basis) == first_panel)) return
        if (.not. any(edge_triangles(:, second_basis) == second_panel)) return
        call triangle_duffy_quadrature( &
            quadrature_degree, xi, eta, weights, status)
        if (status /= 0) return
        reference_triangle = reshape( &
            [0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, 1.0_dp], [2, 3])
        call integrate_adaptive_torus_curved_rwg_pair( &
            vertices, triangles, parameters, edge_vertices, edge_triangles, &
            first_basis, first_panel, second_basis, second_panel, major_radius, &
            minor_radius, wave_number, reference_triangle, reference_triangle, &
            xi, eta, weights, tolerance, present(scalar_value), &
            use_decaying_kernel, 0, max_depth, value, scalar_integral, status)
        if (present(scalar_value)) scalar_value = scalar_integral
    end subroutine integrate_maxwell_torus_curved_adjacent_rwg_pair_3d

    recursive subroutine integrate_adaptive_torus_curved_rwg_pair( &
            vertices, triangles, parameters, edge_vertices, edge_triangles, &
            first_basis, first_panel, second_basis, second_panel, major_radius, &
            minor_radius, wave_number, first_reference, second_reference, xi, &
            eta, weights, tolerance, need_scalar, decaying_kernel, depth, &
            max_depth, value, scalar_value, status)
        real(dp), intent(in) :: vertices(:, :), parameters(:, :)
        real(dp), intent(in) :: major_radius, minor_radius, wave_number
        integer, intent(in) :: triangles(:, :), edge_vertices(:, :)
        integer, intent(in) :: edge_triangles(:, :)
        integer, intent(in) :: first_basis, first_panel, second_basis
        integer, intent(in) :: second_panel, depth, max_depth
        real(dp), intent(in) :: first_reference(2, 3)
        real(dp), intent(in) :: second_reference(2, 3)
        real(dp), intent(in) :: xi(:), eta(:), weights(:), tolerance
        logical, intent(in) :: need_scalar, decaying_kernel
        complex(dp), intent(out) :: value, scalar_value
        integer, intent(out) :: status

        real(dp) :: first_children(2, 3, 4), second_children(2, 3, 4)
        complex(dp) :: coarse, coarse_scalar, contribution
        complex(dp) :: contribution_scalar, refined, refined_scalar
        integer :: first_child, second_child

        call integrate_regular_torus_curved_rwg_pair( &
            vertices, triangles, parameters, edge_vertices, edge_triangles, &
            first_basis, first_panel, second_basis, second_panel, major_radius, &
            minor_radius, wave_number, first_reference, second_reference, xi, &
            eta, weights, coarse, coarse_scalar, decaying_kernel, status)
        if (status /= 0) return
        call subdivide_reference_triangle(first_reference, first_children)
        call subdivide_reference_triangle(second_reference, second_children)
        refined = cmplx(0.0_dp, 0.0_dp, dp)
        refined_scalar = cmplx(0.0_dp, 0.0_dp, dp)
        do first_child = 1, 4
            do second_child = 1, 4
                call integrate_regular_torus_curved_rwg_pair( &
                    vertices, triangles, parameters, edge_vertices, &
                    edge_triangles, first_basis, first_panel, second_basis, &
                    second_panel, major_radius, minor_radius, wave_number, &
                    first_children(:, :, first_child), &
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
                if (torus_reference_children_touch( &
                    parameters, triangles, first_panel, second_panel, &
                    major_radius, minor_radius, &
                    first_children(:, :, first_child), &
                    second_children(:, :, second_child))) then
                    call integrate_adaptive_torus_curved_rwg_pair( &
                        vertices, triangles, parameters, edge_vertices, &
                        edge_triangles, first_basis, first_panel, second_basis, &
                        second_panel, major_radius, minor_radius, wave_number, &
                        first_children(:, :, first_child), &
                        second_children(:, :, second_child), xi, eta, weights, &
                        tolerance, need_scalar, decaying_kernel, depth + 1, &
                        max_depth, contribution, contribution_scalar, status)
                else
                    call integrate_regular_torus_curved_rwg_pair( &
                        vertices, triangles, parameters, edge_vertices, &
                        edge_triangles, first_basis, first_panel, second_basis, &
                        second_panel, major_radius, minor_radius, wave_number, &
                        first_children(:, :, first_child), &
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
    end subroutine integrate_adaptive_torus_curved_rwg_pair

    subroutine integrate_regular_torus_curved_rwg_pair( &
            vertices, triangles, parameters, edge_vertices, edge_triangles, &
            first_basis, first_panel, second_basis, second_panel, major_radius, &
            minor_radius, wave_number, first_reference, second_reference, xi, &
            eta, weights, value, scalar_value, decaying_kernel, status)
        real(dp), intent(in) :: vertices(:, :), parameters(:, :)
        real(dp), intent(in) :: major_radius, minor_radius, wave_number
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
            call evaluate_maxwell_torus_curved_rwg_basis( &
                vertices, triangles, parameters, edge_vertices, &
                edge_triangles, first_basis, first_panel, major_radius, &
                minor_radius, first_xi_eta(1), first_xi_eta(2), first_point, &
                first_value, first_divergence, first_jacobian, status)
            if (status /= 0) return
            do second_node = 1, size(weights)
                second_xi_eta = reference_point( &
                    second_reference, xi(second_node), eta(second_node))
                call evaluate_maxwell_torus_curved_rwg_basis( &
                    vertices, triangles, parameters, edge_vertices, &
                    edge_triangles, second_basis, second_panel, major_radius, &
                    minor_radius, second_xi_eta(1), second_xi_eta(2), &
                    second_point, second_value, second_divergence, &
                    second_jacobian, status)
                if (status /= 0) return
                physical_distance = norm2(first_point - second_point)
                green = torus_boundary_green( &
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
    end subroutine integrate_regular_torus_curved_rwg_pair

    subroutine integrate_maxwell_torus_curved_coincident_rwg_pair_3d( &
            vertices, triangles, parameters, edge_vertices, edge_triangles, &
            first_basis, panel, second_basis, major_radius, minor_radius, &
            wave_number, quadrature_degree, value, status, scalar_value, &
            decaying_kernel)
        real(dp), intent(in) :: vertices(:, :), parameters(:, :)
        real(dp), intent(in) :: major_radius, minor_radius, wave_number
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
        if (major_radius <= minor_radius .or. minor_radius <= 0.0_dp) return
        if (wave_number < 0.0_dp) return
        if (.not. any(edge_triangles(:, first_basis) == panel)) return
        if (.not. any(edge_triangles(:, second_basis) == panel)) return
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
            call evaluate_maxwell_torus_curved_rwg_basis( &
                vertices, triangles, parameters, edge_vertices, &
                edge_triangles, first_basis, panel, major_radius, &
                minor_radius, xi(first_node), eta(first_node), first_point, &
                first_value, first_divergence, first_jacobian, status)
            if (status /= 0) return
            do wedge = 1, 3
                wedge_first = &
                    reference_vertices(:, wedge) - &
                    [xi(first_node), eta(first_node)]
                wedge_second = &
                    reference_vertices(:, modulo(wedge, 3) + 1) - &
                    [xi(first_node), eta(first_node)]
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
                        call evaluate_maxwell_torus_curved_rwg_basis( &
                            vertices, triangles, parameters, edge_vertices, &
                            edge_triangles, second_basis, panel, major_radius, &
                            minor_radius, second_xi, second_eta, second_point, &
                            second_value, second_divergence, second_jacobian, &
                            status)
                        if (status /= 0) return
                        physical_distance = norm2(first_point - second_point)
                        green = torus_boundary_green( &
                            wave_number, physical_distance, &
                            use_decaying_kernel)
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
    end subroutine integrate_maxwell_torus_curved_coincident_rwg_pair_3d

    subroutine evaluate_maxwell_torus_curved_far_field_rwg_3d( &
            vertices, triangles, parameters, major_radius, minor_radius, &
            coefficients, direction, wave_number, impedance, &
            quadrature_degree, far_field, status)
        real(dp), intent(in) :: vertices(:, :), parameters(:, :)
        real(dp), intent(in) :: major_radius, minor_radius, direction(3)
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
        if (major_radius <= minor_radius .or. minor_radius <= 0.0_dp) return
        if (wave_number <= 0.0_dp .or. impedance <= 0.0_dp) return
        if (abs(norm2(direction) - 1.0_dp) > &
            128.0_dp*epsilon(1.0_dp)) return
        call build_maxwell_rwg_surface_space( &
            vertices, triangles, edge_vertices, edge_triangles, status)
        if (status /= 0) return
        if (size(coefficients) /= size(edge_vertices, 2)) then
            status = 1
            return
        end if
        call triangle_duffy_quadrature( &
            quadrature_degree, xi, eta, weights, status)
        if (status /= 0) return
        do panel = 1, size(triangles, 2)
            do node = 1, size(weights)
                surface_current = cmplx(0.0_dp, 0.0_dp, dp)
                do basis = 1, size(edge_vertices, 2)
                    if (.not. any(edge_triangles(:, basis) == panel)) cycle
                    call evaluate_maxwell_torus_curved_rwg_basis( &
                        vertices, triangles, parameters, edge_vertices, &
                        edge_triangles, basis, panel, major_radius, &
                        minor_radius, xi(node), eta(node), point, basis_value, &
                        divergence, jacobian, status)
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
    end subroutine evaluate_maxwell_torus_curved_far_field_rwg_3d

    subroutine evaluate_maxwell_torus_curved_far_field_rwg_3d_jvp( &
            vertices, triangles, parameters, major_radius, minor_radius, &
            coefficients, direction, wave_number, impedance, quadrature_degree, &
            vertices_dot, parameters_dot, major_radius_dot, minor_radius_dot, &
            coefficients_dot, direction_dot, wave_number_dot, impedance_dot, &
            far_field, far_field_dot, status)
        real(dp), intent(in) :: vertices(:, :), parameters(:, :)
        real(dp), intent(in) :: major_radius, minor_radius, direction(3)
        complex(dp), intent(in) :: coefficients(:), coefficients_dot(:)
        real(dp), intent(in) :: wave_number, impedance
        integer, intent(in) :: triangles(:, :), quadrature_degree
        real(dp), intent(in) :: vertices_dot(:, :), parameters_dot(:, :)
        real(dp), intent(in) :: major_radius_dot, minor_radius_dot
        real(dp), intent(in) :: direction_dot(3), wave_number_dot, impedance_dot
        complex(dp), intent(out) :: far_field(3), far_field_dot(3)
        integer, intent(out) :: status

        integer, allocatable :: edge_triangles(:, :), edge_vertices(:, :)
        real(dp), allocatable :: eta(:), weights(:), xi(:)
        real(dp) :: basis_value(3), basis_value_dot(3)
        real(dp) :: divergence, divergence_dot, jacobian, jacobian_dot
        real(dp) :: point(3), point_dot(3)
        complex(dp) :: projection, projection_dot
        real(dp) :: phase_argument_dot
        complex(dp) :: phase, phase_dot, surface_current(3)
        complex(dp) :: surface_current_dot(3), transverse_current(3)
        complex(dp) :: transverse_current_dot(3), integral(3), integral_dot(3)
        complex(dp) :: factor, factor_dot
        integer :: basis, node, panel

        call evaluate_maxwell_torus_curved_far_field_rwg_3d( &
            vertices, triangles, parameters, major_radius, minor_radius, &
            coefficients, direction, wave_number, impedance, quadrature_degree, &
            far_field, status)
        if (status /= 0) return
        if (any(shape(vertices_dot) /= shape(vertices))) then
            status = 1
            return
        end if
        if (any(shape(parameters_dot) /= shape(parameters))) then
            status = 1
            return
        end if
        if (size(coefficients_dot) /= size(coefficients)) then
            status = 1
            return
        end if
        call build_maxwell_rwg_surface_space( &
            vertices, triangles, edge_vertices, edge_triangles, status)
        if (status /= 0) return
        call triangle_duffy_quadrature( &
            quadrature_degree, xi, eta, weights, status)
        if (status /= 0) return
        integral = cmplx(0.0_dp, 0.0_dp, dp)
        integral_dot = cmplx(0.0_dp, 0.0_dp, dp)
        do panel = 1, size(triangles, 2)
            do node = 1, size(weights)
                surface_current = cmplx(0.0_dp, 0.0_dp, dp)
                surface_current_dot = cmplx(0.0_dp, 0.0_dp, dp)
                do basis = 1, size(edge_vertices, 2)
                    if (.not. any(edge_triangles(:, basis) == panel)) cycle
                    call evaluate_maxwell_torus_curved_rwg_basis_jvp( &
                        vertices, triangles, parameters, edge_vertices, &
                        edge_triangles, basis, panel, major_radius, minor_radius, &
                        xi(node), eta(node), vertices_dot, parameters_dot, &
                        major_radius_dot, minor_radius_dot, point, basis_value, &
                        divergence, jacobian, point_dot, basis_value_dot, &
                        divergence_dot, jacobian_dot, status)
                    if (status /= 0) return
                    surface_current = surface_current + &
                        coefficients(basis)*basis_value
                    surface_current_dot = surface_current_dot + &
                        coefficients_dot(basis)*basis_value + &
                        coefficients(basis)*basis_value_dot
                end do
                projection = sum(direction*surface_current)
                projection_dot = sum(direction_dot*surface_current) + &
                    sum(direction*surface_current_dot)
                transverse_current = surface_current - direction*projection
                transverse_current_dot = surface_current_dot - &
                    direction_dot*projection - direction*projection_dot
                phase_argument_dot = -wave_number_dot*sum(direction*point) - &
                    wave_number*(sum(direction_dot*point) + sum(direction*point_dot))
                phase = exp(cmplx( &
                    0.0_dp, -wave_number*sum(direction*point), dp))
                phase_dot = cmplx(0.0_dp, phase_argument_dot, dp)*phase
                integral = integral + weights(node)*jacobian*phase*transverse_current
                integral_dot = integral_dot + weights(node)*( &
                    jacobian_dot*phase*transverse_current + jacobian*phase_dot* &
                    transverse_current + jacobian*phase*transverse_current_dot)
            end do
        end do
        factor = cmplx(0.0_dp, wave_number*impedance/ &
            (4.0_dp*acos(-1.0_dp)), dp)
        factor_dot = cmplx(0.0_dp, (wave_number_dot*impedance + &
            wave_number*impedance_dot)/(4.0_dp*acos(-1.0_dp)), dp)
        far_field_dot = factor_dot*integral + factor*integral_dot
        status = 0
    end subroutine evaluate_maxwell_torus_curved_far_field_rwg_3d_jvp

    subroutine evaluate_maxwell_torus_curved_far_field_rwg_3d_vjp( &
            vertices, triangles, parameters, major_radius, minor_radius, &
            coefficients, direction, wave_number, impedance, quadrature_degree, &
            far_field_bar, vertices_bar, parameters_bar, major_radius_bar, &
            minor_radius_bar, coefficients_bar, direction_bar, wave_number_bar, &
            impedance_bar, status)
        real(dp), intent(in) :: vertices(:, :), parameters(:, :)
        real(dp), intent(in) :: major_radius, minor_radius, direction(3)
        complex(dp), intent(in) :: coefficients(:), far_field_bar(3)
        real(dp), intent(in) :: wave_number, impedance
        integer, intent(in) :: triangles(:, :), quadrature_degree
        real(dp), intent(out) :: vertices_bar(:, :), parameters_bar(:, :)
        real(dp), intent(out) :: major_radius_bar, minor_radius_bar
        real(dp), intent(out) :: direction_bar(3), wave_number_bar, impedance_bar
        complex(dp), allocatable, intent(out) :: coefficients_bar(:)
        integer, intent(out) :: status

        integer, allocatable :: edge_triangles(:, :), edge_vertices(:, :)
        real(dp), allocatable :: eta(:), weights(:), xi(:)
        real(dp) :: basis_value(3), divergence, jacobian, point(3)
        real(dp) :: point_bar(3), value_bar(3), jacobian_bar
        real(dp) :: local_point_bar(3), local_jacobian_bar
        complex(dp) :: projection
        real(dp) :: phase_argument_bar
        real(dp) :: local_major_bar, local_minor_bar
        real(dp) :: local_vertices_bar(3, 2), local_parameters_bar(2, 3)
        complex(dp) :: phase, surface_current(3), transverse_current(3)
        complex(dp) :: surface_current_bar(3), transverse_current_bar(3)
        complex(dp) :: integral(3), integral_bar(3), phase_bar, factor
        complex(dp) :: direction_current_bar(3), primal_far_field(3)
        logical :: explicit_geometry_consumed
        integer :: basis, local, node, panel, status_local, vertex

        vertices_bar = 0.0_dp
        parameters_bar = 0.0_dp
        major_radius_bar = 0.0_dp
        minor_radius_bar = 0.0_dp
        direction_bar = 0.0_dp
        wave_number_bar = 0.0_dp
        impedance_bar = 0.0_dp
        if (allocated(coefficients_bar)) deallocate(coefficients_bar)
        call evaluate_maxwell_torus_curved_far_field_rwg_3d( &
            vertices, triangles, parameters, major_radius, minor_radius, &
            coefficients, direction, wave_number, impedance, quadrature_degree, &
            primal_far_field, status)
        if (status /= 0) return
        if (size(vertices_bar, 1) /= size(vertices, 1)) then
            status = 1
            return
        end if
        if (size(vertices_bar, 2) /= size(vertices, 2)) then
            status = 1
            return
        end if
        if (size(parameters_bar, 1) /= size(parameters, 1)) then
            status = 1
            return
        end if
        if (size(parameters_bar, 2) /= size(parameters, 2)) then
            status = 1
            return
        end if
        call build_maxwell_rwg_surface_space( &
            vertices, triangles, edge_vertices, edge_triangles, status)
        if (status /= 0) return
        if (size(coefficients) /= size(edge_vertices, 2)) then
            status = 1
            return
        end if
        call triangle_duffy_quadrature( &
            quadrature_degree, xi, eta, weights, status)
        if (status /= 0) return
        allocate(coefficients_bar(size(edge_vertices, 2)))
        coefficients_bar = cmplx(0.0_dp, 0.0_dp, dp)
        integral = cmplx(0.0_dp, 0.0_dp, dp)
        do panel = 1, size(triangles, 2)
            do node = 1, size(weights)
                surface_current = cmplx(0.0_dp, 0.0_dp, dp)
                do basis = 1, size(edge_vertices, 2)
                    if (.not. any(edge_triangles(:, basis) == panel)) cycle
                    call evaluate_maxwell_torus_curved_rwg_basis( &
                        vertices, triangles, parameters, edge_vertices, &
                        edge_triangles, basis, panel, major_radius, minor_radius, &
                        xi(node), eta(node), point, basis_value, divergence, &
                        jacobian, status_local)
                    if (status_local /= 0) return
                    surface_current = surface_current + &
                        coefficients(basis)*basis_value
                end do
                projection = sum(direction*surface_current)
                transverse_current = surface_current - direction*projection
                phase = exp(cmplx( &
                    0.0_dp, -wave_number*sum(direction*point), dp))
                integral = integral + weights(node)*jacobian*phase*transverse_current
            end do
        end do
        factor = cmplx(0.0_dp, wave_number*impedance/ &
            (4.0_dp*acos(-1.0_dp)), dp)
        integral_bar = conjg(factor)*far_field_bar
        wave_number_bar = real(sum(conjg(far_field_bar)*cmplx( &
            0.0_dp, impedance/(4.0_dp*acos(-1.0_dp)), dp)*integral), dp)
        impedance_bar = real(sum(conjg(far_field_bar)*cmplx( &
            0.0_dp, wave_number/(4.0_dp*acos(-1.0_dp)), dp)*integral), dp)
        do panel = 1, size(triangles, 2)
            do node = 1, size(weights)
                surface_current = cmplx(0.0_dp, 0.0_dp, dp)
                do basis = 1, size(edge_vertices, 2)
                    if (.not. any(edge_triangles(:, basis) == panel)) cycle
                    call evaluate_maxwell_torus_curved_rwg_basis( &
                        vertices, triangles, parameters, edge_vertices, &
                        edge_triangles, basis, panel, major_radius, minor_radius, &
                        xi(node), eta(node), point, basis_value, divergence, &
                        jacobian, status_local)
                    if (status_local /= 0) return
                    surface_current = surface_current + &
                        coefficients(basis)*basis_value
                end do
                projection = sum(direction*surface_current)
                transverse_current = surface_current - direction*projection
                phase = exp(cmplx( &
                    0.0_dp, -wave_number*sum(direction*point), dp))
                transverse_current_bar = weights(node)*jacobian* &
                    conjg(phase)*integral_bar
                phase_bar = weights(node)*jacobian*sum( &
                    conjg(transverse_current)*integral_bar)
                jacobian_bar = real(sum(conjg(integral_bar)*phase* &
                    transverse_current), dp)*weights(node)
                surface_current_bar = transverse_current_bar - &
                    direction*sum(direction*transverse_current_bar)
                direction_current_bar = -conjg(projection)*transverse_current_bar - &
                    surface_current*sum(direction*conjg(transverse_current_bar))
                direction_bar = direction_bar + real(direction_current_bar, dp)
                phase_argument_bar = real(conjg(phase_bar)* &
                    cmplx(0.0_dp, 1.0_dp, dp)*phase, dp)
                wave_number_bar = wave_number_bar - phase_argument_bar* &
                    sum(direction*point)
                direction_bar = direction_bar - phase_argument_bar*wave_number*point
                point_bar = -phase_argument_bar*wave_number*direction
                explicit_geometry_consumed = .false.
                do basis = 1, size(edge_vertices, 2)
                    if (.not. any(edge_triangles(:, basis) == panel)) cycle
                    call evaluate_maxwell_torus_curved_rwg_basis( &
                        vertices, triangles, parameters, edge_vertices, &
                        edge_triangles, basis, panel, major_radius, minor_radius, &
                        xi(node), eta(node), point, basis_value, divergence, &
                        jacobian, status_local)
                    if (status_local /= 0) return
                    value_bar = real(conjg(surface_current_bar)* &
                        coefficients(basis), dp)
                    coefficients_bar(basis) = coefficients_bar(basis) + &
                        sum(surface_current_bar*cmplx(basis_value, 0.0_dp, dp))
                    if (.not. explicit_geometry_consumed) then
                        local_point_bar = point_bar
                        local_jacobian_bar = jacobian_bar
                        explicit_geometry_consumed = .true.
                    else
                        local_point_bar = 0.0_dp
                        local_jacobian_bar = 0.0_dp
                    end if
                    call evaluate_maxwell_torus_curved_rwg_basis_vjp( &
                        vertices, triangles, parameters, edge_vertices, &
                        edge_triangles, basis, panel, major_radius, minor_radius, &
                        xi(node), eta(node), local_point_bar, value_bar, 0.0_dp, &
                        local_jacobian_bar, local_vertices_bar, local_parameters_bar, &
                        local_major_bar, local_minor_bar, status_local)
                    if (status_local /= 0) return
                    do vertex = 1, 2
                        vertices_bar(:, edge_vertices(vertex, basis)) = &
                            vertices_bar(:, edge_vertices(vertex, basis)) + &
                            local_vertices_bar(:, vertex)
                    end do
                    do local = 1, 3
                        parameters_bar(:, triangles(local, panel)) = &
                            parameters_bar(:, triangles(local, panel)) + &
                            local_parameters_bar(:, local)
                    end do
                    major_radius_bar = major_radius_bar + local_major_bar
                    minor_radius_bar = minor_radius_bar + local_minor_bar
                end do
            end do
        end do
        status = 0
    end subroutine evaluate_maxwell_torus_curved_far_field_rwg_3d_vjp

    subroutine assemble_maxwell_torus_curved_plane_wave_rhs_rwg_3d( &
            vertices, triangles, parameters, major_radius, minor_radius, &
            direction, polarization, wave_number, quadrature_degree, &
            right_hand_side, status)
        real(dp), intent(in) :: vertices(:, :), parameters(:, :)
        real(dp), intent(in) :: major_radius, minor_radius, direction(3)
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
        if (allocated(right_hand_side)) deallocate(right_hand_side)
        if (major_radius <= minor_radius .or. minor_radius <= 0.0_dp) return
        if (wave_number < 0.0_dp) return
        if (abs(norm2(direction) - 1.0_dp) > &
            128.0_dp*epsilon(1.0_dp)) return
        if (abs(sum(polarization*direction)) > &
            128.0_dp*epsilon(1.0_dp)*max(1.0_dp, &
            sqrt(sum(abs(polarization)**2)))) return
        if (sqrt(sum(abs(polarization)**2)) <= tiny(1.0_dp)) return
        call build_maxwell_rwg_surface_space( &
            vertices, triangles, edge_vertices, edge_triangles, status)
        if (status /= 0) return
        if (size(edge_vertices, 2) == 0) then
            status = 1
            return
        end if
        call triangle_duffy_quadrature( &
            quadrature_degree, xi, eta, weights, status)
        if (status /= 0) return
        allocate(right_hand_side(size(edge_vertices, 2)))
        right_hand_side = cmplx(0.0_dp, 0.0_dp, dp)
        do basis = 1, size(edge_vertices, 2)
            do panel = 1, 2
                do node = 1, size(weights)
                    call evaluate_maxwell_torus_curved_rwg_basis( &
                        vertices, triangles, parameters, edge_vertices, &
                        edge_triangles, basis, edge_triangles(panel, basis), &
                        major_radius, minor_radius, xi(node), eta(node), point, &
                        basis_value, divergence, jacobian, status)
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
    end subroutine assemble_maxwell_torus_curved_plane_wave_rhs_rwg_3d

    subroutine assemble_maxwell_torus_curved_plane_wave_rhs_rwg_3d_jvp( &
            vertices, triangles, parameters, major_radius, minor_radius, &
            direction, polarization, wave_number, quadrature_degree, &
            vertices_dot, parameters_dot, major_radius_dot, minor_radius_dot, &
            direction_dot, polarization_dot, wave_number_dot, right_hand_side, &
            right_hand_side_dot, status)
        real(dp), intent(in) :: vertices(:, :), parameters(:, :)
        real(dp), intent(in) :: major_radius, minor_radius, direction(3)
        complex(dp), intent(in) :: polarization(3)
        real(dp), intent(in) :: wave_number
        integer, intent(in) :: triangles(:, :), quadrature_degree
        real(dp), intent(in) :: vertices_dot(:, :), parameters_dot(:, :)
        real(dp), intent(in) :: major_radius_dot, minor_radius_dot
        real(dp), intent(in) :: direction_dot(3), wave_number_dot
        complex(dp), intent(in) :: polarization_dot(3)
        complex(dp), allocatable, intent(out) :: right_hand_side(:)
        complex(dp), allocatable, intent(out) :: right_hand_side_dot(:)
        integer, intent(out) :: status

        integer, allocatable :: edge_triangles(:, :), edge_vertices(:, :)
        real(dp), allocatable :: eta(:), weights(:), xi(:)
        real(dp) :: basis_value(3), basis_value_dot(3)
        real(dp) :: divergence, divergence_dot, jacobian, jacobian_dot
        real(dp) :: point(3), point_dot(3), phase_argument_dot
        complex(dp) :: incident_field(3), incident_field_dot(3)
        complex(dp) :: phase, phase_dot, contraction, contraction_dot
        integer :: basis, node, panel

        if (allocated(right_hand_side_dot)) deallocate(right_hand_side_dot)
        call assemble_maxwell_torus_curved_plane_wave_rhs_rwg_3d( &
            vertices, triangles, parameters, major_radius, minor_radius, direction, &
            polarization, wave_number, quadrature_degree, right_hand_side, status)
        if (status /= 0) return
        if (any(shape(vertices_dot) /= shape(vertices))) then
            status = 1
            return
        end if
        if (any(shape(parameters_dot) /= shape(parameters))) then
            status = 1
            return
        end if
        call build_maxwell_rwg_surface_space( &
            vertices, triangles, edge_vertices, edge_triangles, status)
        if (status /= 0) return
        call triangle_duffy_quadrature( &
            quadrature_degree, xi, eta, weights, status)
        if (status /= 0) return
        allocate(right_hand_side_dot(size(edge_vertices, 2)))
        right_hand_side_dot = cmplx(0.0_dp, 0.0_dp, dp)
        do basis = 1, size(edge_vertices, 2)
            do panel = 1, 2
                do node = 1, size(weights)
                    call evaluate_maxwell_torus_curved_rwg_basis_jvp( &
                        vertices, triangles, parameters, edge_vertices, &
                        edge_triangles, basis, edge_triangles(panel, basis), &
                        major_radius, minor_radius, xi(node), eta(node), &
                        vertices_dot, parameters_dot, major_radius_dot, &
                        minor_radius_dot, point, basis_value, divergence, jacobian, &
                        point_dot, basis_value_dot, divergence_dot, jacobian_dot, &
                        status)
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
                    contraction = sum(cmplx(basis_value, 0.0_dp, dp)* &
                        incident_field)
                    contraction_dot = sum(cmplx(basis_value_dot, 0.0_dp, dp)* &
                        incident_field + cmplx(basis_value, 0.0_dp, dp)* &
                        incident_field_dot)
                    right_hand_side_dot(basis) = right_hand_side_dot(basis) - &
                        weights(node)*(jacobian_dot*contraction + jacobian* &
                        contraction_dot)
                end do
            end do
        end do
        status = 0
    end subroutine assemble_maxwell_torus_curved_plane_wave_rhs_rwg_3d_jvp

    subroutine assemble_maxwell_torus_curved_plane_wave_rhs_rwg_3d_vjp( &
            vertices, triangles, parameters, major_radius, minor_radius, &
            direction, polarization, wave_number, quadrature_degree, &
            right_hand_side_bar, vertices_bar, parameters_bar, major_radius_bar, &
            minor_radius_bar, direction_bar, polarization_bar, wave_number_bar, &
            status)
        real(dp), intent(in) :: vertices(:, :), parameters(:, :)
        real(dp), intent(in) :: major_radius, minor_radius, direction(3)
        complex(dp), intent(in) :: polarization(3)
        real(dp), intent(in) :: wave_number
        integer, intent(in) :: triangles(:, :), quadrature_degree
        complex(dp), intent(in) :: right_hand_side_bar(:)
        real(dp), intent(out) :: vertices_bar(:, :), parameters_bar(:, :)
        real(dp), intent(out) :: major_radius_bar, minor_radius_bar
        real(dp), intent(out) :: direction_bar(3), wave_number_bar
        complex(dp), allocatable, intent(out) :: polarization_bar(:)
        integer, intent(out) :: status

        integer, allocatable :: edge_triangles(:, :), edge_vertices(:, :)
        real(dp), allocatable :: eta(:), weights(:), xi(:)
        real(dp) :: basis_value(3), divergence, jacobian, point(3)
        real(dp) :: point_bar(3), value_bar(3), jacobian_bar
        real(dp) :: phase_argument_bar, local_major_bar, local_minor_bar
        real(dp) :: local_vertices_bar(3, 2), local_parameters_bar(2, 3)
        complex(dp) :: contraction, phase, phase_bar, seed
        integer :: basis, local, node, panel, status_local, vertex

        vertices_bar = 0.0_dp
        parameters_bar = 0.0_dp
        major_radius_bar = 0.0_dp
        minor_radius_bar = 0.0_dp
        direction_bar = 0.0_dp
        wave_number_bar = 0.0_dp
        if (allocated(polarization_bar)) deallocate(polarization_bar)
        status = 1
        if (size(vertices, 1) /= 3) return
        if (size(parameters, 1) /= 2) return
        if (size(parameters, 2) /= size(vertices, 2)) return
        if (size(vertices_bar, 1) /= size(vertices, 1)) return
        if (size(vertices_bar, 2) /= size(vertices, 2)) return
        if (size(parameters_bar, 1) /= size(parameters, 1)) return
        if (size(parameters_bar, 2) /= size(parameters, 2)) return
        if (major_radius <= minor_radius) return
        if (minor_radius <= 0.0_dp) return
        if (wave_number < 0.0_dp) return
        if (abs(norm2(direction) - 1.0_dp) > &
            128.0_dp*epsilon(1.0_dp)) return
        if (abs(sum(polarization*direction)) > &
            128.0_dp*epsilon(1.0_dp)*max(1.0_dp, &
            sqrt(sum(abs(polarization)**2)))) return
        if (sqrt(sum(abs(polarization)**2)) <= tiny(1.0_dp)) return
        call build_maxwell_rwg_surface_space( &
            vertices, triangles, edge_vertices, edge_triangles, status)
        if (status /= 0) return
        if (size(right_hand_side_bar) /= size(edge_vertices, 2)) return
        call triangle_duffy_quadrature( &
            quadrature_degree, xi, eta, weights, status)
        if (status /= 0) return
        allocate(polarization_bar(3))
        polarization_bar = cmplx(0.0_dp, 0.0_dp, dp)
        do basis = 1, size(edge_vertices, 2)
            do panel = 1, 2
                do node = 1, size(weights)
                    call evaluate_maxwell_torus_curved_rwg_basis( &
                        vertices, triangles, parameters, edge_vertices, &
                        edge_triangles, basis, edge_triangles(panel, basis), &
                        major_radius, minor_radius, xi(node), eta(node), point, &
                        basis_value, divergence, jacobian, status_local)
                    if (status_local /= 0) return
                    phase = exp(cmplx( &
                        0.0_dp, wave_number*dot_product(direction, point), dp))
                    contraction = sum(cmplx(basis_value, 0.0_dp, dp)*polarization)
                    seed = -weights(node)*jacobian*conjg(phase)* &
                        right_hand_side_bar(basis)
                    phase_bar = -weights(node)*jacobian*conjg(contraction)* &
                        right_hand_side_bar(basis)
                    jacobian_bar = real(conjg(right_hand_side_bar(basis))* &
                        (-weights(node)*phase*contraction), dp)
                    value_bar = real(conjg(seed)*polarization, dp)
                    polarization_bar = polarization_bar + &
                        seed*cmplx(basis_value, 0.0_dp, dp)
                    phase_argument_bar = real(conjg(phase_bar)* &
                        cmplx(0.0_dp, 1.0_dp, dp)*phase, dp)
                    wave_number_bar = wave_number_bar + phase_argument_bar* &
                        dot_product(direction, point)
                    direction_bar = direction_bar + phase_argument_bar* &
                        wave_number*point
                    point_bar = phase_argument_bar*wave_number*direction
                    call evaluate_maxwell_torus_curved_rwg_basis_vjp( &
                        vertices, triangles, parameters, edge_vertices, &
                        edge_triangles, basis, edge_triangles(panel, basis), &
                        major_radius, minor_radius, xi(node), eta(node), point_bar, &
                        value_bar, 0.0_dp, jacobian_bar, local_vertices_bar, &
                        local_parameters_bar, local_major_bar, local_minor_bar, &
                        status_local)
                    if (status_local /= 0) return
                    do vertex = 1, 2
                        vertices_bar(:, edge_vertices(vertex, basis)) = &
                            vertices_bar(:, edge_vertices(vertex, basis)) + &
                            local_vertices_bar(:, vertex)
                    end do
                    do local = 1, 3
                        parameters_bar(:, triangles(local, edge_triangles(panel, basis))) = &
                            parameters_bar(:, triangles(local, edge_triangles(panel, basis))) + &
                            local_parameters_bar(:, local)
                    end do
                    major_radius_bar = major_radius_bar + local_major_bar
                    minor_radius_bar = minor_radius_bar + local_minor_bar
                end do
            end do
        end do
        status = 0
    end subroutine assemble_maxwell_torus_curved_plane_wave_rhs_rwg_3d_vjp

    subroutine assemble_maxwell_torus_curved_rwg_mass_matrix( &
            vertices, triangles, parameters, major_radius, minor_radius, &
            quadrature_degree, matrix, status)
        real(dp), intent(in) :: vertices(:, :), parameters(:, :)
        real(dp), intent(in) :: major_radius, minor_radius
        integer, intent(in) :: triangles(:, :), quadrature_degree
        real(dp), allocatable, intent(out) :: matrix(:, :)
        integer, intent(out) :: status

        integer, allocatable :: edge_triangles(:, :), edge_vertices(:, :)
        real(dp), allocatable :: eta(:), values(:, :), weights(:), xi(:)
        real(dp) :: divergence, jacobian, point(3)
        integer :: basis, node, panel, test_basis

        status = 1
        if (allocated(matrix)) deallocate(matrix)
        if (size(vertices, 1) /= 3) return
        if (size(parameters, 1) /= 2) return
        if (size(parameters, 2) /= size(vertices, 2)) return
        if (major_radius <= minor_radius .or. minor_radius <= 0.0_dp) return
        call build_maxwell_rwg_surface_space( &
            vertices, triangles, edge_vertices, edge_triangles, status)
        if (status /= 0) return
        if (size(edge_vertices, 2) == 0) then
            status = 1
            return
        end if
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
                    call evaluate_maxwell_torus_curved_rwg_basis( &
                        vertices, triangles, parameters, edge_vertices, &
                        edge_triangles, basis, panel, major_radius, &
                        minor_radius, xi(node), eta(node), point, &
                        values(:, basis), divergence, jacobian, status)
                    if (status /= 0) return
                end do
                do basis = 1, size(edge_vertices, 2)
                    if (.not. any(edge_triangles(:, basis) == panel)) cycle
                    do test_basis = 1, size(edge_vertices, 2)
                        if (.not. any( &
                            edge_triangles(:, test_basis) == panel)) cycle
                        matrix(test_basis, basis) = &
                            matrix(test_basis, basis) + &
                            weights(node)*jacobian*dot_product( &
                            values(:, test_basis), values(:, basis))
                    end do
                end do
            end do
        end do
        status = 0
    end subroutine assemble_maxwell_torus_curved_rwg_mass_matrix

    subroutine assemble_maxwell_torus_curved_rwg_mass_matrix_jvp( &
            vertices, triangles, parameters, major_radius, minor_radius, &
            quadrature_degree, vertices_dot, parameters_dot, major_radius_dot, &
            minor_radius_dot, matrix, matrix_dot, status)
        real(dp), intent(in) :: vertices(:, :), parameters(:, :)
        real(dp), intent(in) :: major_radius, minor_radius
        integer, intent(in) :: triangles(:, :), quadrature_degree
        real(dp), intent(in) :: vertices_dot(:, :), parameters_dot(:, :)
        real(dp), intent(in) :: major_radius_dot, minor_radius_dot
        real(dp), allocatable, intent(out) :: matrix(:, :), matrix_dot(:, :)
        integer, intent(out) :: status

        integer, allocatable :: edge_triangles(:, :), edge_vertices(:, :)
        real(dp), allocatable :: eta(:), values(:, :), values_dot(:, :)
        real(dp), allocatable :: weights(:), xi(:)
        real(dp) :: divergence, divergence_dot, jacobian, jacobian_dot
        real(dp) :: point(3), point_dot(3), value_dot(3)
        integer :: basis, node, panel, test_basis

        status = 1
        if (allocated(matrix_dot)) deallocate(matrix_dot)
        call assemble_maxwell_torus_curved_rwg_mass_matrix( &
            vertices, triangles, parameters, major_radius, minor_radius, &
            quadrature_degree, matrix, status)
        if (status /= 0) return
        if (any(shape(vertices_dot) /= shape(vertices))) return
        if (any(shape(parameters_dot) /= shape(parameters))) return
        call build_maxwell_rwg_surface_space( &
            vertices, triangles, edge_vertices, edge_triangles, status)
        if (status /= 0) return
        call triangle_duffy_quadrature( &
            quadrature_degree, xi, eta, weights, status)
        if (status /= 0) return
        allocate( &
            matrix_dot(size(edge_vertices, 2), size(edge_vertices, 2)), &
            values(3, size(edge_vertices, 2)), &
            values_dot(3, size(edge_vertices, 2)))
        matrix_dot = 0.0_dp
        do panel = 1, size(triangles, 2)
            do node = 1, size(weights)
                values = 0.0_dp
                values_dot = 0.0_dp
                do basis = 1, size(edge_vertices, 2)
                    if (.not. any(edge_triangles(:, basis) == panel)) cycle
                    call evaluate_maxwell_torus_curved_rwg_basis_jvp( &
                        vertices, triangles, parameters, edge_vertices, &
                        edge_triangles, basis, panel, major_radius, minor_radius, &
                        xi(node), eta(node), vertices_dot, parameters_dot, &
                        major_radius_dot, minor_radius_dot, point, values(:, basis), &
                        divergence, jacobian, point_dot, value_dot, &
                        divergence_dot, jacobian_dot, status)
                    if (status /= 0) return
                    values_dot(:, basis) = value_dot
                end do
                do basis = 1, size(edge_vertices, 2)
                    if (.not. any( &
                        edge_triangles(:, basis) == panel)) cycle
                    do test_basis = 1, size(edge_vertices, 2)
                        if (.not. any( &
                            edge_triangles(:, test_basis) == panel)) cycle
                        matrix_dot(test_basis, basis) = &
                            matrix_dot(test_basis, basis) + weights(node)*( &
                            jacobian_dot*dot_product( &
                            values(:, test_basis), values(:, basis)) + jacobian*( &
                            dot_product(values_dot(:, test_basis), values(:, basis)) + &
                            dot_product(values(:, test_basis), values_dot(:, basis))))
                    end do
                end do
            end do
        end do
        status = 0
    end subroutine assemble_maxwell_torus_curved_rwg_mass_matrix_jvp

    subroutine assemble_maxwell_torus_curved_rwg_mass_matrix_vjp( &
            vertices, triangles, parameters, major_radius, minor_radius, &
            quadrature_degree, matrix_bar, vertices_bar, parameters_bar, &
            major_radius_bar, minor_radius_bar, status)
        real(dp), intent(in) :: vertices(:, :), parameters(:, :)
        real(dp), intent(in) :: major_radius, minor_radius
        integer, intent(in) :: triangles(:, :), quadrature_degree
        real(dp), intent(in) :: matrix_bar(:, :)
        real(dp), intent(out) :: vertices_bar(:, :), parameters_bar(:, :)
        real(dp), intent(out) :: major_radius_bar, minor_radius_bar
        integer, intent(out) :: status

        integer, allocatable :: edge_triangles(:, :), edge_vertices(:, :)
        real(dp), allocatable :: eta(:), values(:, :), values_bar(:, :), weights(:)
        real(dp), allocatable :: xi(:)
        real(dp) :: divergence, jacobian, jacobian_bar, point(3)
        real(dp) :: value(3), point_bar(3), surface_divergence_bar
        real(dp) :: local_major_bar, local_minor_bar
        real(dp) :: local_vertices_bar(3, 2), local_parameters_bar(2, 3)
        real(dp) :: local_jacobian_bar
        logical :: jacobian_consumed
        integer :: basis, node, panel, test_basis, vertex, local

        vertices_bar = 0.0_dp
        parameters_bar = 0.0_dp
        major_radius_bar = 0.0_dp
        minor_radius_bar = 0.0_dp
        point_bar = 0.0_dp
        surface_divergence_bar = 0.0_dp
        status = 1
        if (size(vertices, 1) /= 3 .or. size(parameters, 1) /= 2) return
        if (size(parameters, 2) /= size(vertices, 2)) return
        if (major_radius <= minor_radius .or. minor_radius <= 0.0_dp) return
        call build_maxwell_rwg_surface_space( &
            vertices, triangles, edge_vertices, edge_triangles, status)
        if (status /= 0) return
        if (size(matrix_bar, 1) /= size(edge_vertices, 2) .or. &
            size(matrix_bar, 2) /= size(edge_vertices, 2)) return
        call triangle_duffy_quadrature( &
            quadrature_degree, xi, eta, weights, status)
        if (status /= 0) return
        allocate( &
            values(3, size(edge_vertices, 2)), &
            values_bar(3, size(edge_vertices, 2)))
        do panel = 1, size(triangles, 2)
            do node = 1, size(weights)
                values = 0.0_dp
                values_bar = 0.0_dp
                jacobian = 0.0_dp
                jacobian_consumed = .false.
                do basis = 1, size(edge_vertices, 2)
                    if (.not. any(edge_triangles(:, basis) == panel)) cycle
                    call evaluate_maxwell_torus_curved_rwg_basis( &
                        vertices, triangles, parameters, edge_vertices, &
                        edge_triangles, basis, panel, major_radius, minor_radius, &
                        xi(node), eta(node), point, values(:, basis), divergence, &
                        jacobian, status)
                    if (status /= 0) return
                end do
                jacobian_bar = 0.0_dp
                do basis = 1, size(edge_vertices, 2)
                    if (.not. any( &
                        edge_triangles(:, basis) == panel)) cycle
                    do test_basis = 1, size(edge_vertices, 2)
                        if (.not. any( &
                            edge_triangles(:, test_basis) == panel)) cycle
                        values_bar(:, test_basis) = values_bar(:, test_basis) + &
                            weights(node)*jacobian*matrix_bar(test_basis, basis)* &
                            values(:, basis)
                        values_bar(:, basis) = values_bar(:, basis) + &
                            weights(node)*jacobian*matrix_bar(test_basis, basis)* &
                            values(:, test_basis)
                        jacobian_bar = jacobian_bar + weights(node)* &
                            matrix_bar(test_basis, basis)*dot_product( &
                            values(:, test_basis), values(:, basis))
                    end do
                end do
                do basis = 1, size(edge_vertices, 2)
                    if (.not. any( &
                        edge_triangles(:, basis) == panel)) cycle
                    if (.not. jacobian_consumed) then
                        local_jacobian_bar = jacobian_bar
                        jacobian_consumed = .true.
                    else
                        local_jacobian_bar = 0.0_dp
                    end if
                    call evaluate_maxwell_torus_curved_rwg_basis_vjp( &
                        vertices, triangles, parameters, edge_vertices, &
                        edge_triangles, basis, panel, major_radius, minor_radius, &
                        xi(node), eta(node), point_bar, values_bar(:, basis), &
                        surface_divergence_bar, local_jacobian_bar, local_vertices_bar, &
                        local_parameters_bar, local_major_bar, local_minor_bar, &
                        status)
                    if (status /= 0) return
                    do vertex = 1, 2
                        vertices_bar(:, edge_vertices(vertex, basis)) = &
                            vertices_bar(:, edge_vertices(vertex, basis)) + &
                            local_vertices_bar(:, vertex)
                    end do
                    do local = 1, 3
                        parameters_bar(:, triangles(local, panel)) = &
                            parameters_bar(:, triangles(local, panel)) + &
                            local_parameters_bar(:, local)
                    end do
                    major_radius_bar = major_radius_bar + local_major_bar
                    minor_radius_bar = minor_radius_bar + local_minor_bar
                end do
            end do
        end do
        status = 0
    end subroutine assemble_maxwell_torus_curved_rwg_mass_matrix_vjp

    pure subroutine evaluate_maxwell_torus_curved_rwg_basis( &
            vertices, triangles, parameters, edge_vertices, edge_triangles, &
            basis, panel, major_radius, minor_radius, xi, eta, point, value, &
            surface_divergence, surface_jacobian, status)
        real(dp), intent(in) :: vertices(:, :), parameters(:, :)
        real(dp), intent(in) :: major_radius, minor_radius, xi, eta
        integer, intent(in) :: triangles(:, :), edge_vertices(:, :)
        integer, intent(in) :: edge_triangles(:, :), basis, panel
        real(dp), intent(out) :: point(3), value(3), surface_divergence
        real(dp), intent(out) :: surface_jacobian
        integer, intent(out) :: status

        real(dp) :: edge_length, opposite_coordinates(2)
        real(dp) :: panel_parameters(2, 3), tangent_eta(3), tangent_xi(3)
        integer :: local, next, opposite, orientation

        point = 0.0_dp
        value = 0.0_dp
        surface_divergence = 0.0_dp
        surface_jacobian = 0.0_dp
        status = 1
        if (size(vertices, 1) /= 3 .or. size(parameters, 1) /= 2) return
        if (size(parameters, 2) /= size(vertices, 2)) return
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
            panel_parameters(:, local) = &
                parameters(:, triangles(local, panel))
        end do
        call evaluate_torus_curved_panel( &
            panel_parameters, major_radius, minor_radius, xi, eta, point, &
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
    end subroutine evaluate_maxwell_torus_curved_rwg_basis

    pure subroutine evaluate_maxwell_torus_curved_rwg_basis_jvp( &
            vertices, triangles, parameters, edge_vertices, edge_triangles, &
            basis, panel, major_radius, minor_radius, xi, eta, vertices_dot, &
            parameters_dot, major_radius_dot, minor_radius_dot, point, value, &
            surface_divergence, surface_jacobian, point_dot, value_dot, &
            surface_divergence_dot, surface_jacobian_dot, status, xi_dot, eta_dot)
        real(dp), intent(in) :: vertices(:, :), parameters(:, :)
        real(dp), intent(in) :: major_radius, minor_radius, xi, eta
        integer, intent(in) :: triangles(:, :), edge_vertices(:, :)
        integer, intent(in) :: edge_triangles(:, :), basis, panel
        real(dp), intent(in) :: vertices_dot(:, :), parameters_dot(:, :)
        real(dp), intent(in) :: major_radius_dot, minor_radius_dot
        real(dp), intent(out) :: point(3), value(3), surface_divergence
        real(dp), intent(out) :: surface_jacobian, point_dot(3), value_dot(3)
        real(dp), intent(out) :: surface_divergence_dot, surface_jacobian_dot
        integer, intent(out) :: status
        real(dp), optional, intent(in) :: xi_dot, eta_dot

        real(dp) :: edge(3), edge_dot(3), edge_length, edge_length_dot
        real(dp) :: panel_parameters(2, 3), panel_parameters_dot(2, 3)
        real(dp) :: tangent_eta(3), tangent_eta_dot(3)
        real(dp) :: tangent_xi(3), tangent_xi_dot(3)
        real(dp) :: vector(3), vector_dot(3), coefficient, coefficient_dot
        real(dp) :: opposite_coordinates(2), jacobian
        real(dp) :: local_xi_dot, local_eta_dot
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
        call evaluate_maxwell_torus_curved_rwg_basis( &
            vertices, triangles, parameters, edge_vertices, edge_triangles, &
            basis, panel, major_radius, minor_radius, xi, eta, point, value, &
            surface_divergence, surface_jacobian, status)
        if (status /= 0) return
        if (any(shape(vertices_dot) /= shape(vertices))) return
        if (any(shape(parameters_dot) /= shape(parameters))) return
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
            panel_parameters(:, local) = parameters(:, triangles(local, panel))
            panel_parameters_dot(:, local) = &
                parameters_dot(:, triangles(local, panel))
        end do
        call evaluate_torus_curved_panel( &
            panel_parameters, major_radius, minor_radius, xi, eta, point, &
            tangent_xi, tangent_eta, surface_jacobian, status)
        if (status /= 0) return
        call evaluate_torus_curved_panel_jvp( &
            panel_parameters, major_radius, minor_radius, xi, eta, &
            panel_parameters_dot, major_radius_dot, minor_radius_dot, &
            local_xi_dot, local_eta_dot, point_dot, tangent_xi_dot, tangent_eta_dot, &
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
        jacobian = surface_jacobian
        coefficient = real(orientation, dp)*edge_length/jacobian
        coefficient_dot = real(orientation, dp)*( &
            edge_length_dot/jacobian - &
            edge_length*surface_jacobian_dot/jacobian**2)
        value_dot = coefficient_dot*vector + coefficient*vector_dot
        surface_divergence_dot = 2.0_dp*real(orientation, dp)*( &
            edge_length_dot/jacobian - &
            edge_length*surface_jacobian_dot/jacobian**2)
        status = 0
    end subroutine evaluate_maxwell_torus_curved_rwg_basis_jvp

    pure subroutine evaluate_maxwell_torus_curved_rwg_basis_vjp( &
            vertices, triangles, parameters, edge_vertices, edge_triangles, &
            basis, panel, major_radius, minor_radius, xi, eta, point_bar, &
            value_bar, surface_divergence_bar, surface_jacobian_bar, &
            edge_vertices_bar, panel_parameters_bar, major_radius_bar, &
            minor_radius_bar, status, xi_bar, eta_bar)
        real(dp), intent(in) :: vertices(:, :), parameters(:, :)
        real(dp), intent(in) :: major_radius, minor_radius, xi, eta
        integer, intent(in) :: triangles(:, :), edge_vertices(:, :)
        integer, intent(in) :: edge_triangles(:, :), basis, panel
        real(dp), intent(in) :: point_bar(3), value_bar(3)
        real(dp), intent(in) :: surface_divergence_bar, surface_jacobian_bar
        real(dp), intent(out) :: edge_vertices_bar(3, 2)
        real(dp), intent(out) :: panel_parameters_bar(2, 3)
        real(dp), intent(out) :: major_radius_bar, minor_radius_bar
        integer, intent(out) :: status
        real(dp), optional, intent(out) :: xi_bar, eta_bar

        real(dp) :: edge(3), edge_bar(3), edge_length, edge_length_bar
        real(dp) :: dummy_point(3), dummy_value(3), dummy_divergence
        real(dp) :: panel_parameters(2, 3), tangent_eta(3), tangent_eta_bar(3)
        real(dp) :: tangent_xi(3), tangent_xi_bar(3)
        real(dp) :: vector(3), vector_bar(3), coefficient, coefficient_bar
        real(dp) :: jacobian, jacobian_bar, local_major_bar, local_minor_bar
        real(dp) :: opposite_coordinates(2)
        real(dp) :: local_xi_bar, local_eta_bar
        real(dp) :: panel_xi_bar, panel_eta_bar
        integer :: local, next, opposite, orientation

        edge_vertices_bar = 0.0_dp
        panel_parameters_bar = 0.0_dp
        major_radius_bar = 0.0_dp
        minor_radius_bar = 0.0_dp
        local_xi_bar = 0.0_dp
        local_eta_bar = 0.0_dp
        if (present(xi_bar)) xi_bar = 0.0_dp
        if (present(eta_bar)) eta_bar = 0.0_dp
        status = 1
        call evaluate_maxwell_torus_curved_rwg_basis( &
            vertices, triangles, parameters, edge_vertices, edge_triangles, &
            basis, panel, major_radius, minor_radius, xi, eta, dummy_point, &
            dummy_value, dummy_divergence, jacobian, status)
        if (status /= 0) return
        ! The call above only supplies the public basis outputs in their
        ! declared order; evaluate the panel tangents explicitly below.
        do local = 1, 3
            panel_parameters(:, local) = parameters(:, triangles(local, panel))
        end do
        call evaluate_torus_curved_panel( &
            panel_parameters, major_radius, minor_radius, xi, eta, vector, &
            tangent_xi, tangent_eta, jacobian, status)
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
        vector = (xi - opposite_coordinates(1))*tangent_xi + &
            (eta - opposite_coordinates(2))*tangent_eta
        coefficient = real(orientation, dp)*edge_length/jacobian
        coefficient_bar = dot_product(value_bar, vector)
        vector_bar = coefficient*value_bar
        edge_length_bar = real(orientation, dp)*coefficient_bar/jacobian + &
            2.0_dp*real(orientation, dp)*surface_divergence_bar/jacobian
        jacobian_bar = surface_jacobian_bar - &
            real(orientation, dp)*edge_length*coefficient_bar/jacobian**2 - &
            2.0_dp*real(orientation, dp)*edge_length* &
            surface_divergence_bar/jacobian**2
        edge_bar = edge_length_bar*edge/edge_length
        edge_vertices_bar(:, 1) = -edge_bar
        edge_vertices_bar(:, 2) = edge_bar
        tangent_xi_bar = (xi - opposite_coordinates(1))*vector_bar
        tangent_eta_bar = (eta - opposite_coordinates(2))*vector_bar
        local_xi_bar = dot_product(vector_bar, tangent_xi)
        local_eta_bar = dot_product(vector_bar, tangent_eta)
        call evaluate_torus_curved_panel_vjp( &
            panel_parameters, major_radius, minor_radius, xi, eta, point_bar, &
            tangent_xi_bar, tangent_eta_bar, jacobian_bar, &
            panel_parameters_bar, local_major_bar, local_minor_bar, &
            panel_xi_bar, panel_eta_bar, status)
        if (status /= 0) return
        local_xi_bar = local_xi_bar + panel_xi_bar
        local_eta_bar = local_eta_bar + panel_eta_bar
        major_radius_bar = local_major_bar
        minor_radius_bar = local_minor_bar
        if (present(xi_bar)) xi_bar = local_xi_bar
        if (present(eta_bar)) eta_bar = local_eta_bar
        status = 0
    end subroutine evaluate_maxwell_torus_curved_rwg_basis_vjp

    pure logical function torus_reference_children_touch( &
            parameters, triangles, first_panel, second_panel, major_radius, &
            minor_radius, first_reference, second_reference) result(touch)
        real(dp), intent(in) :: parameters(:, :)
        integer, intent(in) :: triangles(:, :), first_panel, second_panel
        real(dp), intent(in) :: major_radius, minor_radius
        real(dp), intent(in) :: first_reference(2, 3)
        real(dp), intent(in) :: second_reference(2, 3)

        real(dp) :: first_panel_parameters(2, 3), first_point(3), jacobian
        real(dp) :: scale, second_panel_parameters(2, 3), second_point(3)
        real(dp) :: tangent_eta(3), tangent_xi(3)
        integer :: first_vertex, second_vertex, status

        scale = max(1.0_dp, major_radius + minor_radius)
        touch = .false.
        do first_vertex = 1, 3
            first_panel_parameters(:, first_vertex) = &
                parameters(:, triangles(first_vertex, first_panel))
            second_panel_parameters(:, first_vertex) = &
                parameters(:, triangles(first_vertex, second_panel))
        end do
        do first_vertex = 1, 3
            call evaluate_torus_curved_panel( &
                first_panel_parameters, major_radius, minor_radius, &
                first_reference(1, first_vertex), &
                first_reference(2, first_vertex), first_point, tangent_xi, &
                tangent_eta, jacobian, status)
            if (status /= 0) return
            do second_vertex = 1, 3
                call evaluate_torus_curved_panel( &
                    second_panel_parameters, major_radius, minor_radius, &
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
    end function torus_reference_children_touch

    pure subroutine map_refined_torus_point_to_parent( &
            refined_panel_parameters, parent_panel_parameters, xi, eta, &
            parent_xi, parent_eta, status)
        real(dp), intent(in) :: refined_panel_parameters(2, 3)
        real(dp), intent(in) :: parent_panel_parameters(2, 3), xi, eta
        real(dp), intent(out) :: parent_xi, parent_eta
        integer, intent(out) :: status

        real(dp) :: determinant, parent_parameters(2, 3)
        real(dp) :: refined_parameters(2, 3), right_hand_side(2), target(2)

        parent_xi = 0.0_dp
        parent_eta = 0.0_dp
        status = 1
        parent_parameters = parent_panel_parameters
        refined_parameters = refined_panel_parameters
        call unwrap_torus_parameter( &
            parent_parameters(:, 1), parent_parameters(:, 2))
        call unwrap_torus_parameter( &
            parent_parameters(:, 1), parent_parameters(:, 3))
        call unwrap_torus_parameter( &
            refined_parameters(:, 1), refined_parameters(:, 2))
        call unwrap_torus_parameter( &
            refined_parameters(:, 1), refined_parameters(:, 3))
        target = refined_parameters(:, 1) + &
            xi*(refined_parameters(:, 2) - refined_parameters(:, 1)) + &
            eta*(refined_parameters(:, 3) - refined_parameters(:, 1))
        call unwrap_torus_parameter(parent_parameters(:, 1), target)
        right_hand_side = target - parent_parameters(:, 1)
        determinant = &
            (parent_parameters(1, 2) - parent_parameters(1, 1))* &
            (parent_parameters(2, 3) - parent_parameters(2, 1)) - &
            (parent_parameters(1, 3) - parent_parameters(1, 1))* &
            (parent_parameters(2, 2) - parent_parameters(2, 1))
        if (abs(determinant) <= tiny(1.0_dp)) return
        parent_xi = ( &
            right_hand_side(1)* &
            (parent_parameters(2, 3) - parent_parameters(2, 1)) - &
            right_hand_side(2)* &
            (parent_parameters(1, 3) - parent_parameters(1, 1)))/determinant
        parent_eta = ( &
            (parent_parameters(1, 2) - parent_parameters(1, 1))* &
            right_hand_side(2) - &
            (parent_parameters(2, 2) - parent_parameters(2, 1))* &
            right_hand_side(1))/determinant
        if (parent_xi < -256.0_dp*epsilon(1.0_dp)) return
        if (parent_eta < -256.0_dp*epsilon(1.0_dp)) return
        if (parent_xi + parent_eta > 1.0_dp + &
            256.0_dp*epsilon(1.0_dp)) return
        parent_xi = max(0.0_dp, parent_xi)
        parent_eta = max(0.0_dp, parent_eta)
        status = 0
    end subroutine map_refined_torus_point_to_parent

    pure subroutine map_refined_torus_point_to_parent_jvp( &
            refined_panel_parameters, parent_panel_parameters, xi, eta, &
            refined_panel_parameters_dot, parent_panel_parameters_dot, &
            xi_dot, eta_dot, parent_xi, parent_eta, parent_xi_dot, &
            parent_eta_dot, status)
        real(dp), intent(in) :: refined_panel_parameters(2, 3)
        real(dp), intent(in) :: parent_panel_parameters(2, 3), xi, eta
        real(dp), intent(in) :: refined_panel_parameters_dot(2, 3)
        real(dp), intent(in) :: parent_panel_parameters_dot(2, 3)
        real(dp), intent(in) :: xi_dot, eta_dot
        real(dp), intent(out) :: parent_xi, parent_eta
        real(dp), intent(out) :: parent_xi_dot, parent_eta_dot
        integer, intent(out) :: status

        real(dp) :: determinant, determinant_dot
        real(dp) :: parent_parameters(2, 3), parent_parameters_dot(2, 3)
        real(dp) :: refined_parameters(2, 3)
        real(dp) :: refined_parameters_dot_local(2, 3)
        real(dp) :: target(2), target_dot(2), right_hand_side(2)
        real(dp) :: right_hand_side_dot(2), first_edge(2), second_edge(2)
        real(dp) :: first_edge_dot(2), second_edge_dot(2)
        real(dp) :: numerator_xi, numerator_eta
        real(dp) :: numerator_xi_dot, numerator_eta_dot

        parent_xi = 0.0_dp
        parent_eta = 0.0_dp
        parent_xi_dot = 0.0_dp
        parent_eta_dot = 0.0_dp
        status = 1
        parent_parameters = parent_panel_parameters
        parent_parameters_dot = parent_panel_parameters_dot
        refined_parameters = refined_panel_parameters
        refined_parameters_dot_local = refined_panel_parameters_dot
        call unwrap_torus_parameter( &
            parent_parameters(:, 1), parent_parameters(:, 2))
        call unwrap_torus_parameter( &
            parent_parameters(:, 1), parent_parameters(:, 3))
        call unwrap_torus_parameter( &
            refined_parameters(:, 1), refined_parameters(:, 2))
        call unwrap_torus_parameter( &
            refined_parameters(:, 1), refined_parameters(:, 3))
        target = refined_parameters(:, 1) + &
            xi*(refined_parameters(:, 2) - refined_parameters(:, 1)) + &
            eta*(refined_parameters(:, 3) - refined_parameters(:, 1))
        target_dot = refined_parameters_dot_local(:, 1) + &
            xi_dot*(refined_parameters(:, 2) - refined_parameters(:, 1)) + &
            eta_dot*(refined_parameters(:, 3) - refined_parameters(:, 1)) + &
            xi*(refined_parameters_dot_local(:, 2) - &
            refined_parameters_dot_local(:, 1)) + &
            eta*(refined_parameters_dot_local(:, 3) - &
            refined_parameters_dot_local(:, 1))
        call unwrap_torus_parameter(parent_parameters(:, 1), target)
        right_hand_side = target - parent_parameters(:, 1)
        right_hand_side_dot = target_dot - parent_parameters_dot(:, 1)
        first_edge = parent_parameters(:, 2) - parent_parameters(:, 1)
        second_edge = parent_parameters(:, 3) - parent_parameters(:, 1)
        first_edge_dot = parent_parameters_dot(:, 2) - &
            parent_parameters_dot(:, 1)
        second_edge_dot = parent_parameters_dot(:, 3) - &
            parent_parameters_dot(:, 1)
        determinant = first_edge(1)*second_edge(2) - &
            first_edge(2)*second_edge(1)
        determinant_dot = first_edge_dot(1)*second_edge(2) + &
            first_edge(1)*second_edge_dot(2) - &
            first_edge_dot(2)*second_edge(1) - &
            first_edge(2)*second_edge_dot(1)
        if (abs(determinant) <= tiny(1.0_dp)) return
        numerator_xi = right_hand_side(1)*second_edge(2) - &
            right_hand_side(2)*second_edge(1)
        numerator_eta = first_edge(1)*right_hand_side(2) - &
            first_edge(2)*right_hand_side(1)
        numerator_xi_dot = right_hand_side_dot(1)*second_edge(2) + &
            right_hand_side(1)*second_edge_dot(2) - &
            right_hand_side_dot(2)*second_edge(1) - &
            right_hand_side(2)*second_edge_dot(1)
        numerator_eta_dot = first_edge_dot(1)*right_hand_side(2) + &
            first_edge(1)*right_hand_side_dot(2) - &
            first_edge_dot(2)*right_hand_side(1) - &
            first_edge(2)*right_hand_side_dot(1)
        parent_xi = numerator_xi/determinant
        parent_eta = numerator_eta/determinant
        parent_xi_dot = (numerator_xi_dot - parent_xi*determinant_dot)/ &
            determinant
        parent_eta_dot = (numerator_eta_dot - parent_eta*determinant_dot)/ &
            determinant
        if (parent_xi < -256.0_dp*epsilon(1.0_dp)) return
        if (parent_eta < -256.0_dp*epsilon(1.0_dp)) return
        if (parent_xi + parent_eta > 1.0_dp + &
            256.0_dp*epsilon(1.0_dp)) return
        parent_xi = max(0.0_dp, parent_xi)
        parent_eta = max(0.0_dp, parent_eta)
        status = 0
    end subroutine map_refined_torus_point_to_parent_jvp

    pure subroutine map_refined_torus_point_to_parent_vjp( &
            refined_panel_parameters, parent_panel_parameters, xi, eta, &
            parent_xi_bar, parent_eta_bar, refined_panel_parameters_bar, &
            parent_panel_parameters_bar, status)
        real(dp), intent(in) :: refined_panel_parameters(2, 3)
        real(dp), intent(in) :: parent_panel_parameters(2, 3), xi, eta
        real(dp), intent(in) :: parent_xi_bar, parent_eta_bar
        real(dp), intent(out) :: refined_panel_parameters_bar(2, 3)
        real(dp), intent(out) :: parent_panel_parameters_bar(2, 3)
        integer, intent(out) :: status

        real(dp) :: determinant, parent_parameters(2, 3)
        real(dp) :: refined_parameters(2, 3), right_hand_side(2)
        real(dp) :: target(2), first_edge(2), second_edge(2)
        real(dp) :: right_hand_side_bar(2), first_edge_bar(2)
        real(dp) :: second_edge_bar(2), target_bar(2)
        real(dp) :: parent_xi, parent_eta

        refined_panel_parameters_bar = 0.0_dp
        parent_panel_parameters_bar = 0.0_dp
        status = 1
        call map_refined_torus_point_to_parent( &
            refined_panel_parameters, parent_panel_parameters, xi, eta, &
            parent_xi, parent_eta, status)
        if (status /= 0) return
        parent_parameters = parent_panel_parameters
        refined_parameters = refined_panel_parameters
        call unwrap_torus_parameter( &
            parent_parameters(:, 1), parent_parameters(:, 2))
        call unwrap_torus_parameter( &
            parent_parameters(:, 1), parent_parameters(:, 3))
        call unwrap_torus_parameter( &
            refined_parameters(:, 1), refined_parameters(:, 2))
        call unwrap_torus_parameter( &
            refined_parameters(:, 1), refined_parameters(:, 3))
        target = refined_parameters(:, 1) + &
            xi*(refined_parameters(:, 2) - refined_parameters(:, 1)) + &
            eta*(refined_parameters(:, 3) - refined_parameters(:, 1))
        call unwrap_torus_parameter(parent_parameters(:, 1), target)
        right_hand_side = target - parent_parameters(:, 1)
        first_edge = parent_parameters(:, 2) - parent_parameters(:, 1)
        second_edge = parent_parameters(:, 3) - parent_parameters(:, 1)
        determinant = first_edge(1)*second_edge(2) - &
            first_edge(2)*second_edge(1)
        if (abs(determinant) <= tiny(1.0_dp)) return
        right_hand_side_bar = [ &
            (second_edge(2)*parent_xi_bar - second_edge(1)*parent_eta_bar)/ &
            determinant, &
            (-first_edge(2)*parent_xi_bar + first_edge(1)*parent_eta_bar)/ &
            determinant]
        first_edge_bar = -right_hand_side_bar*parent_xi
        second_edge_bar = -right_hand_side_bar*parent_eta
        target_bar = right_hand_side_bar
        parent_panel_parameters_bar(:, 1) = -right_hand_side_bar - &
            first_edge_bar - second_edge_bar
        parent_panel_parameters_bar(:, 2) = first_edge_bar
        parent_panel_parameters_bar(:, 3) = second_edge_bar
        refined_panel_parameters_bar(:, 1) = &
            (1.0_dp - xi - eta)*target_bar
        refined_panel_parameters_bar(:, 2) = xi*target_bar
        refined_panel_parameters_bar(:, 3) = eta*target_bar
        status = 0
    end subroutine map_refined_torus_point_to_parent_vjp

    pure subroutine unwrap_torus_parameter(reference, value)
        real(dp), intent(in) :: reference(2)
        real(dp), intent(inout) :: value(2)

        integer :: coordinate

        do coordinate = 1, 2
            do while (value(coordinate) - reference(coordinate) > &
                    acos(-1.0_dp))
                value(coordinate) = value(coordinate) - 2.0_dp*acos(-1.0_dp)
            end do
            do while (value(coordinate) - reference(coordinate) < &
                    -acos(-1.0_dp))
                value(coordinate) = value(coordinate) + 2.0_dp*acos(-1.0_dp)
            end do
        end do
    end subroutine unwrap_torus_parameter

    pure function torus_unit_normal(point, major_radius) result(normal)
        real(dp), intent(in) :: point(3), major_radius
        real(dp) :: normal(3)
        real(dp) :: cylindrical_radius

        cylindrical_radius = sqrt(point(1)**2 + point(2)**2)
        normal = [ &
            (cylindrical_radius - major_radius)*point(1)/cylindrical_radius, &
            (cylindrical_radius - major_radius)*point(2)/cylindrical_radius, &
            point(3)]
        normal = normal/norm2(normal)
    end function torus_unit_normal

    pure subroutine torus_unit_normal_jvp( &
            point, major_radius, point_dot, major_radius_dot, normal_dot)
        real(dp), intent(in) :: point(3), major_radius, point_dot(3)
        real(dp), intent(in) :: major_radius_dot
        real(dp), intent(out) :: normal_dot(3)

        real(dp) :: cylindrical_radius, cylindrical_radius_dot
        real(dp) :: radial_difference, radial_difference_dot
        real(dp) :: raw(3), raw_dot(3), normal(3), raw_norm

        cylindrical_radius = sqrt(point(1)**2 + point(2)**2)
        if (cylindrical_radius <= tiny(1.0_dp)) then
            normal_dot = 0.0_dp
            return
        end if
        cylindrical_radius_dot = (point(1)*point_dot(1) + &
            point(2)*point_dot(2))/cylindrical_radius
        radial_difference = cylindrical_radius - major_radius
        radial_difference_dot = cylindrical_radius_dot - major_radius_dot
        raw = [ &
            radial_difference*point(1)/cylindrical_radius, &
            radial_difference*point(2)/cylindrical_radius, point(3)]
        raw_dot = [ &
            radial_difference_dot*point(1)/cylindrical_radius + &
            radial_difference*point_dot(1)/cylindrical_radius - &
            radial_difference*point(1)*cylindrical_radius_dot/ &
            cylindrical_radius**2, &
            radial_difference_dot*point(2)/cylindrical_radius + &
            radial_difference*point_dot(2)/cylindrical_radius - &
            radial_difference*point(2)*cylindrical_radius_dot/ &
            cylindrical_radius**2, point_dot(3)]
        raw_norm = norm2(raw)
        if (raw_norm <= tiny(1.0_dp)) then
            normal_dot = 0.0_dp
            return
        end if
        normal = raw/raw_norm
        normal_dot = (raw_dot - normal*dot_product(normal, raw_dot))/raw_norm
    end subroutine torus_unit_normal_jvp

    pure subroutine torus_unit_normal_vjp( &
            point, major_radius, normal_bar, point_bar, major_radius_bar)
        real(dp), intent(in) :: point(3), major_radius, normal_bar(3)
        real(dp), intent(out) :: point_bar(3), major_radius_bar

        real(dp) :: cylindrical_radius, radial_difference, raw_norm
        real(dp) :: raw(3), normal(3), raw_bar(3)
        real(dp) :: radial_difference_bar, cylindrical_radius_bar

        point_bar = 0.0_dp
        major_radius_bar = 0.0_dp
        cylindrical_radius = sqrt(point(1)**2 + point(2)**2)
        if (cylindrical_radius <= tiny(1.0_dp)) return
        radial_difference = cylindrical_radius - major_radius
        raw = [ &
            radial_difference*point(1)/cylindrical_radius, &
            radial_difference*point(2)/cylindrical_radius, point(3)]
        raw_norm = norm2(raw)
        if (raw_norm <= tiny(1.0_dp)) return
        normal = raw/raw_norm
        raw_bar = (normal_bar - normal*dot_product(normal, normal_bar))/raw_norm
        point_bar(3) = raw_bar(3)
        radial_difference_bar = raw_bar(1)*point(1)/cylindrical_radius + &
            raw_bar(2)*point(2)/cylindrical_radius
        point_bar(1) = point_bar(1) + &
            raw_bar(1)*radial_difference/cylindrical_radius
        point_bar(2) = point_bar(2) + &
            raw_bar(2)*radial_difference/cylindrical_radius
        cylindrical_radius_bar = -radial_difference*raw_bar(1)*point(1)/ &
            cylindrical_radius**2 - &
            radial_difference*raw_bar(2)*point(2)/cylindrical_radius**2
        point_bar(1) = point_bar(1) + &
            cylindrical_radius_bar*point(1)/cylindrical_radius
        point_bar(2) = point_bar(2) + &
            cylindrical_radius_bar*point(2)/cylindrical_radius
        major_radius_bar = -radial_difference_bar
    end subroutine torus_unit_normal_vjp

    pure function real_cross_product(first, second) result(product)
        real(dp), intent(in) :: first(3), second(3)
        real(dp) :: product(3)

        product = [ &
            first(2)*second(3) - first(3)*second(2), &
            first(3)*second(1) - first(1)*second(3), &
            first(1)*second(2) - first(2)*second(1)]
    end function real_cross_product

    pure function complex_cross_product(first, second) result(product)
        complex(dp), intent(in) :: first(3), second(3)
        complex(dp) :: product(3)

        product = [ &
            first(2)*second(3) - first(3)*second(2), &
            first(3)*second(1) - first(1)*second(3), &
            first(1)*second(2) - first(2)*second(1)]
    end function complex_cross_product

    pure function torus_boundary_green( &
            wave_number, distance, decaying_kernel) result(value)
        real(dp), intent(in) :: wave_number, distance
        logical, intent(in) :: decaying_kernel
        complex(dp) :: value

        if (decaying_kernel) then
            value = cmplx(exp(-wave_number*distance)/ &
                (4.0_dp*acos(-1.0_dp)*distance), 0.0_dp, dp)
        else
            value = exp(cmplx(0.0_dp, wave_number*distance, dp))/ &
                (4.0_dp*acos(-1.0_dp)*distance)
        end if
    end function torus_boundary_green

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

end module fortfem_maxwell_torus_curved_rwg
