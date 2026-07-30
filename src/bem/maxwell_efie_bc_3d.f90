module fortfem_maxwell_efie_bc_3d
    !! Electric-field operators with Buffa--Christiansen trial and test traces.
    use fortfem_helmholtz_galerkin_3d, only: &
        assemble_helmholtz_single_layer_p0_adaptive_3d
    use fortfem_kinds, only: dp
    use fortfem_maxwell_bc_surface, only: build_maxwell_bc_transformation
    use fortfem_maxwell_efie_rwg_3d, only: &
        assemble_maxwell_rwg_potential_operators_3d
    use fortfem_maxwell_localized_rwg_surface, only: &
        evaluate_maxwell_localized_rwg_basis
    use fortfem_maxwell_rwg_surface, only: &
        build_maxwell_rwg_surface_space, evaluate_maxwell_rwg_basis
    use fortfem_triangle_duffy_quadrature, only: triangle_duffy_quadrature
    implicit none
    private

    public :: assemble_maxwell_bc_scalar_potential_3d
    public :: assemble_maxwell_bc_potential_operators_3d
    public :: assemble_maxwell_efie_bc_3d
    public :: assemble_maxwell_efie_bc_imaginary_3d
    public :: build_maxwell_bc_panel_divergence
    public :: build_maxwell_bc_to_refined_rwg

contains

    subroutine assemble_maxwell_efie_bc_imaginary_3d( &
            vertices, triangles, decay_rate, impedance, quadrature_degree, &
            tolerance, max_depth, matrix, status)
        real(dp), intent(in) :: vertices(:, :), decay_rate, impedance
        integer, intent(in) :: triangles(:, :), quadrature_degree, max_depth
        real(dp), intent(in) :: tolerance
        complex(dp), allocatable, intent(out) :: matrix(:, :)
        integer, intent(out) :: status

        complex(dp), allocatable :: refined_scalar(:, :), refined_vector(:, :)
        integer, allocatable :: refined_triangles(:, :)
        real(dp), allocatable :: correction_scalar(:, :)
        real(dp), allocatable :: correction_vector(:, :)
        real(dp), allocatable :: refined_vertices(:, :), transformation(:, :)
        complex(dp), allocatable :: complex_transformation(:, :)
        complex(dp), allocatable :: scalar_potential(:, :)
        complex(dp), allocatable :: vector_potential(:, :)

        status = 1
        if (decay_rate <= 0.0_dp .or. impedance <= 0.0_dp) return
        call build_maxwell_bc_to_refined_rwg( &
            vertices, triangles, refined_vertices, refined_triangles, &
            transformation, status)
        if (status /= 0) return
        call assemble_maxwell_rwg_potential_operators_3d( &
            refined_vertices, refined_triangles, 0.0_dp, quadrature_degree, &
            tolerance, max_depth, refined_vector, refined_scalar, status)
        if (status /= 0) return
        call assemble_modified_kernel_correction( &
            refined_vertices, refined_triangles, decay_rate, quadrature_degree, &
            correction_vector, correction_scalar, status)
        if (status /= 0) return
        refined_vector = refined_vector + &
            cmplx(correction_vector, 0.0_dp, dp)
        refined_scalar = refined_scalar + &
            cmplx(correction_scalar, 0.0_dp, dp)
        complex_transformation = cmplx(transformation, 0.0_dp, dp)
        vector_potential = matmul( &
            transpose(complex_transformation), &
            matmul(refined_vector, complex_transformation))
        scalar_potential = matmul( &
            transpose(complex_transformation), &
            matmul(refined_scalar, complex_transformation))
        matrix = -impedance*( &
            decay_rate*vector_potential + scalar_potential/decay_rate)
        status = 0
    end subroutine assemble_maxwell_efie_bc_imaginary_3d

    subroutine assemble_modified_kernel_correction( &
            vertices, triangles, decay_rate, quadrature_degree, vector_matrix, &
            scalar_matrix, status)
        real(dp), intent(in) :: vertices(:, :), decay_rate
        integer, intent(in) :: triangles(:, :), quadrature_degree
        real(dp), allocatable, intent(out) :: vector_matrix(:, :)
        real(dp), allocatable, intent(out) :: scalar_matrix(:, :)
        integer, intent(out) :: status

        integer, allocatable :: edge_panels(:, :), edge_vertices(:, :)
        real(dp), allocatable :: eta(:), weights(:), xi(:)
        real(dp) :: first_divergence, first_jacobian, first_point(3)
        real(dp) :: first_value(3), kernel, radius, second_divergence
        real(dp) :: second_jacobian, second_point(3), second_value(3)
        integer :: first, first_node, first_panel, first_slot
        integer :: second, second_node, second_panel, second_slot

        call build_maxwell_rwg_surface_space( &
            vertices, triangles, edge_vertices, edge_panels, status)
        if (status /= 0) return
        call triangle_duffy_quadrature( &
            quadrature_degree, xi, eta, weights, status)
        if (status /= 0) return
        allocate( &
            vector_matrix(size(edge_vertices, 2), size(edge_vertices, 2)), &
            scalar_matrix(size(edge_vertices, 2), size(edge_vertices, 2)))
        vector_matrix = 0.0_dp
        scalar_matrix = 0.0_dp
        do first = 1, size(edge_vertices, 2)
            do second = 1, first
                do first_slot = 1, 2
                    first_panel = edge_panels(first_slot, first)
                    first_jacobian = 2.0_dp*triangle_area( &
                        vertices(:, triangles(:, first_panel)))
                    do second_slot = 1, 2
                        second_panel = edge_panels(second_slot, second)
                        second_jacobian = 2.0_dp*triangle_area( &
                            vertices(:, triangles(:, second_panel)))
                        do first_node = 1, size(weights)
                            first_point = triangle_point( &
                                vertices(:, triangles(:, first_panel)), &
                                xi(first_node), eta(first_node))
                            call evaluate_maxwell_rwg_basis( &
                                vertices, triangles, edge_vertices, edge_panels, &
                                first, first_panel, first_point, first_value, &
                                first_divergence, status)
                            if (status /= 0) return
                            do second_node = 1, size(weights)
                                second_point = triangle_point( &
                                    vertices(:, triangles(:, second_panel)), &
                                    xi(second_node), eta(second_node))
                                call evaluate_maxwell_rwg_basis( &
                                    vertices, triangles, edge_vertices, &
                                    edge_panels, second, second_panel, &
                                    second_point, second_value, &
                                    second_divergence, status)
                                if (status /= 0) return
                                radius = norm2(first_point - second_point)
                                kernel = modified_kernel_correction( &
                                    decay_rate, radius)
                                vector_matrix(first, second) = &
                                    vector_matrix(first, second) + &
                                    first_jacobian*second_jacobian* &
                                    weights(first_node)*weights(second_node)* &
                                    kernel*dot_product(first_value, second_value)
                                scalar_matrix(first, second) = &
                                    scalar_matrix(first, second) + &
                                    first_jacobian*second_jacobian* &
                                    weights(first_node)*weights(second_node)* &
                                    kernel*first_divergence*second_divergence
                            end do
                        end do
                    end do
                end do
                vector_matrix(second, first) = vector_matrix(first, second)
                scalar_matrix(second, first) = scalar_matrix(first, second)
            end do
        end do
        status = 0
    end subroutine assemble_modified_kernel_correction

    subroutine assemble_maxwell_efie_bc_3d( &
            vertices, triangles, wave_number, impedance, quadrature_degree, &
            tolerance, max_depth, matrix, status)
        real(dp), intent(in) :: vertices(:, :), wave_number, impedance
        integer, intent(in) :: triangles(:, :), quadrature_degree, max_depth
        real(dp), intent(in) :: tolerance
        complex(dp), allocatable, intent(out) :: matrix(:, :)
        integer, intent(out) :: status

        complex(dp), allocatable :: scalar_potential(:, :)
        complex(dp), allocatable :: vector_potential(:, :)

        status = 1
        if (wave_number <= 0.0_dp .or. impedance <= 0.0_dp) return
        call assemble_maxwell_bc_potential_operators_3d( &
            vertices, triangles, wave_number, quadrature_degree, tolerance, &
            max_depth, vector_potential, scalar_potential, status)
        if (status /= 0) return
        matrix = cmplx(0.0_dp, wave_number*impedance, dp)*vector_potential - &
            cmplx(0.0_dp, impedance/wave_number, dp)*scalar_potential
        status = 0
    end subroutine assemble_maxwell_efie_bc_3d

    subroutine assemble_maxwell_bc_potential_operators_3d( &
            vertices, triangles, wave_number, quadrature_degree, tolerance, &
            max_depth, vector_potential, scalar_potential, status)
        real(dp), intent(in) :: vertices(:, :), wave_number, tolerance
        integer, intent(in) :: triangles(:, :), quadrature_degree, max_depth
        complex(dp), allocatable, intent(out) :: vector_potential(:, :)
        complex(dp), allocatable, intent(out) :: scalar_potential(:, :)
        integer, intent(out) :: status

        complex(dp), allocatable :: refined_scalar(:, :)
        complex(dp), allocatable :: refined_vector(:, :)
        integer, allocatable :: refined_triangles(:, :)
        real(dp), allocatable :: refined_vertices(:, :), transformation(:, :)
        complex(dp), allocatable :: complex_transformation(:, :)

        call build_maxwell_bc_to_refined_rwg( &
            vertices, triangles, refined_vertices, refined_triangles, &
            transformation, status)
        if (status /= 0) return
        call assemble_maxwell_rwg_potential_operators_3d( &
            refined_vertices, refined_triangles, wave_number, &
            quadrature_degree, tolerance, max_depth, refined_vector, &
            refined_scalar, status)
        if (status /= 0) return
        complex_transformation = cmplx(transformation, 0.0_dp, dp)
        vector_potential = matmul( &
            transpose(complex_transformation), &
            matmul(refined_vector, complex_transformation))
        scalar_potential = matmul( &
            transpose(complex_transformation), &
            matmul(refined_scalar, complex_transformation))
        status = 0
    end subroutine assemble_maxwell_bc_potential_operators_3d

    subroutine build_maxwell_bc_to_refined_rwg( &
            vertices, triangles, refined_vertices, refined_triangles, &
            transformation, status, sphere_radius)
        real(dp), intent(in) :: vertices(:, :)
        integer, intent(in) :: triangles(:, :)
        real(dp), allocatable, intent(out) :: refined_vertices(:, :)
        integer, allocatable, intent(out) :: refined_triangles(:, :)
        real(dp), allocatable, intent(out) :: transformation(:, :)
        integer, intent(out) :: status
        real(dp), optional, intent(in) :: sphere_radius

        integer, allocatable :: edge_panels(:, :), edge_vertices(:, :)
        real(dp), allocatable :: localized_transformation(:, :)
        real(dp) :: divergence, global_value(3), local_value(3), point(3)
        real(dp) :: orientation, candidate
        integer :: basis, edge, local_edge, panel, row, slot

        if (present(sphere_radius)) then
            call build_maxwell_bc_transformation( &
                vertices, triangles, refined_vertices, refined_triangles, &
                localized_transformation, status, sphere_radius=sphere_radius)
        else
            call build_maxwell_bc_transformation( &
                vertices, triangles, refined_vertices, refined_triangles, &
                localized_transformation, status)
        end if
        if (status /= 0) return
        call build_maxwell_rwg_surface_space( &
            refined_vertices, refined_triangles, edge_vertices, edge_panels, &
            status)
        if (status /= 0) return
        allocate(transformation( &
            size(edge_vertices, 2), size(localized_transformation, 2)))
        transformation = 0.0_dp
        do edge = 1, size(edge_vertices, 2)
            do slot = 1, 2
                panel = edge_panels(slot, edge)
                local_edge = find_local_edge( &
                    refined_triangles(:, panel), edge_vertices(:, edge))
                if (local_edge == 0) then
                    status = 2
                    return
                end if
                point = sum( &
                    refined_vertices(:, refined_triangles(:, panel)), &
                    dim=2)/3.0_dp
                call evaluate_maxwell_rwg_basis( &
                    refined_vertices, refined_triangles, edge_vertices, &
                    edge_panels, edge, panel, point, global_value, divergence, &
                    status)
                if (status /= 0) return
                call evaluate_maxwell_localized_rwg_basis( &
                    refined_vertices, refined_triangles, panel, local_edge, &
                    point, local_value, divergence, status)
                if (status /= 0) return
                orientation = sign(1.0_dp, dot_product(global_value, local_value))
                row = 3*(panel - 1) + local_edge
                do basis = 1, size(transformation, 2)
                    candidate = orientation* &
                        localized_transformation(row, basis)
                    if (slot == 1) then
                        transformation(edge, basis) = candidate
                    else if (abs(candidate - transformation(edge, basis)) > &
                            256.0_dp*epsilon(1.0_dp)*max( &
                            1.0_dp, abs(candidate))) then
                        status = 3
                        return
                    end if
                end do
            end do
        end do
        status = 0
    end subroutine build_maxwell_bc_to_refined_rwg

    subroutine assemble_maxwell_bc_scalar_potential_3d( &
            vertices, triangles, wave_number, quadrature_degree, tolerance, &
            max_depth, matrix, status)
        real(dp), intent(in) :: vertices(:, :), wave_number, tolerance
        integer, intent(in) :: triangles(:, :), quadrature_degree, max_depth
        complex(dp), allocatable, intent(out) :: matrix(:, :)
        integer, intent(out) :: status

        complex(dp), allocatable :: panel_operator(:, :)
        real(dp), allocatable :: divergence(:, :)
        real(dp), allocatable :: refined_vertices(:, :), transformation(:, :)
        integer, allocatable :: refined_triangles(:, :)

        status = 1
        if (wave_number < 0.0_dp) return
        call build_maxwell_bc_transformation( &
            vertices, triangles, refined_vertices, refined_triangles, &
            transformation, status)
        if (status /= 0) return
        call build_maxwell_bc_panel_divergence( &
            vertices, triangles, divergence, status)
        if (status /= 0) return
        if (size(divergence, 2) /= size(transformation, 2)) return
        call assemble_helmholtz_single_layer_p0_adaptive_3d( &
            refined_vertices, refined_triangles, wave_number, &
            quadrature_degree, tolerance, max_depth, panel_operator, status)
        if (status /= 0) return
        matrix = matmul( &
            transpose(cmplx(divergence, 0.0_dp, dp)), &
            matmul(panel_operator, cmplx(divergence, 0.0_dp, dp)))
        status = 0
    end subroutine assemble_maxwell_bc_scalar_potential_3d

    subroutine build_maxwell_bc_panel_divergence( &
            vertices, triangles, divergence, status)
        real(dp), intent(in) :: vertices(:, :)
        integer, intent(in) :: triangles(:, :)
        real(dp), allocatable, intent(out) :: divergence(:, :)
        integer, intent(out) :: status

        integer, allocatable :: refined_triangles(:, :)
        real(dp), allocatable :: refined_vertices(:, :)
        real(dp), allocatable :: transformation(:, :)
        real(dp) :: local_divergence, local_value(3), point(3)
        integer :: basis, local_edge, refined_panel, row

        call build_maxwell_bc_transformation( &
            vertices, triangles, refined_vertices, refined_triangles, &
            transformation, status)
        if (status /= 0) return
        allocate(divergence( &
            size(refined_triangles, 2), size(transformation, 2)))
        divergence = 0.0_dp
        do refined_panel = 1, size(refined_triangles, 2)
            point = sum( &
                refined_vertices(:, refined_triangles(:, refined_panel)), &
                dim=2)/3.0_dp
            do local_edge = 1, 3
                call evaluate_maxwell_localized_rwg_basis( &
                    refined_vertices, refined_triangles, refined_panel, &
                    local_edge, point, local_value, local_divergence, status)
                if (status /= 0) return
                row = 3*(refined_panel - 1) + local_edge
                do basis = 1, size(transformation, 2)
                    divergence(refined_panel, basis) = &
                        divergence(refined_panel, basis) + &
                        transformation(row, basis)*local_divergence
                end do
            end do
        end do
        status = 0
    end subroutine build_maxwell_bc_panel_divergence

    pure integer function find_local_edge(element, edge) result(local_edge)
        integer, intent(in) :: element(3), edge(2)

        integer :: local_vertices(2), first, second

        do local_edge = 1, 3
            select case (local_edge)
            case (1)
                local_vertices = [1, 2]
            case (2)
                local_vertices = [3, 1]
            case default
                local_vertices = [2, 3]
            end select
            first = element(local_vertices(1))
            second = element(local_vertices(2))
            if (min(first, second) == edge(1) .and. &
                max(first, second) == edge(2)) return
        end do
        local_edge = 0
    end function find_local_edge

    pure real(dp) function modified_kernel_correction( &
            decay_rate, radius) result(value)
        real(dp), intent(in) :: decay_rate, radius

        if (radius <= 64.0_dp*epsilon(1.0_dp)/decay_rate) then
            value = -decay_rate/(4.0_dp*acos(-1.0_dp))
        else
            value = expm1_stable(-decay_rate*radius)/ &
                (4.0_dp*acos(-1.0_dp)*radius)
        end if
    end function modified_kernel_correction

    pure real(dp) function expm1_stable(value) result(result)
        real(dp), intent(in) :: value

        if (abs(value) < 1.0e-5_dp) then
            result = value*(1.0_dp + value*(0.5_dp + &
                value*(1.0_dp/6.0_dp + value/24.0_dp)))
        else
            result = exp(value) - 1.0_dp
        end if
    end function expm1_stable

    pure function triangle_point(vertices, xi, eta) result(point)
        real(dp), intent(in) :: vertices(3, 3), xi, eta
        real(dp) :: point(3)

        point = vertices(:, 1) + xi*(vertices(:, 2) - vertices(:, 1)) + &
            eta*(vertices(:, 3) - vertices(:, 1))
    end function triangle_point

    pure function triangle_area(vertices) result(area)
        real(dp), intent(in) :: vertices(3, 3)
        real(dp) :: area

        area = 0.5_dp*norm2(cross_product( &
            vertices(:, 2) - vertices(:, 1), &
            vertices(:, 3) - vertices(:, 1)))
    end function triangle_area

    pure function cross_product(first, second) result(product)
        real(dp), intent(in) :: first(3), second(3)
        real(dp) :: product(3)

        product = [ &
            first(2)*second(3) - first(3)*second(2), &
            first(3)*second(1) - first(1)*second(3), &
            first(1)*second(2) - first(2)*second(1)]
    end function cross_product

end module fortfem_maxwell_efie_bc_3d
