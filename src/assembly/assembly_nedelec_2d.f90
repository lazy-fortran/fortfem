module fortfem_assembly_nedelec_2d
    use fortfem_kinds, only: dp
    use fortfem_mesh_2d, only: mesh_2d_t
    use fortfem_basis_edge_2d, only: evaluate_edge_basis_curl_2d, &
        evaluate_edge_basis_2d_piola
    use fortfem_gauss_quadrature_2d, only: &
        gauss_quadrature_triangle_t, get_gauss_quadrature_triangle
    implicit none
    private

    public :: assemble_nedelec_curl_mass_element
    public :: assemble_nedelec_curl_mass
    public :: assemble_nedelec_weighted_element
    public :: assemble_nedelec_weighted
    public :: assemble_nedelec_axisymmetric_fourier

    abstract interface
        pure function scalar_coefficient_2d(x, y) result(value)
            import :: dp
            real(dp), intent(in) :: x, y
            real(dp) :: value
        end function scalar_coefficient_2d
    end interface

contains

    subroutine assemble_nedelec_curl_mass_element( &
            mesh, triangle_idx, element_matrix)
        type(mesh_2d_t), intent(in) :: mesh
        integer, intent(in) :: triangle_idx
        real(dp), intent(out) :: element_matrix(3, 3)

        call assemble_nedelec_weighted_element(mesh, triangle_idx, &
            unit_coefficient, unit_coefficient, 2, element_matrix)
    end subroutine assemble_nedelec_curl_mass_element

    subroutine assemble_nedelec_weighted_element(mesh, triangle_idx, &
            curl_coefficient, mass_coefficient, quadrature_degree, &
            element_matrix)
        type(mesh_2d_t), intent(in) :: mesh
        integer, intent(in) :: triangle_idx, quadrature_degree
        procedure(scalar_coefficient_2d) :: curl_coefficient
        procedure(scalar_coefficient_2d) :: mass_coefficient
        real(dp), intent(out) :: element_matrix(3, 3)

        type(gauss_quadrature_triangle_t) :: quadrature
        real(dp) :: basis_values(2, 3), basis_curls(3)
        real(dp) :: area, det_jacobian, physical_weight
        real(dp) :: curl_weight, mass_weight, x, y
        real(dp) :: x1, y1, x2, y2, x3, y3, xi, eta
        integer :: i, j, q

        x1 = mesh%vertices(1, mesh%triangles(1, triangle_idx))
        y1 = mesh%vertices(2, mesh%triangles(1, triangle_idx))
        x2 = mesh%vertices(1, mesh%triangles(2, triangle_idx))
        y2 = mesh%vertices(2, mesh%triangles(2, triangle_idx))
        x3 = mesh%vertices(1, mesh%triangles(3, triangle_idx))
        y3 = mesh%vertices(2, mesh%triangles(3, triangle_idx))
        det_jacobian = (x2 - x1) * (y3 - y1) - &
            (x3 - x1) * (y2 - y1)
        if (det_jacobian <= 0.0_dp) then
            error stop "Nedelec assembly requires counter-clockwise triangles"
        end if
        area = 0.5_dp * det_jacobian

        quadrature = get_gauss_quadrature_triangle(quadrature_degree)
        element_matrix = 0.0_dp
        do q = 1, quadrature%n_points
            xi = quadrature%xi(q)
            eta = quadrature%eta(q)
            x = x1 + (x2 - x1) * xi + (x3 - x1) * eta
            y = y1 + (y2 - y1) * xi + (y3 - y1) * eta
            curl_weight = curl_coefficient(x, y)
            mass_weight = mass_coefficient(x, y)
            call evaluate_edge_basis_2d_piola(mesh, triangle_idx, &
                xi, eta, basis_values)
            call evaluate_edge_basis_curl_2d( &
                xi, eta, area, basis_curls)
            physical_weight = det_jacobian * quadrature%weights(q)

            do j = 1, 3
                do i = 1, 3
                    element_matrix(i, j) = element_matrix(i, j) + &
                        physical_weight * ( &
                        curl_weight * basis_curls(i) * basis_curls(j) + &
                        mass_weight * dot_product( &
                        basis_values(:, i), basis_values(:, j)))
                end do
            end do
        end do
        call quadrature%destroy()
    end subroutine assemble_nedelec_weighted_element

    subroutine assemble_nedelec_curl_mass(mesh, matrix)
        type(mesh_2d_t), intent(inout) :: mesh
        real(dp), intent(out) :: matrix(:, :)

        call assemble_nedelec_weighted(mesh, unit_coefficient, &
            unit_coefficient, 2, matrix)
    end subroutine assemble_nedelec_curl_mass

    subroutine assemble_nedelec_weighted(mesh, curl_coefficient, &
            mass_coefficient, quadrature_degree, matrix)
        type(mesh_2d_t), intent(inout) :: mesh
        procedure(scalar_coefficient_2d) :: curl_coefficient
        procedure(scalar_coefficient_2d) :: mass_coefficient
        integer, intent(in) :: quadrature_degree
        real(dp), intent(out) :: matrix(:, :)

        real(dp) :: element_matrix(3, 3)
        integer :: edge_dofs(3), edge_orientations(3)
        integer :: triangle_idx, i, j, global_i, global_j

        if (size(matrix, 1) /= mesh%n_edges .or. &
            size(matrix, 2) /= mesh%n_edges) then
            error stop "Nedelec matrix shape must equal the edge DOF count"
        end if
        if (.not. allocated(mesh%edge_to_dof)) then
            call mesh%build_edge_dof_numbering()
        end if

        matrix = 0.0_dp
        do triangle_idx = 1, mesh%n_triangles
            call assemble_nedelec_weighted_element(mesh, triangle_idx, &
                curl_coefficient, mass_coefficient, quadrature_degree, &
                element_matrix)
            call mesh%get_triangle_edge_dofs( &
                triangle_idx, edge_dofs, edge_orientations)

            do j = 1, 3
                global_j = edge_dofs(j) + 1
                do i = 1, 3
                    global_i = edge_dofs(i) + 1
                    matrix(global_i, global_j) = &
                        matrix(global_i, global_j) + &
                        real(edge_orientations(i) * edge_orientations(j), &
                        dp) * element_matrix(i, j)
                end do
            end do
        end do
    end subroutine assemble_nedelec_weighted

    subroutine assemble_nedelec_axisymmetric_fourier( &
            mesh, fourier_mode, quadrature_degree, matrix)
        type(mesh_2d_t), intent(inout) :: mesh
        integer, intent(in) :: fourier_mode, quadrature_degree
        real(dp), intent(out) :: matrix(:, :)

        call assemble_nedelec_weighted(mesh, radial_curl_weight, &
            fourier_mass_weight, quadrature_degree, matrix)

    contains

        pure real(dp) function radial_curl_weight(x, y) result(value)
            real(dp), intent(in) :: x, y

            if (x <= 0.0_dp) then
                error stop "Axisymmetric assembly requires positive R"
            end if
            value = x
            associate (unused_y => [y])
                if (size(unused_y) /= 1) error stop
            end associate
        end function radial_curl_weight

        pure real(dp) function fourier_mass_weight(x, y) result(value)
            real(dp), intent(in) :: x, y

            if (fourier_mode == 0) then
                value = 0.0_dp
            else
                if (x <= 0.0_dp) then
                    error stop "Axisymmetric assembly requires positive R"
                end if
                value = real(fourier_mode, dp)**2 / x
            end if
            associate (unused_y => [y])
                if (size(unused_y) /= 1) error stop
            end associate
        end function fourier_mass_weight

    end subroutine assemble_nedelec_axisymmetric_fourier

    pure real(dp) function unit_coefficient(x, y) result(value)
        real(dp), intent(in) :: x, y

        value = 1.0_dp
        associate (unused_x => x, unused_y => y)
            if (kind(unused_x) /= dp .or. kind(unused_y) /= dp) error stop
        end associate
    end function unit_coefficient

end module fortfem_assembly_nedelec_2d
