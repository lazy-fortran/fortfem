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

contains

    subroutine assemble_nedelec_curl_mass_element( &
            mesh, triangle_idx, element_matrix)
        type(mesh_2d_t), intent(in) :: mesh
        integer, intent(in) :: triangle_idx
        real(dp), intent(out) :: element_matrix(3, 3)

        type(gauss_quadrature_triangle_t) :: quadrature
        real(dp) :: basis_values(2, 3), basis_curls(3)
        real(dp) :: area, det_jacobian, physical_weight
        real(dp) :: x1, y1, x2, y2, x3, y3
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

        quadrature = get_gauss_quadrature_triangle(2)
        element_matrix = 0.0_dp
        do q = 1, quadrature%n_points
            call evaluate_edge_basis_2d_piola(mesh, triangle_idx, &
                quadrature%xi(q), quadrature%eta(q), basis_values)
            call evaluate_edge_basis_curl_2d(quadrature%xi(q), &
                quadrature%eta(q), area, basis_curls)
            physical_weight = det_jacobian * quadrature%weights(q)

            do j = 1, 3
                do i = 1, 3
                    element_matrix(i, j) = element_matrix(i, j) + &
                        physical_weight * (basis_curls(i) * basis_curls(j) + &
                        dot_product(basis_values(:, i), basis_values(:, j)))
                end do
            end do
        end do
        call quadrature%destroy()
    end subroutine assemble_nedelec_curl_mass_element

end module fortfem_assembly_nedelec_2d
