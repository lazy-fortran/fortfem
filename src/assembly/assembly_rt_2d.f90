module fortfem_assembly_rt_2d
    use fortfem_kinds, only: dp
    use fortfem_mesh_2d, only: mesh_2d_t
    use fortfem_basis_rt_2d, only: evaluate_rt_basis_div_2d, &
        evaluate_rt_basis_2d_piola
    use fortfem_gauss_quadrature_2d, only: &
        gauss_quadrature_triangle_t, get_gauss_quadrature_triangle
    implicit none
    private

    public :: assemble_rt_div_mass_element
    public :: assemble_rt_div_mass

contains

    subroutine assemble_rt_div_mass_element( &
            mesh, triangle_idx, element_matrix)
        type(mesh_2d_t), intent(in) :: mesh
        integer, intent(in) :: triangle_idx
        real(dp), intent(out) :: element_matrix(3, 3)

        type(gauss_quadrature_triangle_t) :: quadrature
        real(dp) :: basis_values(2, 3), basis_divergences(3)
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
            error stop "RT assembly requires counter-clockwise triangles"
        end if
        area = 0.5_dp * det_jacobian

        quadrature = get_gauss_quadrature_triangle(2)
        element_matrix = 0.0_dp
        do q = 1, quadrature%n_points
            call evaluate_rt_basis_2d_piola(mesh, triangle_idx, &
                quadrature%xi(q), quadrature%eta(q), basis_values)
            call evaluate_rt_basis_div_2d(quadrature%xi(q), &
                quadrature%eta(q), area, basis_divergences)
            physical_weight = det_jacobian * quadrature%weights(q)

            do j = 1, 3
                do i = 1, 3
                    element_matrix(i, j) = element_matrix(i, j) + &
                        physical_weight * ( &
                        basis_divergences(i) * basis_divergences(j) + &
                        dot_product(basis_values(:, i), basis_values(:, j)))
                end do
            end do
        end do
        call quadrature%destroy()
    end subroutine assemble_rt_div_mass_element

    subroutine assemble_rt_div_mass(mesh, matrix)
        type(mesh_2d_t), intent(inout) :: mesh
        real(dp), intent(out) :: matrix(:, :)

        real(dp) :: element_matrix(3, 3)
        integer :: edge_dofs(3), edge_orientations(3)
        integer :: triangle_idx, i, j, global_i, global_j

        if (size(matrix, 1) /= mesh%n_edges .or. &
            size(matrix, 2) /= mesh%n_edges) then
            error stop "RT matrix shape must equal the edge DOF count"
        end if
        if (.not. allocated(mesh%edge_to_dof)) then
            call mesh%build_edge_dof_numbering()
        end if

        matrix = 0.0_dp
        do triangle_idx = 1, mesh%n_triangles
            call assemble_rt_div_mass_element( &
                mesh, triangle_idx, element_matrix)
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
    end subroutine assemble_rt_div_mass

end module fortfem_assembly_rt_2d
