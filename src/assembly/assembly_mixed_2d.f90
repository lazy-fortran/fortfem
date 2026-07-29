module fortfem_assembly_mixed_2d
    use fortfem_basis_edge_2d, only: evaluate_edge_basis_2d_piola
    use fortfem_basis_rt_2d, only: evaluate_rt_basis_2d_piola
    use fortfem_gauss_quadrature_2d, only: &
        gauss_quadrature_triangle_t, get_gauss_quadrature_triangle
    use fortfem_kinds, only: dp
    use fortfem_mesh_2d, only: mesh_2d_t
    use fortsparse, only: csc_from_triplet, csc_t, fortsparse_status_t
    implicit none
    private

    public :: assemble_nedelec_rt_mass_element
    public :: assemble_nedelec_rt_mass_csc

contains

    subroutine assemble_nedelec_rt_mass_element( &
            mesh, triangle_idx, quadrature_degree, element_matrix)
        type(mesh_2d_t), intent(in) :: mesh
        integer, intent(in) :: triangle_idx, quadrature_degree
        real(dp), intent(out) :: element_matrix(3, 3)

        type(gauss_quadrature_triangle_t) :: quadrature
        real(dp) :: nedelec_values(2, 3), rt_values(2, 3)
        real(dp) :: det_jacobian, physical_weight
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
            error stop "Mixed assembly requires counter-clockwise triangles"
        end if

        quadrature = get_gauss_quadrature_triangle(quadrature_degree)
        element_matrix = 0.0_dp
        do q = 1, quadrature%n_points
            xi = quadrature%xi(q)
            eta = quadrature%eta(q)
            call evaluate_edge_basis_2d_piola( &
                mesh, triangle_idx, xi, eta, nedelec_values)
            call evaluate_rt_basis_2d_piola( &
                mesh, triangle_idx, xi, eta, rt_values)
            physical_weight = det_jacobian * quadrature%weights(q)
            do j = 1, 3
                do i = 1, 3
                    element_matrix(i, j) = element_matrix(i, j) + &
                        physical_weight * dot_product( &
                        nedelec_values(:, i), rt_values(:, j))
                end do
            end do
        end do
        call quadrature%destroy()
    end subroutine assemble_nedelec_rt_mass_element

    subroutine assemble_nedelec_rt_mass_csc( &
            mesh, quadrature_degree, matrix, status)
        type(mesh_2d_t), intent(inout) :: mesh
        integer, intent(in) :: quadrature_degree
        type(csc_t), intent(out) :: matrix
        type(fortsparse_status_t), intent(out) :: status

        integer, allocatable :: columns(:), rows(:)
        real(dp), allocatable :: values(:)
        real(dp) :: element_matrix(3, 3)
        integer :: edge_dofs(3), edge_orientations(3)
        integer :: entry, i, j, triangle_idx

        if (.not. allocated(mesh%edge_to_dof)) then
            call mesh%build_edge_dof_numbering()
        end if
        allocate(rows(9 * mesh%n_triangles))
        allocate(columns(9 * mesh%n_triangles))
        allocate(values(9 * mesh%n_triangles))

        entry = 0
        do triangle_idx = 1, mesh%n_triangles
            call assemble_nedelec_rt_mass_element( &
                mesh, triangle_idx, quadrature_degree, element_matrix)
            call mesh%get_triangle_edge_dofs( &
                triangle_idx, edge_dofs, edge_orientations)
            do j = 1, 3
                do i = 1, 3
                    entry = entry + 1
                    rows(entry) = edge_dofs(i) + 1
                    columns(entry) = edge_dofs(j) + 1
                    values(entry) = &
                        real(edge_orientations(i) * edge_orientations(j), &
                        dp) * element_matrix(i, j)
                end do
            end do
        end do
        call csc_from_triplet(mesh%n_edges, mesh%n_edges, &
            rows, columns, values, matrix, status)
    end subroutine assemble_nedelec_rt_mass_csc

end module fortfem_assembly_mixed_2d
