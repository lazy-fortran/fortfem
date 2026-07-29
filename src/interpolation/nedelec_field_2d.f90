module fortfem_nedelec_field_2d
    use fortfem_basis_edge_2d, only: evaluate_edge_basis_2d_piola, &
        evaluate_edge_basis_curl_2d
    use fortfem_kinds, only: dp
    use fortfem_mesh_2d, only: mesh_2d_t
    implicit none
    private

    public :: evaluate_nedelec_field_2d

contains

    subroutine evaluate_nedelec_field_2d( &
            mesh, triangle_idx, x, y, dofs, value, curl)
        type(mesh_2d_t), intent(in) :: mesh
        integer, intent(in) :: triangle_idx
        real(dp), intent(in) :: x, y
        complex(dp), intent(in) :: dofs(:)
        complex(dp), intent(out) :: value(2)
        complex(dp), intent(out), optional :: curl

        real(dp) :: basis_curls(3), basis_values(2, 3)
        real(dp) :: area, eta, xi
        complex(dp) :: local_dofs(3)

        call validate_nedelec_field(mesh, dofs)
        call physical_to_reference(mesh, triangle_idx, x, y, xi, eta, area)
        call local_nedelec_dofs(mesh, triangle_idx, dofs, local_dofs)
        call evaluate_edge_basis_2d_piola( &
            mesh, triangle_idx, xi, eta, basis_values)
        value = matmul(basis_values, local_dofs)

        if (present(curl)) then
            call evaluate_edge_basis_curl_2d(xi, eta, area, basis_curls)
            curl = dot_product(basis_curls, local_dofs)
        end if
    end subroutine evaluate_nedelec_field_2d

    subroutine local_nedelec_dofs(mesh, triangle_idx, dofs, local_dofs)
        type(mesh_2d_t), intent(in) :: mesh
        integer, intent(in) :: triangle_idx
        complex(dp), intent(in) :: dofs(:)
        complex(dp), intent(out) :: local_dofs(3)

        integer :: edge_dofs(3), edge_orientations(3), edge

        call mesh%get_triangle_edge_dofs( &
            triangle_idx, edge_dofs, edge_orientations)
        do edge = 1, 3
            local_dofs(edge) = real(edge_orientations(edge), dp) * &
                dofs(edge_dofs(edge) + 1)
        end do
    end subroutine local_nedelec_dofs

    subroutine physical_to_reference( &
            mesh, triangle_idx, x, y, xi, eta, area)
        type(mesh_2d_t), intent(in) :: mesh
        integer, intent(in) :: triangle_idx
        real(dp), intent(in) :: x, y
        real(dp), intent(out) :: xi, eta, area

        real(dp) :: determinant, jacobian(2, 2), offset(2), vertex_a(2)

        vertex_a = mesh%vertices(:, mesh%triangles(1, triangle_idx))
        jacobian(:, 1) = mesh%vertices(:, &
            mesh%triangles(2, triangle_idx)) - vertex_a
        jacobian(:, 2) = mesh%vertices(:, &
            mesh%triangles(3, triangle_idx)) - vertex_a
        determinant = jacobian(1, 1) * jacobian(2, 2) - &
            jacobian(1, 2) * jacobian(2, 1)
        if (determinant <= 0.0_dp) then
            error stop &
                "Nedelec reconstruction requires counter-clockwise triangles"
        end if

        offset = [x, y] - vertex_a
        xi = (jacobian(2, 2) * offset(1) - &
            jacobian(1, 2) * offset(2)) / determinant
        eta = (-jacobian(2, 1) * offset(1) + &
            jacobian(1, 1) * offset(2)) / determinant
        area = 0.5_dp * determinant
    end subroutine physical_to_reference

    subroutine validate_nedelec_field(mesh, dofs)
        type(mesh_2d_t), intent(in) :: mesh
        complex(dp), intent(in) :: dofs(:)

        if (.not. allocated(mesh%edge_to_dof)) then
            error stop "Nedelec field evaluation requires edge DOF numbering"
        end if
        if (size(dofs) /= mesh%n_edges) then
            error stop "Nedelec field coefficient array has the wrong size"
        end if
    end subroutine validate_nedelec_field

end module fortfem_nedelec_field_2d
