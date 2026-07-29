module fortfem_rt_field_2d
    use fortfem_kinds, only: dp
    use fortfem_mesh_2d, only: mesh_2d_t
    use fortfem_basis_rt_2d, only: evaluate_rt_basis_2d_piola
    use fortfem_gauss_quadrature_2d, only: &
        gauss_quadrature_triangle_t, get_gauss_quadrature_triangle
    implicit none
    private

    public :: evaluate_rt_field_2d
    public :: reconstruct_axisymmetric_fourier_toroidal
    public :: rt_l2_norm

contains

    subroutine evaluate_rt_field_2d( &
            mesh, triangle_idx, x, y, dofs, value, divergence)
        type(mesh_2d_t), intent(in) :: mesh
        integer, intent(in) :: triangle_idx
        real(dp), intent(in) :: x, y
        complex(dp), intent(in) :: dofs(:)
        complex(dp), intent(out) :: value(2)
        complex(dp), intent(out), optional :: divergence

        real(dp) :: basis_values(2, 3), xi, eta, area
        complex(dp) :: local_dofs(3)

        call validate_rt_field(mesh, dofs)
        call physical_to_reference(mesh, triangle_idx, x, y, xi, eta, area)
        call local_rt_dofs(mesh, triangle_idx, dofs, local_dofs)
        call evaluate_rt_basis_2d_piola( &
            mesh, triangle_idx, xi, eta, basis_values)

        value = matmul(basis_values, local_dofs)
        if (present(divergence)) divergence = sum(local_dofs) / area
    end subroutine evaluate_rt_field_2d

    real(dp) function rt_l2_norm(mesh, dofs) result(norm)
        type(mesh_2d_t), intent(in) :: mesh
        complex(dp), intent(in) :: dofs(:)

        type(gauss_quadrature_triangle_t) :: quadrature
        real(dp) :: point(2), vertex_a(2), jacobian(2, 2), det_jacobian
        complex(dp) :: value(2)
        integer :: triangle_idx, q

        call validate_rt_field(mesh, dofs)
        quadrature = get_gauss_quadrature_triangle(2)
        norm = 0.0_dp

        do triangle_idx = 1, mesh%n_triangles
            vertex_a = mesh%vertices(:, mesh%triangles(1, triangle_idx))
            jacobian(:, 1) = mesh%vertices(:, &
                mesh%triangles(2, triangle_idx)) - vertex_a
            jacobian(:, 2) = mesh%vertices(:, &
                mesh%triangles(3, triangle_idx)) - vertex_a
            det_jacobian = jacobian(1, 1) * jacobian(2, 2) - &
                jacobian(1, 2) * jacobian(2, 1)

            do q = 1, quadrature%n_points
                point = vertex_a + &
                    quadrature%xi(q) * jacobian(:, 1) + &
                    quadrature%eta(q) * jacobian(:, 2)
                call evaluate_rt_field_2d(mesh, triangle_idx, &
                    point(1), point(2), dofs, value)
                norm = norm + det_jacobian * quadrature%weights(q) * &
                    real(dot_product(value, value), dp)
            end do
        end do

        norm = sqrt(norm)
        call quadrature%destroy()
    end function rt_l2_norm

    subroutine reconstruct_axisymmetric_fourier_toroidal( &
            mesh, mode, dofs, toroidal)
        type(mesh_2d_t), intent(in) :: mesh
        integer, intent(in) :: mode
        complex(dp), intent(in) :: dofs(:)
        complex(dp), intent(out) :: toroidal(:)

        complex(dp) :: local_dofs(3)
        real(dp) :: area
        integer :: triangle_idx

        call validate_rt_field(mesh, dofs)
        if (mode == 0) then
            error stop "Toroidal reconstruction requires a nonzero Fourier mode"
        end if
        if (size(toroidal) /= mesh%n_triangles) then
            error stop "Toroidal reconstruction output has the wrong size"
        end if

        do triangle_idx = 1, mesh%n_triangles
            call local_rt_dofs(mesh, triangle_idx, dofs, local_dofs)
            area = triangle_area(mesh, triangle_idx)
            toroidal(triangle_idx) = cmplx(0.0_dp, 1.0_dp, dp) * &
                sum(local_dofs) / (real(mode, dp) * area)
        end do
    end subroutine reconstruct_axisymmetric_fourier_toroidal

    subroutine local_rt_dofs(mesh, triangle_idx, dofs, local_dofs)
        type(mesh_2d_t), intent(in) :: mesh
        integer, intent(in) :: triangle_idx
        complex(dp), intent(in) :: dofs(:)
        complex(dp), intent(out) :: local_dofs(3)

        integer :: edge_dofs(3), edge_orientations(3), i

        call mesh%get_triangle_edge_dofs( &
            triangle_idx, edge_dofs, edge_orientations)
        do i = 1, 3
            local_dofs(i) = real(edge_orientations(i), dp) * &
                dofs(edge_dofs(i) + 1)
        end do
    end subroutine local_rt_dofs

    subroutine physical_to_reference( &
            mesh, triangle_idx, x, y, xi, eta, area)
        type(mesh_2d_t), intent(in) :: mesh
        integer, intent(in) :: triangle_idx
        real(dp), intent(in) :: x, y
        real(dp), intent(out) :: xi, eta, area

        real(dp) :: vertex_a(2), jacobian(2, 2), offset(2), determinant

        vertex_a = mesh%vertices(:, mesh%triangles(1, triangle_idx))
        jacobian(:, 1) = mesh%vertices(:, &
            mesh%triangles(2, triangle_idx)) - vertex_a
        jacobian(:, 2) = mesh%vertices(:, &
            mesh%triangles(3, triangle_idx)) - vertex_a
        determinant = jacobian(1, 1) * jacobian(2, 2) - &
            jacobian(1, 2) * jacobian(2, 1)
        if (determinant <= 0.0_dp) then
            error stop "RT reconstruction requires counter-clockwise triangles"
        end if

        offset = [x, y] - vertex_a
        xi = (jacobian(2, 2) * offset(1) - &
            jacobian(1, 2) * offset(2)) / determinant
        eta = (-jacobian(2, 1) * offset(1) + &
            jacobian(1, 1) * offset(2)) / determinant
        area = 0.5_dp * determinant
    end subroutine physical_to_reference

    real(dp) function triangle_area(mesh, triangle_idx) result(area)
        type(mesh_2d_t), intent(in) :: mesh
        integer, intent(in) :: triangle_idx

        real(dp) :: a(2), b(2), c(2)

        a = mesh%vertices(:, mesh%triangles(1, triangle_idx))
        b = mesh%vertices(:, mesh%triangles(2, triangle_idx))
        c = mesh%vertices(:, mesh%triangles(3, triangle_idx))
        area = 0.5_dp * ((b(1) - a(1)) * (c(2) - a(2)) - &
            (b(2) - a(2)) * (c(1) - a(1)))
        if (area <= 0.0_dp) then
            error stop "RT reconstruction requires counter-clockwise triangles"
        end if
    end function triangle_area

    subroutine validate_rt_field(mesh, dofs)
        type(mesh_2d_t), intent(in) :: mesh
        complex(dp), intent(in) :: dofs(:)

        if (.not. allocated(mesh%edge_to_dof)) then
            error stop "RT field evaluation requires edge DOF numbering"
        end if
        if (size(dofs) /= mesh%n_edges) then
            error stop "RT field coefficient array has the wrong size"
        end if
    end subroutine validate_rt_field

end module fortfem_rt_field_2d
