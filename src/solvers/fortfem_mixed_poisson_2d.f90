module fortfem_mixed_poisson_2d
    use fortfem_assembly_rt_arbitrary_order_2d, only: &
        assemble_triangle_rt_div_mass_csc, &
        assemble_triangle_rt_divergence_csc
    use fortfem_kinds, only: dp
    use fortfem_mesh_2d, only: mesh_2d_t
    use fortfem_sparse_direct, only: sparse_direct_solve_csc
    use fortsparse, only: csc_from_triplet, csc_t, &
        FORTSPARSE_INTERNAL_ERROR, FORTSPARSE_INVALID_MATRIX, &
        fortsparse_status_t, status_set
    implicit none

    private

    public :: solve_mixed_poisson_rt0

contains

    !> Solve q + grad(u) = 0 and div(q) = f with RT0 fluxes, DG0 pressure,
    !> and homogeneous Dirichlet pressure data.
    subroutine solve_mixed_poisson_rt0( &
            mesh, cell_source, flux_dofs, cell_pressure, status)
        type(mesh_2d_t), intent(inout) :: mesh
        real(dp), intent(in) :: cell_source(:)
        real(dp), allocatable, intent(out) :: flux_dofs(:)
        real(dp), allocatable, intent(out) :: cell_pressure(:)
        type(fortsparse_status_t), intent(out) :: status

        type(csc_t) :: divergence, flux_mass, system
        integer, allocatable :: columns(:), rows(:)
        real(dp), allocatable :: right_hand_side(:), solution(:), values(:)
        real(dp) :: determinant, vertex_a(2), vertex_b(2), vertex_c(2)
        integer :: column, entry, flux_count, matrix_entry
        integer :: pressure_count, solve_status, system_size, triangle

        call status_set( &
            status, FORTSPARSE_INVALID_MATRIX, &
            "RT0-DG0 mixed Poisson solve failed")
        if (mesh%n_triangles < 1) return
        if (size(cell_source) /= mesh%n_triangles) return

        call assemble_triangle_rt_div_mass_csc( &
            mesh, 0, 2, flux_mass, status, 0.0_dp, 1.0_dp)
        if (status%code /= 0) return
        call assemble_triangle_rt_divergence_csc( &
            mesh, 0, 2, divergence, status)
        if (status%code /= 0) return

        flux_count = flux_mass%nrow
        pressure_count = divergence%nrow
        if (flux_mass%ncol /= flux_count) return
        if (divergence%ncol /= flux_count) return
        if (pressure_count /= mesh%n_triangles) return
        system_size = flux_count + pressure_count
        allocate(rows(flux_mass%nnz + 2 * divergence%nnz))
        allocate(columns(size(rows)), values(size(rows)))

        entry = 0
        do column = 1, flux_mass%ncol
            do matrix_entry = flux_mass%col_ptr(column), &
                    flux_mass%col_ptr(column + 1) - 1
                entry = entry + 1
                rows(entry) = flux_mass%row_idx(matrix_entry)
                columns(entry) = column
                values(entry) = flux_mass%val(matrix_entry)
            end do
        end do
        do column = 1, divergence%ncol
            do matrix_entry = divergence%col_ptr(column), &
                    divergence%col_ptr(column + 1) - 1
                entry = entry + 1
                rows(entry) = column
                columns(entry) = &
                    flux_count + divergence%row_idx(matrix_entry)
                values(entry) = -divergence%val(matrix_entry)

                entry = entry + 1
                rows(entry) = &
                    flux_count + divergence%row_idx(matrix_entry)
                columns(entry) = column
                values(entry) = divergence%val(matrix_entry)
            end do
        end do
        call csc_from_triplet( &
            system_size, system_size, rows, columns, values, system, status)
        if (status%code /= 0) return

        allocate(right_hand_side(system_size), solution(system_size))
        right_hand_side = 0.0_dp
        do triangle = 1, mesh%n_triangles
            vertex_a = mesh%vertices(:, mesh%triangles(1, triangle))
            vertex_b = mesh%vertices(:, mesh%triangles(2, triangle))
            vertex_c = mesh%vertices(:, mesh%triangles(3, triangle))
            determinant = &
                (vertex_b(1) - vertex_a(1)) * &
                (vertex_c(2) - vertex_a(2)) - &
                (vertex_b(2) - vertex_a(2)) * &
                (vertex_c(1) - vertex_a(1))
            if (determinant <= 0.0_dp) return
            right_hand_side(flux_count + triangle) = &
                0.5_dp * determinant * cell_source(triangle)
        end do

        call sparse_direct_solve_csc( &
            system_size, system%col_ptr, system%row_idx, system%val, &
            right_hand_side, solution, solve_status)
        if (solve_status /= 0) then
            call status_set( &
                status, FORTSPARSE_INTERNAL_ERROR, &
                "RT0-DG0 mixed Poisson sparse solve failed")
            return
        end if
        allocate(flux_dofs(flux_count), cell_pressure(pressure_count))
        flux_dofs = solution(:flux_count)
        cell_pressure = solution(flux_count + 1:)
        call status_set(status, 0, "")
    end subroutine solve_mixed_poisson_rt0

end module fortfem_mixed_poisson_2d
