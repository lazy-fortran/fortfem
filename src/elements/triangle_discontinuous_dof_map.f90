module fortfem_triangle_discontinuous_dof_map
    use fortfem_mesh_2d, only: mesh_2d_t
    implicit none

    private

    public :: build_triangle_discontinuous_dof_map

contains

    subroutine build_triangle_discontinuous_dof_map( &
            mesh, degree, global_dofs, global_dof_count, status)
        type(mesh_2d_t), intent(in) :: mesh
        integer, intent(in) :: degree
        integer, allocatable, intent(out) :: global_dofs(:, :)
        integer, intent(out) :: global_dof_count, status

        integer :: local_dof, local_dof_count, triangle

        global_dof_count = 0
        status = 1
        if (degree < 0 .or. mesh%n_triangles < 1) return

        local_dof_count = (degree + 1) * (degree + 2) / 2
        global_dof_count = local_dof_count * mesh%n_triangles
        allocate(global_dofs(local_dof_count, mesh%n_triangles))
        do triangle = 1, mesh%n_triangles
            do local_dof = 1, local_dof_count
                global_dofs(local_dof, triangle) = &
                    (triangle - 1) * local_dof_count + local_dof
            end do
        end do
        status = 0
    end subroutine build_triangle_discontinuous_dof_map

end module fortfem_triangle_discontinuous_dof_map
