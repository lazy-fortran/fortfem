module fortfem_edge_interpolation_2d
    use fortfem_kinds, only: dp
    use fortfem_mesh_2d, only: mesh_2d_t
    use fortnum_quadrature, only: gauss_legendre_ab
    implicit none
    private

    public :: interpolate_nedelec_edge_dofs
    public :: interpolate_rt_edge_dofs

    abstract interface
        pure subroutine vector_field_2d(x, y, value)
            import dp
            real(dp), intent(in) :: x, y
            real(dp), intent(out) :: value(2)
        end subroutine vector_field_2d
    end interface

contains

    subroutine interpolate_nedelec_edge_dofs( &
            mesh, field, quadrature_points, dofs)
        type(mesh_2d_t), intent(inout) :: mesh
        procedure(vector_field_2d) :: field
        integer, intent(in) :: quadrature_points
        real(dp), intent(out) :: dofs(:)

        call interpolate_edge_dofs( &
            mesh, field, quadrature_points, .false., dofs)
    end subroutine interpolate_nedelec_edge_dofs

    subroutine interpolate_rt_edge_dofs(mesh, field, quadrature_points, dofs)
        type(mesh_2d_t), intent(inout) :: mesh
        procedure(vector_field_2d) :: field
        integer, intent(in) :: quadrature_points
        real(dp), intent(out) :: dofs(:)

        call interpolate_edge_dofs( &
            mesh, field, quadrature_points, .true., dofs)
    end subroutine interpolate_rt_edge_dofs

    subroutine interpolate_edge_dofs( &
            mesh, field, quadrature_points, use_normal, dofs)
        type(mesh_2d_t), intent(inout) :: mesh
        procedure(vector_field_2d) :: field
        integer, intent(in) :: quadrature_points
        logical, intent(in) :: use_normal
        real(dp), intent(out) :: dofs(:)

        real(dp), allocatable :: nodes(:), weights(:)
        real(dp) :: point_a(2), edge_vector(2), direction(2)
        real(dp) :: point(2), value(2), moment
        integer :: edge_id, quadrature_id, dof_id

        if (.not. allocated(mesh%edges)) then
            error stop "Edge connectivity must be built before interpolation"
        end if
        if (quadrature_points < 1) then
            error stop "Edge interpolation needs at least one quadrature point"
        end if
        if (size(dofs) /= mesh%n_edges) then
            error stop "Edge interpolation output has the wrong size"
        end if
        if (.not. allocated(mesh%edge_to_dof)) then
            call mesh%build_edge_dof_numbering()
        end if

        allocate(nodes(quadrature_points), weights(quadrature_points))
        call gauss_legendre_ab( &
            quadrature_points, 0.0_dp, 1.0_dp, nodes, weights)

        do edge_id = 1, mesh%n_edges
            point_a = mesh%vertices(:, mesh%edges(1, edge_id))
            edge_vector = mesh%vertices(:, mesh%edges(2, edge_id)) - point_a
            if (use_normal) then
                direction = [edge_vector(2), -edge_vector(1)]
            else
                direction = edge_vector
            end if

            moment = 0.0_dp
            do quadrature_id = 1, quadrature_points
                point = point_a + nodes(quadrature_id) * edge_vector
                call field(point(1), point(2), value)
                moment = moment + weights(quadrature_id) * dot_product( &
                    value, direction)
            end do

            dof_id = mesh%edge_to_dof(edge_id) + 1
            dofs(dof_id) = moment
        end do
    end subroutine interpolate_edge_dofs

end module fortfem_edge_interpolation_2d
