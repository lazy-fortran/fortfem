module fortfem_api_plot_mesh
    use, intrinsic :: ieee_arithmetic, only: ieee_quiet_nan, ieee_value
    use fortfem_kinds, only: dp
    use fortfem_api_types, only: mesh_t
    use fortplot_figure, only: figure_t
    implicit none

    private

    public :: prepare_mesh_plot
    public :: save_mesh_figure
    public :: build_mesh_edge_path

contains

    subroutine prepare_mesh_plot(mesh, fig, labels, title_text)
        type(mesh_t), intent(inout) :: mesh
        type(figure_t), intent(inout) :: fig
        logical, intent(in) :: labels
        character(len=*), intent(in) :: title_text

        real(dp), allocatable :: x_edges(:), y_edges(:)
        real(dp) :: x_min, x_max, y_min, y_max, margin

        x_min = minval(mesh%data%vertices(1, :))
        x_max = maxval(mesh%data%vertices(1, :))
        y_min = minval(mesh%data%vertices(2, :))
        y_max = maxval(mesh%data%vertices(2, :))
        margin = 0.1_dp*max(x_max - x_min, y_max - y_min)

        call fig%initialize()

        if (.not. allocated(mesh%data%edges)) then
            call mesh%data%build_connectivity()
        end if

        call build_mesh_edge_path( &
            mesh%data%vertices, mesh%data%edges, x_edges, y_edges)
        call fig%add_plot( &
            x_edges, y_edges, color=[0.12_dp, 0.31_dp, 0.45_dp])

        call fig%set_xlabel("x")
        call fig%set_ylabel("y")
        call fig%set_title(trim(title_text))

        call fig%set_xlim(x_min - margin, x_max + margin)
        call fig%set_ylim(y_min - margin, y_max + margin)

        deallocate (x_edges, y_edges)
    end subroutine prepare_mesh_plot

    pure subroutine build_mesh_edge_path(vertices, edges, x_path, y_path)
        real(dp), intent(in) :: vertices(:, :)
        integer, intent(in) :: edges(:, :)
        real(dp), allocatable, intent(out) :: x_path(:), y_path(:)

        real(dp) :: separator
        integer :: edge, offset

        separator = ieee_value(0.0_dp, ieee_quiet_nan)
        allocate (x_path(3*size(edges, 2)), y_path(3*size(edges, 2)))

        do edge = 1, size(edges, 2)
            offset = 3*(edge - 1)
            x_path(offset + 1:offset + 2) = vertices(1, edges(:, edge))
            y_path(offset + 1:offset + 2) = vertices(2, edges(:, edge))
            x_path(offset + 3) = separator
            y_path(offset + 3) = separator
        end do
    end subroutine build_mesh_edge_path

    subroutine save_mesh_figure(fig, mesh, output_filename)
        type(figure_t), intent(inout) :: fig
        type(mesh_t), intent(in) :: mesh
        character(len=*), intent(in) :: output_filename

        call fig%savefig(trim(output_filename))

        write (*, *) "Mesh plot saved to: ", trim(output_filename)
        write (*, *) "Mesh info:"
        write (*, *) "  Vertices: ", mesh%data%n_vertices
        write (*, *) "  Triangles: ", mesh%data%n_triangles
        write (*, *) "  Edges: ", mesh%data%n_edges
    end subroutine save_mesh_figure

end module fortfem_api_plot_mesh
