module fortfem_api_plot_mesh
   use fortfem_kinds, only: dp
   use fortfem_api_types, only: mesh_t
   use fortplot_figure, only: figure_t
   implicit none

   private

   public :: prepare_mesh_plot
   public :: save_mesh_figure

contains

   subroutine prepare_mesh_plot(mesh, fig, labels, title_text)
      type(mesh_t), intent(inout) :: mesh
      type(figure_t), intent(inout) :: fig
      logical, intent(in) :: labels
      character(len=*), intent(in) :: title_text

      real(dp), allocatable :: x_tri(:), y_tri(:)
      integer :: ntri_plot
      real(dp) :: x_min, x_max, y_min, y_max, margin

      x_min = minval(mesh%data%vertices(1, :))
      x_max = maxval(mesh%data%vertices(1, :))
      y_min = minval(mesh%data%vertices(2, :))
      y_max = maxval(mesh%data%vertices(2, :))
      margin = 0.1_dp*max(x_max - x_min, y_max - y_min)

      call fig%initialize()

      if (.not. allocated(mesh%data%triangles)) then
         call mesh%data%build_connectivity()
      end if

      allocate (x_tri(4), y_tri(4))

      ntri_plot = min(mesh%data%n_triangles, fig%state%max_plots)
      call add_mesh_triangles_to_figure(mesh, fig, ntri_plot, x_tri, y_tri)

      call fig%set_xlabel("x")
      call fig%set_ylabel("y")
      call fig%set_title(trim(title_text))

      call fig%set_xlim(x_min - margin, x_max + margin)
      call fig%set_ylim(y_min - margin, y_max + margin)

      deallocate (x_tri, y_tri)
   end subroutine prepare_mesh_plot

   subroutine add_mesh_triangles_to_figure(mesh, fig, ntri_plot, x_tri, y_tri)
      type(mesh_t), intent(in) :: mesh
      type(figure_t), intent(inout) :: fig
      integer, intent(in) :: ntri_plot
      real(dp), intent(inout) :: x_tri(:), y_tri(:)

      integer :: t, v1, v2, v3

      do t = 1, ntri_plot
         v1 = mesh%data%triangles(1, t)
         v2 = mesh%data%triangles(2, t)
         v3 = mesh%data%triangles(3, t)

         x_tri(1) = mesh%data%vertices(1, v1)
         y_tri(1) = mesh%data%vertices(2, v1)
         x_tri(2) = mesh%data%vertices(1, v2)
         y_tri(2) = mesh%data%vertices(2, v2)
         x_tri(3) = mesh%data%vertices(1, v3)
         y_tri(3) = mesh%data%vertices(2, v3)
         x_tri(4) = x_tri(1)
         y_tri(4) = y_tri(1)

         call fig%add_plot(x_tri, y_tri)
      end do
   end subroutine add_mesh_triangles_to_figure

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
