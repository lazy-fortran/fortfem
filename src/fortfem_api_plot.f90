module fortfem_api_plot
   use fortfem_kinds
   use fortfem_api_types, only: function_t, vector_function_t, mesh_t
   use fortfem_api_plot_interpolation, only: compute_scalar_plot_grid, &
                                             compute_vector_plot_grid
   use fortfem_api_plot_mesh, only: prepare_mesh_plot, save_mesh_figure
   implicit none

   private

   public :: plot
   public :: plot_function_scalar
   public :: plot_vector_function
   public :: plot_mesh

   interface plot
      module procedure plot_function_scalar
      module procedure plot_vector_function
      module procedure plot_mesh
   end interface plot

contains

   subroutine plot_function_scalar(uh, filename, title, colormap)
      use fortplot, only: figure, contour_filled, xlabel, ylabel, &
                          plot_title => title, savefig, pcolormesh, add_plot
      type(function_t), intent(in) :: uh
      character(len=*), intent(in), optional :: filename
      character(len=*), intent(in), optional :: title
      character(len=*), intent(in), optional :: colormap

      integer, parameter :: nx = 40, ny = 40
      real(dp) :: x_grid(nx + 1)
      real(dp) :: y_grid(ny + 1)
      real(dp) :: z_grid(nx, ny)
      character(len=64) :: output_filename
      character(len=128) :: title_text
      character(len=32) :: cmap

      call resolve_plot_filename_and_title(filename, title, "solution.png", &
                                           "FEM Solution", output_filename, &
                                           title_text)

      if (present(colormap)) then
         cmap = colormap
      else
         cmap = "viridis"
      end if

      call compute_scalar_plot_grid(uh, nx, ny, x_grid, y_grid, z_grid)

      call render_scalar_solution(uh, x_grid, y_grid, z_grid, title_text, &
                                  cmap, output_filename)
   end subroutine plot_function_scalar

   subroutine plot_vector_function(Eh, filename, title, plot_type)
      use fortplot, only: figure, streamplot, xlabel, ylabel, &
                          plot_title => title, savefig
      type(vector_function_t), intent(in) :: Eh
      character(len=*), intent(in), optional :: filename
      character(len=*), intent(in), optional :: title
      character(len=*), intent(in), optional :: plot_type

      integer, parameter :: nx = 20, ny = 20
      real(dp) :: x_grid(nx)
      real(dp) :: y_grid(ny)
      real(dp) :: u_grid(nx, ny)
      real(dp) :: v_grid(nx, ny)
      character(len=64) :: output_filename
      character(len=128) :: title_text
      character(len=32) :: ptype

      call resolve_plot_filename_and_title(filename, title, &
                                           "vector_solution.png", &
                                           "Vector FEM Solution", &
                                           output_filename, title_text)

      if (present(plot_type)) then
         ptype = plot_type
      else
         ptype = "streamplot"
      end if

      call compute_vector_plot_grid(Eh, nx, ny, x_grid, y_grid, u_grid, &
                                    v_grid)

      call render_vector_solution(x_grid, y_grid, u_grid, v_grid, &
                                  title_text, ptype, output_filename)
   end subroutine plot_vector_function

   subroutine plot_mesh(mesh, filename, title, show_labels)
      use fortplot_figure, only: figure_t
      type(mesh_t), intent(inout) :: mesh
      character(len=*), intent(in), optional :: filename
      character(len=*), intent(in), optional :: title
      logical, intent(in), optional :: show_labels

      type(figure_t) :: fig
      character(len=64) :: output_filename
      character(len=128) :: title_text
      logical :: labels

      call resolve_plot_filename_and_title(filename, title, "mesh.png", &
                                           "FEM Mesh", output_filename, title_text)

      if (present(show_labels)) then
         labels = show_labels
      else
         labels = .false.
      end if

      call prepare_mesh_plot(mesh, fig, labels, title_text)

      call save_mesh_figure(fig, mesh, output_filename)
   end subroutine plot_mesh

   subroutine resolve_plot_filename_and_title(filename, title, default_filename, &
                                              default_title, output_filename, &
                                              title_text)
      character(len=*), intent(in), optional :: filename
      character(len=*), intent(in), optional :: title
      character(len=*), intent(in) :: default_filename
      character(len=*), intent(in) :: default_title
      character(len=*), intent(out) :: output_filename
      character(len=*), intent(out) :: title_text

      if (present(filename)) then
         output_filename = filename
      else
         output_filename = default_filename
      end if

      if (present(title)) then
         title_text = title
      else
         title_text = default_title
      end if
   end subroutine resolve_plot_filename_and_title

   subroutine render_scalar_solution(uh, x_grid, y_grid, z_grid, &
                                     title_text, cmap, output_filename)
      use fortplot, only: figure, contour_filled, xlabel, ylabel, &
                          plot_title => title, savefig, pcolormesh, add_plot
      type(function_t), intent(in) :: uh
      real(dp), intent(in) :: x_grid(:)
      real(dp), intent(in) :: y_grid(:)
      real(dp), intent(in) :: z_grid(:, :)
      character(len=*), intent(in) :: title_text
      character(len=*), intent(in) :: cmap
      character(len=*), intent(in) :: output_filename

      real(dp) :: x_edges(2), y_edges(2)
      real(dp), parameter :: black(3) = [0.0_dp, 0.0_dp, 0.0_dp]
      integer :: v1, v2, e

      call figure()
      call plot_title(trim(title_text))
      call xlabel("x")
      call ylabel("y")
      call pcolormesh(x_grid, y_grid, z_grid, colormap=trim(cmap))

      if (.not. allocated(uh%space%mesh%data%edges)) then
         call uh%space%mesh%data%build_connectivity()
      end if

      do e = 1, uh%space%mesh%data%n_edges
         v1 = uh%space%mesh%data%edges(1, e)
         v2 = uh%space%mesh%data%edges(2, e)

         x_edges(1) = uh%space%mesh%data%vertices(1, v1)
         x_edges(2) = uh%space%mesh%data%vertices(1, v2)
         y_edges(1) = uh%space%mesh%data%vertices(2, v1)
         y_edges(2) = uh%space%mesh%data%vertices(2, v2)

         call add_plot(x_edges, y_edges, color=black)
      end do

      call savefig(trim(output_filename))

      write (*, *) "Plot saved to: ", trim(output_filename)
      write (*, *) "Solution range: [", minval(uh%values), ",", &
         maxval(uh%values), "]"
   end subroutine render_scalar_solution

   subroutine render_vector_solution(x_grid, y_grid, u_grid, v_grid, &
                                     title_text, ptype, output_filename)
      use fortplot, only: figure, streamplot, xlabel, ylabel, &
                          plot_title => title, savefig
      real(dp), intent(in) :: x_grid(:)
      real(dp), intent(in) :: y_grid(:)
      real(dp), intent(in) :: u_grid(:, :)
      real(dp), intent(in) :: v_grid(:, :)
      character(len=*), intent(in) :: title_text
      character(len=*), intent(in) :: ptype
      character(len=*), intent(in) :: output_filename

      call figure()
      call plot_title(trim(title_text))
      call xlabel("x")
      call ylabel("y")

      select case (trim(ptype))
      case ("streamplot")
         call streamplot(x_grid, y_grid, u_grid, v_grid)
      case default
         call streamplot(x_grid, y_grid, u_grid, v_grid)
      end select

      call savefig(trim(output_filename))

      write (*, *) "Vector plot saved to: ", trim(output_filename)
      write (*, *) "Vector magnitude range: [", &
         minval(sqrt(u_grid**2 + v_grid**2)), ",", &
         maxval(sqrt(u_grid**2 + v_grid**2)), "]"
   end subroutine render_vector_solution

end module fortfem_api_plot
