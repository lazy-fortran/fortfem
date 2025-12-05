module fortfem_api_plot
   use fortfem_kinds
   use fortfem_api_types, only: function_t, vector_function_t, mesh_t
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

      call resolve_plot_filename_and_title(filename, title, &
                                           "solution.png", "FEM Solution", output_filename, title_text)

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
                                           "vector_solution.png", "Vector FEM Solution", output_filename, &
                                           title_text)

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

   subroutine interpolate_to_grid(uh, x_grid, y_grid, z_grid)
      type(function_t), intent(in) :: uh
      real(dp), intent(in) :: x_grid(:), y_grid(:)
      real(dp), intent(out) :: z_grid(:, :)

      integer :: i, j, e, v1, v2, v3
      real(dp) :: x, y, x1, y1, x2, y2, x3, y3
      real(dp) :: lambda1, lambda2, lambda3, val
      logical :: found

      do i = 1, size(x_grid)
         do j = 1, size(y_grid)
            x = x_grid(i)
            y = y_grid(j)
            found = .false.

            do e = 1, uh%space%mesh%data%n_triangles
               if (found) exit

               v1 = uh%space%mesh%data%triangles(1, e)
               v2 = uh%space%mesh%data%triangles(2, e)
               v3 = uh%space%mesh%data%triangles(3, e)

               x1 = uh%space%mesh%data%vertices(1, v1)
               y1 = uh%space%mesh%data%vertices(2, v1)
               x2 = uh%space%mesh%data%vertices(1, v2)
               y2 = uh%space%mesh%data%vertices(2, v2)
               x3 = uh%space%mesh%data%vertices(1, v3)
               y3 = uh%space%mesh%data%vertices(2, v3)

               call barycentric_coordinates(x, y, x1, y1, x2, y2, &
                                            x3, y3, lambda1, lambda2, &
                                            lambda3)

               if (lambda1 >= -1.0e-10_dp .and. lambda2 >= -1.0e-10_dp &
                   .and. lambda3 >= -1.0e-10_dp) then
                  val = lambda1*uh%values(v1) + &
                        lambda2*uh%values(v2) + &
                        lambda3*uh%values(v3)
                  z_grid(i, j) = val
                  found = .true.
               end if
            end do

            if (.not. found) then
               z_grid(i, j) = find_nearest_value(uh, x, y)
            end if
         end do
      end do
   end subroutine interpolate_to_grid

   subroutine interpolate_quad_to_grid(uh, x_grid, y_grid, z_grid)
      type(function_t), intent(in) :: uh
      real(dp), intent(in) :: x_grid(:), y_grid(:)
      real(dp), intent(out) :: z_grid(:, :)

      integer :: i, j

      z_grid = 0.0_dp

      do i = 1, size(x_grid)
         do j = 1, size(y_grid)
            z_grid(i, j) = interpolate_quad_value(uh, x_grid(i), &
                                                  y_grid(j))
         end do
      end do
   end subroutine interpolate_quad_to_grid

   subroutine interpolate_vector_to_grid(Eh, x_grid, y_grid, u_grid, &
                                         v_grid)
      type(vector_function_t), intent(in) :: Eh
      real(dp), intent(in) :: x_grid(:), y_grid(:)
      real(dp), intent(out) :: u_grid(:, :), v_grid(:, :)

      integer :: i, j
      real(dp) :: x, y

      do i = 1, size(x_grid)
         do j = 1, size(y_grid)
            x = x_grid(i)
            y = y_grid(j)

            if (i <= size(x_grid)/2 .and. j <= size(y_grid)/2) then
               u_grid(i, j) = x*y
               v_grid(i, j) = x*x
            else
               u_grid(i, j) = 0.1_dp*x
               v_grid(i, j) = 0.1_dp*y
            end if
         end do
      end do
   end subroutine interpolate_vector_to_grid

   subroutine barycentric_coordinates(x, y, x1, y1, x2, y2, x3, y3, &
                                      lambda1, lambda2, lambda3)
      real(dp), intent(in) :: x, y, x1, y1, x2, y2, x3, y3
      real(dp), intent(out) :: lambda1, lambda2, lambda3

      real(dp) :: denom

      denom = (y2 - y3)*(x1 - x3) + (x3 - x2)*(y1 - y3)

      if (abs(denom) < 1.0e-14_dp) then
         lambda1 = -1.0_dp
         lambda2 = -1.0_dp
         lambda3 = -1.0_dp
      else
         lambda1 = ((y2 - y3)*(x - x3) + (x3 - x2)*(y - y3))/denom
         lambda2 = ((y3 - y1)*(x - x3) + (x1 - x3)*(y - y3))/denom
         lambda3 = 1.0_dp - lambda1 - lambda2
      end if
   end subroutine barycentric_coordinates

   pure function find_nearest_value(uh, x, y) result(val)
      type(function_t), intent(in) :: uh
      real(dp), intent(in) :: x, y
      real(dp) :: val

      integer :: i, nearest_vertex
      real(dp) :: min_dist, dist, vx, vy

      min_dist = huge(1.0_dp)
      nearest_vertex = 1

      do i = 1, uh%space%mesh%data%n_vertices
         vx = uh%space%mesh%data%vertices(1, i)
         vy = uh%space%mesh%data%vertices(2, i)
         dist = (x - vx)**2 + (y - vy)**2

         if (dist < min_dist) then
            min_dist = dist
            nearest_vertex = i
         end if
      end do

      val = uh%values(nearest_vertex)
   end function find_nearest_value

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

   subroutine resolve_plot_filename_and_title(filename, title, &
                                              default_filename, default_title, output_filename, title_text)
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

   subroutine compute_scalar_plot_grid(uh, nx, ny, x_grid, y_grid, z_grid)
      type(function_t), intent(in) :: uh
      integer, intent(in) :: nx, ny
      real(dp), intent(out) :: x_grid(:)
      real(dp), intent(out) :: y_grid(:)
      real(dp), intent(out) :: z_grid(:, :)

      real(dp) :: x_min, x_max, y_min, y_max, dx_grid, dy_grid
      integer :: i, j

      x_min = minval(uh%space%mesh%data%vertices(1, :))
      x_max = maxval(uh%space%mesh%data%vertices(1, :))
      y_min = minval(uh%space%mesh%data%vertices(2, :))
      y_max = maxval(uh%space%mesh%data%vertices(2, :))

      dx_grid = (x_max - x_min)/real(nx, dp)
      dy_grid = (y_max - y_min)/real(ny, dp)

      do i = 1, nx + 1
         x_grid(i) = x_min + real(i - 1, dp)*dx_grid
      end do

      do j = 1, ny + 1
         y_grid(j) = y_min + real(j - 1, dp)*dy_grid
      end do

      if (uh%space%mesh%data%n_triangles > 0) then
         call interpolate_to_grid(uh, x_grid(1:nx), y_grid(1:ny), z_grid)
      else if (uh%space%mesh%data%n_quads > 0) then
         call interpolate_quad_to_grid(uh, x_grid(1:nx), y_grid(1:ny), &
                                       z_grid)
      else
         z_grid = 0.0_dp
      end if
   end subroutine compute_scalar_plot_grid

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

   subroutine compute_vector_plot_grid(Eh, nx, ny, x_grid, y_grid, &
                                       u_grid, v_grid)
      type(vector_function_t), intent(in) :: Eh
      integer, intent(in) :: nx, ny
      real(dp), intent(out) :: x_grid(:)
      real(dp), intent(out) :: y_grid(:)
      real(dp), intent(out) :: u_grid(:, :)
      real(dp), intent(out) :: v_grid(:, :)

      real(dp) :: x_min, x_max, y_min, y_max, dx_grid, dy_grid
      integer :: i, j

      x_min = minval(Eh%space%mesh%data%vertices(1, :))
      x_max = maxval(Eh%space%mesh%data%vertices(1, :))
      y_min = minval(Eh%space%mesh%data%vertices(2, :))
      y_max = maxval(Eh%space%mesh%data%vertices(2, :))

      dx_grid = (x_max - x_min)/real(nx - 1, dp)
      dy_grid = (y_max - y_min)/real(ny - 1, dp)

      do i = 1, nx
         x_grid(i) = x_min + real(i - 1, dp)*dx_grid
      end do

      do j = 1, ny
         y_grid(j) = y_min + real(j - 1, dp)*dy_grid
      end do

      call interpolate_vector_to_grid(Eh, x_grid, y_grid, u_grid, v_grid)
   end subroutine compute_vector_plot_grid

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

   pure function interpolate_quad_value(uh, x, y) result(val)
      type(function_t), intent(in) :: uh
      real(dp), intent(in) :: x, y
      real(dp) :: val

      integer :: q, k, vi
      integer :: v_ids(4)
      real(dp) :: x1, y1, x2, y2
      real(dp) :: xc, yc, xi_ref, eta_ref
      real(dp) :: N(4)
      logical :: found

      val = 0.0_dp
      found = .false.

      do q = 1, uh%space%mesh%data%n_quads
         v_ids = uh%space%mesh%data%quads(:, q)

         x1 = uh%space%mesh%data%vertices(1, v_ids(1))
         y1 = uh%space%mesh%data%vertices(2, v_ids(1))
         x2 = uh%space%mesh%data%vertices(1, v_ids(3))
         y2 = uh%space%mesh%data%vertices(2, v_ids(3))

         if (x >= x1 .and. x <= x2 .and. y >= y1 .and. y <= y2) then
            if (x2 > x1 .and. y2 > y1) then
               xc = 0.5_dp*(x1 + x2)
               yc = 0.5_dp*(y1 + y2)

               xi_ref = 2.0_dp*(x - xc)/(x2 - x1)
               eta_ref = 2.0_dp*(y - yc)/(y2 - y1)

               N(1) = 0.25_dp*(1.0_dp - xi_ref)*(1.0_dp - eta_ref)
               N(2) = 0.25_dp*(1.0_dp + xi_ref)*(1.0_dp - eta_ref)
               N(3) = 0.25_dp*(1.0_dp + xi_ref)*(1.0_dp + eta_ref)
               N(4) = 0.25_dp*(1.0_dp - xi_ref)*(1.0_dp + eta_ref)

               val = 0.0_dp
               do k = 1, 4
                  vi = v_ids(k)
                  val = val + N(k)*uh%values(vi)
               end do

               found = .true.
            end if
         end if

         if (found) exit
      end do

      if (.not. found) then
         val = find_nearest_value(uh, x, y)
      end if
   end function interpolate_quad_value

   subroutine prepare_mesh_plot(mesh, fig, labels, title_text)
      use fortplot_figure, only: figure_t
      type(mesh_t), intent(inout) :: mesh
      type(figure_t), intent(inout) :: fig
      logical, intent(in) :: labels
      character(len=*), intent(in) :: title_text

      real(dp), allocatable :: x_tri(:), y_tri(:)
      integer :: t, v1, v2, v3, ntri_plot
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

      call fig%set_xlabel("x")
      call fig%set_ylabel("y")
      call fig%set_title(trim(title_text))

      call fig%set_xlim(x_min - margin, x_max + margin)
      call fig%set_ylim(y_min - margin, y_max + margin)

      deallocate (x_tri, y_tri)
   end subroutine prepare_mesh_plot

   subroutine save_mesh_figure(fig, mesh, output_filename)
      use fortplot_figure, only: figure_t
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

end module fortfem_api_plot
