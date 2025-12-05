module fortfem_api_mesh_boundaries
   use fortfem_kinds, only: dp, pi
   use fortfem_boundary, only: boundary_t
   implicit none

   private

   public :: circle_boundary
   public :: rectangle_boundary
   public :: line_segment
   public :: arc_segment
   public :: l_shape_boundary

contains

   function circle_boundary(center, radius, n) result(boundary)
      real(dp), intent(in) :: center(2), radius
      integer, intent(in) :: n
      type(boundary_t) :: boundary
      integer :: i
      real(dp) :: theta

      boundary%n_points = n
      allocate (boundary%points(2, n))
      allocate (boundary%labels(n - 1))

      do i = 1, n
         theta = 2.0_dp*pi*real(i - 1, dp)/real(n, dp)
         boundary%points(1, i) = center(1) + radius*cos(theta)
         boundary%points(2, i) = center(2) + radius*sin(theta)
      end do

      boundary%labels = 1
      boundary%is_closed = .true.
   end function circle_boundary

   function rectangle_boundary(domain, n) result(boundary)
      real(dp), intent(in) :: domain(4)
      integer, intent(in) :: n
      type(boundary_t) :: boundary
      integer :: i, idx
      real(dp) :: t

      boundary%n_points = 4*n
      allocate (boundary%points(2, 4*n))
      allocate (boundary%labels(4*n - 1))

      idx = 0

      do i = 1, n
         idx = idx + 1
         t = real(i - 1, dp)/real(n - 1, dp)
         boundary%points(1, idx) = domain(1) + t*(domain(2) - domain(1))
         boundary%points(2, idx) = domain(3)
      end do

      do i = 1, n
         idx = idx + 1
         t = real(i - 1, dp)/real(n - 1, dp)
         boundary%points(1, idx) = domain(2)
         boundary%points(2, idx) = domain(3) + t*(domain(4) - domain(3))
      end do

      do i = 1, n
         idx = idx + 1
         t = real(i - 1, dp)/real(n - 1, dp)
         boundary%points(1, idx) = domain(2) - t*(domain(2) - domain(1))
         boundary%points(2, idx) = domain(4)
      end do

      do i = 1, n
         idx = idx + 1
         t = real(i - 1, dp)/real(n - 1, dp)
         boundary%points(1, idx) = domain(1)
         boundary%points(2, idx) = domain(4) - t*(domain(4) - domain(3))
      end do

      call set_rectangle_boundary_labels(boundary, n)

      boundary%is_closed = .true.
   end function rectangle_boundary

   pure subroutine set_rectangle_boundary_labels(boundary, n)
      type(boundary_t), intent(inout) :: boundary
      integer, intent(in) :: n
      integer :: i

      do i = 1, n - 1
         boundary%labels(i) = 1
      end do
      do i = n, 2*n - 2
         boundary%labels(i) = 2
      end do
      do i = 2*n - 1, 3*n - 3
         boundary%labels(i) = 3
      end do
      do i = 3*n - 2, 4*n - 1
         boundary%labels(i) = 4
      end do
   end subroutine set_rectangle_boundary_labels

   function line_segment(p1, p2, n) result(boundary)
      real(dp), intent(in) :: p1(2), p2(2)
      integer, intent(in) :: n
      type(boundary_t) :: boundary
      integer :: i
      real(dp) :: t

      boundary%n_points = n
      allocate (boundary%points(2, n))
      allocate (boundary%labels(n - 1))

      do i = 1, n
         t = real(i - 1, dp)/real(n - 1, dp)
         boundary%points(:, i) = p1 + t*(p2 - p1)
      end do

      boundary%labels = 1
      boundary%is_closed = .false.
   end function line_segment

   function arc_segment(p1, p2, center, n) result(boundary)
      real(dp), intent(in) :: p1(2), p2(2), center(2)
      integer, intent(in) :: n
      type(boundary_t) :: boundary
      integer :: i
      real(dp) :: radius1, radius2, radius
      real(dp) :: theta1, theta2, dtheta, angle

      boundary%n_points = n
      allocate (boundary%points(2, n))
      allocate (boundary%labels(n - 1))

      radius1 = sqrt((p1(1) - center(1))**2 + (p1(2) - center(2))**2)
      radius2 = sqrt((p2(1) - center(1))**2 + (p2(2) - center(2))**2)
      radius = 0.5_dp*(radius1 + radius2)

      theta1 = atan2(p1(2) - center(2), p1(1) - center(1))
      theta2 = atan2(p2(2) - center(2), p2(1) - center(1))
      dtheta = theta2 - theta1
      if (dtheta <= 0.0_dp) dtheta = dtheta + 2.0_dp*pi

      do i = 1, n
         angle = theta1 + dtheta*real(i - 1, dp)/ &
                 real(max(n - 1, 1), dp)
         boundary%points(1, i) = center(1) + radius*cos(angle)
         boundary%points(2, i) = center(2) + radius*sin(angle)
      end do

      boundary%labels = 1
      boundary%is_closed = .false.
   end function arc_segment

   function l_shape_boundary(size, n) result(boundary)
      real(dp), intent(in) :: size
      integer, intent(in) :: n
      type(boundary_t) :: boundary
      integer :: idx, n_per_segment
      real(dp) :: s

      s = size

      n_per_segment = max(n - 1, 1)
      boundary%n_points = 6*n_per_segment
      allocate (boundary%points(2, boundary%n_points))
      allocate (boundary%labels(boundary%n_points))

      idx = 0

      call add_l_shape_segment(boundary%points, idx, n_per_segment, &
                               0.0_dp, 0.0_dp, s, 0.0_dp)
      call add_l_shape_segment(boundary%points, idx, n_per_segment, &
                               s, 0.0_dp, 0.0_dp, s)
      call add_l_shape_segment(boundary%points, idx, n_per_segment, &
                               s, s, s, 0.0_dp)
      call add_l_shape_segment(boundary%points, idx, n_per_segment, &
                               2.0_dp*s, s, 0.0_dp, s)
      call add_l_shape_segment(boundary%points, idx, n_per_segment, &
                               2.0_dp*s, 2.0_dp*s, -2.0_dp*s, 0.0_dp)
      call add_l_shape_segment(boundary%points, idx, n_per_segment, &
                               0.0_dp, 2.0_dp*s, 0.0_dp, -2.0_dp*s)

      boundary%labels = 1
      boundary%is_closed = .true.
   end function l_shape_boundary

   pure subroutine add_l_shape_segment(points, idx, n_per_segment, &
                                       start_x, start_y, dx, dy)
      real(dp), intent(inout) :: points(:, :)
      integer, intent(inout) :: idx
      integer, intent(in) :: n_per_segment
      real(dp), intent(in) :: start_x, start_y, dx, dy
      integer :: i
      real(dp) :: t

      do i = 0, n_per_segment - 1
         idx = idx + 1
         t = real(i, dp)/real(n_per_segment, dp)
         points(1, idx) = start_x + dx*t
         points(2, idx) = start_y + dy*t
      end do
   end subroutine add_l_shape_segment

end module fortfem_api_mesh_boundaries

