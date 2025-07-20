module mesh_1d
    ! Simple 1D mesh generation
    use fortfem_kinds, only: dp
    implicit none
    
    private
    public :: mesh_1d_t, create_interval, create_segments
    
    type :: mesh_1d_t
        integer :: n_points = 0
        real(dp), allocatable :: points(:)      ! (n_points)
        integer, allocatable :: segments(:,:)   ! (2, n_segments)  
        integer :: n_segments = 0
    end type mesh_1d_t
    
contains

    subroutine create_interval(mesh, a, b, n)
        !> Create uniform 1D mesh on interval [a,b] with n points
        type(mesh_1d_t), intent(out) :: mesh
        real(dp), intent(in) :: a, b
        integer, intent(in) :: n
        
        integer :: i
        real(dp) :: dx
        
        mesh%n_points = n
        mesh%n_segments = n - 1
        
        allocate(mesh%points(n))
        allocate(mesh%segments(2, n-1))
        
        ! Generate uniform points
        dx = (b - a) / real(n - 1, dp)
        do i = 1, n
            mesh%points(i) = a + (i - 1) * dx
        end do
        
        ! Generate segments
        do i = 1, n - 1
            mesh%segments(1, i) = i
            mesh%segments(2, i) = i + 1
        end do
        
    end subroutine create_interval

    subroutine create_segments(mesh, segment_points, densities)
        !> Create 1D mesh from multiple segments with different densities
        type(mesh_1d_t), intent(out) :: mesh
        real(dp), intent(in) :: segment_points(:)   ! endpoints of segments
        integer, intent(in) :: densities(:)         ! points per segment
        
        integer :: n_segments, total_points, i, j, idx
        real(dp) :: dx
        
        n_segments = size(segment_points) - 1
        total_points = sum(densities) - n_segments + 1  ! Account for shared endpoints
        
        mesh%n_points = total_points
        mesh%n_segments = total_points - 1
        
        allocate(mesh%points(total_points))
        allocate(mesh%segments(2, mesh%n_segments))
        
        ! Generate points
        idx = 1
        do i = 1, n_segments
            dx = (segment_points(i+1) - segment_points(i)) / real(densities(i) - 1, dp)
            
            do j = 1, densities(i)
                if (i > 1 .and. j == 1) cycle  ! Skip duplicate endpoint
                mesh%points(idx) = segment_points(i) + (j - 1) * dx
                idx = idx + 1
            end do
        end do
        
        ! Generate segments
        do i = 1, mesh%n_segments
            mesh%segments(1, i) = i
            mesh%segments(2, i) = i + 1
        end do
        
    end subroutine create_segments

end module mesh_1d