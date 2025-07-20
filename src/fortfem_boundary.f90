module fortfem_boundary
    ! Boundary representation for mesh generation
    use fortfem_kinds
    implicit none
    
    private
    public :: boundary_t
    
    ! Boundary type for defining domains
    type :: boundary_t
        integer :: n_points = 0
        real(dp), allocatable :: points(:,:)     ! (2, n_points)
        integer, allocatable :: labels(:)       ! (n_points-1) segment labels
        logical :: is_closed = .false.
    end type boundary_t
    
end module fortfem_boundary