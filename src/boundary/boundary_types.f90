module fortfem_boundary_types
    use fortfem_kinds, only: dp
    implicit none
    private

    public :: boundary_t

    type :: boundary_t
        integer :: n_points = 0
        real(dp), allocatable :: points(:, :)
        integer, allocatable :: labels(:)
        logical :: is_closed = .false.
    end type boundary_t

end module fortfem_boundary_types
