module fortfem_edge_moment_orientation
    use fortfem_kinds, only: dp
    implicit none

    private

    public :: apply_edge_moment_orientation

contains

    pure subroutine apply_edge_moment_orientation( &
            order, orientation, local_dofs, oriented_dofs, status)
        integer, intent(in) :: order, orientation
        real(dp), intent(in) :: local_dofs(:)
        real(dp), intent(out) :: oriented_dofs(:)
        integer, intent(out) :: status

        integer :: moment

        oriented_dofs = 0.0_dp
        status = 1
        if (order < 1 .or. size(local_dofs) /= order) return
        if (size(oriented_dofs) /= order) return
        if (orientation /= -1 .and. orientation /= 1) return

        oriented_dofs = local_dofs
        if (orientation == -1) then
            do moment = 1, order
                if (mod(moment, 2) == 1) then
                    oriented_dofs(moment) = -oriented_dofs(moment)
                end if
            end do
        end if
        status = 0
    end subroutine apply_edge_moment_orientation

end module fortfem_edge_moment_orientation
