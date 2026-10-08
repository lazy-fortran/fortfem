module fortfem_tetra_nedelec_first_order
    use fortfem_kinds, only: dp
    use fortfem_generated_tetra_whitney, only: generated_tetra_whitney
    implicit none

    private

    public :: evaluate_tetra_nedelec_first_order

contains

    pure subroutine evaluate_tetra_nedelec_first_order( &
            point, values, curls, status)
        real(dp), intent(in) :: point(3)
        real(dp), intent(out) :: values(3, 6), curls(3, 6)
        integer, intent(out) :: status

        real(dp) :: lambda(4), tolerance

        values = 0.0_dp
        curls = 0.0_dp
        status = 1
        tolerance = 64.0_dp * epsilon(1.0_dp)
        lambda = [ &
            1.0_dp - point(1) - point(2) - point(3), &
            point(1), point(2), point(3)]
        if (any(lambda < -tolerance)) return
        if (any(lambda > 1.0_dp + tolerance)) return

        call generated_tetra_whitney(point(1), point(2), point(3), values, curls)
        status = 0
    end subroutine evaluate_tetra_nedelec_first_order

end module fortfem_tetra_nedelec_first_order
