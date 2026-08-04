! Extracted from src/generated/fortfem_fci_polygon_edge_area_products.f90 by
! tools/fortad/extract.py. It is the primal fortsym generated
! its products from, with the offload directives removed:
! they are the consumer's choice, not part of what the
! function means.

pure subroutine generated_fci_polygon_edge_area(x_start, y_start, x_end, y_end, edge_area)
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    real(dp), intent(in) :: x_start, y_start, x_end, y_end
    real(dp), intent(out) :: edge_area

    edge_area = (-x_end*y_start + x_start*y_end)/2.0000000000000000E+000_dp

end subroutine generated_fci_polygon_edge_area
