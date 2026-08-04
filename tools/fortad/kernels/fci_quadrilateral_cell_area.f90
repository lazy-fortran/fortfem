! Extracted from src/generated/fortfem_fci_quadrilateral_area_products.f90 by
! tools/fortad/extract.py. It is the primal fortsym generated
! its products from, with the offload directives removed:
! they are the consumer's choice, not part of what the
! function means.

pure subroutine generated_fci_quadrilateral_cell_area(x_1, y_1, x_2, y_2, x_3, y_3, x_4, y_4, area)
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    real(dp), intent(in) :: x_1, y_1, x_2, y_2, x_3, y_3, x_4, y_4
    real(dp), intent(out) :: area

    area = (x_1*y_2 - x_1*y_4 - x_2*y_1 + x_2*y_3 - x_3*y_2 + x_3*y_4 + x_4*y_1 - x_4*y_3)/ &
        2.0000000000000000E+000_dp

end subroutine generated_fci_quadrilateral_cell_area
