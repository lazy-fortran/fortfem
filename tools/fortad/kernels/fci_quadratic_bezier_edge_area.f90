! Extracted from src/generated/fortfem_fci_quadratic_bezier_edge_area_products.f90 by
! tools/fortad/extract.py. It is the primal fortsym generated
! its products from, with the offload directives removed:
! they are the consumer's choice, not part of what the
! function means.

pure subroutine generated_fci_quadratic_bezier_edge_area(x_start, y_start, control_x, control_y, &
        x_end, y_end, edge_area)
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    real(dp), intent(in) :: x_start, y_start, control_x, control_y, x_end, y_end
    real(dp), intent(out) :: edge_area
    real(dp) :: t1, t2, t3, t4

    t1 = control_y - y_start
    t2 = y_end + y_start - control_y*2.0000000000000000E+000_dp
    t3 = control_x - x_start
    t4 = x_end + x_start - control_x*2.0000000000000000E+000_dp
    edge_area = (x_start*t1*2.0000000000000000E+000_dp + x_start*t2 - y_start*t3* &
        2.0000000000000000E+000_dp - y_start*t4 + (t3*t2*2.0000000000000000E+000_dp - t1*t4* &
        2.0000000000000000E+000_dp)/3.0000000000000000E+000_dp)*5.0000000000000000E-001_dp

end subroutine generated_fci_quadratic_bezier_edge_area
