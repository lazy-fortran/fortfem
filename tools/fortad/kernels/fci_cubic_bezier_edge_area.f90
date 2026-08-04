! Extracted from src/generated/fortfem_fci_cubic_bezier_edge_area_products.f90 by
! tools/fortad/extract.py. It is the primal fortsym generated
! its products from, with the offload directives removed:
! they are the consumer's choice, not part of what the
! function means.

pure subroutine generated_fci_cubic_bezier_edge_area(x_start, y_start, control_1_x, control_1_y, &
        control_2_x, control_2_y, x_end, y_end, edge_area)
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    real(dp), intent(in) :: x_start, y_start, control_1_x, control_1_y, control_2_x, control_2_y, &
        x_end, y_end
    real(dp), intent(out) :: edge_area
    real(dp) :: t1, t2, t3, t4, t5, t6, t7, t8, t9

    t1 = 3.0000000000000000E+000_dp**2
    t2 = -x_start
    t3 = control_1_x - x_start
    t4 = control_2_y + y_start - control_1_y*2.0000000000000000E+000_dp
    t5 = -y_start
    t6 = control_1_y - y_start
    t7 = control_2_x + x_start - control_1_x*2.0000000000000000E+000_dp
    t8 = y_end - y_start + (control_1_y - control_2_y)*3.0000000000000000E+000_dp
    t9 = x_end - x_start + (control_1_x - control_2_x)*3.0000000000000000E+000_dp
    edge_area = (t1*t3*t4 - t1*t6*t7)*5.0000000000000000E-001_dp*2.0000000000000000E+000_dp/ &
        3.0000000000000000E+000_dp + (-t1*t3*t4 + t1*t6*t7)*5.0000000000000000E-001_dp* &
        1.0000000000000000E+000_dp/3.0000000000000000E+000_dp + (-t3*t8*3.0000000000000000E+000_dp + t6*t9* &
        3.0000000000000000E+000_dp)*5.0000000000000000E-001_dp*1.0000000000000000E+000_dp/ &
        4.0000000000000000E+000_dp + (t3*t8*3.0000000000000000E+000_dp - t6*t9*3.0000000000000000E+000_dp)* &
        5.0000000000000000E-001_dp*3.0000000000000000E+000_dp/4.0000000000000000E+000_dp + (-t7*t8* &
        3.0000000000000000E+000_dp + t4*t9*3.0000000000000000E+000_dp)*5.0000000000000000E-001_dp* &
        2.0000000000000000E+000_dp/5.0000000000000000E+000_dp + (t7*t8*3.0000000000000000E+000_dp - t4*t9* &
        3.0000000000000000E+000_dp)*5.0000000000000000E-001_dp*3.0000000000000000E+000_dp/ &
        5.0000000000000000E+000_dp + (x_start*t6*3.0000000000000000E+000_dp - y_start*t3* &
        3.0000000000000000E+000_dp)*5.0000000000000000E-001_dp*1.0000000000000000E+000_dp + (x_start*t4* &
        3.0000000000000000E+000_dp - y_start*t7*3.0000000000000000E+000_dp)*5.0000000000000000E-001_dp + &
        (x_start*t8 - y_start*t9)*5.0000000000000000E-001_dp + 0.0000000000000000E+000_dp

end subroutine generated_fci_cubic_bezier_edge_area
