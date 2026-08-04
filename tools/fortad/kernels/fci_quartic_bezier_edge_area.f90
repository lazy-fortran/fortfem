! Extracted from src/generated/fortfem_fci_quartic_bezier_edge_area_products.f90 by
! tools/fortad/extract.py. It is the primal fortsym generated
! its products from, with the offload directives removed:
! they are the consumer's choice, not part of what the
! function means.

pure subroutine generated_fci_quartic_bezier_edge_area(x_start, y_start, control_1_x, control_1_y, &
        control_2_x, control_2_y, control_3_x, control_3_y, x_end, y_end, edge_area)
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    real(dp), intent(in) :: x_start, y_start, control_1_x, control_1_y, control_2_x, control_2_y, &
        control_3_x, control_3_y, x_end, y_end
    real(dp), intent(out) :: edge_area
    real(dp) :: t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11

    t1 = -x_start
    t2 = control_1_x - x_start
    t3 = control_2_y + y_start - control_1_y*2.0000000000000000E+000_dp
    t4 = -y_start
    t5 = control_1_y - y_start
    t6 = control_2_x + x_start - control_1_x*2.0000000000000000E+000_dp
    t7 = 4.0000000000000000E+000_dp**2
    t8 = control_3_y + control_1_y*3.0000000000000000E+000_dp - control_2_y* &
        3.0000000000000000E+000_dp - y_start
    t9 = control_3_x + control_1_x*3.0000000000000000E+000_dp - control_2_x* &
        3.0000000000000000E+000_dp - x_start
    t10 = y_end + y_start - control_1_y*4.0000000000000000E+000_dp + control_2_y* &
        6.0000000000000000E+000_dp - control_3_y*4.0000000000000000E+000_dp
    t11 = x_end + x_start - control_1_x*4.0000000000000000E+000_dp + control_2_x* &
        6.0000000000000000E+000_dp - control_3_x*4.0000000000000000E+000_dp
    edge_area = (-t2*t3*4.0000000000000000E+000_dp*6.0000000000000000E+000_dp + t5*t6* &
        4.0000000000000000E+000_dp*6.0000000000000000E+000_dp)*5.0000000000000000E-001_dp* &
        1.0000000000000000E+000_dp/3.0000000000000000E+000_dp + (t2*t3*4.0000000000000000E+000_dp* &
        6.0000000000000000E+000_dp - t5*t6*4.0000000000000000E+000_dp*6.0000000000000000E+000_dp)* &
        5.0000000000000000E-001_dp*2.0000000000000000E+000_dp/3.0000000000000000E+000_dp + (t7*t2*t8 - t7* &
        t5*t9)*5.0000000000000000E-001_dp*3.0000000000000000E+000_dp/4.0000000000000000E+000_dp + (-t7*t2* &
        t8 + t7*t5*t9)*5.0000000000000000E-001_dp*1.0000000000000000E+000_dp/4.0000000000000000E+000_dp + (- &
        t2*t10*4.0000000000000000E+000_dp + t5*t11*4.0000000000000000E+000_dp)*5.0000000000000000E-001_dp* &
        1.0000000000000000E+000_dp/5.0000000000000000E+000_dp + (t2*t10*4.0000000000000000E+000_dp - t5*t11* &
        4.0000000000000000E+000_dp)*5.0000000000000000E-001_dp*4.0000000000000000E+000_dp/ &
        5.0000000000000000E+000_dp + (-t6*t8*4.0000000000000000E+000_dp*6.0000000000000000E+000_dp + t3*t9* &
        4.0000000000000000E+000_dp*6.0000000000000000E+000_dp)*5.0000000000000000E-001_dp* &
        2.0000000000000000E+000_dp/5.0000000000000000E+000_dp + (t6*t8*4.0000000000000000E+000_dp* &
        6.0000000000000000E+000_dp - t3*t9*4.0000000000000000E+000_dp*6.0000000000000000E+000_dp)* &
        5.0000000000000000E-001_dp*3.0000000000000000E+000_dp/5.0000000000000000E+000_dp + (-t6*t10* &
        6.0000000000000000E+000_dp + t3*t11*6.0000000000000000E+000_dp)*5.0000000000000000E-001_dp* &
        2.0000000000000000E+000_dp/6.0000000000000000E+000_dp + (t6*t10*6.0000000000000000E+000_dp - t3*t11* &
        6.0000000000000000E+000_dp)*5.0000000000000000E-001_dp*4.0000000000000000E+000_dp/ &
        6.0000000000000000E+000_dp + (-t9*t10*4.0000000000000000E+000_dp + t8*t11* &
        4.0000000000000000E+000_dp)*5.0000000000000000E-001_dp*3.0000000000000000E+000_dp/ &
        7.0000000000000000E+000_dp + (t9*t10*4.0000000000000000E+000_dp - t8*t11*4.0000000000000000E+000_dp) &
        *5.0000000000000000E-001_dp*4.0000000000000000E+000_dp/7.0000000000000000E+000_dp + (x_start*t5* &
        4.0000000000000000E+000_dp - y_start*t2*4.0000000000000000E+000_dp)*5.0000000000000000E-001_dp* &
        1.0000000000000000E+000_dp + (x_start*t3*6.0000000000000000E+000_dp - y_start*t6* &
        6.0000000000000000E+000_dp)*5.0000000000000000E-001_dp + (x_start*t8*4.0000000000000000E+000_dp - &
        y_start*t9*4.0000000000000000E+000_dp)*5.0000000000000000E-001_dp + (x_start*t10 - y_start*t11)* &
        5.0000000000000000E-001_dp + 0.0000000000000000E+000_dp

end subroutine generated_fci_quartic_bezier_edge_area
