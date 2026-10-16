! Extracted from src/generated/fortfem_fci_curved_quadrilateral_area_products.f90 by
! tools/fortad/extract.py. It is the primal fortsym generated
! its products from, with the offload directives removed:
! they are the consumer's choice, not part of what the
! function means.

pure subroutine generated_fci_curved_quadrilateral_cell_area(x_1, y_1, x_2, y_2, x_3, y_3, x_4, y_4, &
        control_x_1, control_y_1, control_x_2, control_y_2, control_x_3, control_y_3, control_x_4, &
        control_y_4, area)
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    real(dp), intent(in) :: x_1, y_1, x_2, y_2, x_3, y_3, x_4, y_4, control_x_1, control_y_1, &
        control_x_2, control_y_2, control_x_3, control_y_3, control_x_4, control_y_4
    real(dp), intent(out) :: area
    real(dp) :: t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13, t14, t15, t16

    t1 = control_y_1 - y_1
    t2 = y_1 + y_2 - control_y_1*2.0000000000000000E+000_dp
    t3 = control_x_1 - x_1
    t4 = x_1 + x_2 - control_x_1*2.0000000000000000E+000_dp
    t5 = control_y_2 - y_2
    t6 = y_2 + y_3 - control_y_2*2.0000000000000000E+000_dp
    t7 = control_x_2 - x_2
    t8 = x_2 + x_3 - control_x_2*2.0000000000000000E+000_dp
    t9 = control_y_3 - y_3
    t10 = y_3 + y_4 - control_y_3*2.0000000000000000E+000_dp
    t11 = control_x_3 - x_3
    t12 = x_3 + x_4 - control_x_3*2.0000000000000000E+000_dp
    t13 = control_y_4 - y_4
    t14 = y_1 + y_4 - control_y_4*2.0000000000000000E+000_dp
    t15 = control_x_4 - x_4
    t16 = x_1 + x_4 - control_x_4*2.0000000000000000E+000_dp
    area = (x_1*t1*2.0000000000000000E+000_dp + x_1*t2 - y_1*t3*2.0000000000000000E+000_dp - y_1* &
        t4 + (t3*t2*2.0000000000000000E+000_dp - t1*t4*2.0000000000000000E+000_dp)/ &
        3.0000000000000000E+000_dp)*5.0000000000000000E-001_dp + (x_2*t5*2.0000000000000000E+000_dp + x_2* &
        t6 - y_2*t7*2.0000000000000000E+000_dp - y_2*t8 + (t7*t6*2.0000000000000000E+000_dp - t5*t8* &
        2.0000000000000000E+000_dp)/3.0000000000000000E+000_dp)*5.0000000000000000E-001_dp + (x_3*t9* &
        2.0000000000000000E+000_dp + x_3*t10 - y_3*t11*2.0000000000000000E+000_dp - y_3*t12 + (t11*t10* &
        2.0000000000000000E+000_dp - t9*t12*2.0000000000000000E+000_dp)/3.0000000000000000E+000_dp)* &
        5.0000000000000000E-001_dp + (x_4*t13*2.0000000000000000E+000_dp + x_4*t14 - y_4*t15* &
        2.0000000000000000E+000_dp - y_4*t16 + (t15*t14*2.0000000000000000E+000_dp - t13*t16* &
        2.0000000000000000E+000_dp)/3.0000000000000000E+000_dp)*5.0000000000000000E-001_dp + &
        0.0000000000000000E+000_dp

end subroutine generated_fci_curved_quadrilateral_cell_area
