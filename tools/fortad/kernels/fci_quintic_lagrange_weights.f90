! Extracted from src/generated/fortfem_fci_quintic_lagrange.f90 by
! tools/fortad/extract.py. It is the primal fortsym generated
! its products from, with the offload directives removed:
! they are the consumer's choice, not part of what the
! function means.

pure subroutine generated_fci_quintic_lagrange_weights(target, node_1, node_2, node_3, node_4, &
        node_5, node_6, weight_1, weight_2, weight_3, weight_4, weight_5, weight_6)
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    real(dp), intent(in) :: target, node_1, node_2, node_3, node_4, node_5, node_6
    real(dp), intent(out) :: weight_1, weight_2, weight_3, weight_4, weight_5, weight_6
    real(dp) :: t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12

    t1 = -node_2
    t2 = -node_3
    t3 = -node_4
    t4 = -node_5
    t5 = -node_6
    t6 = target - node_2
    t7 = target - node_3
    t8 = target - node_4
    t9 = target - node_5
    t10 = target - node_6
    t11 = -node_1
    t12 = target - node_1
    weight_1 = t6*t7*t8*t9*t10/(node_1 - node_2)/(node_1 - node_3)/(node_1 - node_4)/(node_1 - &
        node_5)/(node_1 - node_6)
    weight_2 = t12*t7*t8*t9*t10/(node_2 - node_1)/(node_2 - node_3)/(node_2 - node_4)/(node_2 - &
        node_5)/(node_2 - node_6)
    weight_3 = t12*t6*t8*t9*t10/(node_3 - node_1)/(node_3 - node_2)/(node_3 - node_4)/(node_3 - &
        node_5)/(node_3 - node_6)
    weight_4 = t12*t6*t7*t9*t10/(node_4 - node_1)/(node_4 - node_2)/(node_4 - node_3)/(node_4 - &
        node_5)/(node_4 - node_6)
    weight_5 = t12*t6*t7*t8*t10/(node_5 - node_1)/(node_5 - node_2)/(node_5 - node_3)/(node_5 - &
        node_4)/(node_5 - node_6)
    weight_6 = t12*t6*t7*t8*t9/(node_6 - node_1)/(node_6 - node_2)/(node_6 - node_3)/(node_6 - &
        node_4)/(node_6 - node_5)

end subroutine generated_fci_quintic_lagrange_weights
