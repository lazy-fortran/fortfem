! Extracted from src/generated/fortfem_fci_quartic_lagrange.f90 by
! tools/fortad/extract.py. It is the primal fortsym generated
! its products from, with the offload directives removed:
! they are the consumer's choice, not part of what the
! function means.

pure subroutine generated_fci_quartic_lagrange_weights(target, node_1, node_2, node_3, node_4, &
        node_5, weight_1, weight_2, weight_3, weight_4, weight_5)
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    real(dp), intent(in) :: target, node_1, node_2, node_3, node_4, node_5
    real(dp), intent(out) :: weight_1, weight_2, weight_3, weight_4, weight_5
    real(dp) :: t1, t2, t3, t4, t5, t6, t7, t8, t9, t10

    t1 = -node_2
    t2 = -node_3
    t3 = -node_4
    t4 = -node_5
    t5 = target - node_2
    t6 = target - node_3
    t7 = target - node_4
    t8 = target - node_5
    t9 = -node_1
    t10 = target - node_1
    weight_1 = t5*t6*t7*t8/(node_1 - node_2)/(node_1 - node_3)/(node_1 - node_4)/(node_1 - node_5)
    weight_2 = t10*t6*t7*t8/(node_2 - node_1)/(node_2 - node_3)/(node_2 - node_4)/(node_2 - node_5)
    weight_3 = t10*t5*t7*t8/(node_3 - node_1)/(node_3 - node_2)/(node_3 - node_4)/(node_3 - node_5)
    weight_4 = t10*t5*t6*t8/(node_4 - node_1)/(node_4 - node_2)/(node_4 - node_3)/(node_4 - node_5)
    weight_5 = t10*t5*t6*t7/(node_5 - node_1)/(node_5 - node_2)/(node_5 - node_3)/(node_5 - node_4)

end subroutine generated_fci_quartic_lagrange_weights
