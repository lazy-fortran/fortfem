! Extracted from src/generated/fortfem_fci_parallel_products.f90 by
! tools/fortad/extract.py. It is the primal fortsym generated
! its products from, with the offload directives removed:
! they are the consumer's choice, not part of what the
! function means.

pure subroutine generated_fci_cubic_lagrange_weights(target, node_1, node_2, node_3, node_4, &
        weight_1, weight_2, weight_3, weight_4)
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    real(dp), intent(in) :: target, node_1, node_2, node_3, node_4
    real(dp), intent(out) :: weight_1, weight_2, weight_3, weight_4
    real(dp) :: t1, t2, t3, t4, t5, t6, t7, t8

    t1 = -node_2
    t2 = -node_3
    t3 = -node_4
    t4 = target - node_2
    t5 = target - node_3
    t6 = target - node_4
    t7 = -node_1
    t8 = target - node_1
    weight_1 = t4*t5*t6/(node_1 - node_2)/(node_1 - node_3)/(node_1 - node_4)
    weight_2 = t8*t5*t6/(node_2 - node_1)/(node_2 - node_3)/(node_2 - node_4)
    weight_3 = t8*t4*t6/(node_3 - node_1)/(node_3 - node_2)/(node_3 - node_4)
    weight_4 = t8*t4*t5/(node_4 - node_1)/(node_4 - node_2)/(node_4 - node_3)

end subroutine generated_fci_cubic_lagrange_weights
