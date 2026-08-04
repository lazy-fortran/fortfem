! Extracted from src/generated/fortfem_fci_parallel_products.f90 by
! tools/fortad/extract.py. It is the primal fortsym generated
! its products from, with the offload directives removed:
! they are the consumer's choice, not part of what the
! function means.

pure subroutine generated_fci_quadratic_lagrange_weights(target, node_1, node_2, node_3, weight_1, &
        weight_2, weight_3)
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    real(dp), intent(in) :: target, node_1, node_2, node_3
    real(dp), intent(out) :: weight_1, weight_2, weight_3
    real(dp) :: t1, t2, t3, t4, t5, t6

    t1 = -node_2
    t2 = -node_3
    t3 = target - node_2
    t4 = target - node_3
    t5 = -node_1
    t6 = target - node_1
    weight_1 = t3*t4/(node_1 - node_2)/(node_1 - node_3)
    weight_2 = t6*t4/(node_2 - node_1)/(node_2 - node_3)
    weight_3 = t6*t3/(node_3 - node_1)/(node_3 - node_2)

end subroutine generated_fci_quadratic_lagrange_weights
