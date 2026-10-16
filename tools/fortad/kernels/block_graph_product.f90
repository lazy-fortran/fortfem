! Extracted from src/generated/fortfem_block_graph_products.f90 by
! tools/fortad/extract.py. It is the primal fortsym generated
! its products from, with the offload directives removed:
! they are the consumer's choice, not part of what the
! function means.

pure subroutine generated_block_graph_product(block_value, state_value, contribution)
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    real(dp), intent(in) :: block_value, state_value
    real(dp), intent(out) :: contribution

    contribution = block_value*state_value

end subroutine generated_block_graph_product
