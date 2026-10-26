module fortfem_triangle_feec_operators
    !! Discrete gradient P_k -> N_k on the reference triangle: the Nedelec
    !! degrees of freedom of the Lagrange gradients. The table is computed in
    !! exact rational arithmetic by FortSym and rounded once
    !! (tools/codegen/app/gen_triangle_feec_tables.f90).
    use fortfem_kinds, only: dp
    use fortfem_generated_triangle_feec_tables, only: &
        triangle_feec_operator_tables
    implicit none

    private

    public :: build_triangle_discrete_gradient

contains

    subroutine build_triangle_discrete_gradient(order, matrix, status)
        !! status = 1 for order < 1 or order above the generated table range.
        integer, intent(in) :: order
        real(dp), allocatable, intent(out) :: matrix(:, :)
        integer, intent(out) :: status

        real(dp), allocatable :: curl(:, :), divergence(:, :)

        status = 1
        if (order < 1) return
        call triangle_feec_operator_tables(order, matrix, curl, divergence, &
            status)
    end subroutine build_triangle_discrete_gradient

end module fortfem_triangle_feec_operators
