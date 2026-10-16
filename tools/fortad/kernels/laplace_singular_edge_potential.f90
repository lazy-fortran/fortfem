! Extracted from src/generated/fortfem_laplace_singular_edge_products.f90 by
! tools/fortad/extract.py. It is the primal fortsym generated
! its products from, with the offload directives removed:
! they are the consumer's choice, not part of what the
! function means.

pure subroutine generated_laplace_singular_edge_potential(point_1, point_2, point_3, first_vertex_1, &
        first_vertex_2, first_vertex_3, second_vertex_1, second_vertex_2, second_vertex_3, value)
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    real(dp), intent(in) :: point_1, point_2, point_3, first_vertex_1, first_vertex_2, &
        first_vertex_3, second_vertex_1, second_vertex_2, second_vertex_3
    real(dp), intent(out) :: value
    real(dp) :: t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13, t14, t15, t16, t17, t18

    t1 = -point_1
    t2 = first_vertex_1 - point_1
    t3 = -point_2
    t4 = first_vertex_2 - point_2
    t5 = -point_3
    t6 = first_vertex_3 - point_3
    t7 = sqrt(t2*t2 + t4*t4 + t6*t6)
    t8 = second_vertex_1 - point_1
    t9 = second_vertex_2 - point_2
    t10 = second_vertex_3 - point_3
    t11 = sqrt(t8*t8 + t9*t9 + t10*t10)
    t12 = second_vertex_1 - point_1 - t2
    t13 = second_vertex_2 - point_2 - t4
    t14 = second_vertex_3 - point_3 - t6
    t15 = sqrt(t12*t12 + t13*t13 + t14*t14)
    t16 = t2*t9 - t4*t8
    t17 = -t2*t10 + t6*t8
    t18 = t4*t10 - t6*t9
    value = log((t7 + t11 + t15)/(t7 + t11 - t15))*sqrt(t16*t16 + t17*t17 + t18*t18)/t15

end subroutine generated_laplace_singular_edge_potential
