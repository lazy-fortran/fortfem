! Extracted from src/generated/fortfem_surface_triangle_geometry_products.f90 by
! tools/fortad/extract.py. It is the primal fortsym generated
! its products from, with the offload directives removed:
! they are the consumer's choice, not part of what the
! function means.

pure subroutine generated_surface_triangle_geometry_3d(vertex_11, vertex_21, vertex_31, vertex_12, &
        vertex_22, vertex_32, vertex_13, vertex_23, vertex_33, xi, eta, point_1, point_2, point_3, jacobian, &
        normal_1, normal_2, normal_3)
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    real(dp), intent(in) :: vertex_11, vertex_21, vertex_31, vertex_12, vertex_22, vertex_32, &
        vertex_13, vertex_23, vertex_33, xi, eta
    real(dp), intent(out) :: point_1, point_2, point_3, jacobian, normal_1, normal_2, normal_3
    real(dp) :: t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13

    t1 = -vertex_11
    t2 = vertex_13 - vertex_11
    t3 = vertex_12 - vertex_11
    t4 = -vertex_21
    t5 = vertex_23 - vertex_21
    t6 = vertex_22 - vertex_21
    t7 = -vertex_31
    t8 = vertex_33 - vertex_31
    t9 = vertex_32 - vertex_31
    t10 = t3*t5 - t2*t6
    t11 = -t3*t8 + t2*t9
    t12 = t6*t8 - t5*t9
    t13 = sqrt(t10*t10 + t11*t11 + t12*t12)
    point_1 = vertex_11 + eta*t2 + xi*t3
    point_2 = vertex_21 + eta*t5 + xi*t6
    point_3 = vertex_31 + eta*t8 + xi*t9
    jacobian = t13
    normal_1 = t12/t13
    normal_2 = t11/t13
    normal_3 = t10/t13

end subroutine generated_surface_triangle_geometry_3d
