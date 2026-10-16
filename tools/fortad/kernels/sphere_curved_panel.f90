! Extracted from src/generated/fortfem_sphere_curved_panel_products.f90 by
! tools/fortad/extract.py. It is the primal fortsym generated
! its products from, with the offload directives removed:
! they are the consumer's choice, not part of what the
! function means.

pure subroutine generated_sphere_curved_panel(vertex_11, vertex_21, vertex_31, vertex_12, vertex_22, &
        vertex_32, vertex_13, vertex_23, vertex_33, radius, xi, eta, point_1, point_2, point_3, &
        tangent_xi_1, tangent_xi_2, tangent_xi_3, tangent_eta_1, tangent_eta_2, tangent_eta_3, &
        surface_jacobian)
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    real(dp), intent(in) :: vertex_11, vertex_21, vertex_31, vertex_12, vertex_22, vertex_32, &
        vertex_13, vertex_23, vertex_33, radius, xi, eta
    real(dp), intent(out) :: point_1, point_2, point_3, tangent_xi_1, tangent_xi_2, tangent_xi_3, &
        tangent_eta_1, tangent_eta_2, tangent_eta_3, surface_jacobian
    real(dp) :: t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13, t14, t15, t16, t17, t18, &
        t19, t20, t21, t22, t23, t24

    t1 = -vertex_11
    t2 = vertex_13 - vertex_11
    t3 = vertex_12 - vertex_11
    t4 = vertex_11 + eta*t2 + xi*t3
    t5 = -vertex_21
    t6 = vertex_23 - vertex_21
    t7 = vertex_22 - vertex_21
    t8 = vertex_21 + eta*t6 + xi*t7
    t9 = -vertex_31
    t10 = vertex_33 - vertex_31
    t11 = vertex_32 - vertex_31
    t12 = vertex_31 + eta*t10 + xi*t11
    t13 = sqrt(t4*t4 + t8*t8 + t12*t12)
    t14 = t4*t3 + t8*t7 + t12*t11
    t15 = t3/t13 - t4*t14/(t13*t13*t13)
    t16 = t7/t13 - t8*t14/(t13*t13*t13)
    t17 = t11/t13 - t12*t14/(t13*t13*t13)
    t18 = t4*t2 + t8*t6 + t12*t10
    t19 = t2/t13 - t4*t18/(t13*t13*t13)
    t20 = t6/t13 - t8*t18/(t13*t13*t13)
    t21 = t10/t13 - t12*t18/(t13*t13*t13)
    t22 = radius*radius*t15*t20 - radius*radius*t19*t16
    t23 = -radius*radius*t15*t21 + radius*radius*t19*t17
    t24 = radius*radius*t16*t21 - radius*radius*t20*t17
    point_1 = radius*t4/t13
    point_2 = radius*t8/t13
    point_3 = radius*t12/t13
    tangent_xi_1 = radius*t15
    tangent_xi_2 = radius*t16
    tangent_xi_3 = radius*t17
    tangent_eta_1 = radius*t19
    tangent_eta_2 = radius*t20
    tangent_eta_3 = radius*t21
    surface_jacobian = sqrt(t22*t22 + t23*t23 + t24*t24)

end subroutine generated_sphere_curved_panel
