! Extracted from src/generated/fortfem_torus_curved_panel_products.f90 by
! tools/fortad/extract.py. It is the primal fortsym generated
! its products from, with the offload directives removed:
! they are the consumer's choice, not part of what the
! function means.

pure subroutine generated_torus_curved_panel(parameter_11, parameter_21, parameter_12, parameter_22, &
        parameter_13, parameter_23, major_radius, minor_radius, xi, eta, point_1, point_2, point_3, &
        tangent_xi_1, tangent_xi_2, tangent_xi_3, tangent_eta_1, tangent_eta_2, tangent_eta_3, &
        surface_jacobian)
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    real(dp), intent(in) :: parameter_11, parameter_21, parameter_12, parameter_22, parameter_13, &
        parameter_23, major_radius, minor_radius, xi, eta
    real(dp), intent(out) :: point_1, point_2, point_3, tangent_xi_1, tangent_xi_2, tangent_xi_3, &
        tangent_eta_1, tangent_eta_2, tangent_eta_3, surface_jacobian
    real(dp) :: t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13, t14, t15, t16, t17, t18, &
        t19, t20, t21, t22

    t1 = -parameter_21
    t2 = parameter_23 - parameter_21
    t3 = parameter_22 - parameter_21
    t4 = parameter_21 + eta*t2 + xi*t3
    t5 = cos(t4)
    t6 = -parameter_11
    t7 = parameter_13 - parameter_11
    t8 = parameter_12 - parameter_11
    t9 = parameter_11 + eta*t7 + xi*t8
    t10 = cos(t9)
    t11 = major_radius + minor_radius*t10
    t12 = sin(t4)
    t13 = sin(t9)
    t14 = -minor_radius*t5*t13*t8 - t12*t11*t3
    t15 = -minor_radius*t13*t12*t8 + t5*t11*t3
    t16 = minor_radius*t10*t8 + t3*0.0000000000000000E+000_dp
    t17 = -minor_radius*t5*t13*t7 - t12*t11*t2
    t18 = -minor_radius*t13*t12*t7 + t5*t11*t2
    t19 = minor_radius*t10*t7 + t2*0.0000000000000000E+000_dp
    t20 = t16*t17 - t19*t14
    t21 = -t16*t18 + t19*t15
    t22 = t14*t18 - t17*t15
    point_1 = t5*t11
    point_2 = t12*t11
    point_3 = minor_radius*t13
    tangent_xi_1 = t14
    tangent_xi_2 = t15
    tangent_xi_3 = t16
    tangent_eta_1 = t17
    tangent_eta_2 = t18
    tangent_eta_3 = t19
    surface_jacobian = sqrt(t20*t20 + t21*t21 + t22*t22)

end subroutine generated_torus_curved_panel
