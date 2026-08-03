program test_triangle_piola_ad
    use check, only: check_condition, check_summary
    use fortfem_feec, only: &
        map_triangle_nedelec_covariant, &
        map_triangle_nedelec_covariant_jvp, &
        map_triangle_nedelec_covariant_vjp, &
        map_triangle_rt_contravariant, &
        map_triangle_rt_contravariant_jvp, &
        map_triangle_rt_contravariant_vjp
    use fortfem_kinds, only: dp
    implicit none

    integer, parameter :: count = 4
    real(dp), parameter :: step = 2.0e-7_dp
    real(dp) :: jacobian(2, 2), jacobian_bar(2, 2), jacobian_dot(2, 2)
    real(dp) :: reference_values(2, count), reference_values_bar(2, count)
    real(dp) :: reference_values_dot(2, count)
    real(dp) :: reference_scalars(count), reference_scalars_bar(count)
    real(dp) :: reference_scalars_dot(count)
    real(dp) :: values(2, count), values_bar(2, count), values_dot(2, count)
    real(dp) :: values_minus(2, count), values_plus(2, count)
    real(dp) :: scalars(count), scalars_bar(count), scalars_dot(count)
    real(dp) :: scalars_minus(count), scalars_plus(count)
    real(dp) :: lhs, rhs
    integer :: status, status_minus, status_plus

    jacobian = reshape([1.2_dp, 0.15_dp, -0.25_dp, 0.9_dp], [2, 2])
    jacobian_dot = reshape([0.03_dp, -0.02_dp, 0.015_dp, 0.025_dp], [2, 2])
    reference_values = reshape([ &
        0.2_dp, -0.4_dp, 0.7_dp, 0.1_dp, &
        -0.3_dp, 0.6_dp, 0.25_dp, -0.5_dp], [2, count])
    reference_values_dot = reshape([ &
        0.01_dp, -0.03_dp, 0.02_dp, 0.04_dp, &
        -0.015_dp, 0.025_dp, -0.02_dp, 0.03_dp], [2, count])
    reference_scalars = [0.3_dp, -0.6_dp, 0.8_dp, 0.15_dp]
    reference_scalars_dot = [-0.02_dp, 0.04_dp, 0.01_dp, -0.03_dp]
    values_bar = reshape([ &
        0.05_dp, -0.02_dp, 0.03_dp, 0.07_dp, &
        -0.04_dp, 0.06_dp, -0.01_dp, 0.02_dp], [2, count])
    scalars_bar = [-0.1_dp, 0.08_dp, 0.03_dp, -0.05_dp]

    call check_rt()
    call check_nedelec()
    call check_summary("Triangle Piola derivatives")

contains

    subroutine check_rt()
        call map_triangle_rt_contravariant( &
            jacobian, reference_values, reference_scalars, values, scalars, &
            status)
        call map_triangle_rt_contravariant_jvp( &
            jacobian, reference_values, reference_scalars, jacobian_dot, &
            reference_values_dot, reference_scalars_dot, values_dot, &
            scalars_dot, status)
        call map_triangle_rt_contravariant( &
            jacobian + step*jacobian_dot, &
            reference_values + step*reference_values_dot, &
            reference_scalars + step*reference_scalars_dot, values_plus, &
            scalars_plus, status_plus)
        call map_triangle_rt_contravariant( &
            jacobian - step*jacobian_dot, &
            reference_values - step*reference_values_dot, &
            reference_scalars - step*reference_scalars_dot, values_minus, &
            scalars_minus, status_minus)
        call check_products("RT")
        call map_triangle_rt_contravariant_vjp( &
            jacobian, reference_values, reference_scalars, values_bar, &
            scalars_bar, jacobian_bar, reference_values_bar, &
            reference_scalars_bar, status)
        call check_adjoint("RT")
    end subroutine check_rt

    subroutine check_nedelec()
        call map_triangle_nedelec_covariant( &
            jacobian, reference_values, reference_scalars, values, scalars, &
            status)
        call map_triangle_nedelec_covariant_jvp( &
            jacobian, reference_values, reference_scalars, jacobian_dot, &
            reference_values_dot, reference_scalars_dot, values_dot, &
            scalars_dot, status)
        call map_triangle_nedelec_covariant( &
            jacobian + step*jacobian_dot, &
            reference_values + step*reference_values_dot, &
            reference_scalars + step*reference_scalars_dot, values_plus, &
            scalars_plus, status_plus)
        call map_triangle_nedelec_covariant( &
            jacobian - step*jacobian_dot, &
            reference_values - step*reference_values_dot, &
            reference_scalars - step*reference_scalars_dot, values_minus, &
            scalars_minus, status_minus)
        call check_products("Nedelec")
        call map_triangle_nedelec_covariant_vjp( &
            jacobian, reference_values, reference_scalars, values_bar, &
            scalars_bar, jacobian_bar, reference_values_bar, &
            reference_scalars_bar, status)
        call check_adjoint("Nedelec")
    end subroutine check_nedelec

    subroutine check_products(label)
        character(*), intent(in) :: label

        call check_condition( &
            status == 0 .and. status_plus == 0 .and. status_minus == 0, &
            label//" Piola JVP accepts a valid direction")
        call check_condition(maxval(abs( &
            values_dot - (values_plus - values_minus)/(2.0_dp*step))) < &
            2.0e-9_dp, label//" Piola value JVP matches finite differences")
        call check_condition(maxval(abs( &
            scalars_dot - (scalars_plus - scalars_minus)/(2.0_dp*step))) < &
            2.0e-9_dp, label//" Piola scalar JVP matches finite differences")
    end subroutine check_products

    subroutine check_adjoint(label)
        character(*), intent(in) :: label

        lhs = sum(values_bar*values_dot) + dot_product(scalars_bar, scalars_dot)
        rhs = sum(jacobian_bar*jacobian_dot) + &
            sum(reference_values_bar*reference_values_dot) + &
            dot_product(reference_scalars_bar, reference_scalars_dot)
        call check_condition(status == 0, label//" Piola VJP succeeds")
        call check_condition( &
            abs(lhs - rhs) < 2.0e-12_dp*max(1.0_dp, abs(lhs), abs(rhs)), &
            label//" Piola products obey the adjoint identity")
    end subroutine check_adjoint

end program test_triangle_piola_ad
