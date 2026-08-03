program test_tetra_rt_piola_ad
    use check, only: check_condition, check_summary
    use fortfem_feec, only: map_tetra_rt_contravariant, &
        map_tetra_rt_contravariant_jvp, map_tetra_rt_contravariant_vjp
    use fortfem_kinds, only: dp
    implicit none

    integer, parameter :: count = 5
    real(dp), parameter :: step = 2.0e-7_dp
    real(dp) :: jacobian(3, 3), jacobian_dot(3, 3), jacobian_bar(3, 3)
    real(dp) :: reference_values(3, count), reference_values_dot(3, count)
    real(dp) :: reference_values_bar(3, count)
    real(dp) :: reference_divergences(count), reference_divergences_dot(count)
    real(dp) :: reference_divergences_bar(count)
    real(dp) :: physical_values(3, count), physical_values_dot(3, count)
    real(dp) :: physical_values_bar(3, count)
    real(dp) :: physical_divergences(count), physical_divergences_dot(count)
    real(dp) :: physical_divergences_bar(count)
    real(dp) :: values_plus(3, count), values_minus(3, count)
    real(dp) :: divergences_plus(count), divergences_minus(count)
    real(dp) :: lhs, rhs, relative_error
    integer :: basis, component, status, status_plus, status_minus

    jacobian = reshape([ &
        1.2_dp, 0.1_dp, -0.05_dp, &
        0.2_dp, 0.95_dp, 0.08_dp, &
        -0.1_dp, 0.15_dp, 1.1_dp], [3, 3])
    jacobian_dot = reshape([ &
        0.03_dp, -0.02_dp, 0.01_dp, &
        -0.01_dp, 0.04_dp, -0.03_dp, &
        0.02_dp, 0.01_dp, -0.02_dp], [3, 3])
    do basis = 1, count
        reference_divergences(basis) = 0.08_dp*basis - 0.03_dp
        reference_divergences_dot(basis) = -0.02_dp*basis + 0.01_dp
        physical_divergences_bar(basis) = 0.03_dp*basis - 0.015_dp
        do component = 1, 3
            reference_values(component, basis) = &
                0.07_dp*component + 0.03_dp*basis
            reference_values_dot(component, basis) = &
                -0.02_dp*component + 0.01_dp*basis
            physical_values_bar(component, basis) = &
                0.03_dp*component - 0.02_dp*basis
        end do
    end do

    call map_tetra_rt_contravariant( &
        jacobian, reference_values, reference_divergences, physical_values, &
        physical_divergences, status)
    call map_tetra_rt_contravariant_jvp( &
        jacobian, reference_values, reference_divergences, jacobian_dot, &
        reference_values_dot, reference_divergences_dot, physical_values_dot, &
        physical_divergences_dot, status)
    call map_tetra_rt_contravariant( &
        jacobian + step*jacobian_dot, &
        reference_values + step*reference_values_dot, &
        reference_divergences + step*reference_divergences_dot, values_plus, &
        divergences_plus, status_plus)
    call map_tetra_rt_contravariant( &
        jacobian - step*jacobian_dot, &
        reference_values - step*reference_values_dot, &
        reference_divergences - step*reference_divergences_dot, values_minus, &
        divergences_minus, status_minus)
    relative_error = max( &
        maxval(abs(physical_values_dot - &
        (values_plus - values_minus)/(2.0_dp*step))), &
        maxval(abs(physical_divergences_dot - &
        (divergences_plus - divergences_minus)/(2.0_dp*step))))/ &
        max(1.0_dp, maxval(abs(physical_values_dot)), &
        maxval(abs(physical_divergences_dot)))
    call check_condition( &
        status == 0 .and. status_plus == 0 .and. status_minus == 0, &
        "RT Piola JVP accepts a valid affine map")
    call check_condition( &
        relative_error < 2.0e-9_dp, &
        "RT Piola JVP matches a complete central difference")

    call map_tetra_rt_contravariant_vjp( &
        jacobian, reference_values, reference_divergences, &
        physical_values_bar, physical_divergences_bar, jacobian_bar, &
        reference_values_bar, reference_divergences_bar, status)
    lhs = sum(physical_values_bar*physical_values_dot) + &
        sum(physical_divergences_bar*physical_divergences_dot)
    rhs = sum(jacobian_bar*jacobian_dot) + &
        sum(reference_values_bar*reference_values_dot) + &
        sum(reference_divergences_bar*reference_divergences_dot)
    call check_condition(status == 0, "RT Piola VJP succeeds")
    call check_condition( &
        abs(lhs - rhs) < 2.0e-12_dp*max(1.0_dp, abs(lhs), abs(rhs)), &
        "RT Piola products satisfy the adjoint identity")

    call check_summary("Tetrahedral RT Piola derivatives")
end program test_tetra_rt_piola_ad
