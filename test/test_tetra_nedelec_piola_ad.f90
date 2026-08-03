program test_tetra_nedelec_piola_ad
    use check, only: check_condition, check_summary
    use fortfem_tetra_piola_maps, only: map_tetra_nedelec_covariant, &
        map_tetra_nedelec_covariant_jvp, map_tetra_nedelec_covariant_vjp
    use fortfem_kinds, only: dp
    implicit none

    integer, parameter :: count = 5
    real(dp), parameter :: step = 2.0e-7_dp
    real(dp) :: jacobian(3, 3), jacobian_dot(3, 3), jacobian_bar(3, 3)
    real(dp) :: reference_values(3, count), reference_values_dot(3, count)
    real(dp) :: reference_values_bar(3, count)
    real(dp) :: reference_curls(3, count), reference_curls_dot(3, count)
    real(dp) :: reference_curls_bar(3, count)
    real(dp) :: physical_values(3, count), physical_values_dot(3, count)
    real(dp) :: physical_values_bar(3, count)
    real(dp) :: physical_curls(3, count), physical_curls_dot(3, count)
    real(dp) :: physical_curls_bar(3, count)
    real(dp) :: values_plus(3, count), values_minus(3, count)
    real(dp) :: curls_plus(3, count), curls_minus(3, count)
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
        do component = 1, 3
            reference_values(component, basis) = &
                0.07_dp*component + 0.03_dp*basis
            reference_values_dot(component, basis) = &
                -0.02_dp*component + 0.01_dp*basis
            reference_curls(component, basis) = &
                -0.04_dp*component + 0.05_dp*basis
            reference_curls_dot(component, basis) = &
                0.015_dp*component - 0.008_dp*basis
            physical_values_bar(component, basis) = &
                0.03_dp*component - 0.02_dp*basis
            physical_curls_bar(component, basis) = &
                -0.01_dp*component + 0.025_dp*basis
        end do
    end do

    call map_tetra_nedelec_covariant( &
        jacobian, reference_values, reference_curls, physical_values, &
        physical_curls, status)
    call map_tetra_nedelec_covariant_jvp( &
        jacobian, reference_values, reference_curls, jacobian_dot, &
        reference_values_dot, reference_curls_dot, physical_values_dot, &
        physical_curls_dot, status)
    call map_tetra_nedelec_covariant( &
        jacobian + step*jacobian_dot, &
        reference_values + step*reference_values_dot, &
        reference_curls + step*reference_curls_dot, values_plus, curls_plus, &
        status_plus)
    call map_tetra_nedelec_covariant( &
        jacobian - step*jacobian_dot, &
        reference_values - step*reference_values_dot, &
        reference_curls - step*reference_curls_dot, values_minus, curls_minus, &
        status_minus)
    relative_error = max( &
        maxval(abs(physical_values_dot - &
        (values_plus - values_minus)/(2.0_dp*step))), &
        maxval(abs(physical_curls_dot - &
        (curls_plus - curls_minus)/(2.0_dp*step))))/ &
        max(1.0_dp, maxval(abs(physical_values_dot)), &
        maxval(abs(physical_curls_dot)))
    call check_condition( &
        status == 0 .and. status_plus == 0 .and. status_minus == 0, &
        "Nedelec Piola JVP accepts a valid affine map")
    call check_condition( &
        relative_error < 2.0e-9_dp, &
        "Nedelec Piola JVP matches a complete central difference")

    call map_tetra_nedelec_covariant_vjp( &
        jacobian, reference_values, reference_curls, physical_values_bar, &
        physical_curls_bar, jacobian_bar, reference_values_bar, &
        reference_curls_bar, status)
    lhs = sum(physical_values_bar*physical_values_dot) + &
        sum(physical_curls_bar*physical_curls_dot)
    rhs = sum(jacobian_bar*jacobian_dot) + &
        sum(reference_values_bar*reference_values_dot) + &
        sum(reference_curls_bar*reference_curls_dot)
    call check_condition(status == 0, "Nedelec Piola VJP succeeds")
    call check_condition( &
        abs(lhs - rhs) < 2.0e-12_dp*max(1.0_dp, abs(lhs), abs(rhs)), &
        "Nedelec Piola products satisfy the adjoint identity")

    call check_summary("Tetrahedral Nedelec Piola derivatives")
end program test_tetra_nedelec_piola_ad
