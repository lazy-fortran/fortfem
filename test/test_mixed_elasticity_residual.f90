program test_mixed_elasticity_residual
    use, intrinsic :: iso_fortran_env, only: real64
    use check, only: check_condition, check_summary
    use fortfem_api, only: &
        assemble_mixed_elasticity_residual, &
        assemble_mixed_elasticity_residual_jvp, &
        assemble_mixed_elasticity_residual_vjp
    implicit none

    integer, parameter :: dp = real64
    integer, parameter :: stress_count = 3, displacement_count = 2
    integer, parameter :: equilibrium_count = 2
    real(dp) :: compliance(stress_count, stress_count)
    real(dp) :: strain_map(stress_count, displacement_count)
    real(dp) :: divergence_map(equilibrium_count, stress_count)
    real(dp) :: stress(stress_count), displacement(displacement_count)
    real(dp) :: load(equilibrium_count)
    real(dp) :: compliance_dot(stress_count, stress_count)
    real(dp) :: strain_map_dot(stress_count, displacement_count)
    real(dp) :: divergence_map_dot(equilibrium_count, stress_count)
    real(dp) :: stress_dot(stress_count), displacement_dot(displacement_count)
    real(dp) :: load_dot(equilibrium_count)
    real(dp) :: constitutive_residual(stress_count)
    real(dp) :: equilibrium_residual(equilibrium_count)
    real(dp) :: constitutive_residual_dot(stress_count)
    real(dp) :: equilibrium_residual_dot(equilibrium_count)
    real(dp) :: constitutive_residual_plus(stress_count)
    real(dp) :: equilibrium_residual_plus(equilibrium_count)
    real(dp) :: constitutive_residual_bar(stress_count)
    real(dp) :: equilibrium_residual_bar(equilibrium_count)
    real(dp) :: compliance_bar(stress_count, stress_count)
    real(dp) :: strain_map_bar(stress_count, displacement_count)
    real(dp) :: divergence_map_bar(equilibrium_count, stress_count)
    real(dp) :: stress_bar(stress_count), displacement_bar(displacement_count)
    real(dp) :: load_bar(equilibrium_count)
    real(dp) :: epsilon, finite_difference_error, lhs, rhs
    integer :: status
    logical :: all_passed

    all_passed = .true.
    call build_system( &
        compliance, strain_map, divergence_map, stress, displacement, load)
    call assemble_mixed_elasticity_residual( &
        compliance, strain_map, divergence_map, stress, displacement, load, &
        constitutive_residual, equilibrium_residual, status)
    call record_condition(status == 0, "mixed elasticity residual assembles")
    call record_condition(maxval(abs(constitutive_residual - &
        (matmul(compliance, stress) - matmul(strain_map, displacement)))) &
        < 1.0e-14_dp, &
        "constitutive block matches its matrix oracle")
    call record_condition(maxval(abs(equilibrium_residual - &
        (matmul(divergence_map, stress) - load))) < 1.0e-14_dp, &
        "equilibrium block matches its matrix oracle")

    compliance_dot = reshape([ &
        0.01_dp, -0.02_dp, 0.03_dp, 0.00_dp, 0.02_dp, -0.01_dp, &
        -0.03_dp, 0.01_dp, 0.02_dp], shape(compliance_dot))
    strain_map_dot = reshape([0.02_dp, -0.01_dp, 0.03_dp, 0.01_dp, -0.02_dp, 0.02_dp], &
        shape(strain_map_dot))
    divergence_map_dot = reshape([ &
        0.01_dp, -0.02_dp, 0.02_dp, 0.01_dp, -0.01_dp, 0.03_dp], &
        shape(divergence_map_dot))
    stress_dot = [0.03_dp, -0.02_dp, 0.01_dp]
    displacement_dot = [-0.01_dp, 0.04_dp]
    load_dot = [0.02_dp, -0.03_dp]
    call assemble_mixed_elasticity_residual_jvp( &
        compliance, strain_map, divergence_map, stress, displacement, load, &
        compliance_dot, strain_map_dot, divergence_map_dot, stress_dot, &
        displacement_dot, load_dot, constitutive_residual_dot, &
        equilibrium_residual_dot, status)
    call record_condition(status == 0, "mixed elasticity residual JVP assembles")
    epsilon = 1.0e-7_dp
    call assemble_mixed_elasticity_residual( &
        compliance + epsilon*compliance_dot, strain_map + epsilon*strain_map_dot, &
        divergence_map + epsilon*divergence_map_dot, stress + epsilon*stress_dot, &
        displacement + epsilon*displacement_dot, load + epsilon*load_dot, &
        constitutive_residual_plus, equilibrium_residual_plus, status)
    finite_difference_error = max( &
        maxval(abs(constitutive_residual_dot - &
        (constitutive_residual_plus - constitutive_residual)/epsilon)), &
        maxval(abs(equilibrium_residual_dot - &
        (equilibrium_residual_plus - equilibrium_residual)/epsilon)))
    call record_condition(finite_difference_error < 2.0e-8_dp, &
        "mixed elasticity residual JVP matches a forward difference")

    constitutive_residual_bar = [0.2_dp, -0.3_dp, 0.4_dp]
    equilibrium_residual_bar = [-0.1_dp, 0.5_dp]
    call assemble_mixed_elasticity_residual_vjp( &
        compliance, strain_map, divergence_map, stress, displacement, load, &
        constitutive_residual_bar, equilibrium_residual_bar, compliance_bar, &
        strain_map_bar, divergence_map_bar, stress_bar, displacement_bar, &
        load_bar, status)
    call record_condition(status == 0, "mixed elasticity residual VJP assembles")
    lhs = dot_product(constitutive_residual_bar, constitutive_residual_dot) + &
        dot_product(equilibrium_residual_bar, equilibrium_residual_dot)
    rhs = sum(compliance_bar*compliance_dot) + sum(strain_map_bar*strain_map_dot) + &
        sum(divergence_map_bar*divergence_map_dot) + &
        dot_product(stress_bar, stress_dot) + &
        dot_product(displacement_bar, displacement_dot) + &
        dot_product(load_bar, load_dot)
    call record_condition(abs(lhs - rhs) < 2.0e-13_dp, &
        "mixed elasticity residual VJP satisfies the adjoint identity")

    call assemble_mixed_elasticity_residual( &
        compliance(:stress_count - 1, :), strain_map, divergence_map, stress, &
        displacement, load, constitutive_residual, equilibrium_residual, status)
    call record_condition(status /= 0, &
        "incompatible mixed elasticity dimensions are rejected")

    call check_summary("mixed elasticity residual")
    if (.not. all_passed) error stop 1

contains

    subroutine build_system( &
            compliance, strain_map, divergence_map, stress, displacement, load)
        real(dp), intent(out) :: compliance(:, :), strain_map(:, :)
        real(dp), intent(out) :: divergence_map(:, :)
        real(dp), intent(out) :: stress(:), displacement(:), load(:)

        compliance = reshape([ &
            2.0_dp, 0.1_dp, -0.2_dp, 0.0_dp, 1.5_dp, 0.3_dp, &
            0.4_dp, -0.1_dp, 1.2_dp], shape(compliance))
        strain_map = reshape([0.5_dp, -0.2_dp, 0.1_dp, 0.7_dp, -0.3_dp, 0.4_dp], &
            shape(strain_map))
        divergence_map = reshape([1.0_dp, -0.4_dp, 0.2_dp, -0.3_dp, 0.8_dp, 0.5_dp], &
            shape(divergence_map))
        stress = [0.3_dp, -0.4_dp, 0.2_dp]
        displacement = [0.6_dp, -0.1_dp]
        load = [0.2_dp, -0.5_dp]
    end subroutine build_system

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_mixed_elasticity_residual
