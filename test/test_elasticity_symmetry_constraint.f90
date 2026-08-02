program test_elasticity_symmetry_constraint
    use, intrinsic :: iso_fortran_env, only: real64
    use check, only: check_condition, check_summary
    use fortfem_feec, only: &
        assemble_elasticity_symmetry_constraint, &
        assemble_elasticity_symmetry_constraint_jvp, &
        assemble_elasticity_symmetry_constraint_vjp
    implicit none

    integer, parameter :: dp = real64, constraint_count = 2, stress_count = 3
    real(dp) :: symmetry_map(constraint_count, stress_count), stress(stress_count)
    real(dp) :: target(constraint_count)
    real(dp) :: symmetry_map_dot(constraint_count, stress_count)
    real(dp) :: stress_dot(stress_count), target_dot(constraint_count)
    real(dp) :: residual(constraint_count), residual_dot(constraint_count)
    real(dp) :: residual_plus(constraint_count), residual_bar(constraint_count)
    real(dp) :: symmetry_map_bar(constraint_count, stress_count)
    real(dp) :: stress_bar(stress_count)
    real(dp) :: target_bar(constraint_count), epsilon, error, lhs, rhs
    integer :: status
    logical :: all_passed

    all_passed = .true.
    symmetry_map = reshape([ &
        1.0_dp, -0.2_dp, 0.3_dp, 0.4_dp, 0.8_dp, -0.5_dp], shape(symmetry_map))
    stress = [0.3_dp, -0.4_dp, 0.2_dp]
    target = [0.1_dp, -0.2_dp]
    call assemble_elasticity_symmetry_constraint( &
        symmetry_map, stress, target, residual, status)
    call record_condition(status == 0, "elasticity symmetry constraint assembles")
    call record_condition(maxval(abs(residual - &
        (matmul(symmetry_map, stress) - target))) < 1.0e-14_dp, &
        "elasticity symmetry constraint matches its matrix oracle")

    symmetry_map_dot = reshape([ &
        0.02_dp, -0.01_dp, 0.03_dp, -0.04_dp, 0.01_dp, 0.02_dp], &
        shape(symmetry_map_dot))
    stress_dot = [-0.02_dp, 0.03_dp, 0.01_dp]
    target_dot = [0.04_dp, -0.01_dp]
    call assemble_elasticity_symmetry_constraint_jvp( &
        symmetry_map, stress, target, symmetry_map_dot, stress_dot, target_dot, &
        residual_dot, status)
    call record_condition(status == 0, "elasticity symmetry constraint JVP assembles")
    epsilon = 1.0e-7_dp
    call assemble_elasticity_symmetry_constraint( &
        symmetry_map + epsilon*symmetry_map_dot, stress + epsilon*stress_dot, &
        target + epsilon*target_dot, residual_plus, status)
    error = maxval(abs(residual_dot - (residual_plus - residual)/epsilon))
    call record_condition(error < 2.0e-8_dp, &
        "elasticity symmetry constraint JVP matches a forward difference")

    residual_bar = [0.2_dp, -0.6_dp]
    call assemble_elasticity_symmetry_constraint_vjp( &
        symmetry_map, stress, target, residual_bar, symmetry_map_bar, stress_bar, &
        target_bar, status)
    call record_condition(status == 0, "elasticity symmetry constraint VJP assembles")
    lhs = dot_product(residual_bar, residual_dot)
    rhs = sum(symmetry_map_bar*symmetry_map_dot) + &
        dot_product(stress_bar, stress_dot) + dot_product(target_bar, target_dot)
    call record_condition(abs(lhs - rhs) < 2.0e-13_dp, &
        "elasticity symmetry constraint VJP satisfies the adjoint identity")

    call assemble_elasticity_symmetry_constraint( &
        symmetry_map(:constraint_count - 1, :), stress, target, residual, status)
    call record_condition(status /= 0, "incompatible symmetry dimensions are rejected")

    call check_summary("elasticity symmetry constraint")
    if (.not. all_passed) error stop 1

contains

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_elasticity_symmetry_constraint
