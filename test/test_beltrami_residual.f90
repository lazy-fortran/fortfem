program test_beltrami_residual
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
    use check, only: check_condition, check_summary
    use fortfem_feec, only: &
        assemble_beltrami_residual, &
        assemble_beltrami_residual_jvp, &
        assemble_beltrami_residual_vjp
    use fortfem_kinds, only: dp
    use fortsparse, only: fortsparse_status_t
    implicit none

    integer, parameter :: nregion = 2, nsample = 3, ncomponent = 2
    integer, parameter :: ndivrow = 2, ndivpoint = 2, nconstraint = 2
    real(dp) :: curl_field(nregion, nsample, ncomponent)
    real(dp) :: magnetic_field(nregion, nsample, ncomponent), lambda(nregion)
    real(dp) :: divergence(ndivrow, ndivpoint), divergence_target(ndivrow, ndivpoint)
    real(dp) :: flux(nconstraint), flux_target(nconstraint)
    real(dp) :: helicity(nconstraint), helicity_target(nconstraint)
    real(dp) :: curl_dot(nregion, nsample, ncomponent)
    real(dp) :: magnetic_dot(nregion, nsample, ncomponent), lambda_dot(nregion)
    real(dp) :: divergence_dot(ndivrow, ndivpoint)
    real(dp) :: divergence_target_dot(ndivrow, ndivpoint)
    real(dp) :: flux_dot(nconstraint), flux_target_dot(nconstraint)
    real(dp) :: helicity_dot(nconstraint), helicity_target_dot(nconstraint)
    real(dp) :: residual(nregion, nsample, ncomponent)
    real(dp) :: divergence_residual(ndivrow, ndivpoint)
    real(dp) :: flux_residual(nconstraint), helicity_residual(nconstraint)
    real(dp) :: residual_dot(nregion, nsample, ncomponent)
    real(dp) :: divergence_residual_dot(ndivrow, ndivpoint)
    real(dp) :: flux_residual_dot(nconstraint), helicity_residual_dot(nconstraint)
    real(dp) :: residual_plus(nregion, nsample, ncomponent)
    real(dp) :: divergence_residual_plus(ndivrow, ndivpoint)
    real(dp) :: flux_residual_plus(nconstraint), helicity_residual_plus(nconstraint)
    real(dp) :: residual_bar(nregion, nsample, ncomponent)
    real(dp) :: divergence_residual_bar(ndivrow, ndivpoint)
    real(dp) :: flux_residual_bar(nconstraint), helicity_residual_bar(nconstraint)
    real(dp) :: curl_bar(nregion, nsample, ncomponent)
    real(dp) :: magnetic_bar(nregion, nsample, ncomponent), lambda_bar(nregion)
    real(dp) :: divergence_bar(ndivrow, ndivpoint)
    real(dp) :: divergence_target_bar(ndivrow, ndivpoint)
    real(dp) :: flux_bar(nconstraint), flux_target_bar(nconstraint)
    real(dp) :: helicity_bar(nconstraint), helicity_target_bar(nconstraint)
    real(dp) :: epsilon, finite_difference_error, lhs, rhs
    real(dp) :: bad_lambda(nregion)
    real(dp) :: curl_bad(nregion - 1, nsample, ncomponent)
    real(dp) :: div_none(0, 0), flux_none(0), helicity_none(0)
    type(fortsparse_status_t) :: status
    logical :: all_passed

    all_passed = .true.
    call build_data( &
        curl_field, magnetic_field, lambda, divergence, divergence_target, flux, &
        flux_target, helicity, helicity_target)
    call assemble_beltrami_residual( &
        curl_field, magnetic_field, lambda, divergence, divergence_target, flux, &
        flux_target, helicity, helicity_target, residual, divergence_residual, &
        flux_residual, helicity_residual, status)
    call record_condition(status%code == 0, "Beltrami residual assembles")
    call record_condition(maxval(abs(residual - &
        (curl_field - spread(spread(lambda, 2, nsample), 3, ncomponent) * &
        magnetic_field))) < &
        1.0e-14_dp, "curl-eigen residual matches its manufactured oracle")
    call record_condition(maxval(abs(divergence_residual - &
        (divergence - divergence_target))) < 1.0e-14_dp, &
        "divergence rows match their supplied target")
    call record_condition(maxval(abs(flux_residual - (flux - flux_target))) < &
        1.0e-14_dp, "flux rows match their supplied target")
    call record_condition(maxval(abs(helicity_residual - &
        (helicity - helicity_target))) < 1.0e-14_dp, &
        "helicity rows match their supplied target")

    call build_directions( &
        curl_dot, magnetic_dot, lambda_dot, divergence_dot, divergence_target_dot, &
        flux_dot, flux_target_dot, helicity_dot, helicity_target_dot)
    call assemble_beltrami_residual_jvp( &
        curl_field, magnetic_field, lambda, divergence, divergence_target, flux, &
        flux_target, helicity, helicity_target, curl_dot, magnetic_dot, lambda_dot, &
        divergence_dot, divergence_target_dot, flux_dot, flux_target_dot, helicity_dot, &
        helicity_target_dot, residual_dot, divergence_residual_dot, flux_residual_dot, &
        helicity_residual_dot, status)
    call record_condition(status%code == 0, "Beltrami residual JVP assembles")
    epsilon = 1.0e-7_dp
    call assemble_beltrami_residual( &
        curl_field + epsilon*curl_dot, magnetic_field + epsilon*magnetic_dot, &
        lambda + epsilon*lambda_dot, divergence + epsilon*divergence_dot, &
        divergence_target + epsilon*divergence_target_dot, flux + epsilon*flux_dot, &
        flux_target + epsilon*flux_target_dot, helicity + epsilon*helicity_dot, &
        helicity_target + epsilon*helicity_target_dot, residual_plus, &
        divergence_residual_plus, flux_residual_plus, helicity_residual_plus, status)
    finite_difference_error = maxval(abs(residual_dot - &
        (residual_plus - residual)/epsilon))
    finite_difference_error = max(finite_difference_error, maxval(abs( &
        divergence_residual_dot - (divergence_residual_plus - &
        divergence_residual)/epsilon)))
    finite_difference_error = max(finite_difference_error, maxval(abs( &
        flux_residual_dot - (flux_residual_plus - flux_residual)/epsilon)))
    finite_difference_error = max(finite_difference_error, maxval(abs( &
        helicity_residual_dot - (helicity_residual_plus - &
        helicity_residual)/epsilon)))
    call record_condition(finite_difference_error < 2.0e-8_dp, &
        "Beltrami residual JVP matches a forward difference")

    residual_bar = reshape([ &
        0.2_dp, -0.3_dp, 0.4_dp, -0.1_dp, 0.5_dp, -0.2_dp, &
        -0.4_dp, 0.1_dp, 0.3_dp, -0.6_dp, 0.7_dp, 0.2_dp], &
        shape(residual_bar))
    divergence_residual_bar = reshape([0.3_dp, -0.2_dp, 0.4_dp, 0.1_dp], &
        shape(divergence_residual_bar))
    flux_residual_bar = [0.6_dp, -0.5_dp]
    helicity_residual_bar = [-0.4_dp, 0.8_dp]
    call assemble_beltrami_residual_vjp( &
        curl_field, magnetic_field, lambda, divergence, divergence_target, flux, &
        flux_target, helicity, helicity_target, residual_bar, &
        divergence_residual_bar, flux_residual_bar, helicity_residual_bar, curl_bar, &
        magnetic_bar, lambda_bar, divergence_bar, divergence_target_bar, flux_bar, &
        flux_target_bar, helicity_bar, helicity_target_bar, status)
    call record_condition(status%code == 0, "Beltrami residual VJP assembles")
    lhs = sum(residual_bar*residual_dot) + &
        sum(divergence_residual_bar*divergence_residual_dot) + &
        dot_product(flux_residual_bar, flux_residual_dot) + &
        dot_product(helicity_residual_bar, helicity_residual_dot)
    rhs = sum(curl_bar*curl_dot) + sum(magnetic_bar*magnetic_dot) + &
        dot_product(lambda_bar, lambda_dot) + sum(divergence_bar*divergence_dot) + &
        sum(divergence_target_bar*divergence_target_dot) + &
        dot_product(flux_bar, flux_dot) + dot_product(flux_target_bar, flux_target_dot) + &
        dot_product(helicity_bar, helicity_dot) + &
        dot_product(helicity_target_bar, helicity_target_dot)
    call record_condition(abs(lhs - rhs) < 2.0e-13_dp, &
        "Beltrami residual VJP satisfies the real adjoint identity")

    call assemble_beltrami_residual( &
        curl_bad, magnetic_field, lambda, divergence, divergence_target, flux, &
        flux_target, helicity, helicity_target, residual, divergence_residual, &
        flux_residual, helicity_residual, status)
    call record_condition(status%code /= 0, &
        "incompatible region dimensions are rejected")
    bad_lambda = lambda
    bad_lambda(1) = ieee_value(0.0_dp, ieee_quiet_nan)
    call assemble_beltrami_residual( &
        curl_field, magnetic_field, bad_lambda, divergence, divergence_target, flux, &
        flux_target, helicity, helicity_target, residual, divergence_residual, &
        flux_residual, helicity_residual, status)
    call record_condition(status%code /= 0, "non-finite Beltrami inputs are rejected")

    call test_zero_sized_constraint_rows(curl_field, magnetic_field, lambda)
    call check_summary("Beltrami residual")
    if (.not. all_passed) error stop 1

contains

    subroutine build_data( &
            curl_field, magnetic_field, lambda, divergence, divergence_target, flux, &
            flux_target, helicity, helicity_target)
        real(dp), intent(out) :: curl_field(:, :, :), magnetic_field(:, :, :), lambda(:)
        real(dp), intent(out) :: divergence(:, :), divergence_target(:, :)
        real(dp), intent(out) :: flux(:), flux_target(:), helicity(:), helicity_target(:)

        curl_field = reshape([ &
            0.8_dp, -0.2_dp, 0.4_dp, 1.1_dp, -0.5_dp, 0.3_dp, &
            0.7_dp, -0.1_dp, 0.6_dp, 1.2_dp, -0.4_dp, 0.2_dp], &
            shape(curl_field))
        magnetic_field = reshape([ &
            0.3_dp, 0.7_dp, -0.4_dp, 0.2_dp, 0.5_dp, -0.6_dp, &
            -0.1_dp, 0.8_dp, 0.9_dp, -0.3_dp, 0.4_dp, 0.2_dp], &
            shape(magnetic_field))
        lambda = [0.6_dp, -0.35_dp]
        divergence = reshape([0.2_dp, -0.4_dp, 0.7_dp, -0.1_dp], shape(divergence))
        divergence_target = reshape([ &
            -0.1_dp, 0.3_dp, 0.5_dp, -0.2_dp], shape(divergence_target))
        flux = [0.9_dp, -0.6_dp]
        flux_target = [0.2_dp, -0.4_dp]
        helicity = [-0.7_dp, 0.8_dp]
        helicity_target = [0.1_dp, 0.3_dp]
    end subroutine build_data

    subroutine build_directions( &
            curl_dot, magnetic_dot, lambda_dot, divergence_dot, divergence_target_dot, &
            flux_dot, flux_target_dot, helicity_dot, helicity_target_dot)
        real(dp), intent(out) :: curl_dot(:, :, :), magnetic_dot(:, :, :), lambda_dot(:)
        real(dp), intent(out) :: divergence_dot(:, :), divergence_target_dot(:, :)
        real(dp), intent(out) :: flux_dot(:), flux_target_dot(:), helicity_dot(:)
        real(dp), intent(out) :: helicity_target_dot(:)

        curl_dot = reshape([ &
            0.02_dp, -0.03_dp, 0.01_dp, 0.04_dp, -0.01_dp, 0.02_dp, &
            0.05_dp, -0.02_dp, 0.03_dp, -0.04_dp, 0.06_dp, 0.01_dp], &
            shape(curl_dot))
        magnetic_dot = reshape([ &
            -0.01_dp, 0.04_dp, 0.02_dp, -0.03_dp, 0.05_dp, 0.01_dp, &
            0.03_dp, -0.02_dp, 0.04_dp, 0.01_dp, -0.05_dp, 0.02_dp], &
            shape(magnetic_dot))
        lambda_dot = [0.03_dp, -0.02_dp]
        divergence_dot = reshape([ &
            0.02_dp, -0.01_dp, 0.04_dp, -0.03_dp], shape(divergence_dot))
        divergence_target_dot = reshape([ &
            -0.03_dp, 0.01_dp, 0.02_dp, 0.05_dp], &
            shape(divergence_target_dot))
        flux_dot = [0.04_dp, -0.02_dp]
        flux_target_dot = [0.01_dp, 0.03_dp]
        helicity_dot = [-0.02_dp, 0.05_dp]
        helicity_target_dot = [0.03_dp, -0.04_dp]
    end subroutine build_directions

    subroutine test_zero_sized_constraint_rows(curl_field, magnetic_field, lambda)
        real(dp), intent(in) :: curl_field(:, :, :), magnetic_field(:, :, :), lambda(:)
        real(dp) :: residual_local(size(curl_field, 1), size(curl_field, 2), &
            size(curl_field, 3))
        real(dp) :: divergence_local(0, 0), flux_local(0), helicity_local(0)
        real(dp) :: divergence_target_local(0, 0), flux_target_local(0)
        real(dp) :: helicity_target_local(0)
        real(dp) :: divergence_residual_local(0, 0), flux_residual_local(0)
        real(dp) :: helicity_residual_local(0)
        type(fortsparse_status_t) :: local_status

        call assemble_beltrami_residual( &
            curl_field, magnetic_field, lambda, divergence_local, &
            divergence_target_local, flux_local, flux_target_local, helicity_local, &
            helicity_target_local, residual_local, divergence_residual_local, &
            flux_residual_local, helicity_residual_local, local_status)
        call record_condition(local_status%code == 0, &
            "zero-sized constraint rows are accepted as optional")
    end subroutine test_zero_sized_constraint_rows

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_beltrami_residual
