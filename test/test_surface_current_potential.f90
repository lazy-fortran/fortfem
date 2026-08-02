program test_surface_current_potential
    use check, only: check_condition, check_summary
    use fortfem_kinds, only: dp
    use fortsparse, only: fortsparse_status_t
    use fortfem_surface_current_potential, only: &
        surface_current_potential_metadata_t, &
        initialize_surface_current_potential_metadata, &
        validate_surface_current_potential_metadata, &
        assemble_surface_current_potential, &
        assemble_surface_current_potential_jvp, &
        assemble_surface_current_potential_vjp
    implicit none

    integer, parameter :: nq = 3, np = 2, nh = 2
    real(dp), parameter :: normals(nq, 3) = reshape([ &
        0.0_dp, 0.0_dp, 1.0_dp, &
        0.0_dp, 1.0_dp, 0.0_dp, &
        1.0_dp, 0.0_dp, 0.0_dp], [nq, 3])
    real(dp), parameter :: weights(nq) = [1.0_dp, 2.0_dp, 1.5_dp]
    real(dp), parameter :: gradients(nq, np, 3) = reshape([ &
        1.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
        0.0_dp, 0.0_dp, 1.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, &
        0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, 1.0_dp], [nq, np, 3])
    real(dp), parameter :: loop_basis(nq, nh, 3) = reshape([ &
        0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, &
        1.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
        0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, 1.0_dp], [nq, nh, 3])
    real(dp), parameter :: scalar_coefficients(np) = [0.4_dp, -0.7_dp]
    real(dp), parameter :: harmonic_coefficients(nh) = [0.8_dp, -0.3_dp]
    real(dp), parameter :: gradients_dot(nq, np, 3) = reshape([ &
        0.1_dp, 0.1_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
        0.0_dp, 0.0_dp, 0.3_dp, -0.2_dp, 0.0_dp, 0.0_dp, &
        0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.2_dp, -0.1_dp], [nq, np, 3])
    real(dp), parameter :: loop_basis_dot(nq, nh, 3) = reshape([ &
        0.0_dp, 0.2_dp, 0.0_dp, 0.1_dp, 0.0_dp, 0.0_dp, &
        0.1_dp, 0.0_dp, -0.1_dp, -0.1_dp, 0.0_dp, 0.0_dp, &
        0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.1_dp, 0.1_dp], [nq, nh, 3])
    real(dp), parameter :: normals_dot(nq, 3) = 0.0_dp
    real(dp), parameter :: weights_dot(nq) = [0.1_dp, -0.2_dp, 0.05_dp]
    real(dp), parameter :: scalar_coefficients_dot(np) = [0.2_dp, 0.1_dp]
    real(dp), parameter :: harmonic_coefficients_dot(nh) = [-0.2_dp, 0.15_dp]
    real(dp), parameter :: current_bar(nq, 3) = reshape([ &
        0.3_dp, -0.4_dp, 0.2_dp, 0.1_dp, 0.2_dp, -0.5_dp, &
        -0.6_dp, 0.4_dp, 0.7_dp], [nq, 3])
    real(dp), parameter :: integrated_bar(3) = [0.4_dp, -0.3_dp, 0.2_dp]
    real(dp), parameter :: periods_bar(nh) = [0.7_dp, -0.5_dp]
    real(dp), parameter :: eps = 1.0e-6_dp

    type(surface_current_potential_metadata_t) :: metadata
    type(fortsparse_status_t) :: status
    real(dp) :: current(nq, 3), integrated(3), periods(nh)
    real(dp) :: current_dot(nq, 3), integrated_dot(3), periods_dot(nh)
    real(dp) :: current_plus(nq, 3), current_minus(nq, 3)
    real(dp) :: integrated_plus(3), integrated_minus(3)
    real(dp) :: periods_plus(nh), periods_minus(nh)
    real(dp) :: gradients_bar(nq, np, 3), scalar_coefficients_bar(np)
    real(dp) :: loop_basis_bar(nq, nh, 3), harmonic_coefficients_bar(nh)
    real(dp) :: normals_bar(nq, 3), weights_bar(nq)
    real(dp) :: lhs, rhs
    real(dp) :: oracle_current(nq, 3), oracle_integrated(3), oracle_periods(nh)
    real(dp) :: base_gradient(3), base_current(3)
    real(dp) :: bad_weights(nq)
    integer :: q, p, h

    call initialize_surface_current_potential_metadata( &
        metadata, nh, nh, 1, -1, .true., .true., status)
    call check_condition(status%code == 0, &
        "surface-current potential metadata initializes")
    call validate_surface_current_potential_metadata(metadata, status)
    call check_condition(status%code == 0 .and. metadata%orientation_sign == -1 .and. &
        metadata%gauge_dof == 1 .and. metadata%topology_fixed, &
        "surface-current potential metadata validates orientation and gauge")

    call assemble_surface_current_potential( &
        gradients, scalar_coefficients, loop_basis, harmonic_coefficients, &
        normals, weights, metadata, current, integrated, periods, status)
    call check_condition(status%code == 0, &
        "surface-current potential accepts tangent gradients and loop basis")

    oracle_current = 0.0_dp
    oracle_integrated = 0.0_dp
    oracle_periods = 0.0_dp
    do q = 1, nq
        base_gradient = 0.0_dp
        do p = 1, np
            base_gradient = base_gradient + gradients(q, p, :)*scalar_coefficients(p)
        end do
        call cross_product(normals(q, :), base_gradient, base_current)
        do h = 1, nh
            base_current = base_current + loop_basis(q, h, :)*harmonic_coefficients(h)
        end do
        oracle_current(q, :) = metadata%orientation_sign*base_current
        oracle_integrated = oracle_integrated + weights(q)*oracle_current(q, :)
        do h = 1, nh
            oracle_periods(h) = oracle_periods(h) + weights(q)* &
                dot_product(loop_basis(q, h, :), oracle_current(q, :))
        end do
    end do
    call check_condition(maxval(abs(current - oracle_current)) < 1.0e-14_dp .and. &
        maxval(abs(integrated - oracle_integrated)) < 1.0e-14_dp .and. &
        maxval(abs(periods - oracle_periods)) < 1.0e-14_dp, &
        "surface-current potential matches independent current and ledger oracle")

    call assemble_surface_current_potential_jvp( &
        gradients, scalar_coefficients, loop_basis, harmonic_coefficients, &
        normals, weights, metadata, gradients_dot, scalar_coefficients_dot, &
        loop_basis_dot, harmonic_coefficients_dot, normals_dot, weights_dot, &
        current_dot, integrated_dot, periods_dot, status)
    call assemble_surface_current_potential( &
        gradients + eps*gradients_dot, &
        scalar_coefficients + eps*scalar_coefficients_dot, &
        loop_basis + eps*loop_basis_dot, &
        harmonic_coefficients + eps*harmonic_coefficients_dot, &
        normals + eps*normals_dot, weights + eps*weights_dot, metadata, &
        current_plus, integrated_plus, periods_plus, status)
    call assemble_surface_current_potential( &
        gradients - eps*gradients_dot, &
        scalar_coefficients - eps*scalar_coefficients_dot, &
        loop_basis - eps*loop_basis_dot, &
        harmonic_coefficients - eps*harmonic_coefficients_dot, &
        normals - eps*normals_dot, weights - eps*weights_dot, metadata, &
        current_minus, integrated_minus, periods_minus, status)
    call check_condition(status%code == 0 .and. &
        maxval(abs((current_plus - current_minus)/(2.0_dp*eps) - current_dot)) < &
            1.0e-8_dp .and. &
        maxval(abs((integrated_plus - integrated_minus)/(2.0_dp*eps) - &
            integrated_dot)) < &
            1.0e-8_dp .and. &
        maxval(abs((periods_plus - periods_minus)/(2.0_dp*eps) - &
            periods_dot)) < 1.0e-8_dp, &
        "surface-current potential JVP matches centered differences")

    call assemble_surface_current_potential_vjp( &
        gradients, scalar_coefficients, loop_basis, harmonic_coefficients, &
        normals, weights, metadata, current, integrated, periods, current_bar, &
        integrated_bar, periods_bar, gradients_bar, scalar_coefficients_bar, &
        loop_basis_bar, harmonic_coefficients_bar, normals_bar, weights_bar, status)
    lhs = sum(current_bar*current_dot) + dot_product(integrated_bar, integrated_dot) + &
        dot_product(periods_bar, periods_dot)
    rhs = sum(gradients_bar*gradients_dot) + dot_product(scalar_coefficients_bar, &
        scalar_coefficients_dot) + sum(loop_basis_bar*loop_basis_dot) + &
        dot_product(harmonic_coefficients_bar, harmonic_coefficients_dot) + &
        sum(normals_bar*normals_dot) + dot_product(weights_bar, weights_dot)
    call check_condition(status%code == 0 .and. abs(lhs - rhs) < 1.0e-12_dp, &
        "surface-current potential VJP satisfies the real dot-product identity")

    bad_weights = weights
    bad_weights(2) = 0.0_dp
    call assemble_surface_current_potential( &
        gradients, scalar_coefficients, loop_basis, harmonic_coefficients, &
        normals, bad_weights, metadata, current, integrated, periods, status)
    call check_condition(status%code /= 0, &
        "surface-current potential rejects non-positive surface weights")

    call check_summary("surface current potential")

contains

    pure subroutine cross_product(first, second, result)
        real(dp), intent(in) :: first(3), second(3)
        real(dp), intent(out) :: result(3)

        result = [first(2)*second(3) - first(3)*second(2), &
            first(3)*second(1) - first(1)*second(3), &
            first(1)*second(2) - first(2)*second(1)]
    end subroutine cross_product

end program test_surface_current_potential
