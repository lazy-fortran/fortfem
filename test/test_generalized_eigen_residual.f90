program test_generalized_eigen_residual
    use, intrinsic :: iso_fortran_env, only: real64
    use check, only: check_condition, check_summary
    use fortfem_api, only: &
        assemble_generalized_eigen_residual, &
        assemble_generalized_eigen_residual_jvp, &
        assemble_generalized_eigen_residual_vjp
    implicit none

    integer, parameter :: dp = real64, n = 3
    complex(dp) :: stiffness(n, n), mass(n, n), state(n), eigenvalue
    complex(dp) :: stiffness_dot(n, n), mass_dot(n, n), state_dot(n)
    complex(dp) :: eigenvalue_dot, residual(n), residual_dot(n)
    complex(dp) :: stiffness_plus(n, n), mass_plus(n, n), state_plus(n)
    complex(dp) :: residual_plus(n), residual_bar(n), stiffness_bar(n, n)
    complex(dp) :: mass_bar(n, n), state_bar(n), eigenvalue_bar
    integer :: status
    real(dp) :: epsilon, finite_difference_error, lhs, rhs
    logical :: all_passed

    all_passed = .true.
    call build_system(stiffness, mass, state, eigenvalue)
    call assemble_generalized_eigen_residual( &
        stiffness, mass, eigenvalue, state, residual, status)
    call record_condition(status == 0, "generalized eigen residual assembles")
    call record_condition(maxval(abs(residual - &
        (matmul(stiffness, state) - eigenvalue*matmul(mass, state)))) < 1.0e-14_dp, &
        "generalized eigen residual matches the matrix oracle")

    stiffness_dot = cmplx(0.0_dp, 0.0_dp, dp)
    mass_dot = cmplx(0.0_dp, 0.0_dp, dp)
    stiffness_dot(1, 2) = cmplx(0.04_dp, -0.02_dp, dp)
    stiffness_dot(3, 1) = cmplx(-0.03_dp, 0.05_dp, dp)
    mass_dot(2, 3) = cmplx(0.02_dp, 0.01_dp, dp)
    state_dot = [ &
        cmplx(0.03_dp, -0.01_dp, dp), cmplx(-0.02_dp, 0.04_dp, dp), &
        cmplx(0.01_dp, 0.02_dp, dp)]
    eigenvalue_dot = cmplx(-0.05_dp, 0.03_dp, dp)
    call assemble_generalized_eigen_residual_jvp( &
        stiffness, mass, eigenvalue, state, stiffness_dot, mass_dot, &
        eigenvalue_dot, state_dot, residual_dot, status)
    call record_condition(status == 0, "generalized eigen residual JVP assembles")
    epsilon = 1.0e-7_dp
    call assemble_generalized_eigen_residual( &
        stiffness + epsilon*stiffness_dot, mass + epsilon*mass_dot, &
        eigenvalue + epsilon*eigenvalue_dot, state + epsilon*state_dot, &
        residual_plus, status)
    finite_difference_error = maxval(abs(residual_dot - &
        (residual_plus - residual)/epsilon))
    call record_condition(finite_difference_error < 2.0e-8_dp, &
        "generalized eigen residual JVP matches a forward difference")

    residual_bar = [ &
        cmplx(0.2_dp, -0.4_dp, dp), cmplx(-0.3_dp, 0.1_dp, dp), &
        cmplx(0.5_dp, 0.2_dp, dp)]
    call assemble_generalized_eigen_residual_vjp( &
        stiffness, mass, eigenvalue, state, residual_bar, stiffness_bar, &
        mass_bar, state_bar, eigenvalue_bar, status)
    call record_condition(status == 0, "generalized eigen residual VJP assembles")
    lhs = real(sum(conjg(residual_bar)*residual_dot), dp)
    rhs = real(sum(conjg(stiffness_bar)*stiffness_dot) + &
        sum(conjg(mass_bar)*mass_dot) + sum(conjg(state_bar)*state_dot) + &
        conjg(eigenvalue_bar)*eigenvalue_dot, dp)
    call record_condition(abs(lhs - rhs) < 2.0e-13_dp, &
        "generalized eigen residual VJP satisfies the real complex adjoint identity")

    call assemble_generalized_eigen_residual( &
        stiffness(:n - 1, :), mass, eigenvalue, state, residual, status)
    call record_condition(status /= 0, &
        "incompatible generalized eigen dimensions are rejected")

    call check_summary("generalized eigen residual")
    if (.not. all_passed) error stop 1

contains

    subroutine build_system(stiffness, mass, state, eigenvalue)
        complex(dp), intent(out) :: stiffness(:, :), mass(:, :), state(:), eigenvalue

        stiffness = reshape([ &
            cmplx(2.0_dp, 0.1_dp, dp), cmplx(0.2_dp, -0.3_dp, dp), &
            cmplx(0.1_dp, 0.4_dp, dp), cmplx(-0.1_dp, 0.2_dp, dp), &
            cmplx(1.5_dp, -0.2_dp, dp), cmplx(0.3_dp, 0.1_dp, dp), &
            cmplx(0.4_dp, -0.1_dp, dp), cmplx(0.2_dp, 0.3_dp, dp), &
            cmplx(1.2_dp, 0.0_dp, dp)], shape(stiffness))
        mass = reshape([ &
            cmplx(1.0_dp, 0.0_dp, dp), cmplx(0.1_dp, 0.02_dp, dp), &
            cmplx(0.0_dp, -0.1_dp, dp), cmplx(-0.1_dp, 0.03_dp, dp), &
            cmplx(1.2_dp, 0.1_dp, dp), cmplx(0.2_dp, 0.0_dp, dp), &
            cmplx(0.05_dp, 0.04_dp, dp), cmplx(-0.2_dp, 0.01_dp, dp), &
            cmplx(0.8_dp, -0.2_dp, dp)], shape(mass))
        state = [ &
            cmplx(0.3_dp, -0.4_dp, dp), cmplx(-0.2_dp, 0.1_dp, dp), &
            cmplx(0.5_dp, 0.2_dp, dp)]
        eigenvalue = cmplx(0.7_dp, -0.2_dp, dp)
    end subroutine build_system

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_generalized_eigen_residual
