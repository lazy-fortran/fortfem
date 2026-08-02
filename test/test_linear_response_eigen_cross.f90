program test_linear_response_eigen_cross
    use, intrinsic :: iso_fortran_env, only: real64
    use check, only: check_condition, check_summary
    use fortfem_interop, only: &
        assemble_linear_response_eigen_cross_residual, &
        assemble_linear_response_eigen_cross_residual_jvp, &
        assemble_linear_response_eigen_cross_residual_vjp, &
        initialize_linear_response_cross_metadata, &
        linear_response_cross_metadata_t, &
        validate_linear_response_cross_metadata
    implicit none

    integer, parameter :: dp = real64, n = 2, p = 2, s = 1
    complex(dp), parameter :: iu = cmplx(0.0_dp, 1.0_dp, dp)
    complex(dp) :: stiffness(n, n), mass(n, n), eigenvalue, eigenstate(n)
    complex(dp) :: response_matrix(p, p), response_coupling(p, n)
    complex(dp) :: response_state(p), response_source(p)
    complex(dp) :: shielding_matrix(s, p), shielding_target(s)
    complex(dp) :: residual(n + p + s), expected(n + p + s)
    complex(dp) :: stiffness_dot(n, n), mass_dot(n, n), eigenvalue_dot
    complex(dp) :: eigenstate_dot(n), response_matrix_dot(p, p)
    complex(dp) :: response_coupling_dot(p, n), response_state_dot(p)
    complex(dp) :: response_source_dot(p), shielding_matrix_dot(s, p)
    complex(dp) :: shielding_target_dot(s), residual_dot(n + p + s)
    complex(dp) :: residual_plus(n + p + s), residual_bar(n + p + s)
    complex(dp) :: stiffness_bar(n, n), mass_bar(n, n), eigenvalue_bar
    complex(dp) :: eigenstate_bar(n), response_matrix_bar(p, p)
    complex(dp) :: response_coupling_bar(p, n), response_state_bar(p)
    complex(dp) :: response_source_bar(p), shielding_matrix_bar(s, p)
    complex(dp) :: shielding_target_bar(s)
    type(linear_response_cross_metadata_t) :: metadata, invalid
    integer :: status
    real(dp) :: epsilon, finite_difference_error, lhs, rhs
    logical :: all_passed

    all_passed = .true.
    call build_inputs(stiffness, mass, eigenvalue, eigenstate, response_matrix, &
        response_coupling, response_state, response_source, shielding_matrix, &
        shielding_target)
    call initialize_linear_response_cross_metadata( &
        metadata, "unit-response", 2.0_dp, "manufactured-cross-oracle", &
        [1, 2], [character(len=64) :: "retained_a", "retained_b"], s, .true., status)
    call record_condition(status == 0, "cross metadata initializes")
    call record_condition(validate_linear_response_cross_metadata(metadata, status), &
        "cross metadata validates")
    call record_condition(metadata%ideal_shielding .and. &
        metadata%retained_response_count == p .and. &
        metadata%shielding_trace_count == s, &
        "retained response and ideal shielding metadata are retained")

    call assemble_linear_response_eigen_cross_residual( &
        stiffness, mass, eigenvalue, eigenstate, response_matrix, response_coupling, &
        response_state, response_source, shielding_matrix, shielding_target, &
        metadata, residual, status)
    call record_condition(status == 0, "cross residual assembles")
    call manual_cross(stiffness, mass, eigenvalue, eigenstate, response_matrix, &
        response_coupling, response_state, response_source, shielding_matrix, &
        shielding_target, expected)
    call record_condition(maxval(abs(residual - expected)) < 2.0e-14_dp, &
        "cross residual matches an independent elementwise matrix oracle")

    call build_tangents(stiffness_dot, mass_dot, eigenvalue_dot, eigenstate_dot, &
        response_matrix_dot, response_coupling_dot, response_state_dot, &
        response_source_dot, shielding_matrix_dot, shielding_target_dot)
    call assemble_linear_response_eigen_cross_residual_jvp( &
        stiffness, mass, eigenvalue, eigenstate, response_matrix, response_coupling, &
        response_state, response_source, shielding_matrix, shielding_target, &
        metadata, stiffness_dot, mass_dot, eigenvalue_dot, eigenstate_dot, &
        response_matrix_dot, response_coupling_dot, response_state_dot, &
        response_source_dot, shielding_matrix_dot, shielding_target_dot, residual_dot, status)
    call record_condition(status == 0, "cross residual JVP assembles")
    epsilon = 1.0e-7_dp
    call assemble_linear_response_eigen_cross_residual( &
        stiffness + epsilon*stiffness_dot, mass + epsilon*mass_dot, &
        eigenvalue + epsilon*eigenvalue_dot, eigenstate + epsilon*eigenstate_dot, &
        response_matrix + epsilon*response_matrix_dot, &
        response_coupling + epsilon*response_coupling_dot, &
        response_state + epsilon*response_state_dot, &
        response_source + epsilon*response_source_dot, &
        shielding_matrix + epsilon*shielding_matrix_dot, &
        shielding_target + epsilon*shielding_target_dot, metadata, residual_plus, status)
    finite_difference_error = maxval(abs(residual_dot - &
        (residual_plus - residual)/epsilon))
    call record_condition(finite_difference_error < 3.0e-8_dp, &
        "cross residual JVP matches a forward difference")

    residual_bar = [ &
        cmplx(0.2_dp, -0.3_dp, dp), cmplx(-0.1_dp, 0.4_dp, dp), &
        cmplx(0.3_dp, 0.2_dp, dp), cmplx(-0.2_dp, -0.1_dp, dp), &
        cmplx(0.5_dp, -0.2_dp, dp)]
    call assemble_linear_response_eigen_cross_residual_vjp( &
        stiffness, mass, eigenvalue, eigenstate, response_matrix, response_coupling, &
        response_state, response_source, shielding_matrix, shielding_target, &
        metadata, residual_bar, stiffness_bar, mass_bar, eigenvalue_bar, &
        eigenstate_bar, response_matrix_bar, response_coupling_bar, response_state_bar, &
        response_source_bar, shielding_matrix_bar, shielding_target_bar, status)
    call record_condition(status == 0, "cross residual VJP assembles")
    lhs = real(sum(conjg(residual_bar)*residual_dot), dp)
    rhs = real(sum(conjg(stiffness_bar)*stiffness_dot) + &
        sum(conjg(mass_bar)*mass_dot) + conjg(eigenvalue_bar)*eigenvalue_dot + &
        sum(conjg(eigenstate_bar)*eigenstate_dot) + &
        sum(conjg(response_matrix_bar)*response_matrix_dot) + &
        sum(conjg(response_coupling_bar)*response_coupling_dot) + &
        sum(conjg(response_state_bar)*response_state_dot) + &
        sum(conjg(response_source_bar)*response_source_dot) + &
        sum(conjg(shielding_matrix_bar)*shielding_matrix_dot) + &
        sum(conjg(shielding_target_bar)*shielding_target_dot), dp)
    call record_condition(abs(lhs - rhs) < 3.0e-13_dp, &
        "cross residual VJP satisfies the real complex adjoint identity")

    invalid = metadata
    invalid%retained_indices(2) = invalid%retained_indices(1)
    call record_condition(.not. validate_linear_response_cross_metadata(invalid, status), &
        "duplicate retained response labels are rejected")

    call check_summary("linear response eigen cross composition")
    if (.not. all_passed) error stop 1

contains

    subroutine build_inputs(stiffness, mass, eigenvalue, eigenstate, response_matrix, &
            response_coupling, response_state, response_source, shielding_matrix, &
            shielding_target)
        complex(dp), intent(out) :: stiffness(:, :), mass(:, :), eigenvalue
        complex(dp), intent(out) :: eigenstate(:), response_matrix(:, :)
        complex(dp), intent(out) :: response_coupling(:, :), response_state(:)
        complex(dp), intent(out) :: response_source(:), shielding_matrix(:, :)
        complex(dp), intent(out) :: shielding_target(:)

        stiffness = reshape([ &
            cmplx(2.1_dp, 0.1_dp, dp), cmplx(-0.2_dp, 0.3_dp, dp), &
            cmplx(0.4_dp, -0.1_dp, dp), cmplx(1.7_dp, -0.2_dp, dp)], shape(stiffness))
        mass = reshape([ &
            cmplx(1.0_dp, 0.0_dp, dp), cmplx(0.1_dp, -0.02_dp, dp), &
            cmplx(-0.05_dp, 0.04_dp, dp), cmplx(0.8_dp, 0.1_dp, dp)], shape(mass))
        eigenvalue = cmplx(0.72_dp, -0.08_dp, dp)
        eigenstate = [cmplx(0.3_dp, -0.2_dp, dp), cmplx(-0.4_dp, 0.5_dp, dp)]
        response_matrix = reshape([ &
            cmplx(1.2_dp, 0.1_dp, dp), cmplx(0.2_dp, -0.05_dp, dp), &
            cmplx(-0.1_dp, 0.07_dp, dp), cmplx(0.9_dp, 0.2_dp, dp)], shape(response_matrix))
        response_coupling = reshape([ &
            cmplx(0.15_dp, -0.02_dp, dp), cmplx(-0.04_dp, 0.08_dp, dp), &
            cmplx(0.05_dp, 0.03_dp, dp), cmplx(0.11_dp, -0.06_dp, dp)], shape(response_coupling))
        response_state = [cmplx(-0.2_dp, 0.1_dp, dp), cmplx(0.4_dp, -0.3_dp, dp)]
        response_source = [cmplx(0.1_dp, -0.2_dp, dp), cmplx(-0.05_dp, 0.08_dp, dp)]
        shielding_matrix = reshape([cmplx(0.6_dp, 0.03_dp, dp), cmplx(-0.2_dp, 0.1_dp, dp)], shape(shielding_matrix))
        shielding_target = [cmplx(0.04_dp, -0.06_dp, dp)]
    end subroutine build_inputs

    subroutine build_tangents(stiffness_dot, mass_dot, eigenvalue_dot, eigenstate_dot, &
            response_matrix_dot, response_coupling_dot, response_state_dot, &
            response_source_dot, shielding_matrix_dot, shielding_target_dot)
        complex(dp), intent(out) :: stiffness_dot(:, :), mass_dot(:, :), eigenvalue_dot
        complex(dp), intent(out) :: eigenstate_dot(:), response_matrix_dot(:, :)
        complex(dp), intent(out) :: response_coupling_dot(:, :), response_state_dot(:)
        complex(dp), intent(out) :: response_source_dot(:), shielding_matrix_dot(:, :)
        complex(dp), intent(out) :: shielding_target_dot(:)

        stiffness_dot = reshape([ &
            cmplx(0.03_dp, -0.01_dp, dp), cmplx(-0.02_dp, 0.04_dp, dp), &
            cmplx(0.01_dp, 0.02_dp, dp), cmplx(-0.03_dp, 0.01_dp, dp)], shape(stiffness_dot))
        mass_dot = reshape([ &
            cmplx(0.01_dp, 0.02_dp, dp), cmplx(-0.03_dp, 0.01_dp, dp), &
            cmplx(0.02_dp, -0.01_dp, dp), cmplx(0.04_dp, 0.00_dp, dp)], shape(mass_dot))
        eigenvalue_dot = cmplx(-0.05_dp, 0.02_dp, dp)
        eigenstate_dot = [cmplx(0.03_dp, 0.01_dp, dp), cmplx(-0.02_dp, 0.04_dp, dp)]
        response_matrix_dot = reshape([ &
            cmplx(0.02_dp, -0.01_dp, dp), cmplx(0.03_dp, 0.02_dp, dp), &
            cmplx(-0.01_dp, 0.04_dp, dp), cmplx(0.01_dp, -0.03_dp, dp)], shape(response_matrix_dot))
        response_coupling_dot = reshape([ &
            cmplx(-0.01_dp, 0.02_dp, dp), cmplx(0.02_dp, -0.01_dp, dp), &
            cmplx(0.03_dp, 0.01_dp, dp), cmplx(-0.02_dp, 0.03_dp, dp)], shape(response_coupling_dot))
        response_state_dot = [cmplx(-0.03_dp, 0.02_dp, dp), cmplx(0.01_dp, -0.04_dp, dp)]
        response_source_dot = [cmplx(0.02_dp, -0.01_dp, dp), cmplx(-0.03_dp, 0.02_dp, dp)]
        shielding_matrix_dot = reshape([cmplx(0.01_dp, 0.02_dp, dp), cmplx(-0.02_dp, 0.01_dp, dp)], shape(shielding_matrix_dot))
        shielding_target_dot = [cmplx(-0.01_dp, 0.03_dp, dp)]
    end subroutine build_tangents

    subroutine manual_cross(stiffness, mass, eigenvalue, eigenstate, response_matrix, &
            response_coupling, response_state, response_source, shielding_matrix, &
            shielding_target, expected)
        complex(dp), intent(in) :: stiffness(:, :), mass(:, :), eigenvalue
        complex(dp), intent(in) :: eigenstate(:), response_matrix(:, :)
        complex(dp), intent(in) :: response_coupling(:, :), response_state(:)
        complex(dp), intent(in) :: response_source(:), shielding_matrix(:, :)
        complex(dp), intent(in) :: shielding_target(:)
        complex(dp), intent(out) :: expected(:)
        integer :: row, column, n_local, p_local

        n_local = size(eigenstate)
        p_local = size(response_state)
        do row = 1, n_local
            expected(row) = -eigenvalue*sum(mass(row, :)*eigenstate)
            do column = 1, n_local
                expected(row) = expected(row) + stiffness(row, column)*eigenstate(column)
            end do
        end do
        do row = 1, p_local
            expected(n_local + row) = -response_source(row)
            do column = 1, p_local
                expected(n_local + row) = expected(n_local + row) + &
                    response_matrix(row, column)*response_state(column)
            end do
            do column = 1, n_local
                expected(n_local + row) = expected(n_local + row) - &
                    response_coupling(row, column)*eigenstate(column)
            end do
        end do
        do row = 1, size(shielding_target)
            expected(n_local + p_local + row) = -shielding_target(row)
            do column = 1, p_local
                expected(n_local + p_local + row) = expected(n_local + p_local + row) + &
                    shielding_matrix(row, column)*response_state(column)
            end do
        end do
    end subroutine manual_cross

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_linear_response_eigen_cross
