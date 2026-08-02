program test_eulerian_nonnested_residual
    use check, only: check_condition, check_summary
    use fortfem_api, only: &
        assemble_eulerian_nonnested_residual, &
        assemble_eulerian_nonnested_residual_jvp, &
        assemble_eulerian_nonnested_residual_vjp, &
        CONTINUATION_EVENT_NONE, CONTINUATION_EVENT_SIGN_CROSSING
    use fortfem_kinds, only: dp
    use fortsparse, only: FORTSPARSE_INVALID_MATRIX, FORTSPARSE_OK, &
        fortsparse_status_t
    implicit none

    logical :: all_passed

    all_passed = .true.
    call check_contract()
    call check_summary("Eulerian non-nested residual")
    if (.not. all_passed) error stop 1

contains

    subroutine check_contract()
        real(dp), parameter :: epsilon_fd = 1.0e-7_dp
        integer, parameter :: force_count = 2, divergence_count = 3
        integer, parameter :: total_count = force_count + divergence_count
        real(dp) :: force(total_count-divergence_count), divergence(divergence_count)
        real(dp) :: stabilization(total_count), force_dot(force_count)
        real(dp) :: divergence_dot(divergence_count), stabilization_dot(total_count)
        real(dp) :: residual(total_count), residual_dot(total_count)
        real(dp) :: residual_plus(total_count), residual_minus(total_count)
        real(dp) :: residual_bar(total_count), force_bar(force_count)
        real(dp) :: divergence_bar(divergence_count), stabilization_bar(total_count)
        real(dp) :: previous_margin(2), current_margin(2), minimum_margin
        real(dp) :: lhs, rhs
        integer :: event_code, event_index
        type(fortsparse_status_t) :: status

        force = [1.2_dp, -0.7_dp]
        divergence = [0.4_dp, -0.2_dp, 0.8_dp]
        stabilization = [-0.1_dp, 0.3_dp, 0.2_dp, -0.4_dp, 0.5_dp]
        force_dot = [-0.3_dp, 0.6_dp]
        divergence_dot = [0.7_dp, -0.4_dp, 0.1_dp]
        stabilization_dot = [0.2_dp, -0.5_dp, 0.3_dp, 0.4_dp, -0.6_dp]
        previous_margin = [0.5_dp, -0.3_dp]
        current_margin = [0.2_dp, 0.4_dp]

        call assemble_eulerian_nonnested_residual( &
            force, divergence, residual, status, stabilization, previous_margin, &
            current_margin, 0.05_dp, event_code, event_index, minimum_margin)
        call record(status%code == FORTSPARSE_OK .and. &
            maxval(abs(residual - ([force, divergence] + stabilization))) < &
            1.0e-14_dp .and. &
            event_code == CONTINUATION_EVENT_SIGN_CROSSING .and. &
            event_index == 2 .and. &
            abs(minimum_margin - 0.2_dp) < 1.0e-14_dp, &
            "value composition and topology event oracle")

        call assemble_eulerian_nonnested_residual_jvp( &
            force, divergence, force_dot, divergence_dot, residual_dot, status, &
            stabilization, stabilization_dot, previous_margin, current_margin, &
            0.05_dp, &
            event_code, event_index, minimum_margin)
        call assemble_eulerian_nonnested_residual( &
            force + epsilon_fd*force_dot, divergence + epsilon_fd*divergence_dot, &
            residual_plus, status, stabilization + epsilon_fd*stabilization_dot)
        call assemble_eulerian_nonnested_residual( &
            force - epsilon_fd*force_dot, divergence - epsilon_fd*divergence_dot, &
            residual_minus, status, stabilization - epsilon_fd*stabilization_dot)
        call record(status%code == FORTSPARSE_OK .and. &
            maxval(abs(residual_dot - (residual_plus - residual_minus) / &
            (2.0_dp*epsilon_fd))) < 2.0e-8_dp, "JVP finite-difference oracle")

        residual_bar = [0.2_dp, -0.4_dp, 0.6_dp, -0.1_dp, 0.3_dp]
        call assemble_eulerian_nonnested_residual_vjp( &
            force, divergence, residual_bar, force_bar, divergence_bar, status, &
            stabilization, stabilization_bar)
        lhs = dot_product(residual_bar, residual_dot)
        rhs = dot_product(force_bar, force_dot) + &
            dot_product(divergence_bar, divergence_dot) + &
            dot_product(stabilization_bar, stabilization_dot)
        call record(status%code == FORTSPARSE_OK .and. abs(lhs-rhs) < 1.0e-14_dp, &
            "VJP dot-product oracle")

        call assemble_eulerian_nonnested_residual( &
            force, divergence, residual, status, [1.0_dp, 2.0_dp])
        call record(status%code == FORTSPARSE_INVALID_MATRIX, &
            "incompatible stabilization length is rejected")

        call assemble_eulerian_nonnested_residual( &
            force, divergence, residual, status, event_tolerance=-1.0_dp)
        call record(status%code == FORTSPARSE_INVALID_MATRIX, &
            "negative topology tolerance is rejected")

        call assemble_eulerian_nonnested_residual( &
            force, divergence, residual, status, event_code=event_code, &
            event_index=event_index, minimum_margin=minimum_margin)
        call record(status%code == FORTSPARSE_OK .and. &
            event_code == CONTINUATION_EVENT_NONE .and. &
            event_index == 0 .and. minimum_margin == 0.0_dp, &
            "event metadata is optional")
    end subroutine check_contract

    subroutine record(condition, message)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: message

        call check_condition(condition, message)
        all_passed = all_passed .and. condition
    end subroutine record

end program test_eulerian_nonnested_residual
