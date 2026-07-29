program test_helmholtz_cfie_circle
    use check, only: check_condition, check_summary
    use fortfem_api, only: &
        evaluate_helmholtz_combined_potential_constant, &
        solve_helmholtz_cfie_constant
    use fortfem_kinds, only: dp
    implicit none

    integer, parameter :: mesh_sizes(3) = [16, 32, 64]
    complex(dp), parameter :: exact_scattered = &
        (0.6192993117265686_dp, -0.5636220237912836_dp)
    real(dp) :: errors(3)
    integer :: mesh_id
    logical :: all_passed

    all_passed = .true.
    do mesh_id = 1, size(mesh_sizes)
        call circle_scattering_error(mesh_sizes(mesh_id), errors(mesh_id))
    end do
    call record_condition(errors(1) > 3.0_dp * errors(2) .and. &
        errors(2) > 3.0_dp * errors(3), &
        "Combined-field circle scattering converges near second order")
    call record_condition(errors(3) < 5.0e-3_dp, &
        "Combined-field solution matches the outgoing Mie series")

    call check_summary("Helmholtz combined-field circle scattering")
    if (.not. all_passed) error stop 1

contains

    subroutine circle_scattering_error(panel_count, error)
        integer, intent(in) :: panel_count
        real(dp), intent(out) :: error

        complex(dp), allocatable :: density(:), incident_trace(:)
        real(dp), allocatable :: panel_end(:, :), panel_start(:, :)
        complex(dp) :: scattered(1)
        real(dp) :: angle, midpoint(2), pi, point(2, 1), residual
        integer :: panel, status

        allocate(panel_start(2, panel_count), panel_end(2, panel_count))
        allocate(density(panel_count), incident_trace(panel_count))
        pi = acos(-1.0_dp)
        do panel = 1, panel_count
            angle = 2.0_dp * pi * real(panel - 1, dp) / real(panel_count, dp)
            panel_start(1, panel) = cos(angle)
            panel_start(2, panel) = sin(angle)
            angle = 2.0_dp * pi * real(panel, dp) / real(panel_count, dp)
            panel_end(1, panel) = cos(angle)
            panel_end(2, panel) = sin(angle)
            midpoint = 0.5_dp * ( &
                panel_start(:, panel) + panel_end(:, panel))
            incident_trace(panel) = exp(cmplx(0.0_dp, 1.3_dp * midpoint(1), dp))
        end do

        call solve_helmholtz_cfie_constant( &
            panel_start, panel_end, 1.3_dp, 1.3_dp, incident_trace, 16, &
            density, residual, status)
        call record_condition(status == 0 .and. residual < 2.0e-13_dp, &
            "Combined-field dense solve reaches its algebraic residual")

        point(:, 1) = 2.0_dp * [cos(0.3_dp), sin(0.3_dp)]
        call evaluate_helmholtz_combined_potential_constant( &
            panel_start, panel_end, 1.3_dp, 1.3_dp, density, point, 16, &
            scattered, status)
        call record_condition(status == 0, &
            "Combined layer potential evaluates outside the obstacle")
        error = abs(scattered(1) - exact_scattered)
    end subroutine circle_scattering_error

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_helmholtz_cfie_circle
